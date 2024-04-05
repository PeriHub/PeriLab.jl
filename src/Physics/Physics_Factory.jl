# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Physics
include("../Support/helpers.jl")
using Reexport
@reexport using .Helpers: check_inf_or_nan, find_active, get_active_update_nodes
include("./Additive/Additive_Factory.jl")
include("./Corrosion/Corrosion_Factory.jl")
include("./Damage/Damage_Factory.jl")
include("./Material/Material_Factory.jl")
include("./Thermal/Thermal_Factory.jl")
include("./Pre_calculation/Pre_Calculation_Factory.jl")
include("../Support/Parameters/parameter_handling.jl")
@reexport using .Parameter_Handling: get_physics_option, get_model_parameter, get_heat_capacity
# in future FEM will be outside of the Physics_Factory
include("../FEM/FEM_Factory.jl")

using .Additive
using .Corrosion
using .Damage
using .Material
using .Pre_calculation
using .Thermal
# in future FEM will be outside of the Physics_Factory
using .FEM
using TimerOutputs
export compute_models
export init_models
export read_properties
export init_additive_model_fields
export init_corrosion_model_fields
export init_damage_model_fields
export init_material_model_fields
export init_thermal_model_fields

"""
    compute_models(datamanager::Module, block_nodes::Dict{Int64,Vector{Int64}}, dt::Float64, time::Float64, options::Dict, synchronise_field, to::TimerOutput)

Computes the physics models

# Arguments
- `datamanager::Module`: The datamanager
- `block_nodes::Dict{Int64,Vector{Int64}}`: The block nodes
- `dt::Float64`: The time step
- `time::Float64`: The current time
- `options::Dict`: The options
- `synchronise_field`: The synchronise field
- `to::TimerOutput`: The timer output
# Returns
- `datamanager`: The datamanager
"""
function compute_models(datamanager::Module, block_nodes::Dict{Int64,Vector{Int64}}, dt::Float64, time::Float64, options::Dict, synchronise_field, to::TimerOutput)
    fem_option = datamanager.fem_active()
    if fem_option
        fe_nodes = datamanager.get_field("FE Nodes")
    end
    if options["Additive Models"]
        for block in eachindex(block_nodes)
            if datamanager.check_property(block, "Additive Model")
                @timeit to "compute_additive_model_model" datamanager = Additive.compute_additive_model(datamanager, block_nodes[block], datamanager.get_properties(block, "Additive Model"), time, dt)
            end
        end
    end
    active = datamanager.get_field("Active")
    if options["Damage Models"]
        #tbd damage specific pre_calculation-> in damage template
        for block in eachindex(block_nodes)
            nodes = block_nodes[block]

            active_nodes = view(nodes, find_active(active[nodes]))
            if fem_option
                active_nodes = view(nodes, find_active(.~fe_nodes[active_nodes]))
            end
            if datamanager.check_property(block, "Damage Model") && datamanager.check_property(block, "Material Model")
                datamanager = Damage.set_bond_damage(datamanager, active_nodes)
                @timeit to "damage_pre_calculation" datamanager = compute_damage_pre_calculation(datamanager, options, active_nodes, block, synchronise_field, time, dt, to)
                @timeit to "damage" datamanager = Damage.compute_damage(datamanager, active_nodes, datamanager.get_properties(block, "Damage Model"), block, time, dt)
            end
        end
    end
    update_list = datamanager.get_field("Update List")
    # nodes::Vector{Int64} = []
    # active_nodes::Vector{Int64} = []
    # update_nodes::Vector{Int64} = []
    if fem_option
        nelements = datamanager.get_num_elements()
        # in future the FE analysis can be put into the block loop. Right now it is outside the block.
        datamanager = FEM.eval(datamanager, Vector{Int64}(1:nelements), datamanager.get_properties(1, "FEM"), time, dt)
    end
    for block in eachindex(block_nodes)

        active_nodes, update_nodes = get_active_update_nodes(active, update_list, block_nodes, block)

        if fem_option
            update_nodes = block_nodes[block][find_active(Vector{Bool}(.~fe_nodes[update_nodes]))]
        end
        @timeit to "pre_calculation" datamanager = Pre_calculation.compute(datamanager, update_nodes, datamanager.get_physics_options(), time, dt, to)

        if options["Thermal Models"]
            if datamanager.check_property(block, "Thermal Model")
                @timeit to "compute_thermal_model" datamanager = Thermal.compute_thermal_model(datamanager, update_nodes, datamanager.get_properties(block, "Thermal Model"), time, dt)
            end
        end

        if options["Material Models"]
            if datamanager.check_property(block, "Material Model")
                model_param = datamanager.get_properties(block, "Material Model")
                @timeit to "bond_forces" datamanager = Material.compute_forces(datamanager, update_nodes, model_param, time, dt, to)
                #TODO: I think this needs to stay here as we need the active_nodes not the update_nodes
                @timeit to "distribute_force_densities" datamanager = Material.distribute_force_densities(datamanager, active_nodes)
                if haskey(model_param, "von Mises stresses") && model_param["von Mises stresses"]
                    datamanager = Material.calculate_von_mises_stress(datamanager, active_nodes)
                end
            end
        end
    end
    update_list .= true
    return datamanager

end

"""
    compute_damage_pre_calculation(datamanager::Module, options::Dict, nodes::Union{SubArray,Vector{Int64}}, block::Int64, synchronise_field, time::Float64, dt::Float64)

Compute the damage pre calculation

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `block::Int64`: Block number
- `synchronise_field`: Synchronise function to distribute parameter through cores.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
- `to::TimerOutput`: The timer output.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function compute_damage_pre_calculation(datamanager::Module, options::Dict, nodes::Union{SubArray,Vector{Int64}}, block::Int64, synchronise_field, time::Float64, dt::Float64, to::TimerOutput)

    @timeit to "compute" datamanager = Pre_calculation.compute(datamanager, nodes, datamanager.get_physics_options(), time, dt, to)

    if options["Thermal Models"]
        @timeit to "thermal_model" datamanager = Thermal.compute_thermal_model(datamanager, nodes, datamanager.get_properties(block, "Thermal Model"), time, dt)
    end

    if options["Material Models"]
        @timeit to "compute_forces" datamanager = Material.compute_forces(datamanager, nodes, datamanager.get_properties(block, "Material Model"), time, dt, to)
    end
    datamanager = Damage.compute_damage_pre_calculation(datamanager, nodes, block, datamanager.get_properties(block, "Damage Model"), synchronise_field, time, dt)
    update_list = datamanager.get_field("Update List")
    update_list[nodes] .= false
    return datamanager
end

"""
    get_block_model_definition(params::Dict, block_id::Int64, prop_keys::Vector{String}, properties)

Get block model definition

# Arguments
- `params::Dict`: Parameters.
- `block_id::Int64`: Block id.
- `prop_keys::Vector{String}`: Property keys.
- `properties`: Properties function.
# Returns
- `properties`: Properties function.    
"""
function get_block_model_definition(params::Dict, block_id::Int64, prop_keys::Vector{String}, properties)
    # properties function from datamanager
    if haskey(params["Blocks"], "block_" * string(block_id))
        block = params["Blocks"]["block_"*string(block_id)]
        for model in prop_keys
            if haskey(block, model)
                properties(block_id, model, get_model_parameter(params, model, block[model]))
            end
        end
    end
    return properties
end

"""
    init_models(params::Dict, datamanager::Module, block_nodes::Dict{Int64,Vector{Int64}}, solver_options::Dict)

Initialize models

# Arguments
- `params::Dict`: Parameters.
- `datamanager::Module`: Datamanager.
- `block_nodes::Dict{Int64,Vector{Int64}}`: block nodes.
- `solver_options::Dict`: Solver options.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function init_models(params::Dict, datamanager::Module, block_nodes::Dict{Int64,Vector{Int64}}, solver_options::Dict, to::TimerOutput)
    dof = datamanager.get_dof()
    deformed_coorN, deformed_coorNP1 = datamanager.create_node_field("Deformed Coordinates", Float64, dof)
    deformed_coorN = copy(datamanager.get_field("Coordinates"))
    deformed_coorNP1 = copy(datamanager.get_field("Coordinates"))
    datamanager.create_node_field("Displacements", Float64, dof)
    if solver_options["Additive Models"]
        @info "Init additive models"
        @timeit to "additive_model_fields" datamanager = Physics.init_additive_model_fields(datamanager)
        heat_capacity = datamanager.create_constant_node_field("Specific Heat Capacity", Float64, 1)
        heat_capacity = set_heat_capacity(params, block_nodes, heat_capacity) # includes the neighbors
        for block in eachindex(block_nodes)
            if datamanager.check_property(block, "Additive Model")
                @timeit to "init additive model" datamanager = Additive.init_additive_model(datamanager, block_nodes[block], block)
            end
        end
    end
    if solver_options["Corrosion Models"]
        @info "Init corrosion models"
        @timeit to "corrosion_model_fields" datamanager = Corrosion.init_corrosion_model_fields(datamanager)

        for block in eachindex(block_nodes)
            if datamanager.check_property(block, "Corrosion Model")
                @timeit to "init corrosion model" datamanager = Corrosion.init_corrosion_model(datamanager, block_nodes[block], block)
            end
        end
    end
    if solver_options["Damage Models"]
        @info "Init damage models"
        @timeit to "damage_model_fields" datamanager = Damage.init_damage_model_fields(datamanager, params)
        for block in eachindex(block_nodes)
            if datamanager.check_property(block, "Damage Model")
                @timeit to "init damage model" datamanager = Damage.init_damage_model(datamanager, block_nodes[block], block)
            end
        end

    end
    if solver_options["Material Models"]
        @info "Init material models"
        @timeit to "material model fields" datamanager = Material.init_material_model_fields(datamanager)
        for block in eachindex(block_nodes)
            if datamanager.check_property(block, "Material Model")
                @timeit to "init material" datamanager = Material.init_material_model(datamanager, block_nodes[block], block)
            end
        end
    end
    if solver_options["Thermal Models"]
        @info "Init thermal models"
        @timeit to "thermal model fields" datamanager = Thermal.init_thermal_model_fields(datamanager)
        heat_capacity = datamanager.create_constant_node_field("Specific Heat Capacity", Float64, 1)
        heat_capacity = set_heat_capacity(params, block_nodes, heat_capacity) # includes the neighbors
        for block in eachindex(block_nodes)
            if datamanager.check_property(block, "Thermal Model")
                @timeit to "init thermal model" datamanager = Thermal.init_thermal_model(datamanager, block_nodes[block], block)
            end
        end
    end
    @debug "Init pre calculation models"
    return init_pre_calculation(datamanager, datamanager.get_physics_options())
end


"""
    init_pre_calculation(datamanager::Module, options)

Initialize pre-calculation

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `options`: Options.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function init_pre_calculation(datamanager::Module, options)
    return Pre_calculation.init_pre_calculation(datamanager, options)
end

"""
    read_properties(params::Dict, datamanager::Module, material_model::Bool)

Read properties

# Arguments
- `params::Dict`: Parameters.
- `datamanager::Data_manager`: Datamanager.
- `material_model::Bool`: Material model.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function read_properties(params::Dict, datamanager::Module, material_model::Bool)
    datamanager.init_property()
    blocks = datamanager.get_block_list()
    prop_keys = datamanager.init_property()
    physics_options = datamanager.get_physics_options()

    datamanager.set_physics_options(get_physics_option(params, physics_options))

    for block in blocks
        get_block_model_definition(params, block, prop_keys, datamanager.set_properties)
    end
    if material_model
        dof = datamanager.get_dof()
        for block in blocks
            Material.check_material_symmetry(dof, datamanager.get_properties(block, "Material Model"))
            Material.determine_isotropic_parameter(datamanager, datamanager.get_properties(block, "Material Model"))
        end
    end
    return datamanager
end

"""
    set_heat_capacity(params::Dict, block_nodes::Dict, heat_capacity::SubArray)

Sets the heat capacity of the nodes in the dictionary.

# Arguments
- `params::Dict`: The parameters
- `block_nodes::Dict`: The block nodes
- `heat_capacity::SubArray`: The heat capacity array
# Returns
- `heat_capacity::SubArray`: The heat capacity array
"""
function set_heat_capacity(params::Dict, block_nodes::Dict, heat_capacity::SubArray)
    for block in eachindex(block_nodes)
        heat_capacity[block_nodes[block]] .= get_heat_capacity(params, block)
    end
    return heat_capacity
end

end