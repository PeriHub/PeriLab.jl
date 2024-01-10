# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Physics
include("../Support/helpers.jl")
include("./Additive/Additive_Factory.jl")
include("./Damage/Damage_Factory.jl")
include("./Material/Material_Factory.jl")
include("./Thermal/Thermal_Factory.jl")
include("./Pre_calculation/Pre_Calculation_Factory.jl")
include("../Support/Parameters/parameter_handling.jl")
# in future FEM will be outside of the Physics_Factory
include("../FEM/FEM_Factory.jl")

using .Additive
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

            active_nodes = nodes[find_active(active[nodes])]
            if fem_option
                active_nodes = nodes[find_active(.~fe_nodes[active_nodes])]
            end
            if datamanager.check_property(block, "Damage Model") && datamanager.check_property(block, "Material Model")
                datamanager = Damage.set_bond_damage(datamanager, active_nodes)
                @timeit to "damage_pre_calculation" datamanager = compute_damage_pre_calculation(datamanager, options, active_nodes, block, synchronise_field, time, dt, to)
                @timeit to "damage" datamanager = Damage.compute_damage(datamanager, active_nodes, datamanager.get_properties(block, "Damage Model"), block, time, dt)
            end
        end
    end
    update_list = datamanager.get_field("Update List")
    nodes::Vector{Int64} = []
    active_nodes::Vector{Int64} = []
    update_nodes::Vector{Int64} = []
    if fem_option
        nelements = datamanager.get_num_elements()
        # in future the FE analysis can be put into the block loop. Right now it is outside the block.
        datamanager = FEM.eval(datamanager, Vector{Int64}(1:nelements), datamanager.get_properties(1, "FEM"), time, dt)
    end
    for block in eachindex(block_nodes)
        nodes = block_nodes[block]
        active_index = find_active(active[nodes])
        active_nodes = nodes[active_index]
        update_nodes = view(nodes, find_updatable(active_index, update_list))
        if fem_option
            update_nodes = nodes[find_active(Vector{Bool}(.~fe_nodes[update_nodes]))]
        end
        @timeit to "pre_calculation" datamanager = Pre_calculation.compute(datamanager, update_nodes, datamanager.get_physics_options(), time, dt, to)

        if options["Thermal Models"]
            if datamanager.check_property(block, "Thermal Model")
                @timeit to "compute_thermal_model" datamanager = Thermal.compute_thermal_model(datamanager, update_nodes, datamanager.get_properties(block, "Thermal Model"), time, dt)
            end
        end

        if options["Material Models"]
            if datamanager.check_property(block, "Material Model")
                @timeit to "bond_forces" datamanager = Material.compute_forces(datamanager, update_nodes, datamanager.get_properties(block, "Material Model"), time, dt, to)
                @timeit to "distribute_force_densities" datamanager = Material.distribute_force_densities(datamanager, active_nodes)
            end
        end
    end
    if options["Additive Models"]
        for block in eachindex(block_nodes)
            if datamanager.check_property(block, "Additive Model")
                @timeit to "compute_additive_model_model" datamanager = Additive.compute_additive_model(datamanager, block_nodes[block], datamanager.get_properties(block, "Additive Model"), time, dt)
            end
        end

    end

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
                properties(block_id, model, get_model_parameter(params::Dict, model, block[model]))
            end
        end
    end
    return properties
end

"""
    init_material_model_fields(datamanager::Module)

Initialize material model fields

# Arguments
- `datamanager::Data_manager`: Datamanager.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function init_material_model_fields(datamanager::Module)
    dof = datamanager.get_dof()
    datamanager.create_node_field("Forces", Float64, dof) #-> only if it is an output
    # tbd later in the compute class
    datamanager.create_node_field("Forces", Float64, dof)
    datamanager.create_node_field("Force Densities", Float64, dof)
    datamanager.create_constant_node_field("Acceleration", Float64, dof)
    datamanager.create_node_field("Velocity", Float64, dof)
    datamanager.create_constant_bond_field("Bond Forces", Float64, dof)
    datamanager.set_synch("Bond Forces", false, false)
    datamanager.set_synch("Force Densities", true, false)
    datamanager.set_synch("Velocity", false, true)
    datamanager.set_synch("Displacements", false, true)
    datamanager.set_synch("Acceleration", false, true)
    datamanager.set_synch("Deformed Coordinates", false, true)

    return datamanager
end

"""
    init_damage_model_fields(datamanager::Module)

Initialize damage model fields

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `params::Dict`: Parameters.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function init_damage_model_fields(datamanager::Module, params::Dict)
    dof = datamanager.get_dof()
    datamanager.create_node_field("Damage", Float64, 1)
    block_list = datamanager.get_block_list()
    anistropic_damage = false
    for block_id in block_list
        if !haskey(params["Blocks"]["block_$block_id"], "Damage Model")
            continue
        end
        damageName = params["Blocks"]["block_$block_id"]["Damage Model"]
        damage_parameter = params["Physics"]["Damage Models"][damageName]
        anistropic_damage = haskey(damage_parameter, "Anisotropic Damage")
    end
    if anistropic_damage
        datamanager.create_constant_bond_field("Bond Damage Anisotropic", Float64, dof, 1)
    end
    if length(datamanager.get_inverse_nlist()) == 0
        nlist = datamanager.get_field("Neighborhoodlist")
        inverse_nlist = datamanager.set_inverse_nlist(find_inverse_bond_id(nlist))
    end
    return datamanager
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
    deformed_coorN[:] = copy(datamanager.get_field("Coordinates"))
    deformed_coorNP1[:] = copy(datamanager.get_field("Coordinates"))
    datamanager.create_node_field("Displacements", Float64, dof)
    if solver_options["Additive Models"]
        @timeit to "additive_model_fields" datamanager = Physics.init_additive_model_fields(datamanager)
        heat_capacity = datamanager.create_constant_node_field("Specific Heat Capacity", Float64, 1)
        heat_capacity = set_heat_capacity(params, block_nodes, heat_capacity) # includes the neighbors
        for block in datamanager.get_block_list()
            if datamanager.check_property(block, "Additive Model")
                @timeit to "init additive model" datamanager = Additive.init_damage_model(datamanager, block_nodes[block], block)
            end
        end
    end
    if solver_options["Damage Models"]
        @timeit to "damage_model_fields" datamanager = Physics.init_damage_model_fields(datamanager, params)
        for block in datamanager.get_block_list()
            if datamanager.check_property(block, "Damage Model")
                @timeit to "init damage model" datamanager = Damage.init_damage_model(datamanager, block_nodes[block], block)
            end
        end

    end
    if solver_options["Material Models"]
        @timeit to "material model fields" datamanager = Physics.init_material_model_fields(datamanager)
        for block in datamanager.get_block_list()
            @timeit to "init material" datamanager = Material.init_material_model(datamanager, block_nodes[block], block)
        end
    end
    if solver_options["Thermal Models"]
        @timeit to "thermal model fields" datamanager = Physics.init_thermal_model_fields(datamanager)
        heat_capacity = datamanager.create_constant_node_field("Specific Heat Capacity", Float64, 1)
        heat_capacity = set_heat_capacity(params, block_nodes, heat_capacity) # includes the neighbors
        for block in datamanager.get_block_list()
            @timeit to "init thermal model" datamanager = Thermal.init_thermal_model(datamanager, block_nodes[block], block)
        end
    end

    return init_pre_calculation(datamanager, datamanager.get_physics_options())
end

"""
    init_thermal_model_fields(datamanager::Module)

Initialize thermal model fields

# Arguments
- `datamanager::Data_manager`: Datamanager.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function init_thermal_model_fields(datamanager::Module)
    datamanager.create_node_field("Temperature", Float64, 1)
    datamanager.create_node_field("Heat Flow", Float64, 1)
    datamanager.create_node_field("Specific Volume", Float64, 1)
    datamanager.create_constant_bond_field("Bond Heat Flow", Float64, 1)
    # if it is already initialized via mesh file no new field is created here
    datamanager.create_constant_node_field("Surface_Nodes", Bool, 1, true)
    return datamanager
end

"""
    init_additive_model_fields(datamanager::Module)

Initialize additive model fields

# Arguments
- `datamanager::Data_manager`: Datamanager.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function init_additive_model_fields(datamanager::Module)
    if !("Activation_Time" in datamanager.get_all_field_keys())
        @error "'Activation_Time' is missing. Please define an 'Activation_Time' for each point in the mesh file."
    end
    # must be specified, because it might be that no temperature model has been defined
    datamanager.create_node_field("Temperature", Float64, 1)
    datamanager.create_node_field("Heat Flow", Float64, 1)
    bond_damageN = datamanager.get_field("Bond Damage", "N")
    bond_damageNP1 = datamanager.get_field("Bond Damage", "NP1")
    nnodes = datamanager.get_nnodes()
    if !("Active" in datamanager.get_all_field_keys())
        active = datamanager.create_constant_node_field("Active", Bool, 1, false)
        for iID in 1:nnodes
            bond_damageN[iID][:] .= 0
            bond_damageNP1[iID][:] .= 0
        end
    end
    nlist = datamanager.get_field("Neighborhoodlist")
    inverse_nlist = datamanager.set_inverse_nlist(find_inverse_bond_id(nlist))
    return datamanager
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
            Material.determine_isotropic_parameter(datamanager.get_properties(block, "Material Model"))
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