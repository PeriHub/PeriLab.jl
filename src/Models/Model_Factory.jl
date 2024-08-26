# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Model_Factory
include("../Support/helpers.jl")
include("./Material/material_basis.jl")
using .Helpers: check_inf_or_nan, find_active, get_active_update_nodes, invert
include("./Additive/Additive_Factory.jl")
include("./Corrosion/Corrosion_Factory.jl")
include("./Damage/Damage_Factory.jl")
include("./Material/Material_Factory.jl")
include("./Thermal/Thermal_Factory.jl")
include("./Pre_calculation/Pre_Calculation_Factory.jl")
include("../Support/Parameters/parameter_handling.jl")
include("../Compute/compute_field_values.jl")
using .Parameter_Handling: get_models_option, get_model_parameter, get_heat_capacity
# in future FEM will be outside of the Model_Factory
include("../FEM/FEM_Factory.jl")

using .Additive
using .Corrosion
using .Damage
using .Material
using .Pre_calculation
using .Thermal
# in future FEM will be outside of the Model_Factory
using .FEM
using TimerOutputs
export compute_models
export init_models
export read_properties


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
function init_models(
    params::Dict,
    datamanager::Module,
    block_nodes::Dict{Int64,Vector{Int64}},
    solver_options::Dict,
    to::TimerOutput,
)
    # TODO integrate this correctly
    rotation = datamanager.get_rotation()
    #if rotation

    rotN, rotNP1 =
        datamanager.create_node_field("Rotation", Float64, "Matrix", datamanager.get_dof())
    #    for iID in nodes
    #        rotN[iID, :, :] = Geometry.rotation_tensor(angles[iID, :])
    #        rotNP1 = copy(rotN)
    #    end
    #end

    for model_name in solver_options["Models"]
        datamanager = add_model(datamanager, model_name)
    end
    # TODO order of models

    for (active_model_name, active_model) in pairs(datamanager.get_active_models())
        @info "Init $active_model_name "
        @timeit to "$active_model_name model fields" datamanager =
            active_model.init_fields(datamanager)

        for block in eachindex(block_nodes)
            if datamanager.check_property(block, active_model_name)
                @timeit to "init $active_model_name model" datamanager =
                    active_model.init_model(datamanager, block_nodes[block], block)
                # TODO active_model.fields_for_local_synchronization()
                # put it in datamanager
            end
        end
    end

    if "Additive" in solver_options["Models"] || "Thermal" in solver_options["Models"]
        heat_capacity = datamanager.get_field("Specific Heat Capacity")
        heat_capacity = set_heat_capacity(params, block_nodes, heat_capacity) # includes the neighbors
    end
    @info "Finalize Init Models"
    return datamanager
end

"""
    compute_models(datamanager::Module, block_nodes::Dict{Int64,Vector{Int64}}, dt::Float64, time::Float64, options::Dict, synchronise_field, to::TimerOutput)

Computes the models models

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
function compute_models(
    datamanager::Module,
    block_nodes::Dict{Int64,Vector{Int64}},
    dt::Float64,
    time::Float64,
    options::Dict,
    synchronise_field,
    to::TimerOutput,
)
    fem_option = datamanager.fem_active()
    if fem_option
        # no damage occur here. Therefore, it runs one time
        fe_nodes = datamanager.get_field("FE Nodes")
        nelements = datamanager.get_num_elements()
        # in future the FE analysis can be put into the block loop. Right now it is outside the block.
        datamanager = FEM.eval(
            datamanager,
            Vector{Int64}(1:nelements),
            datamanager.get_properties(1, "FEM"),
            time,
            dt,
        )

    end

    active = datamanager.get_field("Active")
    for (active_model_name, active_model) in pairs(datamanager.get_active_models())
        #synchronise_field(datamanager.local_synch_fiels(active_model_name))
        for block in eachindex(block_nodes)
            nodes = block_nodes[block]
            active_nodes = view(nodes, find_active(active[nodes]))
            if fem_option
                active_nodes = view(nodes, find_active(.~fe_nodes[active_nodes]))
            end
            if datamanager.check_property(block, active_model_name)

                @timeit to "compute $active_model_name model" datamanager =
                    active_model.compute_model(
                        datamanager,
                        active_nodes,
                        datamanager.get_properties(block, active_model_name),
                        time,
                        dt,
                        to,
                    )
            end
        end
    end

    update_nodes = datamanager.get_field("Update List")
    for (active_model_name, active_model) in pairs(datamanager.get_active_models())
        if active_model_name == "Additve Model" || active_model_name == "Damage Model"
            continue
        end
        #synchronise_field(datamanager.local_synch_fiels(active_model_name))
        for block in eachindex(block_nodes)
            nodes = block_nodes[block]
            active_nodes = view(nodes, find_active(active[nodes]))
            if fem_option
                active_nodes = view(nodes, find_active(.~fe_nodes[active_nodes]))
            end
            active_nodes, update_nodes =
                get_active_update_nodes(active, update_list, block_nodes, block)
            if datamanager.check_property(block, active_model_name)

                @timeit to "compute $active_model_name model" datamanager =
                    active_model.compute_model(
                        datamanager,
                        active_nodes,
                        datamanager.get_properties(block, active_model_name),
                        time,
                        dt,
                    )
            end
        end
    end

    # active_nodes::Vector{Int64} = []
    # update_nodes::Vector{Int64} = []


    @timeit to "pre_calculation" datamanager =
        Pre_calculation.compute(datamanager, block_nodes)
    @timeit to "pre_synchronize" Pre_calculation.synchronize(
        datamanager,
        datamanager.get_models_options(),
        synchronise_field,
    )
    for block in eachindex(block_nodes)
        active_nodes, update_nodes =
            get_active_update_nodes(active, update_list, block_nodes, block)
        if fem_option
            update_nodes =
                block_nodes[block][find_active(Vector{Bool}(.~fe_nodes[update_nodes]))]
        end
        if options["Thermal Models"]
            if datamanager.check_property(block, "Thermal Model")
                @timeit to "compute_model" datamanager = Thermal.compute_model(
                    datamanager,
                    update_nodes,
                    datamanager.get_properties(block, "Thermal Model"),
                    time,
                    dt,
                )
            end
        end

        if options["Material Models"]
            for material_model in datamanager.get_material_models()
                @timeit to "synch_field" datamanager =
                    Material.synch_field(datamanager, material_model, synchronise_field)
                #TODO: Check that same fields in seperate damage models arent being synched twice
            end
            if datamanager.check_property(block, "Material Model")
                model_param = datamanager.get_properties(block, "Material Model")
                @timeit to "bond_forces" datamanager = Material.compute_model(
                    datamanager,
                    update_nodes,
                    model_param,
                    time,
                    dt,
                    to,
                )
                #TODO: I think this needs to stay here as we need the active_nodes not the update_nodes
                @timeit to "distribute_force_densities" datamanager =
                    Material.distribute_force_densities(datamanager, active_nodes)


                if !occursin("Correspondence", model_param["Material Model"])
                    if options["Calculate Cauchy"] |
                       options["Calculate von Mises"] |
                       options["Calculate Strain"]
                        datamanager = get_partial_stresses(datamanager, active_nodes)
                    end
                    if options["Calculate von Mises"]
                        datamanager =
                            Material.calculate_von_mises_stress(datamanager, active_nodes)
                    end
                    if options["Calculate Strain"]
                        material_parameter =
                            datamanager.get_properties(block, "Material Model")
                        hookeMatrix = get_Hooke_matrix(
                            material_parameter,
                            material_parameter["Symmetry"],
                            datamanager.get_dof(),
                        )
                        datamanager = Material.calculate_strain(
                            datamanager,
                            active_nodes,
                            invert(hookeMatrix, "Hook matrix not invertable"),
                        )
                    end
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
function compute_damage_pre_calculation(
    datamanager::Module,
    options::Dict,
    nodes::Union{SubArray,Vector{Int64}},
    block::Int64,
    time::Float64,
    dt::Float64,
    to::TimerOutput,
)

    if options["Thermal Models"]
        @timeit to "thermal_model" datamanager = Thermal.compute_model(
            datamanager,
            nodes,
            datamanager.get_properties(block, "Thermal Model"),
            time,
            dt,
        )
    end

    if options["Material Models"]
        @timeit to "compute_model" datamanager = Material.compute_model(
            datamanager,
            nodes,
            datamanager.get_properties(block, "Material Model"),
            time,
            dt,
            to,
        )
    end
    datamanager = Damage.compute_damage_pre_calculation(
        datamanager,
        nodes,
        block,
        datamanager.get_properties(block, "Damage Model"),
        time,
        dt,
    )
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
function get_block_model_definition(
    params::Dict,
    block_id::Int64,
    prop_keys::Vector{String},
    properties,
)
    # properties function from datamanager
    if haskey(params["Blocks"], "block_" * string(block_id))
        block = params["Blocks"]["block_"*string(block_id)]
        for model in prop_keys
            if haskey(block, model)
                properties(
                    block_id,
                    model,
                    get_model_parameter(params, model, block[model]),
                )
            end
        end
    end
    return properties
end


"""
    read_properties(params::Dict, datamanager::Module, material_model::Bool)

Read properties of material.

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
    models_options = datamanager.get_models_options()

    datamanager.set_models_options(get_models_option(params, models_options))

    for block in blocks
        get_block_model_definition(params, block, prop_keys, datamanager.set_properties)
    end
    if material_model
        dof = datamanager.get_dof()
        for block in blocks
            Material.check_material_symmetry(
                dof,
                datamanager.get_properties(block, "Material Model"),
            )
            Material.determine_isotropic_parameter(
                datamanager,
                datamanager.get_properties(block, "Material Model"),
            )
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

"""
    add_model(datamanager::Module, model_name::String)

Includes the models in the datamanager and checks if the model definition is correct or not.

# Arguments
- `datamanager::Module`: Datamanager
- `model_name::String`: The block nodes

# Returns
- `datamanager::Module`: Datamanager
"""
function add_model(datamanager::Module, model_name::String)
    # TODO test missing
    try
        datamanager.add_active_model(model_name * " Model", eval(Meta.parse(model_name)))
        return datamanager
    catch
        @error "Model $model_name is not specified and cannot be included."
        return nothing
    end
end
end
