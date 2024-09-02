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
using .Pre_Calculation
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
    @info "Check pre calculation models are initialized for material models."
    datamanager = Pre_Calculation.check_dependencies(datamanager, block_nodes)

    for model_name in solver_options["Models"]
        datamanager = add_model(datamanager, model_name)
    end

    for (active_model_name, active_model) in pairs(datamanager.get_active_models())
        @info "Init $active_model_name "
        @timeit to "$active_model_name model fields" datamanager =
            active_model.init_fields(datamanager)

        for block in eachindex(block_nodes)
            if datamanager.check_property(block, active_model_name)
                @timeit to "init $active_model_name models" datamanager =
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
    compute_models(datamanager::Module, block_nodes::Dict{Int64,Vector{Int64}}, dt::Float64, time::Float64, options::Vector{String}, synchronise_field, to::TimerOutput)

Computes the models models

# Arguments
- `datamanager::Module`: The datamanager
- `block_nodes::Dict{Int64,Vector{Int64}}`: The block nodes
- `dt::Float64`: The time step
- `time::Float64`: The current time of the solver
- `options::Vector{String}`: The options
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
    options::Vector{String},
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
    # TODO check if pre calculation should run block wise. For mixed model applications it makes sense.
    # TODO add for pre calculation a whole model option, to get the neighbors as well, e.g. for bond associated
    for (active_model_name, active_model) in pairs(datamanager.get_active_models())
        #synchronise_field(datamanager.local_synch_fiels(active_model_name))
        if active_model_name == "Damage Model"
            continue
        end
        ## TODO @jan-timo synchronisation is missing

        #for synch_key in active_model.fields_for_local_synchronization(datamanager.get_properties(block, active_model_name))
        #synchronise_field
        #@timeit to "upload_to_cores" datamanager.synch_manager(
        #    synchronise_field,
        #    "upload_to_cores",
        #)
        #end
        for (block, nodes) in pairs(block_nodes)
            #active_nodes = view(nodes, find_active(active[nodes]))
            active_nodes = view(nodes, find_active(active[nodes]))
            if fem_option
                # find all non-FEM nodes
                active_nodes = view(nodes, find_active(.~fe_nodes[active_nodes]))
            end
            if datamanager.check_property(block, active_model_name)
                # synch
                @timeit to "compute $active_model_name model" datamanager =
                    active_model.compute_model(
                        datamanager,
                        active_nodes,
                        datamanager.get_properties(block, active_model_name),
                        block,
                        time,
                        dt,
                        to,
                    )
            end
        end
    end
    # No update is needed, if no damage occur
    update_list = datamanager.get_field("Update List")
    for (block, nodes) in pairs(block_nodes)
        update_list[nodes] .= false
    end

    for (active_model_name, active_model) in pairs(datamanager.get_active_models())
        if active_model_name == "Additve Model"
            continue
        end
        #synchronise_field(datamanager.local_synch_fiels(active_model_name))
        for block in eachindex(block_nodes)
            nodes = block_nodes[block]
            active_nodes = view(nodes, find_active(active[nodes]))
            if fem_option
                # find all non-FEM nodes
                active_nodes = view(nodes, find_active(.~fe_nodes[active_nodes]))
            end
            active_nodes, update_nodes =
                get_active_update_nodes(active, update_list, block_nodes, block)
            # active or all, or does it not matter?
            if !(active_model_name == "Damage Model")
                update_nodes = active_nodes
            end
            if datamanager.check_property(block, active_model_name)
                # TODO synch
                @timeit to "compute $active_model_name model" datamanager =
                    active_model.compute_model(
                        datamanager,
                        update_nodes,
                        datamanager.get_properties(block, active_model_name),
                        block,
                        time,
                        dt,
                        to,
                    )
            end
        end
    end
    # must be here to avoid double distributions
    # distributes ones over all nodes
    if "Material" in options
        @timeit to "distribute_force_densities" datamanager =
            Material.distribute_force_densities(
                datamanager,
                find_active(active[1:datamanager.get_nnodes()]),
            )
    end

    # if !occursin("Correspondence", model_param["Material Model"])
    #     if options["Calculate Cauchy"] |
    #        options["Calculate von Mises"] |
    #        options["Calculate Strain"]
    #         datamanager = get_partial_stresses(datamanager, active_nodes)
    #     end
    #     if options["Calculate von Mises"]
    #         datamanager =
    #             Material.calculate_von_mises_stress(datamanager, active_nodes)
    #     end
    #     if options["Calculate Strain"]
    #         material_parameter =
    #             datamanager.get_properties(block, "Material Model")
    #         hookeMatrix = get_Hooke_matrix(
    #             material_parameter,
    #             material_parameter["Symmetry"],
    #             datamanager.get_dof(),
    #         )
    #         datamanager = Material.calculate_strain(
    #             datamanager,
    #             active_nodes,
    #             invert(hookeMatrix, "Hook matrix not invertable"),
    #         )
    #     end
    # end


    update_list .= true
    return datamanager

end

"""
    get_block_model_definition(params::Dict, block_id::Int64, prop_keys::Vector{String}, properties)

Get block model definition.

Special case for pre calculation. It is set to all blocks, if no block definition is defined, but pre calculation is.

# Arguments
- `params::Dict`: Parameters.
- `blocks::Vector{Int64}`: List of block id's.
- `prop_keys::Vector{String}`: Property keys.
- `properties`: Properties function.
# Returns
- `properties`: Properties function.
"""
function get_block_model_definition(
    params::Dict,
    blocks::Vector{Int64},
    prop_keys::Vector{String},
    properties,
)
    # properties function from datamanager
    pre_calculation_all_blocks::Bool = true
    for block_id in blocks
        if !haskey(params["Blocks"], "block_" * string(block_id))
            continue
        end
        block = params["Blocks"]["block_"*string(block_id)]
        for model in prop_keys
            if haskey(block, model)
                properties(
                    block_id,
                    model,
                    get_model_parameter(params, model, block[model]),
                )
                if model == "Pre Calculation Model"
                    pre_calculation_all_blocks = false
                end
            end
        end
    end
    # makes sure that pre calculation exists everywhere
    # material needs are set elsewhere
    # maybe not important and can be thrown away
    # TODO check if realy needed
    if haskey(params["Models"], "Pre Calculation Models") && pre_calculation_all_blocks
        for block_id in blocks
            properties(
                block_id,
                model,
                get_model_parameter(params, model, "Pre Calculation Model"),
            )
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
    datamanager.init_properties()
    blocks = datamanager.get_block_list()
    prop_keys = datamanager.init_properties()
    get_block_model_definition(params, blocks, prop_keys, datamanager.set_properties)
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
        # to catch "Pre_Calculation"
        datamanager.add_active_model(
            replace(model_name, "_" => " ") * " Model",
            eval(Meta.parse(model_name)),
        )
        return datamanager
    catch
        @error "Model $model_name is not specified and cannot be included."
        return nothing
    end
end


end
