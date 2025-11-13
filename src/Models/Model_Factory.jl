# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Model_Factory

using TimerOutputs: @timeit
using ...Data_Manager
using ...Helpers:
                  check_inf_or_nan, find_active_nodes, get_active_update_nodes, invert,
                  determinant
include("./Pre_calculation/Pre_Calculation_Factory.jl")
include("./Surface_correction/Surface_correction.jl")
include("./Contact/Contact_Factory.jl")
include("./Additive/Additive_Factory.jl")
include("./Degradation/Degradation_Factory.jl")
include("./Damage/Damage_Factory.jl")
include("./Material/Material_Factory.jl")
include("./Thermal/Thermal_Factory.jl")
using ...Parameter_Handling: get_model_parameter, get_heat_capacity
using .Additive
using .Degradation
using .Damage
using .Material
using .Pre_Calculation
using .Surface_Correction: init_surface_correction, compute_surface_correction
using .Thermal
using .Contact
# in future FEM will be outside of the Model_Factory
using ..FEM
export compute_models
export init_models
export read_properties

"""
    init_models(params::Dict, block_nodes::Dict{Int64,Vector{Int64}}, solver_options::Dict)

Initialize models

# Arguments
- `params::Dict`: Parameters.
- `block_nodes::Dict{Int64,Vector{Int64}}`: block nodes.
- `solver_options::Dict`: Solver options.
"""
function init_models(params::Dict,
                     block_nodes::Dict{Int64,Vector{Int64}},
                     solver_options::Dict,
                     synchronise_field)
    if "Pre_Calculation" in solver_options["Models"]
        @info "Check pre calculation models are initialized for material models."
        Pre_Calculation.check_dependencies(block_nodes)
        if haskey(params["Models"], "Material Models")
            for mat in keys(params["Models"]["Material Models"])
                if haskey(params["Models"]["Material Models"][mat], "Accuracy Order")
                    Data_Manager.set_accuracy_order(params["Models"]["Material Models"][mat]["Accuracy Order"])
                end
            end
        end
    end
    for model_name in solver_options["Models"]
        add_model(model_name)
    end
    for model_name in solver_options["All Models"]
        add_model(model_name, true)
    end

    if "Additive" in solver_options["Models"] || "Thermal" in solver_options["Models"]
        heat_capacity = Data_Manager.create_constant_node_scalar_field("Specific Heat Capacity",
                                                                       Float64)
        heat_capacity = set_heat_capacity(params, block_nodes, heat_capacity) # includes the neighbors
    end

    if isnothing(Data_Manager.get_step()) || Data_Manager.get_step() == 1
        for (active_model_name, active_model) in pairs(Data_Manager.get_active_models(true))
            @debug "Init $active_model_name fields"
            @timeit "$active_model_name model fields" active_model.init_fields()
        end
    end
    for (active_model_name, active_model) in pairs(Data_Manager.get_active_models())
        @info "Init $active_model_name"

        for block in eachindex(block_nodes)
            if Data_Manager.check_property(block, active_model_name)
                @timeit "init $active_model_name models" active_model.init_model(block_nodes[block],
                                                                                 block)
                @timeit "init fields_for_local_synchronization $active_model_name models" active_model.fields_for_local_synchronization(active_model_name,
                                                                                                                                        block)
                if active_model_name == "Damage Model" &&
                   haskey(Data_Manager.get_properties(block, active_model_name),
                          "Local Damping")
                    Material.init_local_damping(block_nodes[block],
                                                Data_Manager.get_properties(block,
                                                                            "Material Model"),
                                                Data_Manager.get_properties(block,
                                                                            "Damage Model"))
                end
                # put it in Data_Manager
            end
        end
    end

    init_surface_correction(params, local_synch,
                            synchronise_field)

    if solver_options["Calculation"]["Calculate Cauchy"] |
       solver_options["Calculation"]["Calculate von Mises stress"]
        Data_Manager.create_node_tensor_field("Cauchy Stress",
                                              Float64, Data_Manager.get_dof())
    end
    if solver_options["Calculation"]["Calculate Strain"]
        Data_Manager.create_node_tensor_field("Strain", Float64, Data_Manager.get_dof())
    end
    if solver_options["Calculation"]["Calculate von Mises stress"]
        Data_Manager.create_node_scalar_field("von Mises Stress", Float64)
    end

    check_contact(params)
    @info "Finalize Init Models"
end

function check_contact(params::Dict)
    if haskey(params, "Contact")
        return Contact.init_contact_model(params["Contact"])
    end
end
function check_contact(params::Dict, time::Float64, dt::Float64)
    if length(params) != 0
        return Contact.compute_contact_model(params, time,
                                             dt)
    end
end

"""
    compute_models(block_nodes::Dict{Int64,Vector{Int64}}, dt::Float64, time::Float64, options::Vector{String}, synchronise_field)

Computes the material point models

# Arguments
- `block_nodes::Dict{Int64,Vector{Int64}}`: The block nodes
- `dt::Float64`: The time step
- `time::Float64`: The current time of the solver
- `options::Vector{String}`: The options
- `synchronise_field`: The synchronise field
"""
function compute_models(block_nodes::Dict{Int64,Vector{Int64}},
                        dt::Float64,
                        time::Float64,
                        options::Vector{String},
                        synchronise_field)
    fem_option = Data_Manager.fem_active()
    if fem_option
        fe_nodes = Data_Manager.get_field("FE Nodes")
    end

    active_list = Data_Manager.get_field("Active")

    # TODO check if pre calculation should run block wise. For mixed model applications it makes sense.
    # TODO add for pre calculation a whole model option, to get the neighbors as well, e.g. for bond associated
    # TODO check for loop order?
    for (active_model_name, active_model) in pairs(Data_Manager.get_active_models())
        #synchronise_field(Data_Manager.local_synch_fiels(active_model_name))
        if active_model_name == "Damage Model"
            continue
        end

        local_synch(active_model_name, "upload_to_cores", synchronise_field)
        # maybe not needed?
        local_synch(active_model_name,
                    "download_from_cores",
                    synchronise_field)

        for (block, nodes) in pairs(block_nodes)
            # "delete" the view of active nodes
            active_nodes = Data_Manager.get_field("Active Nodes")

            active_nodes = find_active_nodes(active_list,
                                             active_nodes,
                                             nodes,
                                             active_model_name != "Additive Model")

            if fem_option # TODO might lead to problems in 3D
                # find all non-FEM nodes in active nodes
                active_nodes = Data_Manager.get_field("Active Nodes")
                active_nodes = find_active_nodes(fe_nodes,
                                                 active_nodes,
                                                 1:Data_Manager.get_nnodes(),
                                                 false)
                if active_nodes == []
                    continue
                end
            end
            if Data_Manager.check_property(block, active_model_name)
                # synch
                @timeit "compute $active_model_name" active_model.compute_model(active_nodes,
                                                                                Data_Manager.get_properties(block,
                                                                                                            active_model_name),
                                                                                block,
                                                                                time,
                                                                                dt)
            end
        end
    end

    # Why not update_list.=false? -> avoid neighbors
    update_list = Data_Manager.get_field("Update")
    for (block, nodes) in pairs(block_nodes)
        update_list[nodes] .= false
    end

    for (active_model_name, active_model) in pairs(Data_Manager.get_active_models())
        if active_model_name == "Additive Model"
            continue
        end
        local_synch(active_model_name, "upload_to_cores", synchronise_field)
        # maybe not needed?
        local_synch(active_model_name,
                    "download_from_cores",
                    synchronise_field)
        for (block, nodes) in pairs(block_nodes)
            active_nodes = Data_Manager.get_field("Active Nodes")
            update_nodes = Data_Manager.get_field("Update Nodes")
            active_nodes = find_active_nodes(active_list, active_nodes, nodes)
            if fem_option
                # FEM active means FEM nodes
                active_nodes = Data_Manager.get_field("Active Nodes")
                active_nodes = find_active_nodes(fe_nodes,
                                                 active_nodes,
                                                 1:Data_Manager.get_nnodes(),
                                                 false)
                if active_nodes == []
                    continue
                end
            end
            update_nodes = get_update_nodes(active_list,
                                            update_list,
                                            nodes,
                                            update_nodes,
                                            active_nodes,
                                            active_model_name)

            # active or all, or does it not matter?

            if Data_Manager.check_property(block, active_model_name)
                # TODO synch
                @timeit "compute $active_model_name" active_model.compute_model(update_nodes,
                                                                                Data_Manager.get_properties(block,
                                                                                                            active_model_name),
                                                                                block,
                                                                                time,
                                                                                dt)
            end
        end
    end
    # must be here to avoid double distributions
    # distributes ones over all nodes
    if fem_option
        @timeit "FEM" begin
            nelements = Data_Manager.get_num_elements()

            @timeit "eval" FEM.eval_FEM(Vector{Int64}(1:nelements),
                                        Data_Manager.get_properties(1, "FEM"),
                                        time,
                                        dt)
            active_nodes = Data_Manager.get_field("Active Nodes")

            FEM.force_densities(find_active_nodes(fe_nodes, active_nodes,
                                                  1:Data_Manager.get_nnodes(), true))
        end
    end

    if "Material" in options
        if "Damage" in options
            for (block, nodes) in pairs(block_nodes)
                if haskey(Data_Manager.get_properties(block, "Damage Model"),
                          "Local Damping")
                    active_nodes = Data_Manager.get_field("Active Nodes")
                    if fem_option
                        active_nodes = find_active_nodes(active_list,
                                                         active_nodes,
                                                         find_active_nodes(fe_nodes,
                                                                           active_nodes,
                                                                           nodes))
                    else
                        find_active_nodes(active_list, active_nodes, nodes)
                    end
                    @timeit "local_damping_due_to_damage" Material.compute_local_damping(active_nodes,
                                                                                         Data_Manager.get_properties(block,
                                                                                                                     "Damage Model")["Local Damping"],
                                                                                         dt)
                end
            end
        end
        active_nodes = Data_Manager.get_field("Active Nodes")
        active_nodes = find_active_nodes(active_list, active_nodes,
                                         1:Data_Manager.get_nnodes())

        compute_surface_correction(active_nodes,
                                   local_synch,
                                   synchronise_field)

        @timeit "distribute_force_densities" Material.distribute_force_densities(active_nodes)
    end

    if fem_option
        @timeit "coupling" FEM.Coupling.compute_coupling(Data_Manager.get_properties(1,
                                                                                     "FEM"))
    end
    check_contact(Data_Manager.get_contact_properties(), time, dt)
    #=
    Used for shape tensor or other fixed calculations, to avoid an update if its not needed.
    The damage update is done in the second loop.
    =#
    update_list = Data_Manager.get_field("Update")
    for (block, nodes) in pairs(block_nodes)
        update_list[nodes] .= false
    end
end

"""
	compute_stiff_matrix_compatible_models(block_nodes::Dict{Int64,Vector{Int64}}, dt::Float64, time::Float64, options::Vector{String}, synchronise_field)

Computes the models models that are compatible with the stiffness matrix calculation.

# Arguments
- `block_nodes::Dict{Int64,Vector{Int64}}`: The block nodes
- `dt::Float64`: The time step
- `time::Float64`: The current time of the solver
- `options::Vector{String}`: The options
- `synchronise_field`: The synchronise field
"""
function compute_stiff_matrix_compatible_models(block_nodes::Dict{Int64,Vector{Int64}},
                                                dt::Float64,
                                                time::Float64,
                                                options::Vector{String},
                                                synchronise_field)
    active_list = Data_Manager.get_field("Active")

    for (active_model_name, active_model) in pairs(Data_Manager.get_active_models())

        #local_synch(Data_Manager, active_model_name, "upload_to_cores", synchronise_field)
        # maybe not needed?
        #local_synch(Data_Manager,
        #	active_model_name,
        #	"download_from_cores",
        #	synchronise_field)
        if active_model_name == "Material Model" && !("Thermal" in options)
            # we need here an activation trigger for mixed models in future
            continue
        end

        for (block, nodes) in pairs(block_nodes)
            # "delete" the view of active nodes
            active_nodes = Data_Manager.get_field("Active Nodes")

            active_nodes = find_active_nodes(active_list,
                                             active_nodes,
                                             nodes,
                                             active_model_name != "Additive Model")

            if Data_Manager.check_property(block, active_model_name)
                # synch
                @timeit "compute $active_model_name" active_model.compute_model(active_nodes,
                                                                                Data_Manager.get_properties(block,
                                                                                                            active_model_name),
                                                                                block,
                                                                                time,
                                                                                dt)
            end
        end
    end

    if ("Material" in options) && ("Thermal" in options)
        active_nodes = Data_Manager.get_field("Active Nodes")
        active_nodes = find_active_nodes(active_list, active_nodes,
                                         1:Data_Manager.get_nnodes())
        @timeit "distribute_force_densities" Material.distribute_force_densities(active_nodes)
    end
end

function get_update_nodes(active_list,
                          update_list,
                          nodes,
                          update_nodes,
                          active_nodes,
                          active_model_name)
    if active_model_name == "Damage Model"
        return @view active_nodes[:]
    else
        return get_active_update_nodes(active_list, update_list, nodes, update_nodes)
    end
end

"""
	get_block_model_definition(params::Dict, block_id_list::Int64, prop_keys::Vector{String}, properties)

Get block model definition.

Special case for pre calculation. It is set to all blocks, if no block definition is defined, but pre calculation is.

# Arguments
- `params::Dict`: Parameters.
- `block_id_list::Vector{Int64}`: List of block id's.
- `prop_keys::Vector{String}`: Property keys.
- `properties`: Properties function.
# Returns
- `properties`: Properties function.
"""
function get_block_model_definition(params::Dict,
                                    block_name_list::Vector{String},
                                    block_id_list::Vector{Int64},
                                    prop_keys::Vector{String},
                                    properties,
                                    directory::String = "",
                                    material_model::Bool = true)
    # properties function from Data_Manager

    if haskey(params["Models"], "Pre Calculation Global")
        for block_id in block_id_list
            properties(block_id,
                       "Pre Calculation Model",
                       params["Models"]["Pre Calculation Global"])
        end
    end

    for (block_id, block_name) in zip(block_id_list, block_name_list)
        if !haskey(params["Blocks"], block_name)
            continue
        end
        block = params["Blocks"][block_name]
        for model in prop_keys
            if model == "Material Model" && !material_model
                continue
            end
            if haskey(block, model)
                properties(block_id,
                           model,
                           get_model_parameter(params, model, block[model], directory))
            end
        end
    end
    return properties
end

"""
    read_properties(params::Dict, material_model::Bool)

Read properties of material.

# Arguments
- `params::Dict`: Parameters.
- `material_model::Bool`: Material model.
"""
function read_properties(params::Dict, material_model::Bool)
    Data_Manager.init_properties()
    block_name_list = Data_Manager.get_block_name_list()
    block_id_list = Data_Manager.get_block_id_list()
    prop_keys = Data_Manager.init_properties()
    directory = Data_Manager.get_directory()
    get_block_model_definition(params,
                               block_name_list,
                               block_id_list,
                               prop_keys,
                               Data_Manager.set_properties,
                               directory,
                               material_model)
    if material_model
        dof = Data_Manager.get_dof()
        for block in block_id_list
            Material.check_material_symmetry(dof,
                                             Data_Manager.get_properties(block,
                                                                         "Material Model"))
            Material.determine_isotropic_parameter(Data_Manager.get_properties(block,
                                                                               "Material Model"))
        end
    end
end

"""
    set_heat_capacity(params::Dict, block_nodes::Dict, heat_capacity::NodeScalarField{Float64})

Sets the heat capacity of the nodes in the dictionary.

# Arguments
- `params::Dict`: The parameters
- `block_nodes::Dict`: The block nodes
- `heat_capacity::NodeScalarField{Float64}`: The heat capacity array
# Returns
- `heat_capacity::SubArray`: The heat capacity array
"""
function set_heat_capacity(params::Dict, block_nodes::Dict,
                           heat_capacity::NodeScalarField{Float64})
    for block in eachindex(block_nodes)
        heat_capacity[block_nodes[block]] .= get_heat_capacity(params, block)
    end
    return heat_capacity
end

"""
    add_model(model_name::String)

Includes the models in the Data_Manager and checks if the model definition is correct or not.

# Arguments
- `model_name::String`: The block nodes
"""
function add_model(model_name::String, all::Bool = false)
    try
        # to catch "Pre_Calculation"
        Data_Manager.add_active_model(replace(model_name, "_" => " ") * " Model",
                                      eval(Meta.parse(model_name)),
                                      all)
    catch
        @error "Model $model_name is not specified and cannot be included."
        return nothing
    end
end

function local_synch(model, direction, synchronise_field)
    synch_fields = Data_Manager.get_local_synch_fields(model)
    for synch_field in keys(synch_fields)
        synchronise_field(Data_Manager.get_comm(),
                          synch_fields,
                          Data_Manager.get_overlap_map(),
                          Data_Manager.get_field,
                          synch_field,
                          direction)
    end
end
end
