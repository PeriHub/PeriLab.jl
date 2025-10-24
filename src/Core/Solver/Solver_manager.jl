# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Solver_Manager
using TimerOutputs: @timeit
using ..Parameter_Handling:
                            get_density,
                            get_horizon,
                            get_solver_name,
                            get_model_options,
                            get_fem_block,
                            get_calculation_options,
                            get_angles,
                            get_block_names_and_ids,
                            get_solver_params
using ..Helpers
include("../../Models/Material/Material_Basis.jl")
include("../Module_inclusion/set_Modules.jl")
include("../../FEM/FEM_Factory.jl")
include("../../Models/Model_Factory.jl")
include("../BC_manager.jl")
include("Verlet_solver.jl")
include("Static_solver.jl")
using ..MPI_Communication: synch_responder_to_controller,
                           synch_controller_to_responder,
                           synch_controller_bonds_to_responder,
                           synch_controller_bonds_to_responder_flattened

include("../Influence_function.jl")

using .Model_Factory: init_models, read_properties
using .Boundary_Conditions: init_BCs
using .Verlet_Solver
using .FEM
using .Influence_Function

export init
export solver

"""
    init(params::Dict, datamanager::Module)

Initialize the solver

# Arguments
- `params::Dict`: The parameters
- `datamanager::Module`: Datamanager
# Returns
- `block_nodes::Dict{Int64,Vector{Int64}}`: A dictionary mapping block IDs to collections of nodes.
- `bcs::Dict{Any,Any}`: A dictionary containing boundary conditions.
- `datamanager::Module`: The data manager module that provides access to data fields and properties.
- `solver_options::Dict{String,Any}`: A dictionary containing solver options.
"""
function init(params::Dict,
              datamanager::Module,
              step_id::Union{Nothing,Int64} = nothing)
    solver_options = Dict()
    nnodes = datamanager.get_nnodes()
    num_responder = datamanager.get_num_responder()
    block_ids = datamanager.get_field("Block_Id")
    block_nodes_with_neighbors = get_block_nodes(block_ids, nnodes + num_responder)
    block_nodes = get_block_nodes(block_ids, nnodes)
    block_name_list,
    block_id_list = get_block_names_and_ids(params, block_ids,
                                            datamanager.get_mpi_active())
    datamanager.set_block_name_list(block_name_list)
    datamanager.set_block_id_list(block_id_list)
    density = datamanager.create_constant_node_field("Density", Float64, 1)
    horizon = datamanager.create_constant_node_field("Horizon", Float64, 1)
    if datamanager.fem_active()
        fem_block = datamanager.create_constant_node_field("FEM Block", Bool, 1, false)
        fem_block = set_fem_block(params, block_nodes_with_neighbors, fem_block) # includes the neighbors
    end
    active_nodes = datamanager.create_constant_node_field("Active Nodes", Int64, 1)
    update_nodes = datamanager.create_constant_node_field("Update Nodes", Int64, 1)
    datamanager.create_constant_node_field("Update", Bool, 1, true)
    density = set_density(params, block_nodes_with_neighbors, density) # includes the neighbors
    horizon = set_horizon(params, block_nodes_with_neighbors, horizon) # includes the neighbors
    set_angles(datamanager, params, block_nodes_with_neighbors) # includes the Neighbors
    solver_params = isnothing(step_id) ? params["Solver"] :
                    get_solver_params(params, step_id)
    solver_options["Models"] = get_model_options(solver_params)
    solver_options["All Models"] = get_model_options(solver_params)
    solver_options["Calculation"] = get_calculation_options(solver_params)
    if !isnothing(step_id)
        for step in 1:datamanager.get_max_step()
            step_solver_params = get_solver_params(params, step)
            append!(solver_options["All Models"],
                    get_model_options(step_solver_params))
            calc_options = get_calculation_options(step_solver_params)
            for key in keys(calc_options)
                if calc_options[key]
                    solver_options["Calculation"][key] = true
                end
            end
        end
        solver_options["All Models"] = unique(solver_options["All Models"])
    end
    datamanager.create_constant_bond_field("Influence Function", Float64, 1, 1)
    for iblock in eachindex(block_nodes)
        datamanager = Influence_Function.init_influence_function(block_nodes[iblock],
                                                                 datamanager,
                                                                 params["Discretization"])
    end
    datamanager.create_bond_field("Bond Damage", Float64, 1, 1)
    @debug "Read properties"
    read_properties(params, datamanager, "Material" in solver_options["All Models"])
    @debug "Init models"
    @timeit "init_models" datamanager=init_models(params,
                                                  datamanager,
                                                  block_nodes,
                                                  solver_options,
                                                  synchronise_field)
    @debug "Init Boundary Conditions"
    @timeit "init_BCs" bcs=init_BCs(params, datamanager)
    solver_options["Solver"] = get_solver_name(solver_params)
    if get_solver_name(solver_params) == "Verlet"
        @debug "Init " * get_solver_name(solver_params)
        @timeit "init_solver" Verlet_Solver.init_solver(solver_options,
                                                        solver_params,
                                                        bcs,
                                                        datamanager,
                                                        block_nodes)
    elseif solver_options["Solver"] == "Static"
        @debug "Init " * get_solver_name(solver_params)
        @timeit "init_solver" Static_Solver.init_solver(solver_options,
                                                        solver_params,
                                                        bcs,
                                                        datamanager,
                                                        block_nodes)

    else
        @error get_solver_name(solver_params) * " is no valid solver."
        return nothing
    end

    if datamanager.fem_active()
        datamanager = FEM.init_FEM(params, datamanager)
        datamanager = FEM.Coupling.init_coupling(datamanager,
                                                 1:datamanager.get_nnodes(),
                                                 params)
    end
    if !datamanager.has_key("Active")
        active = datamanager.create_constant_node_field("Active", Bool, 1, true)
    end
    #TODO: sync active with datamanager

    datamanager = remove_models(datamanager, solver_options["Models"])

    @debug "Finished Init Solver"
    return block_nodes, bcs, datamanager, solver_options
end

"""
    set_density(params::Dict, block_nodes::Dict, density::Vector{Float64})

Sets the density of the nodes in the dictionary.

# Arguments
- `params::Dict`: The parameters
- `block_nodes::Dict`: A dictionary mapping block IDs to collections of nodes
- `density::Vector{Float64}`: The density
# Returns
- `density::Vector{Float64}`: The density
"""
function set_density(params::Dict, block_nodes::Dict, density::Vector{Float64})
    for block in eachindex(block_nodes)
        density[block_nodes[block]] .= get_density(params, block)
    end
    return density
end

"""
    set_angles(datamanager::Module, params::Dict, block_nodes::Dict)

Sets the density of the nodes in the dictionary.

# Arguments
- `datamanager::Module`: The data manager
- `params::Dict`: The parameters
- `block_nodes::Dict`: A dictionary mapping block IDs to collections of nodes
"""
function set_angles(datamanager::Module, params::Dict, block_nodes::Dict)
    mesh_angles = false
    if "Angles" in datamanager.get_all_field_keys()
        datamanager.set_rotation(true)
        mesh_angles = true
    end
    if "Element Angles" in datamanager.get_all_field_keys()
        datamanager.set_element_rotation(true)
    end

    block_rotation = false
    dof = datamanager.get_dof()
    for block in eachindex(block_nodes)
        if get_angles(params, block, dof) !== nothing
            block_rotation = true
            break
        end
    end
    if block_rotation
        if mesh_angles
            @warn "Angles defined in mesh will be overwritten by block angles"
        end
        datamanager.set_rotation(true)
        angles = datamanager.create_constant_node_field("Angles", Float64, dof)

        for block in eachindex(block_nodes)
            angles_global = get_angles(params, block, dof)
            if isnothing(angles_global)
                angles_global = 0.0
            end
            for iID in block_nodes[block]
                angles[iID, :] .= angles_global
            end
        end
    end
end

"""
    set_fem_block(params::Dict, block_nodes::Dict, fem_block::Vector{Bool})

Sets the fem_block of the nodes in the dictionary.

# Arguments
- `params::Dict`: The parameters
- `block_nodes::Dict`: A dictionary mapping block IDs to collections of nodes
- `fem_block::Vector{Bool}`: The fem_block
# Returns
- `fem_block::Vector{Bool}`: The fem_block
"""
function set_fem_block(params::Dict, block_nodes::Dict, fem_block::Vector{Bool})
    for block in eachindex(block_nodes)
        fem_block[block_nodes[block]] .= get_fem_block(params, block)
    end
    return fem_block
end

"""
    set_horizon(params::Dict, block_nodes::Dict, horizon::Vector{Float64})

Sets the horizon of the nodes in the dictionary.

# Arguments
- `params::Dict`: The parameters
- `block_nodes::Dict`: A dictionary mapping block IDs to collections of nodes
- `horizon::Vector{Float64}`: The horizon
# Returns
- `horizon::Vector{Float64}`: The horizon
"""
function set_horizon(params::Dict, block_nodes::Dict, horizon::Vector{Float64})
    for block in eachindex(block_nodes)
        horizon[block_nodes[block]] .= get_horizon(params, block)
    end
    return horizon
end

"""
    solver(solver_options::Dict{String,Any}, block_nodes::Dict{Int64,Vector{Int64}}, bcs::Dict{Any,Any}, datamanager::Module, outputs::Dict{Int64,Dict{}}, result_files::Vector{Any}, write_results, silent::Bool)

Runs the solver.

# Arguments
- `solver_options::Dict{String,Any}`: The solver options
- `block_nodes::Dict{Int64,Vector{Int64}}`: A dictionary mapping block IDs to collections of nodes
- `bcs::Dict{Any,Any}`: The boundary conditions
- `datamanager::Module`: The data manager module
- `outputs::Dict{Int64,Dict{}}`: A dictionary for output settings
- `result_files::Vector{Any}`: A vector of result files
- `write_results`: A function to write simulation results
- `silent::Bool`: A boolean flag to suppress progress bars
# Returns
- `result_files`: A vector of updated result files
"""
function solver(solver_options::Dict{Any,Any},
                block_nodes::Dict{Int64,Vector{Int64}},
                bcs::Dict{Any,Any},
                datamanager::Module,
                outputs::Dict{Int64,Dict{}},
                result_files::Vector{Dict},
                write_results,
                silent::Bool)
    if solver_options["Solver"] == "Verlet"
        return Verlet_Solver.run_solver(solver_options,
                                        block_nodes,
                                        bcs,
                                        datamanager,
                                        outputs,
                                        result_files,
                                        synchronise_field,
                                        write_results,
                                        compute_parabolic_problems_before_model_evaluation,
                                        compute_parabolic_problems_after_model_evaluation,
                                        silent)
    elseif solver_options["Solver"] == "Static"
        return Static_Solver.run_solver(solver_options,
                                        block_nodes,
                                        bcs,
                                        datamanager,
                                        outputs,
                                        result_files,
                                        synchronise_field,
                                        write_results,
                                        compute_parabolic_problems_before_model_evaluation,
                                        compute_parabolic_problems_after_model_evaluation,
                                        silent)
    end
end

"""
    synchronise_field(comm, synch_fields::Dict, overlap_map, get_field, synch_field::String, direction::String)

Synchronises field.

# Arguments
- `comm`: The MPI communicator
- `synch_fields::Dict`: A dictionary of fields
- `overlap_map`: The overlap map
- `get_field`: The function to get the field
- `synch_field::String`: The field
- `direction::String`: The direction
# Returns
- `nothing`
"""
function synchronise_field(comm,
                           synch_fields::Dict,
                           overlap_map,
                           get_field,
                           synch_field::String,
                           direction::String)
    # might not needed
    if !haskey(synch_fields, synch_field)
        @error "Field $synch_field does not exist in synch_field dictionary"
        return nothing
    end
    if direction == "download_from_cores"
        if synch_fields[synch_field][direction]
            vector = get_field(synch_field, synch_fields[synch_field]["time"])
            return synch_responder_to_controller(comm,
                                                 overlap_map,
                                                 vector,
                                                 synch_fields[synch_field]["dof"])
        end
        return nothing
    end
    if direction == "upload_to_cores"
        if synch_fields[synch_field][direction]
            vector = get_field(synch_field, synch_fields[synch_field]["time"])
            if occursin("Bond", synch_field)
                return synch_controller_bonds_to_responder_flattened(comm,
                                                                     overlap_map,
                                                                     vector,
                                                                     synch_fields[synch_field]["dof"])
            else
                return synch_controller_to_responder(comm,
                                                     overlap_map,
                                                     vector,
                                                     synch_fields[synch_field]["dof"])
            end
        end
        return nothing
    end
    @error "Wrong direction key word $direction in function synchronise_field; it should be download_from_cores or upload_to_cores"
    return nothing
end

"""
    remove_models(datamanager::Module, solver_options::Vector{String})

Sets the active models to false if they are deactivated in the solver. They can be active, because they are defined as model and in the blocks.

# Arguments
- `datamanager::Module`: The MPI communicator
- `solver_options::Vector{String}`: A dictionary of fields
# Returns
- `datamanager`
"""
function remove_models(datamanager::Module, solver_options::Vector{String})
    check = replace.(solver_options .* " Model", "_" => " ")
    for active_model_name in keys(datamanager.get_active_models())
        if !(active_model_name in check)
            datamanager.remove_active_model(active_model_name)
        end
    end
    return datamanager
end

## TODO generalize this interface

function compute_parabolic_problems_before_model_evaluation(active_nodes, datamanager,
                                                            solver_options)
    if !("Thermal" in solver_options["Models"]) &&
       !("Thermal" in solver_options["All Models"])
        return
    end
    temperatureN = datamanager.get_field("Temperature", "N")
    temperatureNP1 = datamanager.get_field("Temperature", "NP1")
    deltaT = datamanager.get_field("Delta Temperature")
    if "Thermal" in solver_options["Models"]
        temperatureNP1[active_nodes] = temperatureN[active_nodes] + deltaT[active_nodes]
    else
        if "Thermal" in solver_options["All Models"]
            temperatureNP1[active_nodes] = temperatureN[active_nodes]
        end
    end
end
function compute_parabolic_problems_after_model_evaluation(active_nodes, datamanager,
                                                           solver_options, dt)
    if !("Thermal" in solver_options["Models"]) &&
       !("Thermal" in solver_options["All Models"])
        return
    end
    deltaT = datamanager.get_field("Delta Temperature")
    flowNP1 = datamanager.get_field("Heat Flow", "NP1")
    density = datamanager.get_field("Density")
    heat_capacity = datamanager.get_field("Specific Heat Capacity")
    if "Thermal" in solver_options["Models"]
        check_inf_or_nan(flowNP1, "Heat Flow")
        # heat capacity check. if it is zero deltaT = 0
        @views deltaT[active_nodes] = -flowNP1[active_nodes] .* dt ./
                                      (density[active_nodes] .*
                                       heat_capacity[active_nodes])
        # if fem_option && time == 0
        #     @warn "Thermal models are not supported for FEM yet."
        # end
    end
end

end
