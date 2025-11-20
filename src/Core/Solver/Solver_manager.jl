# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Solver_Manager
using TimerOutputs: @timeit
using ..Data_Manager
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
include("../../Models/Material/Material_Models/Zero_Energy_Control/global_control.jl")
include("../../Models/Model_Factory.jl")
include("../BC_manager.jl")
include("Verlet_solver.jl")
include("Matrix_linear_static.jl")
include("Static_solver.jl")
using ..MPI_Communication: synch_responder_to_controller,
                           synch_controller_to_responder,
                           synch_controller_bonds_to_responder,
                           synch_controller_bonds_to_responder_flattened

include("../Influence_function.jl")

using .Model_Factory: init_models, read_properties
using .Boundary_Conditions: init_BCs
using .Verlet_Solver
using .Linear_static_matrix_based
using .FEM
using .Influence_Function

export init
export solver

"""
	init(params::Dict)

Initialize the solver

# Arguments
- `params::Dict`: The parameters
# Returns
- `block_nodes::Dict{Int64,Vector{Int64}}`: A dictionary mapping block IDs to collections of nodes.
- `bcs::Dict{Any,Any}`: A dictionary containing boundary conditions.
- `solver_options::Dict{String,Any}`: A dictionary containing solver options.
"""
function init(params::Dict,
              step_id::Union{Nothing,Int64} = nothing)
    solver_options = Dict()
    nnodes = Data_Manager.get_nnodes()
    num_responder = Data_Manager.get_num_responder()
    block_ids = Data_Manager.get_field("Block_Id")
    block_nodes_with_neighbors = get_block_nodes(block_ids, nnodes + num_responder)
    block_nodes = get_block_nodes(block_ids, nnodes)
    block_name_list,
    block_id_list = get_block_names_and_ids(params, block_ids,
                                            Data_Manager.get_mpi_active())
    Data_Manager.set_block_name_list(block_name_list)
    Data_Manager.set_block_id_list(block_id_list)
    density = Data_Manager.create_constant_node_scalar_field("Density", Float64)
    horizon = Data_Manager.create_constant_node_scalar_field("Horizon", Float64)
    if Data_Manager.fem_active()
        fem_block = Data_Manager.create_constant_node_scalar_field("FEM Block", Bool;
                                                                   default_value = false)
        fem_block = set_fem_block(params, block_nodes_with_neighbors, fem_block) # includes the neighbors
    end
    active_nodes = Data_Manager.create_constant_node_scalar_field("Active Nodes", Int64)
    update_nodes = Data_Manager.create_constant_node_scalar_field("Update Nodes", Int64)
    Data_Manager.create_constant_node_scalar_field("Update", Bool; default_value = true)
    density = set_density(params, block_nodes_with_neighbors, density) # includes the neighbors
    horizon = set_horizon(params, block_nodes_with_neighbors, horizon) # includes the neighbors
    set_angles(params, block_nodes_with_neighbors) # includes the Neighbors
    solver_params = isnothing(step_id) ? params["Solver"] :
                    get_solver_params(params, step_id)
    solver_options["Models"] = get_model_options(solver_params)
    solver_options["All Models"] = get_model_options(solver_params)
    solver_options["Calculation"] = get_calculation_options(solver_params)
    if !isnothing(step_id)
        for step in 1:Data_Manager.get_max_step()
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
    Data_Manager.create_constant_bond_scalar_state("Influence Function", Float64;
                                                   default_value = 1)
    for iblock in eachindex(block_nodes)
        Influence_Function.init_influence_function(block_nodes[iblock],
                                                   params["Discretization"])
    end
    Data_Manager.create_bond_scalar_state("Bond Damage", Float64; default_value = 1)
    @debug "Read properties"
    read_properties(params, "Material" in solver_options["All Models"])
    @debug "Init models"
    @timeit "init_models" init_models(params,
                                      block_nodes,
                                      solver_options,
                                      synchronise_field)
    @debug "Init Boundary Conditions"

    @timeit "init_BCs" bcs=init_BCs(params)
    # get name and checks if it is there
    solver_options["Solver"] = get_solver_name(solver_params)

    @info "Init " * get_solver_name(solver_params)
    if get_solver_name(solver_params) == "Verlet"
        @debug "Init " * get_solver_name(solver_params)
        @timeit "init_solver" Verlet_Solver.init_solver(solver_options,
                                                        solver_params,
                                                        bcs,
                                                        block_nodes)
    elseif solver_options["Solver"] == "Static"
        @timeit "init_solver" Static_Solver.init_solver(solver_options,
                                                        solver_params,
                                                        bcs,
                                                        block_nodes)
    elseif solver_options["Solver"] == "Linear Static Matrix Based"
        @timeit "init_solver" Linear_static_matrix_based.init_solver(solver_options,
                                                                     solver_params,
                                                                     bcs,
                                                                     block_nodes)
    elseif solver_options["Solver"] == "Verlet Matrix Based"
        @timeit "init_solver" Linear_static_matrix_based.init_solver(solver_options,
                                                                     solver_params,
                                                                     bcs,
                                                                     block_nodes)
    end

    if Data_Manager.fem_active()
        FEM.init_FEM(params)
        FEM.Coupling.init_coupling(1:Data_Manager.get_nnodes(),
                                   params)
    end
    if !Data_Manager.has_key("Active")
        active = Data_Manager.create_constant_node_scalar_field("Active", Bool;
                                                                default_value = true)
    end
    #TODO: sync active with Data_Manager

    remove_models(solver_options["Models"])

    @debug "Finished Init Solver"
    return block_nodes, bcs, solver_options
end

"""
	set_density(params::Dict, block_nodes::Dict, density::NodeScalarField{Float64})

Sets the density of the nodes in the dictionary.

# Arguments
- `params::Dict`: The parameters
- `block_nodes::Dict`: A dictionary mapping block IDs to collections of nodes
- `density::NodeScalarField{Float64}`: The density
# Returns
- `density::NodeScalarField{Float64}`: The density
"""
function set_density(params::Dict, block_nodes::Dict, density::NodeScalarField{Float64})
    for block in eachindex(block_nodes)
        density[block_nodes[block]] .= get_density(params, block)
    end
    return density
end

"""
	set_angles(params::Dict, block_nodes::Dict)

Sets the density of the nodes in the dictionary.

# Arguments
- `params::Dict`: The parameters
- `block_nodes::Dict`: A dictionary mapping block IDs to collections of nodes
"""
function set_angles(params::Dict, block_nodes::Dict)
    mesh_angles = false
    if "Angles" in Data_Manager.get_all_field_keys()
        Data_Manager.set_rotation(true)
        mesh_angles = true
    end
    if "Element Angles" in Data_Manager.get_all_field_keys()
        Data_Manager.set_element_rotation(true)
    end

    block_rotation = false
    dof = Data_Manager.get_dof()
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
        Data_Manager.set_rotation(true)
        angles = Data_Manager.create_constant_node_vector_field("Angles", Float64, dof)

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
	set_horizon(params::Dict, block_nodes::Dict, horizon::NodeScalarField{Float64})

Sets the horizon of the nodes in the dictionary.

# Arguments
- `params::Dict`: The parameters
- `block_nodes::Dict`: A dictionary mapping block IDs to collections of nodes
- `horizon::NodeScalarField{Float64}`: The horizon
# Returns
- `horizon::NodeScalarField{Float64}`: The horizon
"""
function set_horizon(params::Dict, block_nodes::Dict, horizon::NodeScalarField{Float64})
    for block in eachindex(block_nodes)
        horizon[block_nodes[block]] .= get_horizon(params, block)
    end
    return horizon
end

"""
	solver(solver_options::Dict{String,Any}, block_nodes::Dict{Int64,Vector{Int64}}, bcs::Dict{Any,Any}, outputs::Dict{Int64,Dict{}}, result_files::Vector{Any}, write_results, silent::Bool)

Runs the solver.

# Arguments
- `solver_options::Dict{String,Any}`: The solver options
- `block_nodes::Dict{Int64,Vector{Int64}}`: A dictionary mapping block IDs to collections of nodes
- `bcs::Dict{Any,Any}`: The boundary conditions
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
                outputs::Dict{Int64,Dict{}},
                result_files::Vector{Dict},
                write_results,
                silent::Bool)
    if solver_options["Solver"] == "Verlet"
        return Verlet_Solver.run_solver(solver_options,
                                        block_nodes,
                                        bcs,
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
                                        outputs,
                                        result_files,
                                        synchronise_field,
                                        write_results,
                                        compute_parabolic_problems_before_model_evaluation,
                                        compute_parabolic_problems_after_model_evaluation,
                                        silent)
    elseif solver_options["Solver"] == "Linear Static Matrix Based"
        return Linear_static_matrix_based.run_solver(solver_options,
                                                     block_nodes,
                                                     bcs,
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
	remove_models(solver_options::Vector{String})

Sets the active models to false if they are deactivated in the solver. They can be active, because they are defined as model and in the blocks.

# Arguments
- `solver_options::Vector{String}`: A dictionary of fields
"""
function remove_models(solver_options::Vector{String})
    check = replace.(solver_options .* " Model", "_" => " ")
    for active_model_name in keys(Data_Manager.get_active_models())
        if !(active_model_name in check)
            Data_Manager.remove_active_model(active_model_name)
        end
    end
end

## TODO generalize this interface

function compute_parabolic_problems_before_model_evaluation(active_nodes,
                                                            solver_options)
    if !("Thermal" in solver_options["Models"]) &&
       !("Thermal" in solver_options["All Models"])
        return
    end
    temperatureN = Data_Manager.get_field("Temperature", "N")
    temperatureNP1 = Data_Manager.get_field("Temperature", "NP1")
    deltaT = Data_Manager.get_field("Delta Temperature")
    if "Thermal" in solver_options["Models"]
        temperatureNP1[active_nodes] = temperatureN[active_nodes] + deltaT[active_nodes]
    else
        if "Thermal" in solver_options["All Models"]
            temperatureNP1[active_nodes] = temperatureN[active_nodes]
        end
    end
end
function compute_parabolic_problems_after_model_evaluation(active_nodes,
                                                           solver_options, dt)
    if !("Thermal" in solver_options["Models"]) &&
       !("Thermal" in solver_options["All Models"])
        return
    end
    deltaT = Data_Manager.get_field("Delta Temperature")
    flowNP1 = Data_Manager.get_field("Heat Flow", "NP1")
    density = Data_Manager.get_field("Density")
    heat_capacity = Data_Manager.get_field("Specific Heat Capacity")
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
