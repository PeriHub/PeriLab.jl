# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Solver_Manager
using TimerOutputs: @timeit
using ..Data_Manager
using ..PeriLabExceptions: @abort
using ..InputDeck: PeriLabInput, SolverParams, CalculationOptions, calculation_options,
                   solver_name, model_options, solver_step,
                   BlockParams, block_by_id, block_angles, block_names_and_ids
using ..Helpers
using ..ModuleLoader: find_registered_modules
using ..ParameterSpec: model_module
include("../../Models/Material/Material_Basis.jl")
include("../../FEM/FEM_Factory.jl")
include("../../Models/Material/Material_Models/Zero_Energy_Control/Zero_Energy_Control.jl")
include("../../Models/Material/Material_Models/Correspondence/Correspondence_matrix_based.jl")
include("../../Models/Model_Factory.jl")
include("../BC_manager.jl")
include("Solver_Common.jl")
using .Solver_Common

for file in find_registered_modules(@__DIR__, "register_solver")
    include(file)
end

include("../Influence_function.jl")

using .Model_Factory: init_models, read_properties
using .Boundary_Conditions: init_BCs, BoundaryCondition
using .Verlet_Solver
using .Linear_static_matrix_based
using .Newmark
using .Matrix_Verlet
using .FEM
using .Influence_Function

export init
export solver
export TimeWindow, SolverRun, SolverContext, with_steps, create_run

"""
	init(input::PeriLabInput, step_id::Int64)

Initialize the solver

# Arguments
- `input::PeriLabInput`: The typed input
- `step_id::Int64`: The solver step (`-1`: single `Solver`)
# Returns
- `block_nodes::Dict{Int64,Vector{Int64}}`: A dictionary mapping block IDs to collections of nodes.
- `bcs::Dict{String,BoundaryCondition}`: A dictionary containing boundary conditions.
- `solver_run::SolverRun`: the initialized solver step.
"""
function init(input::PeriLabInput, step_id::Int64)
    dof = Data_Manager.get_dof()
    nnodes = Data_Manager.get_nnodes()
    num_responder = Data_Manager.get_num_responder()
    block_ids = Data_Manager.get_field("Block_Id")
    block_nodes_with_neighbors = get_block_nodes(block_ids, nnodes + num_responder)
    block_nodes = get_block_nodes(block_ids, nnodes)
    blocks = input.sections.blocks
    block_name_list,
    block_id_list = block_names_and_ids(blocks, block_ids, Data_Manager.get_mpi_active())
    Data_Manager.set_block_name_list(block_name_list)
    Data_Manager.set_block_id_list(block_id_list)
    density = Data_Manager.create_constant_node_scalar_field("Density", Float64)
    horizon = Data_Manager.create_constant_node_scalar_field("Horizon", Float64)
    deformed_coorN,
    deformed_coorNP1 = Data_Manager.create_node_vector_field("Deformed Coordinates",
                                                             Float64, dof)
    deformed_coorN = copy(Data_Manager.get_field("Coordinates"))
    deformed_coorNP1 = copy(Data_Manager.get_field("Coordinates"))
    if Data_Manager.fem_active()
        fem_block = Data_Manager.create_constant_node_scalar_field("FEM Block", Bool;
                                                                   default_value = false)
        fem_block = set_fem_block(blocks, block_nodes_with_neighbors, fem_block) # includes the neighbors
    end
    Data_Manager.create_constant_node_scalar_field("Active Nodes", Int64)
    Data_Manager.create_constant_node_scalar_field("Update Nodes", Int64)
    Data_Manager.create_constant_node_scalar_field("Update", Bool; default_value = true)
    density = set_density(blocks, block_nodes_with_neighbors, density) # includes the neighbors
    horizon = set_horizon(blocks, block_nodes_with_neighbors, horizon) # includes the neighbors
    set_angles(blocks, block_nodes_with_neighbors) # includes the Neighbors
    solver_params = solver_step(input, step_id)
    models = model_options(solver_params)
    all_models = copy(models)
    calculation = calculation_options(solver_params)
    if step_id != -1
        for step in 1:Data_Manager.get_max_step()
            step_solver_params = solver_step(input, step)
            append!(all_models, model_options(step_solver_params))
            calculation = calculation | calculation_options(step_solver_params)
        end
        all_models = unique(all_models)
    end
    Data_Manager.create_constant_bond_scalar_state("Influence Function", Float64;
                                                   default_value = 1)
    for iblock in eachindex(block_nodes)
        Influence_Function.init_influence_function(block_nodes[iblock],
                                                   input.sections.discretization.influence_function)
    end
    Data_Manager.create_bond_scalar_state("Bond Damage", Float64; default_value = 1)
    @debug "Read properties"
    read_properties(input, "Material" in models)
    @debug "Init models"
    @timeit "init_models" init_models(input,
                                      block_nodes,
                                      models,
                                      all_models,
                                      calculation,
                                      synchronise_field)
    @debug "Init Boundary Conditions"

    @timeit "init_BCs" bcs=init_BCs(input.sections.boundary_conditions)
    @info "Init " * solver_name(solver_params)
    if !Data_Manager.has_key("Active")
        active = Data_Manager.create_constant_node_scalar_field("Active", Bool;
                                                                default_value = true)
    end
    solver_run = create_run(solver_params, models, all_models, calculation, bcs, block_nodes)

    if Data_Manager.fem_active()
        @timeit "init_FEM" FEM.init_FEM(input.sections.fem,
                                        input.materials)
        @timeit "init_coupling" FEM.Coupling.init_coupling(1:Data_Manager.get_nnodes(),
                                                           input.sections.fem)
    end

    #TODO: sync active with Data_Manager

    remove_models(models)

    @debug "Finished Init Solver"
    return block_nodes, bcs, solver_run
end

"""
    create_run(params, models, all_models, calculation, bcs, block_nodes)

Initializes the solver selected by `params.options` (its module's
`init_solver`) and returns the `SolverRun` of the step.
"""
function create_run(params::SolverParams, models::Vector{String},
                    all_models::Vector{String}, calculation::CalculationOptions,
                    bcs::Dict{String,BoundaryCondition},
                    block_nodes::Dict{Int64,Vector{Int64}})
    mod = model_module(params.options)
    # invokelatest: the solver may come from a module loaded at runtime (licensed modules)
    time, state = Base.invokelatest(mod.init_solver, params.options, params, models, bcs,
                                    block_nodes)
    return SolverRun(params, models, all_models, calculation, time, state)
end

"""
	set_density(blocks::Dict{String,BlockParams}, block_nodes::Dict, density::NodeScalarField{Float64})

Sets the density of the nodes in the dictionary.

# Arguments
- `blocks::Dict{String,BlockParams}`: The blocks of the input deck
- `block_nodes::Dict`: A dictionary mapping block IDs to collections of nodes
- `density::NodeScalarField{Float64}`: The density
# Returns
- `density::NodeScalarField{Float64}`: The density
"""
function set_density(blocks::Dict{String,BlockParams}, block_nodes::Dict,
                     density::NodeScalarField{Float64})
    for block in eachindex(block_nodes)
        density[block_nodes[block]] .= block_by_id(blocks, block)[2].density
    end
    return density
end

"""
	set_angles(blocks::Dict{String,BlockParams}, block_nodes::Dict)

Sets the density of the nodes in the dictionary.

# Arguments
- `blocks::Dict{String,BlockParams}`: The blocks of the input deck
- `block_nodes::Dict`: A dictionary mapping block IDs to collections of nodes
"""
function set_angles(blocks::Dict{String,BlockParams}, block_nodes::Dict)
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
        if block_angles(block_by_id(blocks, block)..., dof) !== nothing
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
            angles_global = block_angles(block_by_id(blocks, block)..., dof)
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
	set_fem_block(blocks::Dict{String,BlockParams}, block_nodes::Dict, fem_block::Vector{Bool})

Sets the fem_block of the nodes in the dictionary.

# Arguments
- `blocks::Dict{String,BlockParams}`: The blocks of the input deck
- `block_nodes::Dict`: A dictionary mapping block IDs to collections of nodes
- `fem_block::Vector{Bool}`: The fem_block
# Returns
- `fem_block::Vector{Bool}`: The fem_block
"""
function set_fem_block(blocks::Dict{String,BlockParams}, block_nodes::Dict,
                       fem_block::Vector{Bool})
    for block in eachindex(block_nodes)
        fem_block[block_nodes[block]] .= something(block_by_id(blocks, block)[2].fem, false)
    end
    return fem_block
end

"""
	set_horizon(blocks::Dict{String,BlockParams}, block_nodes::Dict, horizon::NodeScalarField{Float64})

Sets the horizon of the nodes in the dictionary.

# Arguments
- `blocks::Dict{String,BlockParams}`: The blocks of the input deck
- `block_nodes::Dict`: A dictionary mapping block IDs to collections of nodes
- `horizon::NodeScalarField{Float64}`: The horizon
# Returns
- `horizon::NodeScalarField{Float64}`: The horizon
"""
function set_horizon(blocks::Dict{String,BlockParams}, block_nodes::Dict,
                     horizon::NodeScalarField{Float64})
    for block in eachindex(block_nodes)
        horizon[block_nodes[block]] .= block_by_id(blocks, block)[2].horizon
    end
    return horizon
end

"""
    solver(run::SolverRun, ctx::SolverContext)

Runs the solver step `run`; returns the updated result files.
"""
function solver(run::SolverRun, ctx::SolverContext)
    return Base.invokelatest(model_module(run.params.options).run_solver, run, ctx)
end

"""
	remove_models(models::Vector{String})

Sets the active models to false if they are deactivated in the solver. They can be active, because they are defined as model and in the blocks.

# Arguments
- `models::Vector{String}`: The model categories of the step
"""
function remove_models(models::Vector{String})
    check = replace.(models .* " Model", "_" => " ")
    for active_model_name in keys(Data_Manager.get_active_models())
        if !(active_model_name in check)
            Data_Manager.remove_active_model(active_model_name)
        end
    end
end


end
