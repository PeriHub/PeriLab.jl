# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Static_Solver
using LinearAlgebra
using NLsolve
using Printf
using LoopVectorization
using Logging
using TimerOutputs: @timeit

using ...Data_Manager
using ...Helpers: check_inf_or_nan, find_active_nodes, progress_bar, matrix_style
using ...ParameterSpec: @params, register_solver
using ...InputDeck: SolverParams, SolverOptions, start_time, end_time
using ..Solver_Common: TimeWindow, SolverRun, SolverContext, StepClock, StepFailed,
                       synchronise_field,
                       compute_parabolic_problems_before_model_evaluation,
                       compute_parabolic_problems_after_model_evaluation

@params "Quasi-static solver: finds the equilibrium of each step with a nonlinear solver (NLsolve, Anderson acceleration)" struct StaticParams <: SolverOptions
    "Not used by this solver"
    safety_factor::Float64 = opt("Safety Factor"; default = 1.0, min = 0)
    "Time step; -1: the step is split into Number of Steps increments"
    fixed_dt::Float64 = opt("Fixed dt"; default = -1.0, quantity = :time)
    "Not used by this solver"
    numerical_damping::Float64 = opt("Numerical Damping"; default = 0.0, min = 0)
    "Convergence limit on the change of the displacements between two iterations"
    solution_tolerance::Float64 = opt("Solution tolerance"; default = 1e-7, min = 0)
    "Convergence limit on the scaled force residual"
    residual_tolerance::Float64 = opt("Residual tolerance"; default = 1e-7, min = 0)
    "Iterations per step before the step counts as failed"
    maximum_number_of_iterations::Int64 = opt("Maximum number of iterations"; default = 100,
                                              min = 1)
    "Prints the iterations of the nonlinear solver"
    show_solver_iteration::Bool = opt("Show solver iteration"; default = false)
    "Divisor of the force residual before it is compared with the residual tolerance"
    residual_scaling::Float64 = opt("Residual scaling"; default = 1e6)
    "Number of previous iterations used by the Anderson acceleration"
    m::Int64 = opt("m"; default = 15, min = 1)
    "Start guess of the displacements: per direction a value at the minimum and one at the maximum coordinate, separated by spaces, interpolated linearly; zero if not given"
    linear_start_value::Union{Nothing,String} = opt("Linear Start Value"; default = nothing)
    "Not used by this solver"
    nlsolve::Union{Nothing,Bool} = opt("NLSolve"; default = nothing)
    "Not used by this solver"
    solver_type::Union{Nothing,String} = opt("Solver Type"; default = nothing)
end
__init__() = register_solver("Static", StaticParams)
# using ...MPI_Communication: barrier
using ..Model_Factory
using ..Boundary_Conditions: apply_bc_dirichlet, apply_bc_neumann, find_bc_free_dof, BoundaryCondition
using ...Logging_Module: print_table

export init_solver
import ..Solver_Common: prepare_solver, time_step!
"""
    init_solver(options::StaticParams, params::SolverParams, models::Vector{String}, bcs::Dict{String,BoundaryCondition}, block_nodes::Dict{Int64,Vector{Int64}})

Initialize the Static solver for a simulation.

This function sets up the Static solver for a simulation by initializing various parameters.

# Arguments
- `options::StaticParams`: The Static options of the step (tolerances, iterations, start values).
- `params::SolverParams`: The solver parameters of the current step.
- `models::Vector{String}`: The model categories of the step.
- `bcs::Dict{String,BoundaryCondition}`: Boundary conditions
- `block_nodes::Dict{Int64,Vector{Int64}}`: A dictionary mapping block IDs to collections of nodes.

# Returns
- `TimeWindow`: initial and final time, time step and number of steps.
- `nothing`: the run reads its settings from `options`.

# Dependencies
This function may depend on the following functions:
- `start_time` and the fields of `params` (final time, safety factor, fixed time step): The simulation parameters.
- `get_integration_steps`: Used to determine the number of integration steps and adjust the time step.
"""
function init_solver(options::StaticParams,
                     params::SolverParams,
                     models::Vector{String},
                     bcs::Dict{String,BoundaryCondition},
                     block_nodes::Dict{Int64,Vector{Int64}})
    # @info "==============================="
    # @info "==== NLsolve Static Solver ===="
    # @info "==============================="

    mechanical = "Material" in models
    thermal = "Thermal" in models
    initial_time = start_time(params, Data_Manager.get_current_time())
    final_time = end_time(params, Data_Manager.get_current_time())
    fixed_dt = options.fixed_dt
    dof = Data_Manager.get_dof()
    if fixed_dt == -1.0
        if params.number_of_steps !== nothing
            nsteps = params.number_of_steps
            dt = (final_time - initial_time) / nsteps
        else
            nsteps = Int64(1)
            dt = final_time - initial_time
        end
    else
        if params.number_of_steps !== nothing
            @warn "''Number of Steps'' and ''Fixed dt'' are defined. ''Fixed dt'' is used and ''Number of Steps'' from yaml is ignored."
        end
        nsteps = Int64(round((final_time - initial_time) / fixed_dt))
        dt = fixed_dt
    end
    comm = Data_Manager.get_comm()

    # not needed here
    static = options
    solver_specifics = Dict("Solution tolerance" => static.solution_tolerance,
                            "Residual tolerance" => static.residual_tolerance,
                            "Maximum number of iterations" => static.maximum_number_of_iterations,
                            "Show solver iteration" => static.show_solver_iteration,
                            "Residual scaling" => static.residual_scaling,
                            "m" => static.m,
                            "Linear Start Value" => static.linear_start_value === nothing ?
                                                    zeros(2 * dof) :
                                                    parse.(Float64,
                                                           split(static.linear_start_value)))

    residual = Data_Manager.create_constant_node_vector_field("Residual", Float64, dof)

    if !("Start_Values" in Data_Manager.get_all_field_keys())
        start_u = Data_Manager.create_constant_node_vector_field("Start_Values", Float64,
                                                                 dof)
        coor = Data_Manager.get_field("Coordinates")
        ls = solver_specifics["Linear Start Value"]
        if ls[1] == ls[2]
            start_u .= ls[1]
        else
            for idof in 1:dof
                m = (ls[2*idof] - ls[2*idof-1]) /
                    (maximum(coor[:, idof]) - minimum(coor[:, idof]))
                n = ls[2*idof] - m * maximum(coor[:, idof])
                start_u[:, idof] = (m .* coor[:, idof] .+ n) ./ nsteps .* final_time
                n = ls[2*idof] - m * maximum(coor[:, idof])
                start_u[:, idof] = (m .* coor[:, idof] .+ n) ./ nsteps
            end
        end
    else
        start_u = Data_Manager.get_field("Start_Values")
    end
    check_inf_or_nan(start_u, "Start_Values")
    find_bc_free_dof(bcs)

    if Data_Manager.get_rank() == 0
        data = ["Static Solver" "" "" ""
                "Initial time" "."^10 initial_time "s"
                "Final time" "."^10 final_time "s"
                "Maximum number of iterations" "."^10 solver_specifics["Maximum number of iterations"] ""
                "Solution tolerance" "."^10 solver_specifics["Solution tolerance"] ""
                "Residual tolerance" "."^10 solver_specifics["Residual tolerance"] ""]
        if fixed_dt != -1.0
            data = vcat(data,
                        ["Fixed time increment" "."^10 fixed_dt "s";])
        else
            data = vcat(data,
                        ["Time increment" "."^10 dt "s"])
        end
        data = vcat(data,
                    ["Number of steps" "."^10 nsteps ""])

        print_table(data)
    end

    return TimeWindow(initial_time, final_time, dt, nsteps), nothing
end

"""
    StaticWorkspace

Fields and nlsolve settings of the static solver; `start_u` is the start value
of the next step.
"""
mutable struct StaticWorkspace
    bc_free_dof::Vector{Int64}
    scaling::Float64
    xtol::Float64
    ftol::Float64
    iterations::Int64
    show_trace::Bool
    m::Int64
    start_u::Matrix{Float64}
    uNP1::Matrix{Float64}
    forces::Matrix{Float64}
    volume::Vector{Float64}
    active_list::Vector{Bool}
end

"""
    prepare_solver(run::SolverRun{StaticParams}, ctx::SolverContext)

Fetches the fields and the nlsolve settings of the static steps.
"""
function prepare_solver(run::SolverRun{StaticParams}, ctx::SolverContext)
    options = run.params.options
    return StaticWorkspace(convert(Vector{Int64}, Data_Manager.get_bc_free_dof()),
                           options.residual_scaling,
                           options.solution_tolerance,
                           options.residual_tolerance,
                           options.maximum_number_of_iterations,
                           options.show_solver_iteration,
                           options.m,
                           Data_Manager.get_field("Start_Values"),
                           Data_Manager.get_field("Displacements", "NP1"),
                           Data_Manager.get_field("Forces", "NP1"),
                           Data_Manager.get_field("Volume"),
                           Data_Manager.get_field("Active"))
end

"""
    time_step!(ws::StaticWorkspace, run::SolverRun{StaticParams}, ctx::SolverContext, clock::StepClock)

One static step: Dirichlet conditions, then the equilibrium solved with nlsolve
(Anderson). Returns `StepFailed` if nlsolve does not converge.
"""
function time_step!(ws::StaticWorkspace, run::SolverRun{StaticParams}, ctx::SolverContext,
                    clock::StepClock)
    block_nodes = ctx.block_nodes
    bcs = ctx.bcs
    time = clock.time
    step_time = clock.step_time
    dt = clock.dt
    apply_bc_dirichlet(["Forces", "Force Densities"], bcs, time, step_time)
    active_nodes = Data_Manager.get_field("Active Nodes")
    active_nodes = find_active_nodes(ws.active_list, active_nodes,
                                     1:Data_Manager.get_nnodes())
    compute_parabolic_problems_before_model_evaluation(active_nodes, run)
    apply_bc_dirichlet(["Displacements", "Temperature"], bcs, time, step_time + dt)
    sol = nlsolve((residual, U) -> residual!(residual, U, ws.bc_free_dof, block_nodes, dt,
                                             time, run.models, synchronise_field,
                                             ws.scaling),
                  ws.start_u;
                  xtol = ws.xtol,
                  ftol = ws.ftol,
                  iterations = ws.iterations,
                  show_trace = ws.show_trace & !ctx.silent,
                  extended_trace = false,
                  method = :anderson,
                  m = ws.m,)
    compute_parabolic_problems_after_model_evaluation(active_nodes, run, dt)
    ws.start_u = copy(ws.uNP1)
    if !sol.x_converged && !sol.f_converged
        return StepFailed("no convergence after $(ws.iterations) iterations")
    end
    active_nodes = find_active_nodes(ws.active_list, active_nodes,
                                     1:Data_Manager.get_nnodes())
    force_densities = Data_Manager.get_field("Force Densities", "NP1")
    @. @views ws.forces[active_nodes, :] = force_densities[active_nodes, :] *
                                           ws.volume[active_nodes]
    return nothing
end

function residual!(residual,
                   U,
                   bc_free_dof,
                   block_nodes,
                   dt,
                   time,
                   models::Vector{String},
                   synchronise_field,
                   scaling)
    active_list = Data_Manager.get_field("Active")

    deformed_coorNP1 = Data_Manager.get_field("Deformed Coordinates", "NP1")

    force_densities = Data_Manager.get_field("Force Densities", "NP1")
    coor = Data_Manager.get_field("Coordinates")
    active_nodes = Data_Manager.get_field("Active Nodes")
    active_nodes = find_active_nodes(active_list, active_nodes, 1:Data_Manager.get_nnodes())
    external_force_densities = Data_Manager.get_field("External Force Densities")
    # one step more, because of init step (time = 0)
    # displacement  from interface
    # Set bond damages to original state
    # Data_Manager.switch_bonds!(
    #     Data_Manager.get_field("Bond Damage", "NP1"),
    #     Data_Manager.get_field("Bond Damage", "N"),
    # )

    uNP1 = Data_Manager.get_field("Displacements", "NP1")
    uNP1[bc_free_dof] .= U[bc_free_dof]

    # bc_dof = setdiff(1:length(uNP1), bc_free_dof)

    @views deformed_coorNP1[active_nodes,
    :] = coor[active_nodes, :] .+
                                               uNP1[active_nodes, :]

    force_densities[:, :] .= 0 # TODO check where to put it for iterative solver
    forces = Data_Manager.get_field("Forces", "NP1")
    forces[:, :] .= 0

    Data_Manager.synch_manager(synchronise_field, "upload_to_cores")
    # synch
    Model_Factory.compute_models(block_nodes,
                                 dt,
                                 time,
                                 models,
                                 synchronise_field)

    Data_Manager.synch_manager(synchronise_field, "download_from_cores")
    # synch

    # @timeit "apply_bc_neumann" apply_bc_neumann(bcs, time) #-> von neumann
    #active_nodes = Data_Manager.get_field("Active Nodes")
    #active_nodes =            find_active_nodes(active_list, active_nodes, 1:Data_Manager.get_nnodes())

    # check_inf_or_nan(force_densities, "Forces")

    # Richtung muss mit rein
    #println("adsads", force_densities[bc_free_dof] + external_force_densities[bc_free_dof])
    residual .= 0
    residual[bc_free_dof] = (force_densities[bc_free_dof] +
                             external_force_densities[bc_free_dof]) ./ scaling
end

end
