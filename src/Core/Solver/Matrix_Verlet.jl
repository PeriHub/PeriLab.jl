# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Matrix_Verlet
using LinearAlgebra
using SparseArrays
using TimerOutputs
using Printf
using PrettyTables
using Logging
using LinearAlgebra: lu
using ...Data_Manager
using ...PeriLabExceptions: @abort

using ...Helpers: check_inf_or_nan, matrix_style

using ...ParameterSpec: @params, register_solver
using ...InputDeck: SolverParams, SolverOptions, CalculationOptions, start_time, end_time
using ..Solver_Common: TimeWindow, SolverRun, SolverContext, StepClock, synchronise_field,
                       compute_parabolic_problems_before_model_evaluation,
                       compute_parabolic_problems_after_model_evaluation
using ...MPI_Communication: find_and_set_core_value_min, find_and_set_core_value_max
using ..Model_Factory: compute_models, compute_crititical_time_step
using ..Boundary_Conditions: apply_bc_dirichlet, apply_bc_neumann, find_bc_free_dof,
                             BoundaryCondition
using ...Logging_Module: print_table
include("../Model_reduction/Model_reduction.jl")
using .Model_reduction: init_reduce_model, ReducedState, setup_reduced_state,
                        n_modal, physical_part, modal_part,
                        pull_from_nodes!, add_from_nodes!, push_to_nodes!,
                        ModelReductionParams

@params "Explicit Verlet time integration with the stiffness matrix of correspondence materials, optionally with model reduction" struct VerletMatrixParams <:
                                                                                                                                         SolverOptions
    "Factor on the stable time step; values below 1 add a safety margin"
    safety_factor::Float64 = opt("Safety Factor"; default = 1.0, min = 0)
    "Time step to use instead of the stable time step; -1: compute it from the stable time step and the safety factor"
    fixed_dt::Float64 = opt("Fixed dt"; default = -1.0, quantity = :time)
    "Damps the velocity: the velocity of the previous step is scaled by one minus this value"
    numerical_damping::Float64 = opt("Numerical Damping"; default = 0.0, min = 0)
    model_reduction::Union{Nothing,ModelReductionParams} = opt("Model Reduction";
                                                               default = nothing)
end
__init__() = register_solver("Verlet Matrix Based", VerletMatrixParams)
include("../../Compute/compute_field_values.jl")
using ..Correspondence_matrix_based: build_mass_matrix, init_model, init_matrix,
                                     compute_model
export init_solver
import ..Solver_Common: prepare_solver, time_step!

"""
    init_solver(options::VerletMatrixParams, params::SolverParams, models::Vector{String}, bcs::Dict{String,BoundaryCondition}, block_nodes::Dict{Int64,Vector{Int64}})

Initializes the matrix based Verlet solver: stiffness and mass matrix, the
optional model reduction and the time discretization. Returns the `TimeWindow`
and `nothing` (the run reads `Model Reduction` from `options`).
"""
function init_solver(options::VerletMatrixParams,
                     params::SolverParams,
                     models::Vector{String},
                     bcs::Dict{String,BoundaryCondition},
                     block_nodes::Dict{Int64,Vector{Int64}})
    horizon = Data_Manager.get_field("Horizon")
    if Data_Manager.get_rank() > 1
        @warn "Implementation might not work for MPI. Especially for coupling. It has to be tested."
    end
    find_bc_free_dof(bcs)
    initial_time = start_time(params, Data_Manager.get_current_time())
    final_time = end_time(params, Data_Manager.get_current_time())
    nsteps = something(params.number_of_steps, 1)

    mechanical = "Material" in models
    thermal = "Thermal" in models

    safety_factor = options.safety_factor
    fixed_dt = options.fixed_dt
    min_dt = compute_crititical_time_step(block_nodes, mechanical, thermal)
    if fixed_dt == -1.0
        nsteps_temp,
        dt = get_integration_steps(initial_time, final_time,
                                   safety_factor * min_dt)
    else
        nsteps_temp, dt = get_integration_steps(initial_time, final_time, fixed_dt)
    end
    if nsteps == 1
        nsteps = nsteps_temp
    end

    dof = Data_Manager.get_dof()
    density = Data_Manager.get_field("Density")

    model_reduction = something(options.model_reduction, false)
    reduce = model_reduction !== false
    @info "Solver parameter"
    @info "Initial Time: $initial_time"
    @info "Final Time: $final_time"
    @info "Number of Steps: $nsteps"
    @info "dt: $dt"
    @info "Numerical Damping: $(options.numerical_damping)"

    K = Data_Manager.get_stiffness_matrix()

    for (block, nodes) in pairs(block_nodes)
        init_model(nodes, Data_Manager.get_block_models(block).material, block)
    end
    @timeit "init_matrix" init_matrix()
    K = Data_Manager.get_stiffness_matrix()

    density_mass = zeros(length(density) * Data_Manager.get_dof())
    for iID in eachindex(density_mass)
        density_mass[iID] = density[Int(ceil(iID / dof))]
    end

    if reduce
        @info "Model reduction enabled."
        @timeit "Model reduction" init_reduce_model(model_reduction, block_nodes,
                                                    density_mass)
        mass = Data_Manager.get_mass_matrix()
        Data_Manager.set_mass_matrix(lu(sparse(mass)))
        return TimeWindow(initial_time, final_time, dt, nsteps), nothing
    end
    # if not reduced the mass has to be created here.
    mass = Diagonal(density_mass)
    Data_Manager.set_mass_matrix(lu(sparse(mass)))
    # No material point region without a reduction -- compute_models must not
    # contribute anywhere, K alone drives the whole domain. See run_solver.
    Data_Manager.set_reduced_model_pd(Int64[])

    return TimeWindow(initial_time, final_time, dt, nsteps), nothing
end

"""
    MatrixVerletWorkspace

Fields, mass matrix factorization, stiffness matrix and the (possibly reduced)
state vectors of the matrix based Verlet solver.
"""
struct MatrixVerletWorkspace{TM,TK,TS}
    coor::Matrix{Float64}
    active_list::Vector{Bool}
    volume::Vector{Float64}
    external_forces::Matrix{Float64}
    external_force_densities::Matrix{Float64}
    numerical_damping::Float64
    M_fact::TM
    K::TK
    master_nodes::Vector{Int64}
    state::TS
    modal_active::Bool
end

"""
    prepare_solver(run::SolverRun{VerletMatrixParams}, ctx::SolverContext)

Fetches the fields and matrices and sets up the state vectors; with a model
reduction only the master nodes and the modal coordinates are integrated.
"""
function prepare_solver(run::SolverRun{VerletMatrixParams}, ctx::SolverContext)
    material = "Material" in run.models
    K = Data_Manager.get_stiffness_matrix()::AbstractMatrix{Float64}
    # Master nodes and state vectors. Without a reduction every node is a master and the
    # state has no modal part, so the step needs no case distinction: the reduced
    # run differs only in which nodes are masters and in the extra modal entries.
    master_nodes,
    state = setup_reduced_state(something(run.params.options.model_reduction, false), K)
    return MatrixVerletWorkspace(Data_Manager.get_field("Coordinates"),
                                 Data_Manager.get_field("Active"),
                                 Data_Manager.get_field("Volume"),
                                 material ? Data_Manager.get_field("External Forces") :
                                 zeros(0, 0),
                                 material ?
                                 Data_Manager.get_field("External Force Densities") :
                                 zeros(0, 0),
                                 run.params.options.numerical_damping,
                                 Data_Manager.get_mass_matrix(),
                                 K,
                                 master_nodes,
                                 state,
                                 n_modal(state) > 0)
end

"""
    time_step!(ws::MatrixVerletWorkspace, run::SolverRun{VerletMatrixParams}, ctx::SolverContext, clock::StepClock)

One explicit Verlet step with the internal forces from the stiffness matrix
(`f = -K q`); the modal coordinates of a reduced model are integrated with the nodes.
"""
function time_step!(ws::MatrixVerletWorkspace, run::SolverRun{VerletMatrixParams},
                    ctx::SolverContext, clock::StepClock)
    block_nodes = ctx.block_nodes
    bcs = ctx.bcs
    time = clock.time
    step_time = clock.step_time
    dt = clock.dt
    coor = ws.coor
    active_list = ws.active_list
    volume = ws.volume
    external_forces = ws.external_forces
    external_force_densities = ws.external_force_densities
    numerical_damping = ws.numerical_damping
    M_fact = ws.M_fact
    K = ws.K
    master_nodes = ws.master_nodes
    state = ws.state
    modal_active = ws.modal_active
    @inbounds @fastmath begin
        @timeit "solver start" begin
            uN::Matrix{Float64} = Data_Manager.get_field("Displacements", "N")
            uNP1::Matrix{Float64} = Data_Manager.get_field("Displacements", "NP1")
            forces::Matrix{Float64} = Data_Manager.get_field("Forces", "NP1")
            deformed_coorNP1::Matrix{Float64} = Data_Manager.get_field("Deformed Coordinates",
                                                                       "NP1")
            vN::Matrix{Float64} = Data_Manager.get_field("Velocity", "N")
            vNP1::Matrix{Float64} = Data_Manager.get_field("Velocity", "NP1")
            aN::Matrix{Float64} = Data_Manager.get_field("Acceleration", "N")
            aNP1::Matrix{Float64} = Data_Manager.get_field("Acceleration", "NP1")
            force_densities_NP1::Matrix{Float64} = Data_Manager.get_field("Force Densities",
                                                                          "NP1")
            # Force computation is split the same way regardless of whether a
            # reduction is active: K covers every node by default, and
            # compute_models covers exactly the material point nodes that have
            # been swapped out of K -- empty without a reduction, so
            # compute_models contributes nothing anywhere and K alone drives
            # the whole domain. master_nodes already covers every node in that
            # case too, so no case distinction is needed here either.
            active_nodes::Vector{Int64} = master_nodes
            active_list .= false
            active_list[Data_Manager.get_reduced_model_pd()] .= true
        end
        @timeit "compute Velocity" begin
            @. @views vNP1[active_nodes,
                           :] = (1 - numerical_damping) *
                                vN[active_nodes, :] +
                                0.5 * dt * aN[active_nodes, :]

            apply_bc_dirichlet(["Velocity"], bcs, time, step_time)

            @. @views uNP1[active_nodes,
                           :] = uN[active_nodes, :] +
                                dt * vNP1[active_nodes, :]

            # The modal coordinates have no node field, so they are integrated here.
            # Without them the condensed degrees of freedom carry no dynamics and
            # the modal reduction degenerates to a static one.
            if modal_active
                modal_part(state.q_dot,
                           state) .= (1 - numerical_damping) .*
                                     modal_part(state.q_dot, state) .+
                                     0.5 .* dt .*
                                     modal_part(state.q_ddot, state)
                modal_part(state.q, state) .+= dt .*
                                               modal_part(state.q_dot, state)
            end
        end

        @timeit "apply BC" apply_bc_dirichlet(["Displacements", "Temperature"],
                                              bcs,
                                              time,
                                              step_time) #-> Dirichlet
        @timeit "deformed coor" @views deformed_coorNP1 .= coor .+ uNP1

        @timeit "compute models" compute_models(block_nodes,
                                                dt,
                                                time,
                                                run.models,
                                                synchronise_field)

        @timeit "Force computations" begin
            # Displacements of the master nodes into the state; the modal part was
            # already advanced above and stays as it is.
            pull_from_nodes!(state.q, uNP1, master_nodes, state)

            # f = -K q over the whole state, modal coordinates included
            @timeit "Force matrix computations" mul!(state.f, K, state.q)
            state.f .*= -1.0

            # External contributions act on nodes only, so they go into the physical
            # part. A load on a condensed degree of freedom would have to be
            # projected onto the modes, which the reduction does not provide.
            @views @. force_densities_NP1[active_nodes,
                                          :] += external_force_densities[active_nodes,
                                                                         :] +
                                                external_forces[active_nodes,
                                                                :] /
                                                volume[active_nodes]
            add_from_nodes!(state.f, force_densities_NP1, master_nodes, state)

            push_to_nodes!(force_densities_NP1, state.f, master_nodes, state)

            @. @views forces[active_nodes,
                             :] = force_densities_NP1[active_nodes, :] *
                                  volume[active_nodes]
            check_inf_or_nan(force_densities_NP1, "Force Densities")
        end
        @timeit "Accelaration computation" begin
            #state.q_ddot .= M_fact \ state.f
            ldiv!(state.q_ddot, M_fact, state.f)
            push_to_nodes!(aNP1, state.q_ddot, master_nodes, state)
        end
    end
    return nothing
end

"""
	get_integration_steps(initial_time::Float64, end_time::Float64, dt::Float64)

Calculate the number of integration steps and the adjusted time step for a numerical integration process.

# Arguments
- `initial_time::Float64`: The initial time for the integration.
- `end_time::Float64`: The final time for the integration.
- `dt::Float64`: The time step size.

# Returns
A tuple `(nsteps, dt)` where:
- `nsteps::Int64`: The number of integration steps required to cover the specified time range.
- `dt::Float64`: The adjusted time step size to evenly divide the time range.

# Errors
- Throws an error if the `dt` is less than or equal to zero.
"""
function get_integration_steps(initial_time::Float64, end_time::Float64, dt::Float64)
    if !(0 < dt < 1e50)
        @abort "Time step $dt [s] is not valid"
        return nothing
    end
    nsteps::Int64 = ceil((end_time - initial_time) / dt)
    dt = (end_time - initial_time) / nsteps
    return nsteps, dt
end

end
