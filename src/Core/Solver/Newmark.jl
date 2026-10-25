# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Newmark
using Printf
using SparseArrays
using LinearAlgebra
using TimerOutputs: @timeit

using ...Data_Manager
using ...Helpers: check_inf_or_nan, find_active_nodes
using ..Boundary_Conditions: apply_bc_dirichlet, apply_bc_neumann, find_bc_free_dof, BoundaryCondition
using ...ParameterSpec: @params, register_solver
using ...InputDeck: SolverParams, SolverOptions, start_time, end_time
using ..Solver_Common: TimeWindow, SolverRun, SolverContext, StepClock, synchronise_field,
                       compute_parabolic_problems_before_model_evaluation,
                       compute_parabolic_problems_after_model_evaluation

@params "Implicit Newmark time integration with the stiffness matrix of correspondence materials" struct NewmarkParams <: SolverOptions
    "Rebuilds the stiffness matrix in every step, e.g. for damage or large deformations; otherwise the initial matrix is reused"
    matrix_update::Bool = opt("Matrix Update"; default = false)
    "Newmark parameter weighting the accelerations in the velocity update (often called gamma); 0.5 adds no numerical damping"
    newmark_delta::Float64 = opt("Newmark Delta"; default = 0.5)
    "Newmark parameter weighting the accelerations in the displacement update (often called beta); 0.25 (0.5 + delta)^2 if not given"
    newmark_alpha::Union{Nothing,Float64} = opt("Newmark Alpha"; default = nothing)
end
__init__() = register_solver("Newmark", NewmarkParams)
using ..Model_Factory: compute_stiff_matrix_compatible_models,
                       compute_matrix_based_bond_forces

using ..Correspondence_matrix_based: build_mass_matrix, init_model, init_matrix,
                                     compute_model
using ..Model_Factory.Pre_Calculation.Bond_Deformation

export init_solver
import ..Solver_Common: prepare_solver, time_step!


# ────────────────────────────────────────────────────────────
# Cache
# ────────────────────────────────────────────────────────────
mutable struct NewmarkCache
    F_eff::Vector{Float64}
    bc_mask::BitVector
    temp::Vector{Float64}
    u_free::Vector{Float64}
    K_eff_lu::Union{Nothing,Any}
    last_non_BCs::Vector{Int}
    bc_dofs::Vector{Int}
    K_bc::Union{Nothing,SparseMatrixCSC{Float64,Int}}
    F_bc::Vector{Float64}

    function NewmarkCache(n::Int)
        new(Float64[], BitVector(undef, n), Float64[], Float64[],
            nothing, Int[],
            Int[], nothing, Float64[])
    end
end

# ────────────────────────────────────────────────────────────
# Init
# ────────────────────────────────────────────────────────────
"""
    init_solver(options::NewmarkParams, params, models, bcs, block_nodes)

Initialise Newmark-beta solver. Returns the `TimeWindow` and the state
`(alpha = …, beta = …)`.

Default: average acceleration (beta=0.5,alpha=0.25, unconditionally stable) taken from Bathe page 512.
With numerical damping alpha: beta >= 0.5 ,alpha >= 0.25*(0.5+beta)^2
"""
function init_solver(options::NewmarkParams,
                     params::SolverParams,
                     models::Vector{String},
                     bcs::Dict{String,BoundaryCondition},
                     block_nodes::Dict{Int64,Vector{Int64}})
    find_bc_free_dof(bcs)

    Data_Manager.create_constant_node_vector_field("Delta Displacements", Float64,
                                                   Data_Manager.get_dof())

    initial_time = start_time(params, Data_Manager.get_current_time())
    final_time = end_time(params, Data_Manager.get_current_time())
    nsteps = something(params.number_of_steps, 1)
    dt = (final_time - initial_time) / nsteps

    beta = options.newmark_delta # name is delta because of beta. Internally it is beta.
    alpha = something(options.newmark_alpha, 0.25 * (0.5 + beta)^2)

    for (block, nodes) in pairs(block_nodes)
        init_model(nodes, Data_Manager.get_block_models(block).material, block)
    end
    @timeit "init matrix" init_matrix()
    density = Data_Manager.get_field("Density")

    Data_Manager.set_mass_matrix(build_mass_matrix(density, Data_Manager.get_nnodes(),
                                                   Data_Manager.get_dof()))
    Data_Manager.create_node_scalar_field("Damage", Float64)
    coor = Data_Manager.get_field("Coordinates")
    Data_Manager.get_field("Deformed Coordinates", "N") .= coor
    Data_Manager.get_field("Deformed Coordinates", "NP1") .= coor
    return TimeWindow(initial_time, final_time, dt, nsteps), (alpha = alpha, beta = beta)
end

# ────────────────────────────────────────────────────────────
# Run
# ────────────────────────────────────────────────────────────
"""
    NewmarkWorkspace

Newmark coefficients, mass and stiffness matrix, external loads and the solver
cache; `K` is replaced when the matrix is updated, `active_nodes` after each step.
"""
mutable struct NewmarkWorkspace{TM}
    volume::Vector{Float64}
    coor::Matrix{Float64}
    external_forces::Matrix{Float64}
    external_force_densities::Matrix{Float64}
    a0::Float64
    a2::Float64
    a3::Float64
    a6::Float64
    a7::Float64
    matrix_update::Bool
    nnodes::Int64
    dof::Int64
    nodes::Vector{Int64}
    M::TM
    K::AbstractMatrix{Float64}
    non_BCs::Vector{Int64}
    active_nodes::AbstractVector{Int64}
    cache::NewmarkCache
end

"""
    prepare_solver(run::SolverRun{NewmarkParams}, ctx::SolverContext)

Computes the Newmark coefficients and fetches matrices and fields.
"""
function prepare_solver(run::SolverRun{NewmarkParams}, ctx::SolverContext)
    dt = run.time.dt
    nsteps = run.time.nsteps
    beta::Float64 = run.state.beta
    alpha::Float64 = run.state.alpha
    matrix_update = run.params.options.matrix_update
    nnodes = Data_Manager.get_nnodes()
    dof = Data_Manager.get_dof()
    nodes = collect(1:nnodes)
    # following Bathe in notation page 512 in his FEM book. using Beta instead of delta, because delta is already the horizon
    # Newmark coefficients
    a0 = 1.0 / (alpha * dt^2)
    #a1 = beta / (alpha * dt) -> not needed, because no damping matrix
    a2 = 1.0 / (alpha * dt)
    a3 = 1.0 / (2.0 * alpha) - 1.0
    #a4 = beta / (alpha) - 1.0 -> not needed, because no damping matrix
    #a5 = dt / 2.0 * (beta / alpha - 2.0)-> not needed, because no damping matrix
    a6 = dt * (1.0 - beta)
    a7 = beta * dt

    @info "Newmark: alpha=$alpha, gamma=$beta, dt=$dt, nsteps=$nsteps, matrix update: $matrix_update"

    active_nodes = Data_Manager.get_field("Active Nodes")
    active_list = Data_Manager.get_field("Active")
    active_nodes = find_active_nodes(active_list, active_nodes, nodes)

    return NewmarkWorkspace(Data_Manager.get_field("Volume"),
                            Data_Manager.get_field("Coordinates"),
                            Data_Manager.get_field("External Forces"),
                            Data_Manager.get_field("External Force Densities"),
                            a0, a2, a3, a6, a7, matrix_update, nnodes, dof, nodes,
                            Data_Manager.get_mass_matrix(),
                            Data_Manager.get_stiffness_matrix(),
                            convert(Vector{Int64}, Data_Manager.get_bc_free_dof()),
                            active_nodes,
                            NewmarkCache(nnodes * dof))
end

"""
    time_step!(ws::NewmarkWorkspace, run::SolverRun{NewmarkParams}, ctx::SolverContext, clock::StepClock)

One step of the Newmark-beta implicit time integration for structural dynamics.

# Method

Solves the semi-discrete equation of motion

``M \\ddot{u} + K u = F_\\mathrm{ext}``

using the Newmark-beta family of time integrators. At each step the effective
system

``K_\\mathrm{eff}\\, u_{n+1} = F_\\mathrm{eff}``

is solved, with

``K_\\mathrm{eff} = K + \\frac{1}{\\beta\\,\\Delta t^2}\\, M``

``F_\\mathrm{eff} = F_\\mathrm{ext} - F_\\mathrm{int}
  + M\\!\\left(\\frac{1}{\\beta\\,\\Delta t^2}\\, u_n
  + \\frac{1}{\\beta\\,\\Delta t}\\, v_n
  + \\left(\\frac{1}{2\\beta} - 1\\right) a_n\\right)``

followed by the acceleration and velocity updates

``a_{n+1} = \\frac{1}{\\beta\\,\\Delta t^2}(u_{n+1} - u_n)
  - \\frac{1}{\\beta\\,\\Delta t}\\, v_n
  - \\left(\\frac{1}{2\\beta} - 1\\right) a_n``

``v_{n+1} = v_n + \\Delta t\\,(1-\\gamma)\\, a_n + \\gamma\\,\\Delta t\\, a_{n+1}``

# Parameters (`run.params.options` and `run.state`)

| Parameter        | Default                | Description                          |
|:-----------------|:-----------------------|:-------------------------------------|
| `Newmark Beta`   | `0.25 (1+alpha)²`         | Newmark beta                            |
| `Newmark Gamma`  | `0.5 + alpha`             | Newmarkgamma                            |
| `Matrix Update`  | `false`                | Reassemble K every step              |

The numerical damping parameter `alpha` (HHT-alpha) is read from the input deck
via `get_numerical_damping`. For `alpha = 0` the method reduces to the
average acceleration scheme (beta = 1/4,gamma = 1/2), which is unconditionally
stable, second-order accurate, and energy-conserving.

# Mass matrix

Lumped diagonal mass: ``M_i = \\rho_i\\, V_i``.

# Arguments
- `ws::NewmarkWorkspace`:  coefficients, matrices and cache from `prepare_solver`
- `run::SolverRun{NewmarkParams}`:  the solver step (time discretization, models, alpha and beta in `run.state`)
- `ctx::SolverContext`:  block nodes and boundary conditions
- `clock::StepClock`:  time and increment of the step

# References
- Newmark, N. M. (1959). "A method of computation for structural dynamics."
  *J. Eng. Mech. Div.*, 85(3), 67–94.
- Hughes, T. J. R. (2000). *The Finite Element Method*, §9.1–9.3.
"""
function time_step!(ws::NewmarkWorkspace, run::SolverRun{NewmarkParams}, ctx::SolverContext,
                    clock::StepClock)
    block_nodes = ctx.block_nodes
    bcs = ctx.bcs
    time = clock.time
    step_time = clock.step_time
    dt = clock.dt
    volume = ws.volume
    coor = ws.coor
    external_forces = ws.external_forces
    external_force_densities = ws.external_force_densities
    a0, a2, a3, a6, a7 = ws.a0, ws.a2, ws.a3, ws.a6, ws.a7
    matrix_update = ws.matrix_update
    nnodes = ws.nnodes
    dof = ws.dof
    nodes = ws.nodes
    M = ws.M
    K = ws.K
    non_BCs = ws.non_BCs
    active_nodes = ws.active_nodes
    cache = ws.cache
    uN = Data_Manager.get_field("Displacements", "N")
    uNP1 = Data_Manager.get_field("Displacements", "NP1")
    velN = Data_Manager.get_field("Velocity", "N")
    velNP1 = Data_Manager.get_field("Velocity", "NP1")
    aN = Data_Manager.get_field("Acceleration", "N")
    aNP1 = Data_Manager.get_field("Acceleration", "NP1")
    deformed_coorNP1 = Data_Manager.get_field("Deformed Coordinates", "NP1")
    forces = Data_Manager.get_field("Forces", "NP1")
    force_densities_NP1 = Data_Manager.get_field("Force Densities", "NP1")

    deformed_coorNP1 .= coor

    compute_parabolic_problems_before_model_evaluation(active_nodes, run)
    apply_bc_dirichlet([
                           "Displacements",
                           "Forces",
                           "Force Densities",
                           "Temperature"
                       ],
                       bcs, time, step_time)
    @. external_force_densities += external_forces / volume
    @. deformed_coorNP1 = coor + uN
    compute_matrix_based_bond_forces(block_nodes, time, dt) #-> für schädigung davor
    @timeit "compute_models" compute_stiff_matrix_compatible_models(block_nodes,
                                                                    dt, time,
                                                                    run.models,
                                                                    synchronise_field)
    # Stiffness
    if matrix_update
        @timeit "update stiffness" begin
            compute_matrix(collect(1:nnodes))
            K = Data_Manager.get_stiffness_matrix()
            cache.K_eff_lu = nothing
        end
    end

    # Solve
    @timeit "newmark_solve" @views newmark_step!(-K, M, non_BCs,
                                                 uN, uNP1, velN, velNP1, aN, aNP1,
                                                 force_densities_NP1,
                                                 external_force_densities, cache,
                                                 a0, a2, a3, a6, a7, nnodes, dof)

    # Post
    deformed_coorNP1 .= coor .+ uNP1
    @inbounds for iID in nodes, d in 1:dof
        forces[iID, d] = force_densities_NP1[iID, d] * volume[iID]
    end

    active_nodes = Data_Manager.get_field("Active Nodes")
    active_list = Data_Manager.get_field("Active")
    active_nodes = find_active_nodes(active_list, active_nodes, nodes)
    compute_parabolic_problems_after_model_evaluation(active_nodes, run, dt)
    ws.K = K
    ws.active_nodes = active_nodes
    return nothing
end

# ────────────────────────────────────────────────────────────
# Core: one Newmark step
# ────────────────────────────────────────────────────────────
function newmark_step!(K::AbstractMatrix{Float64},
                       M::Vector{Float64},
                       non_BCs::Vector{Int},
                       uN::Matrix{Float64},
                       uNP1::Matrix{Float64},
                       vN::Matrix{Float64},
                       vNP1::Matrix{Float64},
                       aN::Matrix{Float64},
                       aNP1::Matrix{Float64},
                       F_int::Matrix{Float64},
                       F_ext::Matrix{Float64},
                       cache::NewmarkCache,
                       a0, a2, a3, a6, a7,
                       nnodes::Int, dof::Int)
    @timeit "initialize Newmark step" begin
        n_total = nnodes * dof
        n_free = length(non_BCs)
        isempty(non_BCs) && return nothing
        u_vec = vec(uNP1)
        F_int_vec = vec(F_int)
        F_ext_vec = vec(F_ext)
    end

    # Detect BC change ONCE at the start
    bcs_changed = cache.last_non_BCs != non_BCs

    # ── Effective force for free DOFs ──
    @timeit "compute F_eff" begin
        resize!(cache.F_eff, n_free)
        @inbounds for (idx, i) in enumerate(non_BCs)
            m_rhs = M[i] * (a0 * uN[i] + a2 * vN[i] + a3 * aN[i])
            cache.F_eff[idx] = -F_ext_vec[i] - F_int_vec[i] + m_rhs
        end
    end

    @timeit "bc" begin
        # ── Prescribed DOF contributions ──
        has_BCs = n_free < n_total
        if has_BCs
            # bc_mask & bc_dofs nur bei BC-Wechsel neu berechnen
            if bcs_changed
                resize!(cache.bc_mask, n_total)
                fill!(cache.bc_mask, true)
                @inbounds for i in non_BCs
                    cache.bc_mask[i] = false
                end
                cache.bc_dofs = findall(cache.bc_mask)   # ← FIX 1: in cache!
            end
            resize!(cache.temp, n_free)
            mul!(cache.temp, K[non_BCs, cache.bc_mask], @view(u_vec[cache.bc_mask]))
            @inbounds for (idx, i) in enumerate(non_BCs)
                cache.F_eff[idx] += cache.temp[idx]
            end
        else
            empty!(cache.bc_dofs)
        end
    end

    @timeit "factorisation" begin
        # ── Factorise K_eff = K_free + a0·diag(M_free) ──
        if cache.K_eff_lu === nothing || bcs_changed
            @timeit "K_eff factorisation" begin
                K_free = -sparse(K[non_BCs, non_BCs])
                @inbounds for (idx, i) in enumerate(non_BCs)
                    K_free[idx, idx] += a0 * M[i]
                end
            end
            @timeit "lu" cache.K_eff_lu=lu(K_free)
        end
    end

    @timeit "rest of newmark step" begin
        # ── Solve ──
        resize!(cache.u_free, n_free)
        @timeit "ldiv" ldiv!(cache.u_free, cache.K_eff_lu, cache.F_eff)

        @timeit "update uNP1" begin
            @inbounds for (idx, i) in enumerate(non_BCs)
                uNP1[i] = cache.u_free[idx]
            end
        end

        if !isempty(cache.bc_dofs)
            @timeit "update forces" begin
                if cache.K_bc === nothing || bcs_changed
                    cache.K_bc = K[cache.bc_dofs, :]
                    resize!(cache.F_bc, length(cache.bc_dofs))
                end
                # Eine sparse Mat-Vec statt n_bc Zeilen-Loops
                mul!(cache.F_bc, cache.K_bc, u_vec)
                @inbounds for (idx, j) in enumerate(cache.bc_dofs)
                    F_int_vec[j] = cache.F_bc[idx] - F_ext_vec[j]
                end
            end
        end

        # ── Update acceleration and velocity ──
        @timeit "update accel and vel" begin
            @inbounds for i in 1:n_total
                aNP1[i] = a0 * (uNP1[i] - uN[i]) - a2 * vN[i] - a3 * aN[i]
                vNP1[i] = vN[i] + a6 * aN[i] + a7 * aNP1[i]
            end
        end
    end

    # ── FIX 2: last_non_BCs erst am ENDE updaten ──
    if bcs_changed
        cache.last_non_BCs = copy(non_BCs)
    end

    return nothing
end

function compute_matrix(nodes::AbstractVector{Int64})
    isempty(nodes) && return Data_Manager.get_stiffness_matrix()
    compute_model(nodes)
end

end # module Newmark
