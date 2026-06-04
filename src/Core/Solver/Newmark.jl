# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Newmark
using Printf
using SparseArrays
using LinearAlgebra
using TimerOutputs: @timeit
using ProgressBars: set_multiline_postfix, set_postfix

using ...Data_Manager
using ...Helpers: check_inf_or_nan, find_active_nodes, progress_bar
using ..Boundary_Conditions: apply_bc_dirichlet, apply_bc_neumann, find_bc_free_dof
using ...Parameter_Handling: get_initial_time,
                             get_fixed_dt,
                             get_final_time,
                             get_numerical_damping,
                             get_safety_factor,
                             get_max_damage,
                             get_nsteps
using ..Model_Factory: compute_stiff_matrix_compatible_models,
                       compute_matrix_based_bond_forces

include("../../Models/Material/Material_Models/Correspondence/Correspondence_matrix_based.jl")
using .Correspondence_matrix_based
using ..Model_Factory.Pre_Calculation.Bond_Deformation

export init_solver
export run_solver

function solver_name()
    return "Newmark"
end

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

    function NewmarkCache(n::Int)
        new(Float64[], BitVector(undef, n), Float64[], Float64[],
            nothing, Int[])
    end
end

# ────────────────────────────────────────────────────────────
# Init
# ────────────────────────────────────────────────────────────
"""
    init_solver(solver_options, params, bcs, block_nodes)

Initialise Newmark-beta solver.

Default: average acceleration (beta=0.5,alpha=0.25, unconditionally stable) taken from Bathe page 512.
With numerical damping alpha: beta >= 0.5 ,alpha >= 0.25*(0.5+beta)^2
"""
function init_solver(solver_options::Dict{Any,Any},
                     params::Dict,
                     bcs::Dict{Any,Any},
                     block_nodes::Dict{Int64,Vector{Int64}})
    find_bc_free_dof(bcs)

    Data_Manager.create_constant_node_vector_field("Delta Displacements", Float64,
                                                   Data_Manager.get_dof())

    solver_options["Initial Time"] = get_initial_time(params)
    solver_options["Final Time"] = get_final_time(params)
    solver_options["Number of Steps"] = get_nsteps(params)
    solver_options["dt"] = (solver_options["Final Time"] - solver_options["Initial Time"]) /
                           solver_options["Number of Steps"]

    solver_options["Newmark Beta"] = get(params["Newmark"], "Newmark Delta", 0.5) # name is delta because of beta. Internally it is beta.
    solver_options["Newmark Alpha"] = get(params["Newmark"], "Newmark Alpha",
                                          0.25 * (0.5 + solver_options["Newmark Beta"])^2)
    solver_options["Matrix Update"] = get(params["Newmark"], "Matrix Update", false)

    for (block, nodes) in pairs(block_nodes)
        model_param = Data_Manager.get_properties(block, "Material Model")
        Correspondence_matrix_based.init_model(nodes, model_param, block)
    end
    @timeit "init matrix" Correspondence_matrix_based.init_matrix()

    Data_Manager.create_node_scalar_field("Damage", Float64)
    coor = Data_Manager.get_field("Coordinates")
    Data_Manager.get_field("Deformed Coordinates", "N") .= coor
    Data_Manager.get_field("Deformed Coordinates", "NP1") .= coor
end

# ────────────────────────────────────────────────────────────
# Lumped mass vector  M_i = ρ_i V_i
# ────────────────────────────────────────────────────────────
function build_lumped_mass(nnodes::Int, dof::Int)
    density = Data_Manager.get_field("Density")
    volume = Data_Manager.get_field("Volume")
    M = zeros(Float64, nnodes * dof)
    @inbounds for d in 1:dof
        off = (d - 1) * nnodes
        for i in 1:nnodes
            M[off + i] = density[i] * volume[i]
        end
    end
    return M
end

# ────────────────────────────────────────────────────────────
# Run
# ────────────────────────────────────────────────────────────
"""
    run_solver(solver_options, block_nodes, bcs, outputs, result_files,
               synchronise_field, write_results,
               compute_parabolic_problems_before_model_evaluation,
               compute_parabolic_problems_after_model_evaluation, silent)

Run the Newmark-beta implicit time integration for structural dynamics.

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

# Parameters (set in `solver_options`)

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
- `solver_options::Dict{Any,Any}`:  solver configuration (dt, nsteps, beta,gamma, …)
- `block_nodes::Dict{Int64,Vector{Int64}}`:  block → node mapping
- `bcs::Dict{Any,Any}`:  boundary conditions (Dirichlet / Neumann)
- `outputs::Dict{Int64,Dict{}}`:  output configuration per block
- `result_files::Vector{Dict}`:  file handles for result output
- `synchronise_field::Function`:  MPI field synchronisation callback
- `write_results::Function`:  result writer callback
- `compute_parabolic_problems_before_model_evaluation::Function`:  thermal pre-step
- `compute_parabolic_problems_after_model_evaluation::Function`:  thermal post-step
- `silent::Bool`:  suppress progress bar

# Returns
- `result_files::Vector{Dict}`:  updated result file handles

# References
- Newmark, N. M. (1959). "A method of computation for structural dynamics."
  *J. Eng. Mech. Div.*, 85(3), 67–94.
- Hughes, T. J. R. (2000). *The Finite Element Method*, §9.1–9.3.
"""

function run_solver(solver_options::Dict{Any,Any},
                    block_nodes::Dict{Int64,Vector{Int64}},
                    bcs::Dict{Any,Any},
                    outputs::Dict{Int64,Dict{}},
                    result_files::Vector{Dict},
                    synchronise_field::Function,
                    write_results::Function,
                    compute_parabolic_problems_before_model_evaluation::Function,
                    compute_parabolic_problems_after_model_evaluation::Function,
                    silent::Bool)
    atexit(() -> Data_Manager.set_cancel(true))
    @info "Run Newmark-beta Solver"

    volume = Data_Manager.get_field("Volume")
    coor = Data_Manager.get_field("Coordinates")
    external_forces = Data_Manager.get_field("External Forces")
    external_force_densities = Data_Manager.get_field("External Force Densities")

    dt::Float64 = solver_options["dt"]
    nsteps::Int64 = solver_options["Number of Steps"]
    beta::Float64 = solver_options["Newmark Beta"]
    alpha::Float64 = solver_options["Newmark Alpha"]
    matrix_update = solver_options["Matrix Update"]
    time::Float64 = solver_options["Initial Time"]
    final_time::Float64 = solver_options["Final Time"]
    step_time = 0.0

    rank = Data_Manager.get_rank()
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

    @info "Newmark: alpha=$alpha, gamma=$beta, dt=$dt"

    M = build_lumped_mass(nnodes, dof)

    active_nodes = Data_Manager.get_field("Active Nodes")
    active_list = Data_Manager.get_field("Active")
    active_nodes = find_active_nodes(active_list, active_nodes, nodes)

    K = Data_Manager.get_stiffness_matrix()
    non_BCs = Data_Manager.get_bc_free_dof()

    cache = NewmarkCache(nnodes * dof)
    iter = progress_bar(rank, nsteps, silent)

    for idt in iter
        Data_Manager.set_iteration(idt)
        @timeit "Newmark step" begin
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

            compute_parabolic_problems_before_model_evaluation(active_nodes, solver_options)
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
                                                                            solver_options["Models"],
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
            @timeit "newmark_solve" @views newmark_step!(K, M, non_BCs,
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
            compute_parabolic_problems_after_model_evaluation(active_nodes, solver_options,
                                                              dt)

            @timeit "write_results" result_files=write_results(result_files, time, 0.0,
                                                               outputs)

            if rank == 0 && !silent && Data_Manager.get_cancel()
                set_multiline_postfix(iter, "Simulation canceled!")
                break
            end

            @timeit "switch_NP1_to_N" Data_Manager.switch_NP1_to_N()

            time += dt
            step_time += dt
            Data_Manager.set_current_time(time)

            if idt % max(1, ceil(Int, nsteps / 100)) == 0
                @info "Step: $idt / $nsteps [$(@sprintf("%.3e", time)) s]"
            end
            rank == 0 && !silent && set_postfix(iter, t = @sprintf("%.3e", time))
        end
    end
    Data_Manager.set_current_time(final_time)
    return result_files
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
    n_total = nnodes * dof
    n_free = length(non_BCs)
    isempty(non_BCs) && return nothing

    u_n_vec = vec(uN)
    u_vec = vec(uNP1)
    F_int_vec = vec(F_int)
    F_ext_vec = vec(F_ext)

    # ── Effective force for free DOFs ──
    resize!(cache.F_eff, n_free)
    @inbounds for (idx, i) in enumerate(non_BCs)
        m_rhs = M[i] * (a0 * uN[i] + a2 * vN[i] + a3 * aN[i])
        cache.F_eff[idx] = F_ext_vec[i] - F_int_vec[i] + m_rhs
    end

    # ── Prescribed DOF contributions ──
    has_BCs = n_free < n_total
    if has_BCs
        resize!(cache.bc_mask, n_total)
        fill!(cache.bc_mask, true)
        @inbounds for i in non_BCs
            cache.bc_mask[i] = false
        end

        resize!(cache.temp, n_free)
        mul!(cache.temp, K[non_BCs, cache.bc_mask], @view(u_vec[cache.bc_mask]))

        # Mass contribution from prescribed DOFs
        bc_dofs = findall(cache.bc_mask)
        @inbounds for (idx, i) in enumerate(non_BCs)
            # K part
            cache.F_eff[idx] -= cache.temp[idx]
            # M part: a0 * M_i * u_bc already in u_vec for prescribed DOFs
            # (prescribed u_{n+1} is set by Dirichlet BC)
        end
    end

    # ── Factorise K_eff = K_free + a0·diag(M_free) ──
    if cache.K_eff_lu === nothing || cache.last_non_BCs != non_BCs
        @timeit "K_eff factorisation" begin
            K_free = sparse(K[non_BCs, non_BCs])
            @inbounds for (idx, i) in enumerate(non_BCs)
                K_free[idx, idx] += a0 * M[i]
            end
            cache.K_eff_lu = lu(K_free)
            cache.last_non_BCs = copy(non_BCs)
        end
    end
    if !isempty(bc_dofs)
        @inbounds for j in bc_dofs
            f_int = 0.0
            for i in 1:n_total
                f_int += K[j, i] * uNP1[i]
            end

            F_int_vec[j] = f_int - F_ext_vec[j]
        end
    end

    # ── Solve ──
    resize!(cache.u_free, n_free)
    @timeit "ldiv" ldiv!(cache.u_free, cache.K_eff_lu, cache.F_eff)

    @inbounds for (idx, i) in enumerate(non_BCs)
        uNP1[i] = cache.u_free[idx]
    end

    # ── Update acceleration and velocity ──
    @inbounds for i in 1:n_total
        aNP1[i] = a0 * (uNP1[i] - uN[i]) - a2 * vN[i] - a3 * aN[i]
        vNP1[i] = vN[i] + a6 * aN[i] - a7 * aNP1[i]
    end

    return nothing
end

function compute_matrix(nodes::AbstractVector{Int64})
    isempty(nodes) && return Data_Manager.get_stiffness_matrix()
    Correspondence_matrix_based.compute_model(nodes)
end

end # module Newmark
