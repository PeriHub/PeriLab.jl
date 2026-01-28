# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Linear_static_matrix_based
using ProgressBars: set_multiline_postfix, set_postfix
using Printf
using SparseArrays
using LinearAlgebra
using TimerOutputs: @timeit

using ...Data_Manager
using ...Helpers: check_inf_or_nan, find_active_nodes, progress_bar, matrix_style,
                  create_permutation
using ...MPI_Communication: barrier
using ..Boundary_Conditions: apply_bc_dirichlet, apply_bc_neumann, find_bc_free_dof
using ...Parameter_Handling:
                             get_initial_time,
                             get_fixed_dt,
                             get_final_time,
                             get_numerical_damping,
                             get_safety_factor,
                             get_max_damage,
                             get_nsteps

using ..Model_Factory: compute_stiff_matrix_compatible_models

include("../../Models/Material/Material_Models/Correspondence/Correspondence_matrix_based.jl")
using .Correspondence_matrix_based
using ..Model_Factory.Pre_Calculation.Bond_Deformation

export init_solver
export run_solver

mutable struct DisplacementSolverCache
    # Workspace
    F_modified::Vector{Float64}
    bc_mask::BitVector
    temp::Vector{Float64}
    u_free::Vector{Float64}

    # LU cache
    K_free_lu::Union{Nothing,Any}
    last_non_BCs::Vector{Int}

    # Active DOFs cache
    K_active_cached::Union{Nothing,AbstractMatrix{Float64}}
    last_active_dofs::Vector{Int}

    function DisplacementSolverCache(n_total::Int)
        new(Float64[], BitVector(undef, n_total), Float64[], Float64[],
            nothing, Int[], nothing, Int[])
    end
end

"""
	compute_thermodynamic_critical_time_step(nodes::AbstractVector{Int64}, lambda::Float64, Cv::Float64)

Calculate the critical time step for a thermodynamic simulation based on  [OterkusS2014](@cite).

This function iterates over a collection of nodes and computes the critical time step for each node using provided input data and parameters.

# Arguments
- `nodes::AbstractVector{Int64}`: The collection of nodes to calculate the critical time step for.
- `lambda::Float64`: The material parameter used in the calculations.
- `Cv::Float64`: The heat capacity at constant volume used in the calculations.

# Returns
- `Float64`: The calculated critical time step for the thermodynamic simulation.

# Dependencies
This function depends on the following data fields from the `datamanager` module:
- `get_nlist()`: Returns the neighbor list.
- `get_field("Density")`: Returns the density field.
- `get_field("Bond Length")`: Returns the bond distance field.
- `get_field("Volume")`: Returns the volume field.
- `get_field("Number of Neighbors")`: Returns the number of neighbors field.
"""
function compute_thermodynamic_critical_time_step(nodes::AbstractVector{Int64},
                                                  lambda::Union{Float64,Int64})
    critical_time_step::Float64 = 1.0e50
    nlist = Data_Manager.get_nlist()
    density = Data_Manager.get_field("Density")
    undeformed_bond_length = Data_Manager.get_field("Bond Length")
    volume = Data_Manager.get_field("Volume")
    Cv = Data_Manager.get_field("Specific Heat Capacity")
    lambda = matrix_style(lambda)
    eigLam = maximum(eigvals(lambda))

    for iID in nodes
        denominator = get_cs_denominator(volume[nlist[iID]], undeformed_bond_length[iID])
        t = density[iID] * Cv[iID] / (eigLam * denominator)
        critical_time_step = test_timestep(t, critical_time_step)
    end
    return sqrt(critical_time_step)
end

"""
	init_solver(params::Dict, bcs::Dict{Any,Any}, block_nodes::Dict{Int64,Vector{Int64}}, mechanical::Bool, thermo::Bool)

Initialize the Verlet solver for a simulation.

This function sets up the Verlet solver for a simulation by initializing various parameters and calculating the time step based on provided parameters or critical time step calculations.

# Arguments
- `params::Dict`: A dictionary containing simulation parameters.
- `bcs::Dict{Any,Any}`: Boundary conditions
- `block_nodes::Dict{Int64,Vector{Int64}}`: A dictionary mapping block IDs to collections of nodes.
- `mechanical::Bool`: If `true`, mechanical properties are considered in the calculation.
- `thermo::Bool`: If `true`, thermodynamic properties are considered in the calculation.

# Returns
A tuple `(initial_time, dt, nsteps, numerical_damping)` where:
- `initial_time::Float64`: The initial time for the simulation.
- `dt::Float64`: The time step for the simulation.
- `nsteps::Int64`: The number of time integration steps.
- `numerical_damping::Float64`: The numerical damping factor.
- `max_damage::Float64`: The maximum damage in the simulation.

# Dependencies
This function may depend on the following functions:
- `get_initial_time`, `get_final_time`, `get_safety_factor`, `get_fixed_dt`: Used to retrieve simulation parameters.
- `get_integration_steps`: Used to determine the number of integration steps and adjust the time step.
- `find_and_set_core_value_min` and `find_and_set_core_value_max`: Used to set core values in a distributed computing environment.
"""
function init_solver(solver_options::Dict{Any,Any},
                     params::Dict,
                     bcs::Dict{Any,Any},
                     block_nodes::Dict{Int64,Vector{Int64}})
    find_bc_free_dof(bcs)
    delta_u = Data_Manager.create_constant_node_vector_field("Delta Displacements", Float64,
                                                             Data_Manager.get_dof())
    solver_options["Initial Time"] = get_initial_time(params)
    solver_options["Final Time"] = get_final_time(params)

    solver_options["Number of Steps"] = get_nsteps(params)
    ## Remark: For static analysis, the number of steps is mandatory
    # dt cannot be defined here; Number of steps define the value,
    # because for static analysis it is a virtual case
    # for thermal analysis this must be adopted and dt and nsteps must be computed
    solver_options["dt"] = (solver_options["Final Time"] - solver_options["Initial Time"]) /
                           solver_options["Number of Steps"]

    for (block, nodes) in pairs(block_nodes)
        model_param = Data_Manager.get_properties(block, "Material Model")
        Correspondence_matrix_based.init_model(nodes, model_param, block)
    end

    @timeit "init matrix" Correspondence_matrix_based.init_matrix()

    solver_options["Matrix update"] = get(params["Linear Static Matrix Based"],
                                          "Matrix update", false)
    deformed_coorN = Data_Manager.get_field("Deformed Coordinates", "N")
    deformed_coorNP1 = Data_Manager.get_field("Deformed Coordinates", "NP1")
    coor = Data_Manager.get_field("Coordinates")
    deformed_coorN .= coor
    deformed_coorNP1 .= coor

    ### for coupled thermal analysis

    #critical_time_step::Float64 = 1.0e50
    #for iblock in eachindex(block_nodes)
    #	if thermal
    #		lambda = datamanager.get_property(iblock, "Thermal Model",
    #			"Thermal Conductivity")
    #		# if Cv and lambda are not defined it is valid, because an analysis can take place, if material is still analysed
    #		if isnothing(lambda)
    #			if !mechanical
    #				@error "No time step can be calculated, because the heat conduction is not defined."
    #				return nothing
    #			end
    #		else
    #			t = compute_thermodynamic_critical_time_step(block_nodes[iblock],
    #				datamanager,
    #				lambda)
    #			critical_time_step = test_timestep(t, critical_time_step)
    #		end
    #	end
    #end
end

function filter_dofs_to_active_nodes(non_BCs::AbstractVector{Int64},
                                     active_nodes::AbstractVector{Int64},
                                     dof::Int64)
    active_set = Set(active_nodes)

    filtered = Int64[]
    for dof_idx in non_BCs
        node = div(dof_idx - 1, dof) + 1
        if node in active_set
            push!(filtered, dof_idx)
        end
    end

    return filtered
end
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

    @info "Run Linear Static Solver"
    max_damage::Float64 = 0
    volume = Data_Manager.get_field("Volume")
    coor = Data_Manager.get_field("Coordinates")
    bond_force = Data_Manager.get_field("Bond Forces")
    comm = Data_Manager.get_comm()

    if "Material" in solver_options["Models"]
        external_forces = Data_Manager.get_field("External Forces")
        external_force_densities = Data_Manager.get_field("External Force Densities")
    end

    dt::Float64 = solver_options["dt"]
    nsteps::Int64 = solver_options["Number of Steps"]
    time::Float64 = solver_options["Initial Time"]
    step_time::Float64 = 0
    rank = Data_Manager.get_rank()
    iter = progress_bar(rank, nsteps, silent)
    nodes::Vector{Int64} = collect(1:Data_Manager.get_nnodes())
    dof = Data_Manager.get_dof()
    nnodes = Data_Manager.get_nnodes()

    delta_u = Data_Manager.get_field("Delta Displacements")

    matrix_update::Bool = solver_options["Matrix update"]

    active_nodes = Data_Manager.get_field("Active Nodes")
    active_list = Data_Manager.get_field("Active")

    active_nodes = find_active_nodes(active_list, active_nodes, nodes)
    K = Data_Manager.get_stiffness_matrix()

    displacement_solver_cache = DisplacementSolverCache(nnodes * dof)

    # Determine if all nodes are active
    all_nodes_active = length(active_nodes) == nnodes

    non_BCs_global = Data_Manager.get_bc_free_dof()

    for idt in iter
        Data_Manager.set_iteration(idt)
        @timeit "Linear Static" begin
            uN::NodeVectorField{Float64} = Data_Manager.get_field("Displacements", "N")
            uNP1::NodeVectorField{Float64} = Data_Manager.get_field("Displacements", "NP1")
            velocities::NodeVectorField{Float64} = Data_Manager.get_field("Velocity", "NP1")
            deformed_coorN::NodeVectorField{Float64} = Data_Manager.get_field("Deformed Coordinates",
                                                                              "N")
            deformed_coorNP1::NodeVectorField{Float64} = Data_Manager.get_field("Deformed Coordinates",
                                                                                "NP1")

            forces::NodeVectorField{Float64} = Data_Manager.get_field("Forces", "NP1")
            force_densities_N::NodeVectorField{Float64} = Data_Manager.get_field("Force Densities",
                                                                                 "N")
            force_densities_NP1::NodeVectorField{Float64} = Data_Manager.get_field("Force Densities",
                                                                                   "NP1")

            deformed_coorNP1 .= coor

            compute_parabolic_problems_before_model_evaluation(active_nodes, solver_options)
            apply_bc_dirichlet([
                                   "Displacements",
                                   "Forces",
                                   "Force Densities",
                                   "Temperature"
                               ],
                               bcs,
                               time,
                               step_time)

            @. external_force_densities = external_forces / volume

            @timeit "compute_models" compute_stiff_matrix_compatible_models(block_nodes,
                                                                            dt,
                                                                            time,
                                                                            solver_options["Models"],
                                                                            synchronise_field)

            active_nodes = Data_Manager.get_field("Active Nodes")
            active_list = Data_Manager.get_field("Active")
            active_nodes = find_active_nodes(active_list, active_nodes, nodes)

            # Test current status
            current_all_active = length(active_nodes) == nnodes

            if !current_all_active
                @timeit "rebuild matrix (partial active)" begin
                    compute_matrix(active_nodes)
                    K = Data_Manager.get_stiffness_matrix()
                end

                # Build global DOF indices for active nodes
                # Using block-style indexing: global_idx = (d - 1) * nnodes + node
                active_dofs_global = Vector{Int}(undef, length(active_nodes) * dof)
                idx = 1
                for d in 1:dof
                    for node in active_nodes
                        active_dofs_global[idx] = (d - 1) * nnodes + node
                        idx += 1
                    end
                end
                sort!(active_dofs_global)

                # Filter non_BCs to active DOFs
                active_non_BCs_global = intersect(non_BCs_global, active_dofs_global)

                @timeit "compute_displacements (partial)" begin
                    if length(active_non_BCs_global) > 0
                        # Create views to active node data in matrix form
                        u_active_nodes = @view uNP1[active_nodes, :]
                        F_int_active_nodes = @view force_densities_NP1[active_nodes, :]
                        F_ext_active_nodes = @view external_force_densities[active_nodes, :]

                        # Build mapping from global to active indices
                        global_to_active = Dict{Int,Int}()
                        sizehint!(global_to_active, length(active_dofs_global))
                        for (active_idx, global_idx) in enumerate(active_dofs_global)
                            global_to_active[global_idx] = active_idx
                        end

                        # Map BC indices to active system
                        active_non_BCs_local = Vector{Int}(undef,
                                                           length(active_non_BCs_global))
                        @inbounds for (i, g) in enumerate(active_non_BCs_global)
                            active_non_BCs_local[i] = global_to_active[g]
                        end

                        # Extract submatrix for active DOFs
                        K_active = K[active_dofs_global, active_dofs_global]

                        # Solve system using same structure as reference
                        # but with mapping for active nodes
                        compute_displacements_active_subset!(K_active,
                                                             active_dofs_global,
                                                             active_non_BCs_local,
                                                             u_active_nodes,
                                                             F_int_active_nodes,
                                                             F_ext_active_nodes,
                                                             volume[active_nodes],
                                                             active_nodes,
                                                             nnodes,
                                                             dof,
                                                             displacement_solver_cache)
                    end

                    # Update forces for active nodes
                    @inbounds for iID in active_nodes
                        for d in 1:dof
                            forces[iID, d] = force_densities_NP1[iID, d] * volume[iID]
                        end
                    end
                end
            else
                # All nodes are active
                if matrix_update
                    @timeit "update matrix (all active)" begin
                        compute_matrix(active_nodes)
                        K = Data_Manager.get_stiffness_matrix()
                    end
                end

                all_dofs = collect(1:(nnodes * dof))

                @timeit "compute_displacements (all active)" begin
                    @views compute_displacements!(K,
                                                  all_dofs,
                                                  non_BCs_global,
                                                  uNP1,
                                                  force_densities_NP1,
                                                  external_force_densities,
                                                  displacement_solver_cache)

                    # Update forces for all nodes
                    @inbounds for iID in nodes
                        for d in 1:dof
                            forces[iID, d] = force_densities_NP1[iID, d] * volume[iID]
                        end
                    end
                end
            end

            all_nodes_active = current_all_active

            active_nodes = Data_Manager.get_field("Active Nodes")
            active_list = Data_Manager.get_field("Active")
            active_nodes = find_active_nodes(active_list, active_nodes, 1:nnodes)

            compute_parabolic_problems_after_model_evaluation(active_nodes, solver_options,
                                                              dt)

            @. velocities = (uNP1 - uN) / dt
            @. @views deformed_coorNP1 = coor + uNP1

            @timeit "write_results" result_files=write_results(result_files, time,
                                                               max_damage, outputs)

            if rank == 0 && !silent && Data_Manager.get_cancel()
                set_multiline_postfix(iter, "Simulation canceled!")
                break
            end

            @timeit "switch_NP1_to_N" Data_Manager.switch_NP1_to_N()

            time += dt
            step_time += dt
            Data_Manager.set_current_time(time)

            if idt % ceil(nsteps / 100) == 0
                @info "Step: $idt / $(nsteps) [$time s]"
            end
            if rank == 0 && !silent
                set_postfix(iter, t = @sprintf("%.4e", time))
            end

            barrier(comm)
        end
    end
    Data_Manager.set_current_time(time - dt)
    return result_files
end

# New function for active subset that mimics reference implementation
function compute_displacements_active_subset!(K_active::AbstractMatrix{Float64},
                                              active_dofs_global::AbstractVector{Int},
                                              active_non_BCs::AbstractVector{Int},
                                              u_active::AbstractMatrix{Float64},
                                              F_int_active::AbstractMatrix{Float64},
                                              F_ext_active::AbstractMatrix{Float64},
                                              volume_active::AbstractVector{Float64},
                                              active_nodes::AbstractVector{Int},
                                              nnodes::Int,
                                              dof::Int,
                                              solver::DisplacementSolverCache)
    isempty(active_non_BCs) && return nothing

    n_active_nodes = length(active_nodes)
    n_active_dofs = n_active_nodes * dof

    # Vectorize matrices (same as reference)
    u_vec = vec(u_active)
    F_int_vec = vec(F_int_active)
    F_ext_vec = vec(F_ext_active)

    n_free = length(active_non_BCs)
    has_BCs = n_free < n_active_dofs

    # Resize buffers
    length(solver.F_modified) != n_free && resize!(solver.F_modified, n_free)
    length(solver.bc_mask) != n_active_dofs && resize!(solver.bc_mask, n_active_dofs)

    # Modify force vector to account for prescribed displacements (same as reference)
    if has_BCs
        fill!(solver.bc_mask, true)
        @inbounds for i in active_non_BCs
            solver.bc_mask[i] = false
        end

        K_sub = K_active[active_non_BCs, solver.bc_mask]
        u_bc = @view u_vec[solver.bc_mask]

        length(solver.temp) != n_free && resize!(solver.temp, n_free)
        mul!(solver.temp, K_sub, u_bc)

        @inbounds for (idx, i) in enumerate(active_non_BCs)
            F_int_vec[i] += solver.temp[idx]
        end
    end

    # Build modified force vector (same as reference)
    @inbounds for (idx, i) in enumerate(active_non_BCs)
        solver.F_modified[idx] = F_ext_vec[i] - F_int_vec[i]
    end

    # LU factorization and solve (same as reference)
    @timeit "LU factorization" begin
        K_free = K_active[active_non_BCs, active_non_BCs]
        try
            solver.K_free_lu = lu(K_free)
        catch e
            @error "LU factorization failed for reduced system"
            @error "Matrix size: $(size(K_free))"
            @error "Rank: $(rank(K_free))"
            rethrow(e)
        end
    end

    length(solver.u_free) != n_free && resize!(solver.u_free, n_free)

    @timeit "ldiv" ldiv!(solver.u_free, solver.K_free_lu, solver.F_modified)

    # Write solution back (same as reference)
    @inbounds for (idx, i) in enumerate(active_non_BCs)
        u_vec[i] = solver.u_free[idx]
    end

    return nothing
end

# Neue Hilfsfunktion für reduziertes System
function compute_displacements_reduced!(K_active::AbstractMatrix{Float64},
                                        active_non_BCS::AbstractVector{Int},
                                        u_active::AbstractVector{Float64},
                                        F_int_active::AbstractVector{Float64},
                                        F_ext_active::AbstractVector{Float64},
                                        solver::DisplacementSolverCache)
    isempty(active_non_BCS) && return nothing

    n_total = length(u_active)
    n_free = length(active_non_BCS)

    length(solver.F_modified) != n_free && resize!(solver.F_modified, n_free)
    length(solver.bc_mask) != n_total && resize!(solver.bc_mask, n_total)

    has_BCs = n_free < n_total

    if has_BCs
        fill!(solver.bc_mask, true)
        @inbounds for i in active_non_BCS
            solver.bc_mask[i] = false
        end

        K_sub = K_active[active_non_BCS, solver.bc_mask]
        u_bc = @view u_active[solver.bc_mask]

        length(solver.temp) != n_free && resize!(solver.temp, n_free)
        mul!(solver.temp, K_sub, u_bc)

        @inbounds for (idx, i) in enumerate(active_non_BCS)
            F_int_active[i] += solver.temp[idx]
        end
    end

    @inbounds for (idx, i) in enumerate(active_non_BCS)
        solver.F_modified[idx] = F_ext_active[i] - F_int_active[i]
    end

    # LU factorization
    @timeit "LU factorization" begin
        K_free = K_active[active_non_BCS, active_non_BCS]
        try
            solver.K_free_lu = lu(K_free)
        catch e
            @error "LU factorization failed: $e"
            @error "Matrix size: $(size(K_free))"
            @error "Condition number: $(cond(K_free))"
            rethrow(e)
        end
    end

    length(solver.u_free) != n_free && resize!(solver.u_free, n_free)

    @timeit "ldiv" ldiv!(solver.u_free, solver.K_free_lu, solver.F_modified)

    @inbounds for (idx, i) in enumerate(active_non_BCS)
        u_active[i] = solver.u_free[idx]
    end

    return nothing
end

function get_active_dof(active_nodes::AbstractVector{Int64}, dof::Int64, nnodes::Int64)
    n_active = length(active_nodes)

    active_dofs_local = collect(1:(n_active * dof))
    local_to_global = Dict{Int,Int}()

    for (iD, global_iD) in enumerate(active_nodes)
        for j in 1:dof
            local_idx = iD + (j - 1) * n_active
            global_idx = global_iD + (j - 1) * nnodes
            local_to_global[local_idx] = global_idx
        end
    end

    return active_dofs_local, local_to_global
end

function filter_and_map_bc(non_BCs_global::AbstractVector{Int64},
                           local_to_global::Dict{Int,Int})
    active_non_BCs_local = Int[]

    for (local_idx, global_idx) in pairs(local_to_global)
        if global_idx in non_BCs_global
            push!(active_non_BCs_local, local_idx)
        end
    end

    return sort(active_non_BCs_local)
end
function compute_matrix(nodes::AbstractVector{Int64})
    if length(nodes) == 0
        return Data_Manager.get_stiffness_matrix()
    end
    Correspondence_matrix_based.compute_model(nodes)
end

"""
	compute_displacements!(K, non_BCs, u, F, F_temp, K_reduced, lu_fact, temp)

Compute displacements with prescribed displacement boundary conditions.
No memory allocations.

Arguments:
- K: Global stiffness matrix
- non_BCs: Free DOF indices
- u: Displacement vector (contains prescribed values at fixed DOFs)
- F: External force vector
- F_temp: Temporary force vector (pre-allocated)

"""
function compute_displacements!(K::AbstractMatrix{Float64},
                                active_dofs::AbstractVector{Int},
                                active_non_BCS::AbstractVector{Int},
                                u_active::AbstractMatrix{Float64},
                                F_int_active::AbstractMatrix{Float64},
                                F_ext_active::AbstractMatrix{Float64},
                                solver::DisplacementSolverCache)
    isempty(active_non_BCS) && return nothing

    # STEP 1: Cache K[active_dofs, active_dofs] if active_dofs changed
    if solver.K_active_cached === nothing || solver.last_active_dofs != active_dofs
        @timeit "extract K_active" begin
            solver.K_active_cached = K
            solver.last_active_dofs = copy(active_dofs)
        end
        solver.K_free_lu = nothing
    end

    K_active = solver.K_active_cached

    # Step 2: Modify force vector to account for prescribed displacements
    u_vec = vec(u_active)
    F_int_vec = vec(F_int_active)
    F_ext_vec = vec(F_ext_active)
    n_total = length(u_vec)
    n_free = length(active_non_BCS)

    length(solver.F_modified) != n_free && resize!(solver.F_modified, n_free)
    length(solver.bc_mask) != n_total && resize!(solver.bc_mask, n_total)

    has_BCs = n_free < n_total

    if has_BCs
        fill!(solver.bc_mask, true)
        @inbounds for i in active_non_BCS
            solver.bc_mask[i] = false
        end

        K_sub = K_active[active_non_BCS, solver.bc_mask]
        u_bc = @view u_vec[solver.bc_mask]

        length(solver.temp) != n_free && resize!(solver.temp, n_free)
        mul!(solver.temp, K_sub, u_bc)

        @inbounds for (idx, i) in enumerate(active_non_BCS)
            F_int_vec[i] += solver.temp[idx]
        end
    end

    @inbounds for (idx, i) in enumerate(active_non_BCS)
        solver.F_modified[idx] = F_ext_vec[i] - F_int_vec[i]
    end

    # Step 3: LU cachen
    if solver.K_free_lu === nothing || solver.last_non_BCs != active_non_BCS
        @timeit "LU factorization" begin
            K_free = K_active[active_non_BCS, active_non_BCS]

            try
                solver.K_free_lu = lu(K_free)
            catch e
                @error "LU factorization failed. The stiffness matrix may be singular, please check the boundary conditions."
            end
            solver.last_non_BCs = copy(active_non_BCS)
        end
    end

    length(solver.u_free) != n_free && resize!(solver.u_free, n_free)

    @timeit "ldiv" ldiv!(solver.u_free, solver.K_free_lu, solver.F_modified)

    @inbounds for (idx, i) in enumerate(active_non_BCS)
        u_vec[i] = solver.u_free[idx]
    end

    return nothing
end

end
