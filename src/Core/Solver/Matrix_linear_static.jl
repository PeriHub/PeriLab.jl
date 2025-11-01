# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Linear_static_matrix_based
using ProgressBars: set_multiline_postfix, set_postfix
using Printf
using SparseArrays
using LinearAlgebra
using TimerOutputs

using ...Helpers: check_inf_or_nan, find_active_nodes, progress_bar
using ...Parameter_Handling: get_initial_time, get_nsteps, get_final_time
using ...MPI_Communication: barrier
using ..Boundary_Conditions: apply_bc_dirichlet, apply_bc_neumann, find_bc_free_dof

using ...Helpers: check_inf_or_nan, find_active_nodes, progress_bar, matrix_style
using ...Parameter_Handling:
                             get_initial_time,
                             get_fixed_dt,
                             get_final_time,
                             get_numerical_damping,
                             get_safety_factor,
                             get_max_damage

using ..Model_Factory: compute_stiff_matrix_compatible_models

include("../../Models/Material/Material_Models/Correspondence/Correspondence_matrix_based.jl")
using .Correspondence_matrix_based
include("../../Models/Pre_calculation/bond_deformation.jl")
using .Bond_Deformation

"""
	compute_thermodynamic_critical_time_step(nodes::AbstractVector{Int64}, datamanager::Module, lambda::Float64, Cv::Float64)

Calculate the critical time step for a thermodynamic simulation based on  [OterkusS2014](@cite).

This function iterates over a collection of nodes and computes the critical time step for each node using provided input data and parameters.

# Arguments
- `nodes::AbstractVector{Int64}`: The collection of nodes to calculate the critical time step for.
- `datamanager::Module`: The data manager module that provides access to required data fields.
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
                                                  datamanager::Module,
                                                  lambda::Union{Float64,Int64})
    critical_time_step::Float64 = 1.0e50
    nlist = datamanager.get_nlist()
    density = datamanager.get_field("Density")
    undeformed_bond_length = datamanager.get_field("Bond Length")
    volume = datamanager.get_field("Volume")
    Cv = datamanager.get_field("Specific Heat Capacity")
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
	init_solver(params::Dict, bcs::Dict{Any,Any}, datamanager::Module, block_nodes::Dict{Int64,Vector{Int64}}, mechanical::Bool, thermo::Bool)

Initialize the Verlet solver for a simulation.

This function sets up the Verlet solver for a simulation by initializing various parameters and calculating the time step based on provided parameters or critical time step calculations.

# Arguments
- `params::Dict`: A dictionary containing simulation parameters.
- `bcs::Dict{Any,Any}`: Boundary conditions
- `datamanager::Module`: The data manager module that provides access to required data fields and properties.
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
- `compute_crititical_time_step`: Used to calculate the critical time step if `dt` is not fixed.
- `get_integration_steps`: Used to determine the number of integration steps and adjust the time step.
- `find_and_set_core_value_min` and `find_and_set_core_value_max`: Used to set core values in a distributed computing environment.
"""
function init_solver(solver_options::Dict{Any,Any},
                     params::Dict,
                     bcs::Dict{Any,Any},
                     datamanager::Module,
                     block_nodes::Dict{Int64,Vector{Int64}},
                     to::TimerOutput)
    find_bc_free_dof(datamanager, bcs)
    delta_u = datamanager.create_constant_node_field("Delta Displacements", Float64,
                                                     datamanager.get_dof())
    solver_options["Initial Time"] = get_initial_time(params, datamanager)
    solver_options["Final Time"] = get_final_time(params, datamanager)

    solver_options["Number of Steps"] = get_nsteps(params)
    ## Remark: For static analysis, the number of steps is mandatory
    # dt cannot be defined here; Number of steps define the value,
    # because for static analysis it is a virtual case
    # for thermal analysis this must be adopted and dt and nsteps must be computed
    solver_options["dt"] = (solver_options["Final Time"] - solver_options["Initial Time"]) /
                           solver_options["Number of Steps"]

    for (block, nodes) in pairs(block_nodes)
        model_param = datamanager.get_properties(block, "Material Model")
        Correspondence_matrix_based.init_model(datamanager, nodes, model_param)
    end
    @timeit to "init_matrix" Correspondence_matrix_based.init_matrix(datamanager)
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

function run_solver(solver_options::Dict{Any,Any},
                    block_nodes::Dict{Int64,Vector{Int64}},
                    bcs::Dict{Any,Any},
                    datamanager::Module,
                    outputs::Dict{Int64,Dict{}},
                    result_files::Vector{Dict},
                    synchronise_field,
                    write_results,
                    compute_parabolic_problems_before_model_evaluation,
                    compute_parabolic_problems_after_model_evaluation,
                    to::TimerOutputs.TimerOutput,
                    silent::Bool)
    atexit(() -> datamanager.set_cancel(true))

    @info "Run Linear Static Solver"
    max_damage::Float64 = 0
    volume = datamanager.get_field("Volume")
    coor = datamanager.get_field("Coordinates")
    bond_force = datamanager.get_field("Bond Forces")
    comm = datamanager.get_comm()
    # needed to bring K_sparse into the consistent style of the transformation matrix -> vector
    perm = create_permutation(datamanager.get_nnodes(), datamanager.get_dof())
    if "Material" in solver_options["Models"]
        external_forces = datamanager.get_field("External Forces")
        external_force_densities = datamanager.get_field("External Force Densities")
        a = datamanager.get_field("Acceleration")
    end

    dt::Float64 = solver_options["dt"]
    nsteps::Int64 = solver_options["Number of Steps"]
    time::Float64 = solver_options["Initial Time"]
    step_time::Float64 = 0
    rank = datamanager.get_rank()
    iter = progress_bar(rank, nsteps, silent)
    nodes::Vector{Int64} = collect(1:datamanager.get_nnodes())
    nlist = datamanager.get_nlist()
    dof = datamanager.get_dof()

    uN = datamanager.get_field("Displacements", "N")
    uNP1 = datamanager.get_field("Displacements", "NP1")
    delta_u = datamanager.get_field("Delta Displacements")
    deformed_coorNP1 = datamanager.get_field("Deformed Coordinates", "NP1")
    forces = datamanager.get_field("Forces", "NP1")
    force_densities_N = datamanager.get_field("Force Densities", "N")
    force_densities_NP1 = datamanager.get_field("Force Densities", "NP1")
    volume = datamanager.get_field("Volume")
    forces = datamanager.get_field("Forces", "NP1")
    velocities = datamanager.get_field("Velocity", "NP1")
    external_force_densities = datamanager.get_field("External Force Densities")

    @inbounds @fastmath for idt in iter
        datamanager.set_iteration(idt)
        @timeit to "Linear Static" begin
            datamanager = apply_bc_dirichlet(["Displacements", "Forces", "Force Densities"],
                                             bcs,
                                             datamanager, time,
                                             step_time)
            external_force_densities .= external_forces ./ volume # it must be delta external forces
            # reshape
            non_BCs = datamanager.get_bc_free_dof()
            delta_u .= uNP1-uN
            K = datamanager.get_stiffness_matrix()

            @timeit to "update stiffness matrix" K = compute_matrix(datamanager)
            @timeit to "compute_displacements" compute_displacements!(K[perm, perm],
                                                                      non_BCs,
                                                                      delta_u,
                                                                      force_densities_NP1,
                                                                      external_force_densities)
            @timeit to "compute bond forces" Correspondence_matrix_based.compute_bond_force(bond_force,
                                                                                            K,
                                                                                            uNP1,
                                                                                            nodes,
                                                                                            nlist,
                                                                                            dof)
            @timeit to "compute_models" datamanager = compute_stiff_matrix_compatible_models(datamanager,
                                                                                             block_nodes,
                                                                                             dt,
                                                                                             time,
                                                                                             solver_options["Models"],
                                                                                             synchronise_field,
                                                                                             to)
            @timeit to "update stiffness matrix" K = compute_matrix(datamanager)
            @timeit to "compute_displacements" compute_displacements!(K[perm, perm],
                                                                      non_BCs,
                                                                      delta_u,
                                                                      force_densities_NP1,
                                                                      external_force_densities)
            force_densities_NP1 .+= force_densities_N
            for iID in nodes
                @views forces[iID, :] .= force_densities_NP1[iID, :] .* volume[iID]
            end
            uNP1.=delta_u .+ uN
            velocities .= delta_u ./ dt
            @views deformed_coorNP1 .= coor .+ uNP1

            # @timeit to "download_from_cores" datamanager.synch_manager(synchronise_field,
            #                                                            "download_from_cores")

            @timeit to "write_results" result_files=write_results(result_files, time,
                                                                  max_damage, outputs,
                                                                  datamanager)
            # for file in result_files
            #     flush(file)
            # end
            if rank == 0 && !silent && datamanager.get_cancel()
                set_multiline_postfix(iter, "Simulation canceled!")
                break
            end

            time += dt
            step_time += dt
            datamanager.set_current_time(time)

            if idt % ceil(nsteps / 100) == 0
                @info "Step: $idt / $(nsteps+1) [$time s]"
            end
            if rank == 0 && !silent
                set_postfix(iter, t = @sprintf("%.4e", time))
            end

            barrier(comm)
        end
    end
    return result_files
end

function compute_matrix(datamanager::Module)
    nodes = collect(1:datamanager.get_nnodes())
    Bond_Deformation.compute(datamanager,
                             nodes,
                             Dict(),
                             0) # not needed here
    return Correspondence_matrix_based.compute_model(datamanager, nodes)
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
function compute_displacements!(K::SparseMatrixCSC{Float64,Int64},
                                non_BCs::Vector{Int64},
                                u::Matrix{Float64},
                                F_int::Matrix{Float64},
                                F_ext::Matrix{Float64})

    # Compute modified force: F_modified = F - K * u_prescribed
    # (u contains prescribed displacements at fixed DOFs)

    #mul!(vec(F_int), K, vec(u))  # F_temp = K * u (in-place, no allocation)
    #F_int .+= F_ext   # F_temp = F - K*u (in-place)
    ## Update displacement vector at free DOFs
    #u[non_BCs] = K[non_BCs, non_BCs] \ F_int[non_BCs]

    F_modified = copy(vec(F_ext))
    BCs = setdiff(1:length(vec(u)), non_BCs)
    F_int[non_BCs] = (K[non_BCs, BCs] * vec(u)[BCs])
    @views F_modified[non_BCs] .-= F_int[non_BCs]
    @views vec(u)[non_BCs] .= K[non_BCs, non_BCs] \ F_modified[non_BCs]

    #return reshape(F_int, 2, :)

end

function create_permutation(nnodes::Int, dof::Int)
    perm = Vector{Int}(undef, nnodes * dof)
    idx = 1
    for d in 1:dof           # Für jeden DOF
        for n in 1:nnodes     # Für jeden Knoten
            old_idx = (n-1)*dof + d  # Row-major: node, dann dof
            perm[idx] = old_idx
            idx += 1
        end
    end
    return perm
end
end
