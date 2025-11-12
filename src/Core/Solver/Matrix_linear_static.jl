# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# Notizen
# die matrix soll über eine nodes liste aufgebau werden. Nachbarn die dort nicht existieren, sollen nicht berücksichtigt werden
# die knoten nummern sollen so bleiben. d.h. die Matrixgröße bleibt
# dann gibt es einen active node filter. der muss die randbedingungen berücksichtigen und die reduzierte matrix
# nutzen. die reduzierte matrix greift nur auf die werte zu welche existieren
# Herausforderung
#   I,J -> müssten auf Null einträge zeigen können, da die indizierung vorher bestimmt wird
#   idx muss richtig berechnet werden, auch wenn schleifen geskippt werden
#   bc ist ein nodes x dof Feld auf globaler indizierung; diese muss auf die lokale aktivierte indizierung gemappt werden

module Linear_static_matrix_based
using ProgressBars: set_multiline_postfix, set_postfix
using Printf
using SparseArrays
using LinearAlgebra
using TimerOutputs

using ...Data_Manager
using ...Helpers: check_inf_or_nan, find_active_nodes, progress_bar, matrix_style
using ...MPI_Communication: barrier
using ..Boundary_Conditions: apply_bc_dirichlet, apply_bc_neumann, find_bc_free_dof
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
                     block_nodes::Dict{Int64,Vector{Int64}},
                     to::TimerOutput)
    find_bc_free_dof(Data_Manager, bcs)
    delta_u = Data_Manager.create_constant_node_field("Delta Displacements", Float64,
                                                      Data_Manager.get_dof())
    solver_options["Initial Time"] = get_initial_time(params, Data_Manager)
    solver_options["Final Time"] = get_final_time(params, Data_Manager)

    solver_options["Number of Steps"] = get_nsteps(params)
    ## Remark: For static analysis, the number of steps is mandatory
    # dt cannot be defined here; Number of steps define the value,
    # because for static analysis it is a virtual case
    # for thermal analysis this must be adopted and dt and nsteps must be computed
    solver_options["dt"] = (solver_options["Final Time"] - solver_options["Initial Time"]) /
                           solver_options["Number of Steps"]

    for (block, nodes) in pairs(block_nodes)
        model_param = Data_Manager.get_properties(block, "Material Model")
        Correspondence_matrix_based.init_model(Data_Manager, nodes, model_param)
    end
    @timeit to "init_matrix" Correspondence_matrix_based.init_matrix(Data_Manager)

    solver_options["Matrix update"] = get(params["Linear Static Matrix Based"],
                                          "Matrix update", false)

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
                    synchronise_field,
                    write_results,
                    compute_parabolic_problems_before_model_evaluation,
                    compute_parabolic_problems_after_model_evaluation,
                    to::TimerOutputs.TimerOutput,
                    silent::Bool)
    atexit(() -> Data_Manager.set_cancel(true))

    @info "Run Linear Static Solver"
    max_damage::Float64 = 0
    volume = Data_Manager.get_field("Volume")
    coor = Data_Manager.get_field("Coordinates")
    bond_force = Data_Manager.get_field("Bond Forces")
    comm = Data_Manager.get_comm()
    # needed to bring K_sparse into the consistent style of the transformation matrix -> vector

    if "Material" in solver_options["Models"]
        external_forces = Data_Manager.get_field("External Forces")
        external_force_densities = Data_Manager.get_field("External Force Densities")
        a = Data_Manager.get_field("Acceleration")
    end

    dt::Float64 = solver_options["dt"]
    nsteps::Int64 = solver_options["Number of Steps"]
    time::Float64 = solver_options["Initial Time"]
    step_time::Float64 = 0
    rank = Data_Manager.get_rank()
    iter = progress_bar(rank, nsteps, silent)
    nodes::Vector{Int64} = collect(1:Data_Manager.get_nnodes())
    nlist = Data_Manager.get_nlist()
    dof = Data_Manager.get_dof()

    volume = Data_Manager.get_field("Volume")
    forces = Data_Manager.get_field("Forces", "NP1")

    delta_u = Data_Manager.get_field("Delta Displacements")
    external_force_densities = Data_Manager.get_field("External Force Densities")

    matrix_update::Bool = solver_options["Matrix update"]

    active_nodes = Data_Manager.get_field("Active Nodes")
    active_list = Data_Manager.get_field("Active")

    active_nodes = find_active_nodes(active_list, active_nodes,
                                     nodes)

    if matrix_update
        @timeit to "update stiffness matrix" K = compute_matrix(Data_Manager,
                                                                active_nodes)
    else
        K = Data_Manager.get_stiffness_matrix()
    end
    for idt in iter
        Data_Manager.set_iteration(idt)
        @timeit to "Linear Static" begin

            # reshape
            # All 2 time steps fields must be in the loop, because switch NP1 to N changes the adresses
            uN = Data_Manager.get_field("Displacements", "N")
            uNP1 = Data_Manager.get_field("Displacements", "NP1")
            velocities = Data_Manager.get_field("Velocity", "NP1")
            deformed_coorNP1 = Data_Manager.get_field("Deformed Coordinates", "NP1")
            forces = Data_Manager.get_field("Forces", "NP1")
            force_densities_N = Data_Manager.get_field("Force Densities", "N")
            force_densities_NP1 = Data_Manager.get_field("Force Densities", "NP1")

            non_BCs = Data_Manager.get_bc_free_dof()
            #force_densities_NP1[active_nodes, :] .= force_densities_N[active_nodes, :]

            # das muss in model evaluation; dann funktioniert die Reihenfolge auch. in material
            @timeit to "compute bond forces" Correspondence_matrix_based.compute_bond_force(bond_force,
                                                                                            K,
                                                                                            uN,
                                                                                            active_nodes,
                                                                                            nlist,
                                                                                            dof)

            compute_parabolic_problems_before_model_evaluation(active_nodes, Data_Manager,
                                                               solver_options)
            Data_Manager = apply_bc_dirichlet([
                                                  "Displacements",
                                                  "Forces",
                                                  "Force Densities",
                                                  "Temperature"
                                              ],
                                              bcs,
                                              Data_Manager, time,
                                              step_time)

            external_force_densities .= external_forces ./ volume # it must be delta external forces

            @timeit to "compute_models" Data_Manager = compute_stiff_matrix_compatible_models(Data_Manager,
                                                                                              block_nodes,
                                                                                              dt,
                                                                                              time,
                                                                                              solver_options["Models"],
                                                                                              synchronise_field,
                                                                                              to)

            active_nodes = Data_Manager.get_field("Active Nodes")
            active_list = Data_Manager.get_field("Active")

            active_nodes = find_active_nodes(active_list, active_nodes,
                                             nodes)

            if matrix_update
                @timeit to "update stiffness matrix" K = compute_matrix(Data_Manager,
                                                                        active_nodes)

            else
                K = Data_Manager.get_stiffness_matrix()
            end

            #@views external_force_densities[active_nodes, :] += force_densities_NP1[active_nodes, :]
            perm = create_permutation(active_nodes, Data_Manager.get_dof()) # only active node dofs are there

            filtered_perm = filter_and_map_bc(non_BCs, active_nodes, dof,
                                              Data_Manager.get_nnodes())
            #@info (filtered_perm)
            @timeit to "compute_displacements" @views compute_displacements!(K[perm, perm],
                                                                             filtered_perm,
                                                                             uNP1[active_nodes,
                                                                                  :],
                                                                             force_densities_NP1[active_nodes,
                                                                                                 :],
                                                                             external_force_densities[active_nodes,
                                                                                                      :])

            active_nodes = Data_Manager.get_field("Active Nodes")
            active_list = Data_Manager.get_field("Active")

            active_nodes = find_active_nodes(active_list, active_nodes,
                                             1:Data_Manager.get_nnodes())

            compute_parabolic_problems_after_model_evaluation(active_nodes, Data_Manager,
                                                              solver_options, dt)

            for iID in nodes
                @views forces[iID, :] .= force_densities_NP1[iID, :] .* volume[iID]
            end

            velocities .= (uNP1 - uN) ./ dt
            @views deformed_coorNP1 .= coor .+ uNP1

            # @timeit to "download_from_cores" Data_Manager.synch_manager(synchronise_field,
            #                                                            "download_from_cores")

            @timeit to "write_results" result_files=write_results(result_files, time,
                                                                  max_damage, outputs,
                                                                  Data_Manager)

            # for file in result_files
            #     flush(file)
            # end
            if rank == 0 && !silent && Data_Manager.get_cancel()
                set_multiline_postfix(iter, "Simulation canceled!")
                break
            end

            @timeit to "switch_NP1_to_N" Data_Manager.switch_NP1_to_N()

            time += dt
            step_time += dt
            Data_Manager.set_current_time(time)

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
function filter_and_map_bc(non_BCs, active_nodes, dof, nnodes::Int64)
    filtered = []
    for i in eachindex(active_nodes)
        for j in 1:dof
            if active_nodes[i] + (j - 1)*nnodes in non_BCs
                push!(filtered, i+(j - 1) * length(active_nodes))
            end
        end
    end

    return filtered
end

function compute_matrix(nodes::AbstractVector{Int64})
    if length(nodes)==0
        return Data_Manager.get_stiffness_matrix()
    end
    Bond_Deformation.compute(Data_Manager,
                             nodes,
                             Dict(),
                             0) # not needed here
    return Correspondence_matrix_based.compute_model(Data_Manager, nodes)
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
                                non_BCs,
                                u::AbstractMatrix{Float64},
                                F_int::AbstractMatrix{Float64},
                                F_ext::AbstractMatrix{Float64})

    # Get BC DOFs
    BCs = setdiff(1:length(vec(u)), non_BCs)

    # 1. Total force on free DOFs (external + internal/thermal)
    @views F_total = vec(F_ext) .- vec(F_int)
    # 2. Force contribution from prescribed displacements
    # TODO must be optimized

    if isempty(non_BCs)
        return nothing
    end
    if !isempty(BCs)
        F_from_BCs = K[non_BCs, BCs] * vec(u)[BCs]
    else
        F_from_BCs = zeros(length(non_BCs))
    end
    # 3. Modified force: F_total - K_fb * u_b
    F_modified = F_total[non_BCs] .- F_from_BCs
    #@info non_BCs
    # 4. Solve for free DOFs
    @views vec(u)[non_BCs] .= K[non_BCs, non_BCs] \ F_modified

    return nothing
end

function create_permutation(nnodes::Int64, dof::Int64)
    perm = Vector{Int}(undef, nnodes * dof)
    idx = 1
    for d in 1:dof
        for n in 1:nnodes
            old_idx = (n-1)*dof + d  # Row-major: node, dann dof
            perm[idx] = old_idx
            idx += 1
        end
    end
    return perm
end

function create_permutation(nodes::AbstractVector{Int64}, dof::Int64)
    perm = Vector{Int}(undef, length(nodes) * dof)
    idx = 1
    for d in 1:dof
        for n in nodes
            old_idx = (n-1)*dof + d  # Row-major: node, dann dof
            perm[idx] = old_idx
            idx += 1
        end
    end
    return perm
end
end
