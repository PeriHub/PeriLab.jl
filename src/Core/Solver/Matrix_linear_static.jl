# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Matrix_linear_static
using SparseArrays
using LinearAlgebra

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
                     block_nodes::Dict{Int64,Vector{Int64}})
    # hier wird erstmal K bestimmt, bevor es woanders hinkommt.
    # indizes prüfen
    # vec und reshape prüfen

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
    volume = datamanager.get_field("Volume")
    coor = datamanager.get_field("Coordinates")
    comm = datamanager.get_comm()

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
    #nodes::Vector{Int64} = Vector{Int64}(1:datamanager.get_nnodes())
    @inbounds @fastmath for idt in iter
        datamanager.set_iteration(idt)
        @timeit to "Linear Static" begin
            uNP1 = datamanager.get_field("Displacements", "NP1")
            deformed_coorNP1 = datamanager.get_field("Deformed Coordinates", "NP1")
            forces = datamanager.get_field("Forces", "NP1")
            force_densities = datamanager.get_field("Force Densities", "NP1")
            volume = datamanager.get_field("Volume")
            forces = datamanager.get_field("Forces", "NP1")
            K = datamanager.get_stiffness_matrix()
            external_force_densities = datamanager.get_field("External Force Densities")

            datamanager = apply_bc_dirichlet(["Displacements", "Forces", "Force Densities"],
                                             bcs,
                                             datamanager, time,
                                             step_time)
            external_force_densities.=external_forces ./ volume
            # reshape
            compute_displacements!(K,
                                   non_BCs,
                                   vec(uNP1),
                                   vec(force_densities),
                                   vec(external_force_densities))

            @views deformed_coorNP1 := coor .+ uNP1

            if "Material" in solver_options["Models"]
                # TODO rename function -> missleading, because strains are also covered. Has to be something like a factory class
                @timeit to "calculate_stresses" datamanager=calculate_stresses(datamanager,
                                                                               block_nodes,
                                                                               solver_options["Calculation"])
            end

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
                                u::Vector{Float64},
                                F_int::Vector{Float64},
                                F_ext::Vector{Float64})

    # Compute modified force: F_modified = F - K * u_prescribed
    # (u contains prescribed displacements at fixed DOFs)
    mul!(F_int, K, u)  # F_temp = K * u (in-place, no allocation)
    F_int .= F_ext .- F_temp  # F_temp = F - K*u (in-place)

    # Update displacement vector at free DOFs
    @views u[non_BCs] .= K[non_BCs, non_BCs] / F_int[non_BCs]'

    return nothing
end

end
