# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Verlet_Solver
using LinearAlgebra
using ProgressBars: set_multiline_postfix, set_postfix
using Printf
using LoopVectorization
using Logging
using TimerOutputs: @timeit

using ...Data_Manager
using ...Helpers#: check_inf_or_nan, find_active_nodes, progress_bar, matrix_style
using ...Parameter_Handling:
                             get_initial_time,
                             get_fixed_dt,
                             get_final_time,
                             get_numerical_damping,
                             get_safety_factor,
                             get_max_damage
using ...MPI_Communication: find_and_set_core_value_min, find_and_set_core_value_max,
                            barrier
using ..Model_Factory: compute_models, compute_crititical_time_step
using ..Boundary_Conditions: apply_bc_dirichlet, apply_bc_neumann
using ...Logging_Module: print_table
include("../../Compute/compute_field_values.jl")
export init_solver
export run_solver

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
- `compute_crititical_time_step`: Used to calculate the critical time step if `dt` is not fixed.
- `get_integration_steps`: Used to determine the number of integration steps and adjust the time step.
- `find_and_set_core_value_min` and `find_and_set_core_value_max`: Used to set core values in a distributed computing environment.
"""
function init_solver(solver_options::Dict{Any,Any},
                     params::Dict,
                     bcs::Dict{Any,Any},
                     block_nodes::Dict{Int64,Vector{Int64}})
    # @info "======================="
    # @info "==== Verlet Solver ===="
    # @info "======================="

    mechanical = "Material" in solver_options["Models"]
    thermal = "Thermal" in solver_options["Models"]
    initial_time = get_initial_time(params)
    final_time = get_final_time(params)
    safety_factor = get_safety_factor(params)
    fixed_dt = get_fixed_dt(params)
    min_dt = compute_crititical_time_step(block_nodes, mechanical, thermal)
    if fixed_dt == -1.0
        nsteps, dt = get_integration_steps(initial_time, final_time, safety_factor * min_dt)
    else
        nsteps, dt = get_integration_steps(initial_time, final_time, fixed_dt)
    end

    comm = Data_Manager.get_comm()
    dt = find_and_set_core_value_min(comm, dt)
    nsteps = find_and_set_core_value_max(comm, nsteps)
    numerical_damping = get_numerical_damping(params)
    max_damage = get_max_damage(params)

    if Data_Manager.get_rank() == 0
        data = ["Verlet Solver" "" "" ""
                "Initial time" "."^10 initial_time "s"
                "Final time" "."^10 final_time "s"]
        if fixed_dt != -1.0
            data = vcat(data,
                        ["Minimal time increment" "."^10 min_dt "s"
                         "Fixed time increment" "."^10 fixed_dt "s"])
        else
            data = vcat(data,
                        ["Minimal time increment" "."^10 min_dt "s"
                         "Safety factor" "."^10 safety_factor ""
                         "Time increment" "."^10 dt "s"])
        end
        data = vcat(data,
                    ["Number of steps" "."^10 nsteps ""
                     "Numerical Damping" "."^10 numerical_damping ""])

        print_table(data)
    end
    if fixed_dt > min_dt
        @warn "Fixed time increment is larger than minimal time increment. This could lead to instability."
    end
    solver_options["Initial Time"] = initial_time
    solver_options["Final Time"] = final_time
    solver_options["dt"] = dt
    solver_options["Number of Steps"] = nsteps
    solver_options["Numerical Damping"] = numerical_damping
    solver_options["Maximum Damage"] = max_damage
    solver_options["Solver specifics"] = Dict()
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
        @error "Time step $dt [s] is not valid"
        return nothing
    end
    nsteps::Int64 = ceil((end_time - initial_time) / dt)
    dt = (end_time - initial_time) / nsteps
    return nsteps, dt
end

"""
	run_solver(
		solver_options::Dict{Any,Any},
		block_nodes::Dict{Int64,Vector{Int64}},
		bcs::Dict{Any,Any},
		outputs::Dict{Int64,Dict{}},
		result_files::Vector{Any},
		synchronise_field,
		write_results,
		silent::Bool
	)

Run the Verlet solver for a simulation based on the strategy provided in [BobaruF2016](@cite) and  [LittlewoodDJ2023](@cite).

This function performs the Verlet solver simulation, updating various data fields and properties over a specified number of time steps.

# Arguments
- `solver_options::Dict{String,Any}`: A dictionary containing solver options and parameters.
- `block_nodes::Dict{Int64,Vector{Int64}}`: A dictionary mapping block IDs to collections of nodes.
- `bcs::Dict{Any,Any}`: A dictionary containing boundary conditions.
- `outputs::Dict{Int64,Dict{}}`: A dictionary for output settings.
- `result_files::Vector{Any}`: A vector of result files.
- `synchronise_field`: A function for synchronization.
- `write_results`: A function to write simulation results.
- `silent::Bool`: A boolean flag to suppress progress bars.

# Returns
- `result_files`: A vector of updated result files.

# Dependencies
This function depends on various data fields and properties from the `Data_Manager` module, as well as several helper functions. It also relies on solver options and boundary conditions provided in the input parameters.

# Function Workflow
1. Initialize simulation parameters and data fields.
2. Perform Verlet integration over a specified number of time steps.
3. Update data fields and properties based on the solver options.
4. Write simulation results using the `write_results` function.
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

    @info "Run Verlet Solver"
    volume = Data_Manager.get_field("Volume")
    active_list = Data_Manager.get_field("Active")
    density = Data_Manager.get_field("Density")
    coor = Data_Manager.get_field("Coordinates")
    comm = Data_Manager.get_comm()

    if "Material" in solver_options["Models"]
        external_forces = Data_Manager.get_field("External Forces")
        external_force_densities = Data_Manager.get_field("External Force Densities")
        a = Data_Manager.get_field("Acceleration")
    end

    fem_option = Data_Manager.fem_active()
    if fem_option
        lumped_mass = Data_Manager.get_field("Lumped Mass Matrix")
        fe_nodes = Data_Manager.get_field("FE Nodes")
    end

    dt::Float64 = solver_options["dt"]
    nsteps::Int64 = solver_options["Number of Steps"]
    time::Float64 = solver_options["Initial Time"]
    step_time::Float64 = 0
    max_cancel_damage::Float64 = solver_options["Maximum Damage"]
    numerical_damping::Float64 = solver_options["Numerical Damping"]
    max_damage::Float64 = 0
    damage_init::Bool = false
    rank = Data_Manager.get_rank()
    iter = progress_bar(rank, nsteps, silent)
    #nodes::Vector{Int64} = Vector{Int64}(1:Data_Manager.get_nnodes())

    @inbounds @fastmath for idt in iter
        Data_Manager.set_iteration(idt)
        @timeit "Verlet" begin
            if "Material" in solver_options["Models"]
                uNP1::NodeVectorField{Float64} = Data_Manager.get_field("Displacements",
                                                                        "NP1")
                deformed_coorNP1::NodeVectorField{Float64} = Data_Manager.get_field("Deformed Coordinates",
                                                                                    "NP1")
                forces::NodeVectorField{Float64} = Data_Manager.get_field("Forces",
                                                                          "NP1")
                force_densities::NodeVectorField{Float64} = Data_Manager.get_field("Force Densities",
                                                                                   "NP1")
                uN::NodeVectorField{Float64} = Data_Manager.get_field("Displacements",
                                                                      "N")
                vN::NodeVectorField{Float64} = Data_Manager.get_field("Velocity", "N")
                vNP1::NodeVectorField{Float64} = Data_Manager.get_field("Velocity",
                                                                        "NP1")
            end

            # if "Degradation" in solver_options["Models"]
            #     concentrationN = Data_Manager.get_field("Concentration", "N")
            #     concentrationNP1 = Data_Manager.get_field("Concentration", "NP1")
            #     concentration_fluxN = Data_Manager.get_field("Concentration Flux", "N")
            #     concentration_fluxNP1 = Data_Manager.get_field("Concentration Flux", "NP1")
            # end
            if "Damage" in solver_options["Models"]
                damage = Data_Manager.get_damage("NP1")
            end
            active_nodes = Data_Manager.get_field("Active Nodes")
            active_nodes::Vector{Int64} = find_active_nodes(active_list,
                                                            active_nodes,
                                                            1:Data_Manager.get_nnodes())

            # one step more, because of init step (time = 0)
            if "Material" in solver_options["Models"]
                c = 0.5 * dt
                @. @views vNP1[active_nodes, :] = (1 - numerical_damping) .*
                                                  vN[active_nodes, :] +
                                                  c * a[active_nodes, :]

                # @views vNP1[active_nodes, :] .= (1 - numerical_damping) .*
                #                                 vN[active_nodes, :] .+
                #                                 0.5 * dt .*
                #                                 a[active_nodes, :]
                apply_bc_dirichlet(["Velocity"], bcs, time,
                                   step_time)
                @. @views uNP1[active_nodes,
                :] = uN[active_nodes, :] +
                                                  dt * vNP1[active_nodes, :]
            end

            compute_parabolic_problems_before_model_evaluation(active_nodes,
                                                               solver_options)
            # if "Degradation" in solver_options["Models"]
            #     concentrationNP1[active_nodes] = concentrationN[active_nodes] +
            #                                      delta_concentration[active_nodes]
            # end
            apply_bc_dirichlet(["Displacements", "Temperature"],
                               bcs,
                               time,
                               step_time) #-> Dirichlet
            #needed because of optional deformation_gradient, Deformed bonds, etc.
            # all points to guarantee that the neighbors have coor as coordinates if they are not active
            if "Material" in solver_options["Models"]
                @. @views deformed_coorNP1[active_nodes,
                :] = coor[active_nodes, :] +
                                                              uNP1[active_nodes, :]
            else
                deformed_coorNP1 = Data_Manager.get_field("Deformed Coordinates", "NP1")
                @. @views deformed_coorNP1[active_nodes, :] = coor[active_nodes, :]
            end
            @timeit "upload_to_cores" Data_Manager.synch_manager(synchronise_field,
                                                                 "upload_to_cores")
            # synch

            @timeit "compute_models" compute_models(block_nodes,
                                                    dt,
                                                    time,
                                                    solver_options["Models"],
                                                    synchronise_field)
            # update the current active nodes; might have been changed by the additive models

            if "Material" in solver_options["Models"]
                # TODO rename function -> missleading, because strains are also covered. Has to be something like a factory class
                @timeit "calculate_stresses" calculate_stresses(block_nodes,
                                                                solver_options["Calculation"])
            end

            @timeit "download_from_cores" Data_Manager.synch_manager(synchronise_field,
                                                                     "download_from_cores")
            # synch
            apply_bc_dirichlet(["Forces", "Force Densities"],
                               bcs,
                               time,
                               step_time) #-> Dirichlet
            # @timeit "apply_bc_neumann" apply_bc_neumann(bcs, time) #-> von neumann
            active_nodes = get_field("Active Nodes")
            active_nodes = find_active_nodes(active_list, active_nodes,
                                             1:Data_Manager.get_nnodes())
            if "Material" in solver_options["Models"]
                check_inf_or_nan(force_densities, "Forces")

                if fem_option
                    # edit external force densities won't work so easy, because the corresponded volume is in detJ
                    # force density is for FEM part force
                    active_nodes = Data_Manager.get_field("Active Nodes")
                    active_nodes = find_active_nodes(fe_nodes,
                                                     active_nodes,
                                                     1:Data_Manager.get_nnodes())

                    @. forces[active_nodes, :] += external_forces[active_nodes, :]
                    @. force_densities[active_nodes,
                    :] += external_force_densities[active_nodes,
                                                                                    :] +
                                                           external_forces[active_nodes,
                                                                           :] /
                                                           volume[active_nodes]
                    @. a[active_nodes,
                    :] = force_densities[active_nodes, :] /
                                            lumped_mass[active_nodes] # element wise

                    active_nodes = Data_Manager.get_field("Active Nodes")
                    active_nodes = find_active_nodes(fe_nodes,
                                                     active_nodes,
                                                     1:Data_Manager.get_nnodes(),
                                                     false)
                end

                @. @views forces[active_nodes, :] = external_forces[active_nodes, :]
                @. @views force_densities[active_nodes, :] += external_force_densities[active_nodes,
                                                                                       :] +
                                                              external_forces[active_nodes,
                                                                              :] /
                                                              volume[active_nodes]
                @. @views a[active_nodes, :] = force_densities[active_nodes, :] /
                                               density[active_nodes] # element wise
                @. @views forces[active_nodes, :] += force_densities[active_nodes, :] *
                                                     volume[active_nodes]
            end

            compute_parabolic_problems_after_model_evaluation(active_nodes, solver_options,
                                                              dt)

            # if "Degradation" in solver_options["Models"]
            #     delta_concentration[active_nodes] = -concentration_fluxNP1[active_nodes] .*
            #                                         dt
            # end
            if "Damage" in solver_options["Models"] #TODO gather value
                max_damage = maximum(damage[active_nodes])
                if max_damage > max_cancel_damage
                    @info "Simulation cancelled, max. damage reached!"
                    Data_Manager.set_cancel(true)
                end
                if !damage_init && max_damage > 0
                    damage_init = true
                    if rank == 0 && !silent
                        set_multiline_postfix(iter,
                                              "Damage initated in step $idt [$time s]!")
                    end
                end
            end
            @timeit "write_results" result_files=write_results(result_files, time,
                                                               max_damage, outputs)
            # for file in result_files
            #     flush(file)
            # end
            if Data_Manager.get_cancel()
                if rank == 0 && !silent
                    set_multiline_postfix(iter, "Simulation canceled!")
                end
                break
            end
            @timeit "switch_NP1_to_N" Data_Manager.switch_NP1_to_N()

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
    Data_Manager.set_current_time(time - dt)
    return result_files
end

end
