# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Verlet
using LinearAlgebra
using TimerOutputs
using ProgressBars: set_multiline_postfix, set_postfix
using Printf
using LoopVectorization
using PrettyTables
using Logging

include("../../Support/Helpers.jl")
using .Helpers: check_inf_or_nan, find_active_nodes, progress_bar
include("../../Support/Parameters/parameter_handling.jl")
using .Parameter_Handling:
                           get_initial_time,
                           get_fixed_dt,
                           get_final_time,
                           get_numerical_damping,
                           get_safety_factor,
                           get_max_damage

include("../../MPI_communication/MPI_communication.jl")
using .MPI_communication: find_and_set_core_value_min, find_and_set_core_value_max, barrier
include("../BC_manager.jl")
include("../../Models/Model_Factory.jl")
include("../../IO/logging.jl")
include("../../Compute/compute_field_values.jl")
using .Model_Factory
using .Boundary_conditions: apply_bc_dirichlet, apply_bc_neumann
using .Helpers: matrix_style
using .Logging_module: print_table
export init_solver
export run_solver

"""
    compute_thermodynamic_critical_time_step(nodes::AbstractVector{Int64}, datamanager::Module, lambda::Float64, Cv::Float64)

Calculate the critical time step for a thermodynamic simulation based on  [OterkusS2014](@cite).

This function iterates over a collection of nodes and computes the critical time step for each node using provided input data and parameters.

# Arguments
- `nodes::Union{SubArray, Vector{Int64}}`: The collection of nodes to calculate the critical time step for.
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
    get_cs_denominator(volume::AbstractVector{Float64}, undeformed_bond::AbstractVector{Float64})

Calculate the denominator for the critical time step calculation.

# Arguments
- `volume::AbstractVector{Float64}`: The volume field.
- `undeformed_bond::Union{SubArray,Vector{Float64},Vector{Int64}}`: The undeformed bond field.
# Returns
- `Float64`: The denominator for the critical time step calculation.
"""
function get_cs_denominator(volume::AbstractVector{Float64},
                            undeformed_bond::AbstractVector{Float64})::Float64
    return sum(volume ./ undeformed_bond)
end

"""
    compute_mechanical_critical_time_step(nodes::AbstractVector{Int64}, datamanager::Module, bulk_modulus::Float64)

Calculate the critical time step for a mechanical simulation using a bond-based approximation [LittlewoodDJ2013](@cite).

This function iterates over a collection of nodes and computes the critical time step for each node based on the given input data and parameters.

# Arguments
- `nodes::Union{SubArray, Vector{Int64}}`: The collection of nodes to calculate the critical time step for.
- `datamanager::Module`: The data manager module that provides access to required data fields.
- `bulk_modulus::Float64`: The bulk modulus used in the calculations.

# Returns
- `Float64`: The calculated critical time step for the mechanical simulation.

# Dependencies
This function depends on the following data fields from the `datamanager` module:
- `get_nlist()`: Returns the neighbor list.
- `get_field("Density")`: Returns the density field.
- `get_field("Bond Length")`: Returns the bond distance field.
- `get_field("Volume")`: Returns the volume field.
- `get_field("Horizon")`: Returns the horizon field.
"""
function compute_mechanical_critical_time_step(nodes::AbstractVector{Int64},
                                               datamanager::Module,
                                               bulk_modulus::Union{Float64,Int64,SubArray,
                                                                   Vector{Float64}})
    critical_time_step::Float64 = 1.0e50
    nlist = datamanager.get_nlist()
    density = datamanager.get_field("Density")
    undeformed_bond_length = datamanager.get_field("Bond Length")
    volume = datamanager.get_field("Volume")
    horizon = datamanager.get_field("Horizon")

    for iID in nodes
        denominator = get_cs_denominator(volume[nlist[iID]], undeformed_bond_length[iID])
        # TODO Adapt to 2D applications
        springConstant = 18.0 * maximum(bulk_modulus) /
                         (pi * horizon[iID] * horizon[iID] * horizon[iID] * horizon[iID])

        t = density[iID] / (denominator * springConstant)
        critical_time_step = test_timestep(t, critical_time_step)
    end
    return sqrt(2 * critical_time_step)
end

"""
    test_timestep(t::Float64, critical_time_step::Float64)

Compare a time step `t` with a critical time step `critical_time_step` and update `critical_time_step` if `t` is smaller.

# Arguments
- `t::Float64`: The time step to compare with `critical_time_step`.
- `critical_time_step::Float64`: The current critical time step.

# Returns
- `critical_time_step::Float64`: The updated critical time step, which is either the original `critical_time_step` or `t`, whichever is smaller.
"""
function test_timestep(t::Float64, critical_time_step::Float64)
    if t < critical_time_step
        critical_time_step = t
    end
    return critical_time_step
end

"""
    compute_crititical_time_step(datamanager::Module, block_nodes::Dict{Int64,Vector{Int64}}, mechanical::Bool, thermo::Bool)

Calculate the critical time step for a simulation considering both mechanical and thermodynamic aspects.

This function computes the critical time step by considering mechanical and thermodynamic properties of different blocks. The resulting critical time step is based on the smallest critical time step found among the blocks.

# Arguments
- `datamanager::Module`: The data manager module that provides access to required data fields and properties.
- `block_nodes::Dict{Int64, Vector{Int64}}`: A dictionary mapping block IDs to collections of nodes.
- `mechanical::Bool`: If `true`, mechanical properties are considered in the calculation.
- `thermo::Bool`: If `true`, thermodynamic properties are considered in the calculation.

# Returns
- `Float64`: The calculated critical time step based on the smallest critical time step found among the blocks.

# Dependencies
This function may depend on the following functions:
- `compute_thermodynamic_critical_time_step`: Used if `thermo` is `true` to calculate thermodynamic critical time steps.
- `compute_mechanical_critical_time_step`: Used if `mechanical` is `true` to calculate mechanical critical time steps.
- The availability of specific properties from the data manager module.

# Errors
- If required properties are not available in the data manager, it may raise an error message.
"""
function compute_crititical_time_step(datamanager::Module,
                                      block_nodes::Dict{Int64,Vector{Int64}},
                                      mechanical::Bool,
                                      thermal::Bool)
    critical_time_step::Float64 = 1.0e50
    for iblock in eachindex(block_nodes)
        if thermal
            lambda = datamanager.get_property(iblock, "Thermal Model",
                                              "Thermal Conductivity")
            # if Cv and lambda are not defined it is valid, because an analysis can take place, if material is still analysed
            if isnothing(lambda)
                if !mechanical
                    @error "No time step can be calculated, because the heat conduction is not defined."
                    return nothing
                end
            else
                t = compute_thermodynamic_critical_time_step(block_nodes[iblock],
                                                             datamanager,
                                                             lambda)
                critical_time_step = test_timestep(t, critical_time_step)
            end
        end
        if mechanical
            bulk_modulus = datamanager.get_property(iblock, "Material Model",
                                                    "Bulk Modulus")
            nu_xy = datamanager.get_property(iblock, "Material Model", "Poisson's Ratio XY")
            nu_yz = datamanager.get_property(iblock, "Material Model", "Poisson's Ratio YZ")
            nu_xz = datamanager.get_property(iblock, "Material Model", "Poisson's Ratio XZ")
            E_x = datamanager.get_property(iblock, "Material Model", "Young's Modulus X")
            E_y = datamanager.get_property(iblock, "Material Model", "Young's Modulus Y")
            E_z = datamanager.get_property(iblock, "Material Model", "Young's Modulus Z")
            c_44 = datamanager.get_property(iblock, "Material Model", "C44")
            c_55 = datamanager.get_property(iblock, "Material Model", "C55")
            c_66 = datamanager.get_property(iblock, "Material Model", "C66")
            if !isnothing(bulk_modulus)
                bulk_modulus = bulk_modulus
            elseif !isnothing(nu_xy) && !isnothing(nu_yz) && !isnothing(nu_xz)
                s11 = 1 / E_x
                s22 = 1 / E_y
                s33 = 1 / E_z
                s12 = -nu_xy / E_x
                s23 = -nu_yz / E_z
                s13 = -nu_xz / E_z
                bulk_modulus = 1 / (s11 + s22 + s33 + 2 * (s12 + s23 + s13))
            elseif !isnothing(c_44) && !isnothing(c_55) && !isnothing(c_66)
                bulk_modulus = maximum([c_44 / 2, c_55 / 2, c_66 / 2])
                #TODO: temporary solution!!!
            else
                @error "No time step for material is determined because of missing properties."
                return nothing
            end
            t = compute_mechanical_critical_time_step(block_nodes[iblock],
                                                      datamanager,
                                                      bulk_modulus)
            critical_time_step = test_timestep(t, critical_time_step)
        end
    end
    return critical_time_step
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
                     block_nodes::Dict{Int64,Vector{Int64}})
    # @info "======================="
    # @info "==== Verlet Solver ===="
    # @info "======================="

    mechanical = "Material" in solver_options["Models"]
    thermal = "Thermal" in solver_options["Models"]
    initial_time = get_initial_time(params, datamanager)
    final_time = get_final_time(params, datamanager)
    safety_factor = get_safety_factor(params)
    fixed_dt = get_fixed_dt(params)
    if fixed_dt == -1.0
        min_dt = compute_crititical_time_step(datamanager, block_nodes, mechanical, thermal)
        nsteps, dt = get_integration_steps(initial_time, final_time, safety_factor * min_dt)
    else
        nsteps, dt = get_integration_steps(initial_time, final_time, fixed_dt)
    end

    comm = datamanager.get_comm()
    dt = find_and_set_core_value_min(comm, dt)
    nsteps = find_and_set_core_value_max(comm, nsteps)
    numerical_damping = get_numerical_damping(params)
    max_damage = get_max_damage(params)

    if datamanager.get_rank() == 0
        data = ["Verlet Solver" "" "" ""
                "Initial time" "."^10 initial_time "s"
                "Final time" "."^10 final_time "s"]
        if fixed_dt != -1.0
            data = vcat(data,
                        ["Fixed time increment" "."^10 fixed_dt "s";])
        else
            data = vcat(data,
                        ["Minimal time increment" "."^10 min_dt "s"
                         "Safety factor" "."^10 safety_factor ""
                         "Time increment" "."^10 dt "s"])
        end
        data = vcat(data,
                    ["Number of steps" "."^10 nsteps ""
                     "Numerical Damping" "."^10 numerical_damping ""])

        print_table(data, datamanager)
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
    if dt <= 0
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
        datamanager::Module,
        outputs::Dict{Int64,Dict{}},
        result_files::Vector{Any},
        synchronise_field,
        write_results,
        to::TimerOutputs.TimerOutput,
        silent::Bool
    )

Run the Verlet solver for a simulation based on the strategy provided in [BobaruF2016](@cite) and  [LittlewoodDJ2023](@cite).

This function performs the Verlet solver simulation, updating various data fields and properties over a specified number of time steps.

# Arguments
- `solver_options::Dict{String,Any}`: A dictionary containing solver options and parameters.
- `block_nodes::Dict{Int64,Vector{Int64}}`: A dictionary mapping block IDs to collections of nodes.
- `bcs::Dict{Any,Any}`: A dictionary containing boundary conditions.
- `datamanager::Module`: The data manager module that provides access to data fields and properties.
- `outputs::Dict{Int64,Dict{}}`: A dictionary for output settings.
- `result_files::Vector{Any}`: A vector of result files.
- `synchronise_field`: A function for synchronization.
- `write_results`: A function to write simulation results.
- `to::TimerOutputs.TimerOutput`: A timer output.
- `silent::Bool`: A boolean flag to suppress progress bars.

# Returns
- `result_files`: A vector of updated result files.

# Dependencies
This function depends on various data fields and properties from the `datamanager` module, as well as several helper functions. It also relies on solver options and boundary conditions provided in the input parameters.

# Function Workflow
1. Initialize simulation parameters and data fields.
2. Perform Verlet integration over a specified number of time steps.
3. Update data fields and properties based on the solver options.
4. Write simulation results using the `write_results` function.
"""
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

    @info "Run Verlet Solver"
    volume = datamanager.get_field("Volume")
    active_list = datamanager.get_field("Active")
    density = datamanager.get_field("Density")
    coor = datamanager.get_field("Coordinates")
    comm = datamanager.get_comm()

    if "Material" in solver_options["Models"]
        external_forces = datamanager.get_field("External Forces")
        external_force_densities = datamanager.get_field("External Force Densities")
        a = datamanager.get_field("Acceleration")
    end

    fem_option = datamanager.fem_active()
    if fem_option
        lumped_mass = datamanager.get_field("Lumped Mass Matrix")
        fe_nodes = datamanager.get_field("FE Nodes")
    end
    active = datamanager.get_field("Active")

    dt::Float64 = solver_options["dt"]
    nsteps::Int64 = solver_options["Number of Steps"]
    time::Float64 = solver_options["Initial Time"]
    step_time::Float64 = 0
    max_cancel_damage::Float64 = solver_options["Maximum Damage"]
    numerical_damping::Float64 = solver_options["Numerical Damping"]
    max_damage::Float64 = 0
    damage_init::Bool = false
    rank = datamanager.get_rank()
    iter = progress_bar(rank, nsteps, silent)
    #nodes::Vector{Int64} = Vector{Int64}(1:datamanager.get_nnodes())
    @inbounds @fastmath for idt in iter
        datamanager.set_iteration(idt)
        @timeit to "Verlet" begin
            if "Material" in solver_options["Models"]
                uNP1 = datamanager.get_field("Displacements", "NP1")
                deformed_coorNP1 = datamanager.get_field("Deformed Coordinates", "NP1")
                forces = datamanager.get_field("Forces", "NP1")
                force_densities = datamanager.get_field("Force Densities", "NP1")
                uN = datamanager.get_field("Displacements", "N")
                vN = datamanager.get_field("Velocity", "N")
                vNP1 = datamanager.get_field("Velocity", "NP1")
                deformed_coorN = datamanager.get_field("Deformed Coordinates", "N")
                forces = datamanager.get_field("Forces", "NP1")
            end

            # if "Degradation" in solver_options["Models"]
            #     concentrationN = datamanager.get_field("Concentration", "N")
            #     concentrationNP1 = datamanager.get_field("Concentration", "NP1")
            #     concentration_fluxN = datamanager.get_field("Concentration Flux", "N")
            #     concentration_fluxNP1 = datamanager.get_field("Concentration Flux", "NP1")
            # end
            if "Damage" in solver_options["Models"]
                damage = datamanager.get_damage("NP1")
            end
            active_nodes = datamanager.get_field("Active Nodes")
            active_nodes = find_active_nodes(active_list, active_nodes,
                                             1:datamanager.get_nnodes())
            # one step more, because of init step (time = 0)
            if "Material" in solver_options["Models"]
                @views vNP1[active_nodes,
                            :] = (1 - numerical_damping) .*
                                 vN[active_nodes, :] .+
                                 0.5 * dt .* a[active_nodes, :]
                datamanager = apply_bc_dirichlet(["Velocity"], bcs, datamanager, time,
                                                 step_time)
                @views uNP1[active_nodes,
                            :] = uN[active_nodes, :] .+
                                 dt .* vNP1[active_nodes, :]
            end

            compute_parabolic_problems_before_model_evaluation(active_nodes, datamanager,
                                                               solver_options)
            # if "Degradation" in solver_options["Models"]
            #     concentrationNP1[active_nodes] = concentrationN[active_nodes] +
            #                                      delta_concentration[active_nodes]
            # end
            datamanager = apply_bc_dirichlet(["Displacements", "Temperature"],
                                             bcs,
                                             datamanager, time,
                                             step_time) #-> Dirichlet
            #needed because of optional deformation_gradient, Deformed bonds, etc.
            # all points to guarantee that the neighbors have coor as coordinates if they are not active
            if "Material" in solver_options["Models"]
                @views deformed_coorNP1[active_nodes,
                                        :] = coor[active_nodes, :] .+
                                             uNP1[active_nodes, :]
            end
            @timeit to "upload_to_cores" datamanager.synch_manager(synchronise_field,
                                                                   "upload_to_cores")
            # synch

            @timeit to "compute_models" datamanager=Model_Factory.compute_models(datamanager,
                                                                                 block_nodes,
                                                                                 dt,
                                                                                 time,
                                                                                 solver_options["Models"],
                                                                                 synchronise_field,
                                                                                 to)
            # update the current active nodes; might have been changed by the additive models

            if "Material" in solver_options["Models"]
                # TODO rename function -> missleading, because strains are also covered. Has to be something like a factory class
                @timeit to "calculate_stresses" datamanager=calculate_stresses(datamanager,
                                                                               block_nodes,
                                                                               solver_options["Calculation"])
            end

            @timeit to "download_from_cores" datamanager.synch_manager(synchronise_field,
                                                                       "download_from_cores")
            # synch
            datamanager = apply_bc_dirichlet(["Forces", "Force Densities"],
                                             bcs,
                                             datamanager, time,
                                             step_time) #-> Dirichlet
            # @timeit to "apply_bc_neumann" datamanager = Boundary_conditions.apply_bc_neumann(bcs, datamanager, time) #-> von neumann
            active_nodes = datamanager.get_field("Active Nodes")
            active_nodes = find_active_nodes(active_list, active_nodes,
                                             1:datamanager.get_nnodes())
            if "Material" in solver_options["Models"]
                check_inf_or_nan(force_densities, "Forces")

                if fem_option
                    # edit external force densities won't work so easy, because the corresponded volume is in detJ
                    # force density is for FEM part force
                    active_nodes = datamanager.get_field("Active Nodes")
                    active_nodes = find_active_nodes(fe_nodes,
                                                     active_nodes,
                                                     1:datamanager.get_nnodes())

                    forces[active_nodes, :] += external_forces[active_nodes, :]
                    force_densities[active_nodes,
                                    :] += external_force_densities[active_nodes,
                                                                   :] .+
                                          external_forces[active_nodes, :] ./
                                          volume[active_nodes]
                    a[active_nodes,
                      :] = force_densities[active_nodes, :] ./
                           lumped_mass[active_nodes] # element wise

                    active_nodes = datamanager.get_field("Active Nodes")
                    active_nodes = find_active_nodes(fe_nodes,
                                                     active_nodes,
                                                     1:datamanager.get_nnodes(),
                                                     false)
                end

                forces[active_nodes, :] += external_forces[active_nodes, :]
                @views force_densities[active_nodes,
                                       :] += external_force_densities[active_nodes,
                                                                      :] .+
                                             external_forces[active_nodes,
                                                             :] ./
                                             volume[active_nodes]
                @views a[active_nodes,
                         :] = force_densities[active_nodes, :] ./
                              density[active_nodes] # element wise
                @views forces[active_nodes,
                              :] = force_densities[active_nodes, :] .*
                                   volume[active_nodes]
            end

            compute_parabolic_problems_after_model_evaluation(active_nodes, datamanager,
                                                              solver_options, dt)

            # if "Degradation" in solver_options["Models"]
            #     delta_concentration[active_nodes] = -concentration_fluxNP1[active_nodes] .*
            #                                         dt
            # end
            if "Damage" in solver_options["Models"] #TODO gather value
                max_damage = maximum(damage[active_nodes])
                if max_damage > max_cancel_damage
                    datamanager.set_cancel(true)
                end
                if !damage_init && max_damage > 0
                    damage_init = true
                    if rank == 0 && !silent
                        set_multiline_postfix(iter,
                                              "Damage initated in step $idt [$time s]!")
                    end
                end
            end
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
            @timeit to "switch_NP1_to_N" datamanager.switch_NP1_to_N()

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

end
