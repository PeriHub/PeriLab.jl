# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Verlet_Solver
using LinearAlgebra
using Printf
using LoopVectorization
using Logging
using TimerOutputs: @timeit

using ...Data_Manager
using ...PeriLabExceptions: @abort
using ...Helpers#: check_inf_or_nan, find_active_nodes, progress_bar, matrix_style
using ...ParameterSpec: @params, register_solver
using ...InputDeck: SolverParams, SolverOptions, CalculationOptions, start_time, end_time
using ..Solver_Common: TimeWindow, SolverRun, SolverContext, StepClock, synchronise_field,
                       compute_parabolic_problems_before_model_evaluation,
                       compute_parabolic_problems_after_model_evaluation

@params "Explicit Verlet time integration with the stable time step; for dynamic problems" struct VerletParams <: SolverOptions
    "Factor on the stable time step; values below 1 add a safety margin"
    safety_factor::Float64 = opt("Safety Factor"; default = 1.0, min = 0)
    "Time step to use instead of the stable time step; -1: compute it from the stable time step and the safety factor"
    fixed_dt::Float64 = opt("Fixed dt"; default = -1.0, quantity = :time)
    "Damps the velocity: the velocity of the previous step is scaled by one minus this value"
    numerical_damping::Float64 = opt("Numerical Damping"; default = 0.0, min = 0)
end
__init__() = register_solver("Verlet", VerletParams)
using ...MPI_Communication: find_and_set_core_value_min, find_and_set_core_value_max
using ..Model_Factory: compute_models, compute_crititical_time_step
using ..Boundary_Conditions: apply_bc_dirichlet, apply_bc_neumann, BoundaryCondition
using ...Logging_Module: print_table
include("../../Compute/compute_field_values.jl")
export init_solver
import ..Solver_Common: prepare_solver, time_step!

"""
	init_solver(options::VerletParams, params::SolverParams, models::Vector{String}, bcs::Dict{String,BoundaryCondition}, block_nodes::Dict{Int64,Vector{Int64}})

Initialize the Verlet solver for a simulation.

This function sets up the Verlet solver for a simulation by initializing various parameters and calculating the time step based on provided parameters or critical time step calculations.

# Arguments
- `options::VerletParams`: The Verlet options of the step.
- `params::SolverParams`: The solver parameters of the current step.
- `models::Vector{String}`: The model categories of the step.
- `bcs::Dict{String,BoundaryCondition}`: Boundary conditions
- `block_nodes::Dict{Int64,Vector{Int64}}`: A dictionary mapping block IDs to collections of nodes.

# Returns
- `TimeWindow`: initial and final time, time step and number of steps.
- `nothing`: Verlet has no derived state.

# Dependencies
This function may depend on the following functions:
- `start_time` and the fields of `params` (final time, safety factor, fixed time step): The simulation parameters.
- `compute_crititical_time_step`: Used to calculate the critical time step if `dt` is not fixed.
- `get_integration_steps`: Used to determine the number of integration steps and adjust the time step.
- `find_and_set_core_value_min` and `find_and_set_core_value_max`: Used to set core values in a distributed computing environment.
"""
function init_solver(options::VerletParams,
                     params::SolverParams,
                     models::Vector{String},
                     bcs::Dict{String,BoundaryCondition},
                     block_nodes::Dict{Int64,Vector{Int64}})
    # @info "======================="
    # @info "==== Verlet Solver ===="
    # @info "======================="

    mechanical = "Material" in models
    thermal = "Thermal" in models
    initial_time = start_time(params, Data_Manager.get_current_time())
    final_time = end_time(params, Data_Manager.get_current_time())
    safety_factor = options.safety_factor
    fixed_dt = options.fixed_dt
    min_dt = compute_crititical_time_step(block_nodes, mechanical, thermal)
    if fixed_dt == -1.0
        nsteps, dt = get_integration_steps(initial_time, final_time, safety_factor * min_dt)
    else
        nsteps, dt = get_integration_steps(initial_time, final_time, fixed_dt)
    end

    comm = Data_Manager.get_comm()
    dt = find_and_set_core_value_min(comm, dt)
    nsteps = find_and_set_core_value_max(comm, nsteps)
    numerical_damping = options.numerical_damping
    max_damage = params.maximum_damage

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
    return TimeWindow(initial_time, final_time, dt, nsteps), nothing
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

"""
    VerletWorkspace

Fields Verlet reads in every step.
"""
struct VerletWorkspace
    volume::Vector{Float64}
    density::Vector{Float64}
    coor::Matrix{Float64}
    active_list::Vector{Bool}
    fem_option::Bool
    fe_nodes::Vector{Bool}
    numerical_damping::Float64
end

"""
    prepare_solver(run::SolverRun{VerletParams}, ctx::SolverContext)

Fetches the fields the Verlet steps read.
"""
function prepare_solver(run::SolverRun{VerletParams}, ctx::SolverContext)
    fem_option = Data_Manager.fem_active()
    return VerletWorkspace(Data_Manager.get_field("Volume"),
                           Data_Manager.get_field("Density"),
                           Data_Manager.get_field("Coordinates"),
                           Data_Manager.get_field("Active"),
                           fem_option,
                           fem_option ? Data_Manager.get_field("FE Nodes") : Bool[],
                           run.params.options.numerical_damping)
end

"""
    time_step!(ws::VerletWorkspace, run::SolverRun{VerletParams}, ctx::SolverContext, clock::StepClock)

One explicit Verlet step following [BobaruF2016](@cite) and [LittlewoodDJ2023](@cite):
velocities and displacements, Dirichlet conditions, models, stresses and forces.
Writing results, the damage limit, `NP1 → N` and the time are done by the shared loop.
"""
function time_step!(ws::VerletWorkspace, run::SolverRun{VerletParams}, ctx::SolverContext,
                    clock::StepClock)
    block_nodes = ctx.block_nodes
    bcs = ctx.bcs
    time = clock.time
    step_time = clock.step_time
    dt = clock.dt
    volume = ws.volume
    density = ws.density
    coor = ws.coor
    active_list = ws.active_list
    fe_nodes = ws.fe_nodes
    fem_option = ws.fem_option
    numerical_damping = ws.numerical_damping
    @inbounds @fastmath begin
        if "Material" in run.models
            uNP1::NodeVectorField{Float64} = Data_Manager.get_field("Displacements",
                                                                    "NP1")

            uN::NodeVectorField{Float64} = Data_Manager.get_field("Displacements",
                                                                  "N")
            vN::NodeVectorField{Float64} = Data_Manager.get_field("Velocity", "N")
            vNP1::NodeVectorField{Float64} = Data_Manager.get_field("Velocity",
                                                                    "NP1")
            aN::NodeVectorField{Float64} = Data_Manager.get_field("Acceleration", "N")
            aNP1::NodeVectorField{Float64} = Data_Manager.get_field("Acceleration",
                                                                    "NP1")
        end
        deformed_coorNP1::NodeVectorField{Float64} = Data_Manager.get_field("Deformed Coordinates",
                                                                            "NP1")
        active_nodes = Data_Manager.get_field("Active Nodes")::Vector{Int64}
        active_nodes = find_active_nodes(active_list,
                                         active_nodes,
                                         1:Data_Manager.get_nnodes())

        # one step more, because of init step (time = 0)

        update_velocity_and_position!(dt,
                                      numerical_damping,
                                      active_nodes,
                                      "Material" in run.models,
                                      bcs,
                                      time,
                                      step_time)

        compute_parabolic_problems_before_model_evaluation(active_nodes, run)
        # if "Degradation" in run.models
        #     concentrationNP1[active_nodes] = concentrationN[active_nodes] +
        #                                      delta_concentration[active_nodes]
        # end
        apply_bc_dirichlet(["Displacements", "Temperature"],
                           bcs,
                           time,
                           step_time) #-> Dirichlet
        #needed because of optional deformation_gradient, Deformed bonds, etc.
        # all points to guarantee that the neighbors have coor as coordinates if they are not active

        update_deformed_coordinates!(deformed_coorNP1, coor, active_nodes,
                                     "Material" in run.models)

        @timeit "upload_to_cores" Data_Manager.synch_manager(synchronise_field,
                                                             "upload_to_cores")
        # synch

        @timeit "compute_models" compute_models(block_nodes,
                                                dt,
                                                time,
                                                run.models,
                                                synchronise_field)
        # update the current active nodes; might have been changed by the additive models

        if "Material" in run.models
            # TODO rename function -> missleading, because strains are also covered. Has to be something like a factory class
            @timeit "calculate_stresses" calculate_stresses(block_nodes,
                                                            run.calculation)
        end

        @timeit "download_from_cores" Data_Manager.synch_manager(synchronise_field,
                                                                 "download_from_cores")
        # synch
        @timeit "apply_bc_dirichlet" apply_bc_dirichlet(["Forces", "Force Densities"],
                                                        bcs,
                                                        time,
                                                        step_time) #-> Dirichlet
        # @timeit "apply_bc_neumann" apply_bc_neumann(bcs, time) #-> von neumann
        active_nodes = get_field("Active Nodes")::Vector{Int64}
        active_nodes = find_active_nodes(active_list, active_nodes,
                                         1:Data_Manager.get_nnodes())
        if "Material" in run.models
            force_densities = Data_Manager.get_field("Force Densities", "NP1")
            check_inf_or_nan(force_densities, "Force Densities")

            if fem_option
                # edit external force densities won't work so easy, because the corresponded volume is in detJ
                # force density is for FEM part force
                active_nodes = Data_Manager.get_field("Active Nodes")::Vector{Int64}
                active_nodes = find_active_nodes(fe_nodes,
                                                 active_nodes,
                                                 1:Data_Manager.get_nnodes())

                compute_forces_add!(volume, density, active_nodes,
                                    "Material" in run.models)

                active_nodes = Data_Manager.get_field("Active Nodes")::Vector{Int64}
                active_nodes = find_active_nodes(fe_nodes,
                                                 active_nodes,
                                                 1:Data_Manager.get_nnodes(),
                                                 false)
            end
            @timeit "compute_forces!" begin
                compute_forces!(volume, density, active_nodes,
                                "Material" in run.models)
            end
        end

        compute_parabolic_problems_after_model_evaluation(active_nodes, run, dt)
    end
    return nothing
end

function update_velocity_and_position!(dt::Float64,
                                       numerical_damping::Float64,
                                       active_nodes::AbstractVector{Int},
                                       has_material::Bool,
                                       bcs,
                                       time::Float64,
                                       step_time::Float64)
    if !has_material
        return nothing
    end

    vN::NodeVectorField{Float64} = Data_Manager.get_field("Velocity", "N")
    vNP1::NodeVectorField{Float64} = Data_Manager.get_field("Velocity", "NP1")
    aN::NodeVectorField{Float64} = Data_Manager.get_field("Acceleration", "N")
    uN::NodeVectorField{Float64} = Data_Manager.get_field("Displacements", "N")
    uNP1::NodeVectorField{Float64} = Data_Manager.get_field("Displacements", "NP1")

    c::Float64 = 0.5 * dt
    damping::Float64 = 1.0 - numerical_damping

    # Prädiktor-Halbschritt für die Geschwindigkeit
    @inbounds for i in active_nodes
        for d in axes(vNP1, 2)
            vNP1[i, d] = damping * vN[i, d] + c * aN[i, d]
        end
    end

    apply_bc_dirichlet(["Velocity"], bcs, time, step_time)

    @inbounds for i in active_nodes
        for d in axes(uNP1, 2)
            uNP1[i, d] = uN[i, d] + dt * vNP1[i, d]
        end
    end

    return nothing
end

function compute_forces!(volume::Vector{Float64}, density::Vector{Float64},
                         active_nodes::AbstractVector{Int}, has_material::Bool)
    if !has_material
        return nothing
    end
    external_forces::NodeVectorField{Float64} = Data_Manager.get_field("External Forces")
    external_force_densities::NodeVectorField{Float64} = Data_Manager.get_field("External Force Densities")
    forces::NodeVectorField{Float64} = Data_Manager.get_field("Forces", "NP1")
    force_densities::NodeVectorField{Float64} = Data_Manager.get_field("Force Densities",
                                                                       "NP1")
    aNP1::NodeVectorField{Float64} = Data_Manager.get_field("Acceleration", "NP1")
    @inbounds for i in active_nodes
        for d in axes(forces, 2)
            forces[i, d] = external_forces[i, d]
            fd = force_densities[i, d] + external_force_densities[i, d] +
                 external_forces[i, d] / volume[i]
            force_densities[i, d] = fd
            aNP1[i, d] = fd / density[i]
            forces[i, d] += fd * volume[i]
        end
    end
    return nothing
end

function compute_forces_add!(volume, mass, active_nodes, has_material::Bool)
    if !has_material
        return nothing
    end
    external_forces = Data_Manager.get_field("External Forces")
    external_force_densities = Data_Manager.get_field("External Force Densities")
    forces = Data_Manager.get_field("Forces", "NP1")
    force_densities = Data_Manager.get_field("Force Densities", "NP1")
    aNP1 = Data_Manager.get_field("Acceleration", "NP1")
    @inbounds for i in active_nodes
        for d in axes(forces, 2)
            forces[i, d] += external_forces[i, d]
            fd = force_densities[i, d] + external_force_densities[i, d] +
                 external_forces[i, d] / volume[i]
            force_densities[i, d] = fd
            aNP1[i, d] = fd / mass[i]
            forces[i, d] += fd * volume[i]
        end
    end
    return nothing
end
function update_deformed_coordinates!(deformed_coorNP1, coor, active_nodes,
                                      has_material::Bool)
    has_material ?
    uNP1::NodeVectorField{Float64} = Data_Manager.get_field("Displacements",
                                                            "NP1") : nothing
    @inbounds for i in active_nodes
        for d in axes(coor, 2)
            if has_material
                deformed_coorNP1[i, d] = coor[i, d] + uNP1[i, d]
            else
                deformed_coorNP1[i, d] = coor[i, d]
            end
        end
    end
    return nothing
end

end
