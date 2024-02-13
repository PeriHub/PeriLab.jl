# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module External_solver
using DifferentialEquations
using LinearAlgebra
using TimerOutputs
using ProgressBars: set_multiline_postfix, set_postfix
using Printf
using Reexport

include("../../Support/helpers.jl")
@reexport using .Helpers: check_inf_or_nan, find_active, progress_bar
include("../../Support/Parameters/parameter_handling.jl")
@reexport using .Parameter_Handling: get_initial_time, get_fixed_dt, get_final_time, get_numerical_damping, get_safety_factor, get_max_damage, get_nsteps

include("../../MPI_communication/MPI_communication.jl")
include("../BC_manager.jl")
include("../../Physics/Physics_Factory.jl")
using .Physics
using .Boundary_conditions: apply_bc
using .Helpers: matrix_style
export init_solver
export run_solver

"""
    init_solver(params::Dict, datamanager::Module, block_nodes::Dict{Int64,Vector{Int64}}, mechanical::Bool, thermo::Bool)

Initialize the Verlet solver for a simulation.

This function sets up the Verlet solver for a simulation by initializing various parameters and calculating the time step based on provided parameters or critical time step calculations.

# Arguments
- `params::Dict`: A dictionary containing simulation parameters.
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
function init_solver(params::Dict, datamanager::Module, block_nodes::Dict{Int64,Vector{Int64}}, mechanical::Bool, thermo::Bool)
    @info "======================="
    @info "==== Static Solver ===="
    @info "======================="
    initial_time = get_initial_time(params)
    final_time = get_final_time(params)

    @info "Initial time: " * string(initial_time) * " [s]"
    @info "Final time: " * string(final_time) * " [s]"
    dt = get_fixed_dt(params)
    if dt == -1.0
        nsteps = get_nsteps(params)
    else
        nsteps = Int64((final_time - initial_time) / dt)
    end
    dt = (final_time - initial_time) / nsteps
    numerical_damping = get_numerical_damping(params)
    max_damage = get_max_damage(params)

    @info "Time increment: " * string(dt) * " [s]"
    @info "Number of steps: " * string(nsteps)
    @info "Numerical Damping " * string(numerical_damping)
    return initial_time, dt, nsteps, numerical_damping, max_damage
end


"""
    run_solver(
        solver_options::Dict{String,Any},
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
5. Return the updated `result_files` vector.
"""
function run_solver(solver_options::Dict{String,Any}, block_nodes::Dict{Int64,Vector{Int64}}, bcs::Dict{Any,Any}, datamanager::Module, outputs::Dict{Int64,Dict{}}, result_files::Vector{Dict}, synchronise_field, write_results, to::TimerOutputs.TimerOutput, silent::Bool)
    @info "Run Static Solver"
    dof = datamanager.get_dof()
    nnodes = datamanager.get_nnodes()
    volume = datamanager.get_field("Volume")
    density = datamanager.get_field("Density")
    coor = datamanager.get_field("Coordinates")
    uNP1 = datamanager.get_field("Displacements", "NP1")

    deformed_coorNP1 = datamanager.get_field("Deformed Coordinates", "NP1")
    if solver_options["Material Models"]
        forces = datamanager.get_field("Forces", "NP1")
        forces_density = datamanager.get_field("Force Densities", "NP1")
        uN = datamanager.get_field("Displacements", "N")
        vN = datamanager.get_field("Velocity", "N")
        vNP1 = datamanager.get_field("Velocity", "NP1")
        a = datamanager.get_field("Acceleration")
        coor = datamanager.get_field("Coordinates")
        deformed_coorN = datamanager.get_field("Deformed Coordinates", "N")
    end

    active = datamanager.get_field("Active")
    update_list = datamanager.get_field("Update List")

    dt::Float64 = solver_options["dt"]
    nsteps::Int64 = solver_options["nsteps"]
    start_time::Float64 = solver_options["Initial Time"]
    max_cancel_damage::Float64 = solver_options["Maximum Damage"]
    step_time::Float64 = 0
    numerical_damping::Float64 = solver_options["Numerical Damping"]
    max_damage::Float64 = 0
    damage_init::Bool = false
    rank = datamanager.get_rank()
    iter = progress_bar(rank, nsteps, silent)
    nodes::Vector{Int64} = []

    for idt in iter
        # synch
        #einfaches Beispiel bauen;
        # -> das Beispiel von der Website mit dem Datamanager; 
        # -> weniger Punkte
        p = (datamanager, block_nodes, dt, step_time, solver_options, bcs, synchronise_field, to)
        tspan = (step_time, step_time + dt)
        #datamanager = Boundary_conditions.apply_bc(bcs, datamanager, step_time + dt)
        uNP1 = datamanager.get_field("Displacements", "NP1")
        vNP1 = datamanager.get_field("Velocity", "NP1")
        #TODO tolerance in solver options
        prob = SecondOrderODEProblem(ode_interface, vN, uNP1, tspan, p)
        res = DifferentialEquations.solve(prob, save_everystep=false)

        # the switch N to NP1 already happend in solve process
        vNP1, uNP1 .= res.u[end]

        @timeit to "write_results" result_files = write_results(result_files, start_time + step_time, max_damage, outputs, datamanager)

        if datamanager.get_cancel()
            set_multiline_postfix(iter, "Simulation canceled!")
            break
        end
        #@timeit to "switch_NP1_to_N" datamanager.switch_NP1_to_N()
        #update_list .= true
        step_time += dt
        if idt < 10 || nsteps - idt < 10 || idt % ceil(nsteps / 10) == 0
            @info "Step: $idt / $(nsteps+1) [$step_time s]"
        end
        if rank == 0 && !silent
            set_postfix(iter, t=@sprintf("%.4e", step_time))
        end

    end
    # end
    return result_files
end


function ode_interface(ddu, u, p, t)
    datamanager, block_nodes, dt, step_time, solver_options, bcs, synchronise_field, to = p
    println(t)
    update_list = datamanager.get_field("Update List")
    uNP1 = datamanager.get_field("Displacements", "NP1")
    uNP1[:, :] = copy(u[:, :])
    datamanager = Boundary_conditions.apply_bc(bcs, datamanager, step_time + dt)
    uNP1 = datamanager.get_field("Displacements", "NP1")
    deformed_coorNP1 = datamanager.get_field("Deformed Coordinates", "NP1")
    coor = datamanager.get_field("Coordinates")
    density = datamanager.get_field("Density")

    uNP1 = datamanager.get_field("Displacements", "NP1")

    deformed_coorNP1[:, :] = coor[:, :] + uNP1[:, :]

    datamanager.synch_manager(synchronise_field, "upload_to_cores")

    datamanager = Physics.compute_models(datamanager, block_nodes, dt, step_time, solver_options, synchronise_field, to)
    datamanager.synch_manager(synchronise_field, "download_from_cores")
    forces_density = datamanager.get_field("Force Densities", "NP1")
    ddu .= forces_density[:, :] ./ density[:] # element wise
    # to avoid that forces accumulate during iteration
    #    forces_density .= 0
    @timeit to "switch_NP1_to_N" datamanager.switch_NP1_to_N()
    update_list .= true

    #println(t)
end

function compute_residual(residual, sol, po, t)
    #datamanager, block_nodes, dt, step_time, solver_options, bcs, synchronise_field, to = p
    #datamanager = Boundary_conditions.apply_bc(bcs, datamanager, step_time + dt)
    # uNP1 = datamanager.get_field("Displacements", "NP1")
    #TODO external force


    #residual = sol(t[end])

    # force_internal is the reaction of force_external
    # by Newtons law it must be zero

    #residual = force_internal + force_extern
    #residual_norm2 = norm(residual)
    # is used to avoid single values which are out of control
    #residual_norm_inf = maximum(eachcol(abs.(residual)))
    #return residual_norm2 + 20.0 * residual_norm_inf
end

end