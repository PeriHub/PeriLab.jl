# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module External
using LinearAlgebra
using TimerOutputs
using ProgressBars: set_multiline_postfix, set_postfix
using Printf
using Reexport
using PrettyTables
using ImplicitAD
using NLsolve
using Logging

include("../../Support/helpers.jl")
@reexport using .Helpers: check_inf_or_nan, find_active, progress_bar
include("../../Support/Parameters/parameter_handling.jl")
@reexport using .Parameter_Handling: get_initial_time, get_fixed_dt, get_final_time, get_numerical_damping, get_safety_factor, get_max_damage

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

Initialize the External solver for a simulation.

This function sets up the External solver for a simulation by initializing various parameters and calculating the time step based on provided parameters or critical time step calculations.

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
    # @info "======================="
    # @info "==== External Solver ===="
    # @info "======================="

    initial_time = get_initial_time(params)
    final_time = get_final_time(params)
    safety_factor = get_safety_factor(params)
    dt = get_fixed_dt(params)
    nsteps, dt = get_integration_steps(initial_time, final_time, safety_factor * dt)
    comm = datamanager.get_comm()
    dt = find_and_set_core_value_min(comm, dt)
    nsteps = find_and_set_core_value_max(comm, nsteps)
    numerical_damping = get_numerical_damping(params)
    max_damage = get_max_damage(params)

    # @info "Safety Factor: " * string(safety_factor)
    # @info "Time increment: " * string(dt) * " [s]"
    # @info "Number of steps: " * string(nsteps)
    # @info "Numerical Damping " * string(numerical_damping)
    data = [
        "External Solver" "" "" "";
        "Initial time" "."^10 initial_time "s";
        "Final time" "."^10 final_time "s";
        "Minimal time increment" "."^10 dt "s";
        "Fixed time increment" "."^10 dt "s";
        "Safety factor" "."^10 safety_factor "";
        "Time increment" "."^10 dt "s";
        "Number of steps" "."^10 nsteps "";
        "Numerical Damping" "."^10 numerical_damping "";
    ]
    if datamanager.get_rank() == 0
        pretty_table(
            data;
            body_hlines        = [1,9],
            body_hlines_format = Tuple('─' for _ = 1:4),
            cell_alignment     = Dict((1, 1) => :l),
            formatters         = ft_printf("%10.1f", 2),
            highlighters       = (
                hl_cell([(1, 1)], crayon"bold"),
                hl_col(2, crayon"dark_gray")
            ),
            show_header        = false,
            tf                 = tf_borderless
        )
        pretty_table(
            current_logger().loggers[2].logger.stream,
            data;
            body_hlines        = [1,9],
            body_hlines_format = Tuple('─' for _ = 1:4),
            cell_alignment     = Dict((1, 1) => :l),
            formatters         = ft_printf("%10.1f", 2),
            highlighters       = (
                hl_cell([(1, 1)], crayon"bold"),
                hl_col(2, crayon"dark_gray")
            ),
            show_header        = false,
            tf                 = tf_borderless
        )
    end
    return initial_time, dt, nsteps, numerical_damping, max_damage
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

Run the External solver for a simulation based on the strategy provided in [BobaruF2016](@cite) and  [LittlewoodDJ2023](@cite).
    
This function performs the External solver simulation, updating various data fields and properties over a specified number of time steps.

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
2. Perform External integration over a specified number of time steps.
3. Update data fields and properties based on the solver options.
4. Write simulation results using the `write_results` function.
"""
function run_solver(solver_options::Dict{String,Any}, block_nodes::Dict{Int64,Vector{Int64}}, bcs::Dict{Any,Any}, datamanager::Module, outputs::Dict{Int64,Dict{}}, result_files::Vector{Dict}, synchronise_field, write_results, to::TimerOutputs.TimerOutput, silent::Bool)

    atexit(() -> datamanager.set_cancel(true))

    @info "Run External Solver"
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
    
    if solver_options["Damage Models"]
        damage = datamanager.get_damage("NP1")
    end
    active = datamanager.get_field("Active")


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
        @timeit to "External" begin
            nodes = find_active(active[1:nnodes])
            # one step more, because of init step (time = 0)

            vNP1[nodes, :] = (1 - numerical_damping) .* vN[nodes, :] + 0.5 * dt .* a[nodes, :]
            uNP1[nodes, :] = uN[nodes, :] + dt .* vNP1[nodes, :]
     
           
            @timeit to "apply_bc" datamanager = Boundary_conditions.apply_bc(bcs, datamanager, step_time)
            #if solver_options["Material Models"]
            #needed because of optional deformation_gradient, Deformed bonds, etc.
            # all points to guarantee that the neighbors have coor as coordinates if they are not active
            deformed_coorNP1[:, :] = coor[:, :] + uNP1[:, :]

            @timeit to "upload_to_cores" datamanager.synch_manager(synchronise_field, "upload_to_cores")
            # synch

            #@timeit to "compute_models" datamanager = Physics.compute_models(datamanager, block_nodes, dt, step_time, solver_options, synchronise_field, to)
            solving(datamanager, block_nodes, dt, step_time, solver_options, synchronise_field, to)
            @timeit to "download_from_cores" datamanager.synch_manager(synchronise_field, "download_from_cores")
            # synch
            @timeit to "second apply_bc" datamanager = Boundary_conditions.apply_bc(bcs, datamanager, step_time)

           
            check_inf_or_nan(forces_density, "Forces")

            a[nodes, :] = forces_density[nodes, :] ./ density[nodes] # element wise
            forces[nodes, :] = forces_density[nodes, :] .* volume[nodes]
           
            @timeit to "write_results" result_files = write_results(result_files, start_time + step_time, max_damage, outputs, datamanager)
            # for file in result_files
            #     flush(file)
            # end
            if datamanager.get_cancel()
                set_multiline_postfix(iter, "Simulation canceled!")
                break
            end
            @timeit to "switch_NP1_to_N" datamanager.switch_NP1_to_N()

            step_time += dt
            if idt < 10 || nsteps - idt < 10 || idt % ceil(nsteps / 10) == 0
                @info "Step: $idt / $(nsteps+1) [$step_time s]"
            end
            if rank == 0 && !silent
                set_postfix(iter, t=@sprintf("%.4e", step_time))
            end

        end
    end
    return result_files
end



function residual(r,y,x,p)
    return external_force-internal_force
end


function solving(datamanager, block_nodes, dt, step_time, solver_options, synchronise_field, to)
    @timeit to "compute_models" datamanager = Physics.compute_models(datamanager, block_nodes, dt, step_time, solver_options, synchronise_field, to)
end
end