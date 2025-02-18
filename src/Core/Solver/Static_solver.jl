# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Static_solver
using LinearAlgebra
using NLsolve
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
include("../BC_manager.jl")
include("../../Models/Model_Factory.jl")
include("../../IO/logging.jl")
include("../../Compute/compute_field_values.jl")
using .Model_Factory
using .Boundary_conditions:
    apply_bc_dirichlet, apply_bc_neumann, find_bc_free_dof, apply_bc_dirichlet_force
using .Helpers: matrix_style
using .Logging_module: print_table
export init_solver
export run_solver


"""
    init_solver(params::Dict, bcs::Dict{Any,Any}, datamanager::Module, block_nodes::Dict{Int64,Vector{Int64}}, mechanical::Bool, thermo::Bool)

Initialize the Static solver for a simulation.

This function sets up the Static solver for a simulation by initializing various parameters.

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
- `xtol::Float64`: solution tolerance; mininum step size value between two iterations
- `ftol::Float64`: residual tolerance; mininum residual value of the function
- `iterations::Int64`: maximum number of iterations per time step
- `show_trace::Bool`: show the iteration steps of the solver


# Dependencies
This function may depend on the following functions:
- `get_initial_time`, `get_final_time`, `get_safety_factor`, `get_fixed_dt`: Used to retrieve simulation parameters.
- `get_integration_steps`: Used to determine the number of integration steps and adjust the time step.
- `find_and_set_core_value_min` and `find_and_set_core_value_max`: Used to set core values in a distributed computing environment.
"""
function init_solver(
    params::Dict,
    bcs::Dict{Any,Any},
    datamanager::Module,
    block_nodes::Dict{Int64,Vector{Int64}},
    mechanical::Bool,
    thermo::Bool,
)
    # @info "==============================="
    # @info "==== NLsolve Static Solver ===="
    # @info "==============================="

    initial_time = get_initial_time(params)
    final_time = get_final_time(params)
    fixed_dt = get_fixed_dt(params)
    if fixed_dt == -1.0
        if haskey(params, "nstep")
            nsteps = params["nstep"]
            dt = (final_time - initial_time) / (nsteps - 1)
        else
            nsteps = Int64(1)
            dt = final_time - initial_time
        end
    else
        if haskey(params, "nstep")
            @warn "nsteps and fixed dt is defined. Fixed dt is used and nsteps from yaml ignored."
        end
        nsteps = Int64(round((final_time - initial_time) / fixed_dt) + 1)
        dt = (final_time - initial_time) / (nsteps - 1)
    end
    comm = datamanager.get_comm()
    # not needed here
    numerical_damping = get_numerical_damping(params)
    max_damage = get_max_damage(params)
    solver_specifics = Dict(
        "Solution tolerance" => 1e-7,
        "Residual tolerance" => 1e-7,
        "Maximum number of iterations" => 1,
        "Show solver iteration" => false,
    )
    if haskey(params["Static"], "Solution tolerance")
        solver_specifics["Solution tolerance"] = params["Static"]["Solution tolerance"]
    end
    if haskey(params["Static"], "Residual tolerance")
        solver_specifics["Residual tolerance"] = params["Static"]["Residual tolerance"]
    end
    if haskey(params["Static"], "Maximum number of iterations")
        solver_specifics["Maximum number of iterations"] =
            params["Static"]["Maximum number of iterations"]
    end
    if haskey(params["Static"], "Show solver iteration")
        solver_specifics["Show solver iteration"] =
            params["Static"]["Show solver iteration"]
    end
    dof = datamanager.get_dof()
    residual = datamanager.create_constant_node_field("Residual", Float64, dof)
    start_u = datamanager.create_constant_node_field("Start Value", Float64, dof)
    Boundary_conditions.find_bc_free_dof(datamanager, bcs)


    if datamanager.get_rank() == 0
        data = [
            "Static Solver" "" "" ""
            "Initial time" "."^10 initial_time "s"
            "Final time" "."^10 final_time "s"
        ]
        if fixed_dt != -1.0
            data = vcat(
                data,
                [
                    "Fixed time increment" "."^10 fixed_dt "s";
                ],
            )
        else
            data = vcat(
                data,
                [
                    "Time increment" "."^10 dt "s"
                ],
            )
        end
        data = vcat(
            data,
            [
                "Number of steps" "."^10 nsteps ""
            ],
        )

        print_table(data, datamanager)
    end

    return initial_time, dt, nsteps, numerical_damping, max_damage, solver_specifics
end

function run_solver(
    solver_options::Dict{Any,Any},
    block_nodes::Dict{Int64,Vector{Int64}},
    bcs::Dict{Any,Any},
    datamanager::Module,
    outputs::Dict{Int64,Dict{}},
    result_files::Vector{Dict},
    synchronise_field,
    write_results,
    to::TimerOutputs.TimerOutput,
    silent::Bool,
)

    atexit(() -> datamanager.set_cancel(true))

    #show_trace = solver_options["Show iterations"]
    show_trace = false

    xtol = 1e-6  # Genauigkeit der Lösung
    ftol = 1e-6  # Genauigkeit der Funktion
    iterations = 50  # Max. Anzahl an Iterationen

    coor = datamanager.get_field("Coordinates")
    comm = datamanager.get_comm()

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
    bc_free_dof = datamanager.get_bc_free_dof()

    xtol = solver_options["Solver specifics"]["Solution tolerance"]
    ftol = solver_options["Solver specifics"]["Residual tolerance"]
    iterations = solver_options["Solver specifics"]["Maximum number of iterations"]
    show_trace = solver_options["Solver specifics"]["Show solver iteration"]
    uN = datamanager.get_field("Displacements", "N")
    resdiual = datamanager.get_field("Residual")
    start_u = datamanager.get_field("Start Value")
    coor = datamanager.get_field("Coordinates")
    external_forces = datamanager.get_field("External Forces")
    external_force_densities = datamanager.get_field("External Force Densities")
    force_densities = datamanager.get_field("Force Densities", "NP1")
    forces = datamanager.get_field("Forces", "NP1")
    volume = datamanager.get_field("Volume")
    active_nodes = datamanager.get_field("Active Nodes")
    active_list = datamanager.get_field("Active")
    for idt in iter

        datamanager.set_iteration(idt)
        #datamanager = apply_bc_dirichlet(bcs, datamanager, step_time) #-> Dirichlet
        datamanager = apply_bc_dirichlet_force(bcs, datamanager, step_time) #-> Dirichlet
        external_force_densities += external_forces ./ volume
        start_u = copy(uN) .+ 0.001
        sol = nlsolve(
            (residual, U) -> residual!(
                residual,
                U,
                datamanager,
                bc_free_dof,
                bcs,
                block_nodes,
                dt,
                step_time,
                solver_options,
                synchronise_field,
                to,
            ),
            start_u;
            xtol = xtol,
            ftol = ftol,
            iterations = iterations,
            show_trace = show_trace,
        )
        active_nodes =
            find_active_nodes(active_list, active_nodes, 1:datamanager.get_nnodes())
        @views forces[active_nodes, :] =
            force_densities[active_nodes, :] .* volume[active_nodes]

        #datamanager = Boundary_conditions.apply_bc_dirichlet(bcs, datamanager, step_time) #-> Dirichlet

        @timeit to "write_results" result_files = write_results(
            result_files,
            start_time + step_time,
            max_damage,
            outputs,
            datamanager,
        )

        if rank == 0 && !silent && datamanager.get_cancel()
            set_multiline_postfix(iter, "Simulation canceled!")
            break
        end
        @timeit to "switch_NP1_to_N" datamanager.switch_NP1_to_N()

        step_time += dt
        if idt % ceil(nsteps / 100) == 0
            @info "Step: $idt / $(nsteps+1) [$step_time s]"
        end
        # if rank == 0 && !silent
        #     set_postfix(iter, t=@sprintf("%.4e", step_time))
        # end
        MPI.Barrier(comm)


    end
    return result_files
end

function residual!(
    residual,
    U,
    datamanager,
    bc_free_dof,
    bcs,
    block_nodes,
    dt,
    step_time,
    solver_options,
    synchronise_field,
    to,
)

    active_list = datamanager.get_field("Active")

    uNP1 = datamanager.get_field("Displacements", "NP1")

    deformed_coorNP1 = datamanager.get_field("Deformed Coordinates", "NP1")

    force_densities = datamanager.get_field("Force Densities", "NP1")
    coor = datamanager.get_field("Coordinates")
    active_nodes = datamanager.get_field("Active Nodes")
    active_nodes = find_active_nodes(active_list, active_nodes, 1:datamanager.get_nnodes())
    external_force_densities = datamanager.get_field("External Force Densities")
    # one step more, because of init step (time = 0)
    # displacement  from interface


    datamanager = apply_bc_dirichlet(bcs, datamanager, step_time) #-> Dirichlet


    foreach(t -> uNP1[t...] = U[t...], bc_free_dof)
    #println(uNP1[5:8, :], U[5:8, :])

    @views deformed_coorNP1[active_nodes, :] =
        coor[active_nodes, :] .+ uNP1[active_nodes, :]

    force_densities .= 0 # TODO check where to put it for iterative solver

    datamanager.synch_manager(synchronise_field, "upload_to_cores")
    # synch
    datamanager = Model_Factory.compute_models(
        datamanager,
        block_nodes,
        dt,
        step_time,
        solver_options["Models"],
        synchronise_field,
        to,
    )


    datamanager.synch_manager(synchronise_field, "download_from_cores")
    # synch

    # @timeit to "apply_bc_neumann" datamanager = Boundary_conditions.apply_bc_neumann(bcs, datamanager, step_time) #-> von neumann
    #active_nodes = datamanager.get_field("Active Nodes")
    #active_nodes =            find_active_nodes(active_list, active_nodes, 1:datamanager.get_nnodes())

    # check_inf_or_nan(force_densities, "Forces")

    # Richtung muss mit rein
    residual .= 0
    foreach(
        t -> residual[t...] = force_densities[t...] + external_force_densities[t...],
        bc_free_dof,
    )


end


end
