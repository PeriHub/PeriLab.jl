# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Static_Solver
using LinearAlgebra
using NLsolve
using TimerOutputs
using ProgressBars: set_multiline_postfix, set_postfix
using Printf
using LoopVectorization
using Logging

using ...Helpers: check_inf_or_nan, find_active_nodes, progress_bar, matrix_style
using ...Parameter_Handling:
                             get_initial_time,
                             get_fixed_dt,
                             get_final_time,
                             get_numerical_damping,
                             get_safety_factor,
                             get_max_damage
using ...MPI_Communication: barrier
using ..Model_Factory
using ..Boundary_Conditions: apply_bc_dirichlet, apply_bc_neumann, find_bc_free_dof
using ...Logging_Module: print_table
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
"""
function init_solver(solver_options::Dict{Any,Any},
                     params::Dict,
                     bcs::Dict{Any,Any},
                     datamanager::Module,
                     block_nodes::Dict{Int64,Vector{Int64}})
    # @info "==============================="
    # @info "==== NLsolve Static Solver ===="
    # @info "==============================="

    mechanical = "Material" in solver_options["Models"]
    thermal = "Thermal" in solver_options["Models"]
    initial_time = get_initial_time(params, datamanager)
    final_time = get_final_time(params, datamanager)
    fixed_dt = get_fixed_dt(params)
    dof = datamanager.get_dof()
    if fixed_dt == -1.0
        if haskey(params, "Number of Steps")
            nsteps = params["Number of Steps"]
            dt = (final_time - initial_time) / nsteps
        else
            nsteps = Int64(1)
            dt = final_time - initial_time
        end
    else
        if haskey(params, "Number of Steps")
            @warn "''Number of Steps'' and ''Fixed dt'' are defined. ''Fixed dt'' is used and ''Number of Steps'' from yaml is ignored."
        end
        nsteps = Int64(round((final_time - initial_time) / fixed_dt))
        dt = fixed_dt
    end
    comm = datamanager.get_comm()

    # not needed here
    numerical_damping = get_numerical_damping(params)
    max_damage = get_max_damage(params)
    solver_specifics = Dict("Solution tolerance" => 1e-7,
                            "Residual tolerance" => 1e-7,
                            "Maximum number of iterations" => 100,
                            "Show solver iteration" => false,
                            "Residual scaling" => 1e6,
                            "m" => 15,
                            "Linear Start Value" => zeros(2 * dof))
    if haskey(params["Static"], "Residual scaling")
        volume = datamanager.get_field("Volume")
        solver_specifics["Residual scaling"] = params["Static"]["Residual scaling"]# / minimum(volume) / minimum(volume)
    end
    if haskey(params["Static"], "Solution tolerance")
        solver_specifics["Solution tolerance"] = params["Static"]["Solution tolerance"]
    end
    if haskey(params["Static"], "Residual tolerance")
        solver_specifics["Residual tolerance"] = params["Static"]["Residual tolerance"]
    end
    if haskey(params["Static"], "Maximum number of iterations")
        solver_specifics["Maximum number of iterations"] = params["Static"]["Maximum number of iterations"]
    end
    if haskey(params["Static"], "Show solver iteration")
        solver_specifics["Show solver iteration"] = params["Static"]["Show solver iteration"]
    end
    if haskey(params["Static"], "m")
        solver_specifics["m"] = params["Static"]["m"]
    end
    if haskey(params["Static"], "Linear Start Value")
        solver_specifics["Linear Start Value"] = parse.(Float64,
                                                        split(params["Static"]["Linear Start Value"]))
    end

    residual = datamanager.create_constant_node_field("Residual", Float64, dof)

    if !("Start_Values" in datamanager.get_all_field_keys())
        start_u = datamanager.create_constant_node_field("Start_Values", Float64, dof)
        coor = datamanager.get_field("Coordinates")
        ls = solver_specifics["Linear Start Value"]
        if ls[1] == ls[2]
            start_u .= ls[1]
        else
            for idof in 1:dof
                m = (ls[2 * idof] - ls[2 * idof - 1]) /
                    (maximum(coor[:, idof]) - minimum(coor[:, idof]))
                n = ls[2 * idof] - m * maximum(coor[:, idof])
                start_u[:, idof] = (m .* coor[:, idof] .+ n) ./ nsteps .* final_time
                n = ls[2 * idof] - m * maximum(coor[:, idof])
                start_u[:, idof] = (m .* coor[:, idof] .+ n) ./ nsteps
            end
        end
    else
        start_u = datamanager.get_field("Start_Values")
    end
    check_inf_or_nan(start_u, "Start_Values")
    find_bc_free_dof(datamanager, bcs)

    if datamanager.get_rank() == 0
        data = ["Static Solver" "" "" ""
                "Initial time" "."^10 initial_time "s"
                "Final time" "."^10 final_time "s"
                "Maximum number of iterations" "."^10 solver_specifics["Maximum number of iterations"] ""
                "Solution tolerance" "."^10 solver_specifics["Solution tolerance"] ""
                "Residual tolerance" "."^10 solver_specifics["Residual tolerance"] ""]
        if fixed_dt != -1.0
            data = vcat(data,
                        ["Fixed time increment" "."^10 fixed_dt "s";])
        else
            data = vcat(data,
                        ["Time increment" "."^10 dt "s"])
        end
        data = vcat(data,
                    ["Number of steps" "."^10 nsteps ""])

        print_table(data, datamanager)
    end

    solver_options["Initial Time"] = initial_time
    solver_options["Final Time"] = final_time
    solver_options["dt"] = dt
    solver_options["Number of Steps"] = nsteps
    solver_options["Numerical Damping"] = numerical_damping
    solver_options["Maximum Damage"] = max_damage
    solver_options["Solver specifics"] = solver_specifics
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

    coor = datamanager.get_field("Coordinates")
    comm = datamanager.get_comm()

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
    bc_free_dof = datamanager.get_bc_free_dof()
    scaling = solver_options["Solver specifics"]["Residual scaling"]
    xtol = solver_options["Solver specifics"]["Solution tolerance"]
    ftol = solver_options["Solver specifics"]["Residual tolerance"]
    iterations = solver_options["Solver specifics"]["Maximum number of iterations"]
    show_trace = solver_options["Solver specifics"]["Show solver iteration"]
    m = solver_options["Solver specifics"]["m"]
    uN = datamanager.get_field("Displacements", "N")
    uNP1 = datamanager.get_field("Displacements", "NP1")
    resdiual = datamanager.get_field("Residual")
    start_u = datamanager.get_field("Start_Values")
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
        #datamanager = apply_bc_dirichlet(bcs, datamanager, time) #-> Dirichlet
        datamanager = apply_bc_dirichlet(["Forces", "Force Densities"],
                                         bcs,
                                         datamanager,
                                         time,
                                         step_time) #-> Dirichlet
        external_force_densities += external_forces ./ volume
        active_nodes = datamanager.get_field("Active Nodes")
        active_nodes = find_active_nodes(active_list, active_nodes,
                                         1:datamanager.get_nnodes())

        compute_parabolic_problems_before_model_evaluation(active_nodes, datamanager,
                                                           solver_options)
        datamanager = apply_bc_dirichlet(["Displacements", "Temperature"], bcs, datamanager,
                                         time,
                                         step_time + dt) #-> Dirichlet

        sol = nlsolve((residual,
                       U) -> residual!(residual,
                                       U,
                                       datamanager,
                                       bc_free_dof,
                                       block_nodes,
                                       dt,
                                       time,
                                       solver_options,
                                       synchronise_field,
                                       to,
                                       scaling),
                      start_u;
                      xtol = xtol,
                      ftol = ftol,
                      iterations = iterations,
                      show_trace = show_trace & !silent,
                      extended_trace = false,
                      method = :anderson,
                      m = m,)

        compute_parabolic_problems_after_model_evaluation(active_nodes, datamanager,
                                                          solver_options, dt)

        start_u = copy(uNP1)
        if !sol.x_converged && !sol.f_converged
            @info "Failed to converge at step $idt: maximum number of iterations reached"
            datamanager.set_cancel(true)
        end
        # method=:newton, linsolve=KrylovJL_GMRES()
        # method=:broyden
        active_nodes = find_active_nodes(active_list, active_nodes,
                                         1:datamanager.get_nnodes())

        if "Damage" in solver_options["Models"]
            damage = datamanager.get_damage("NP1")
            max_damage = maximum(damage[active_nodes])
            if max_damage > max_cancel_damage
                @info "Maximum damage reached at step $idt: $max_damage"
                datamanager.set_cancel(true)
            end
            if !damage_init && max_damage > 0
                damage_init = true
            end
        end

        if datamanager.get_cancel()
            @info "Canceling at step $idt"
            break
        end

        force_densities = datamanager.get_field("Force Densities", "NP1")

        @views forces[active_nodes,
        :] = force_densities[active_nodes, :] .*
                                         volume[active_nodes]

        @timeit to "write_results" result_files=write_results(result_files, time,
                                                              max_damage, outputs,
                                                              datamanager)

        @timeit to "switch_NP1_to_N" datamanager.switch_NP1_to_N()

        time += dt
        step_time += dt
        datamanager.set_current_time(time)
        if idt % ceil(nsteps / 100) == 0
            @info "Step: $idt / $(nsteps+1) [$time s]"
        end
        # if rank == 0 && !silent
        #     set_postfix(iter, t=@sprintf("%.4e", time))
        # end
        barrier(comm)
    end
    return result_files
end

function residual!(residual,
                   U,
                   datamanager,
                   bc_free_dof,
                   block_nodes,
                   dt,
                   time,
                   solver_options,
                   synchronise_field,
                   to,
                   scaling)
    active_list = datamanager.get_field("Active")

    deformed_coorNP1 = datamanager.get_field("Deformed Coordinates", "NP1")

    force_densities = datamanager.get_field("Force Densities", "NP1")
    coor = datamanager.get_field("Coordinates")
    active_nodes = datamanager.get_field("Active Nodes")
    active_nodes = find_active_nodes(active_list, active_nodes, 1:datamanager.get_nnodes())
    external_force_densities = datamanager.get_field("External Force Densities")
    # one step more, because of init step (time = 0)
    # displacement  from interface
    # Set bond damages to original state
    # datamanager.switch_bonds!(
    #     datamanager.get_field("Bond Damage", "NP1"),
    #     datamanager.get_field("Bond Damage", "N"),
    # )

    uNP1 = datamanager.get_field("Displacements", "NP1")
    uNP1[bc_free_dof] .= U[bc_free_dof]

    # bc_dof = setdiff(1:length(uNP1), bc_free_dof)

    @views deformed_coorNP1[active_nodes,
    :] = coor[active_nodes, :] .+
                                               uNP1[active_nodes, :]

    force_densities[:, :] .= 0 # TODO check where to put it for iterative solver
    forces = datamanager.get_field("Forces", "NP1")
    forces[:, :] .= 0

    datamanager.synch_manager(synchronise_field, "upload_to_cores")
    # synch
    datamanager = Model_Factory.compute_models(datamanager,
                                               block_nodes,
                                               dt,
                                               time,
                                               solver_options["Models"],
                                               synchronise_field,
                                               to)

    datamanager.synch_manager(synchronise_field, "download_from_cores")
    # synch

    # @timeit to "apply_bc_neumann" datamanager = apply_bc_neumann(bcs, datamanager, time) #-> von neumann
    #active_nodes = datamanager.get_field("Active Nodes")
    #active_nodes =            find_active_nodes(active_list, active_nodes, 1:datamanager.get_nnodes())

    # check_inf_or_nan(force_densities, "Forces")

    # Richtung muss mit rein
    #println("adsads", force_densities[bc_free_dof] + external_force_densities[bc_free_dof])
    residual .= 0
    residual[bc_free_dof] = (force_densities[bc_free_dof] +
                             external_force_densities[bc_free_dof]) ./ scaling
end

end
