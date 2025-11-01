# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Matrix_verlet
using LinearAlgebra
using TimerOutputs
using ProgressBars: set_multiline_postfix, set_postfix
using Printf
using LoopVectorization
using PrettyTables
using Logging

using ...Helpers: check_inf_or_nan, find_active_nodes, progress_bar, matrix_style
using ...Parameter_Handling:
                             get_initial_time,
                             get_fixed_dt,
                             get_final_time,
                             get_numerical_damping,
                             get_safety_factor,
                             get_max_damage
using ...MPI_Communication: find_and_set_core_value_min, find_and_set_core_value_max,
                            barrier
using ..Model_Factory: compute_models
using ..Boundary_Conditions: apply_bc_dirichlet, apply_bc_neumann
using ...Logging_Module: print_table
include("../../Compute/compute_field_values.jl")

function init_solver(solver_options::Dict{Any,Any},
                     params::Dict,
                     bcs::Dict{Any,Any},
                     datamanager::Module,
                     block_nodes::Dict{Int64,Vector{Int64}})
    find_bc_free_dof(datamanager, bcs)
    solver_options["Initial Time"] = get_initial_time(params, datamanager)
    solver_options["Final Time"] = get_final_time(params, datamanager)

    solver_options["Number of Steps"] = get_nsteps(params)
    solver_options["dt"] = (solver_options["Final Time"] - solver_options["Initial Time"]) /
                           solver_options["Number of Steps"]

    for (block, nodes) in pairs(block_nodes)
        model_param = datamanager.get_properties(block, "Material Model")
        Correspondence_matrix_based.init_model(datamanager, nodes, model_param)
    end
    @timeit to "init_matrix" Correspondence_matrix_based.init_matrix(datamanager)

    density = datamanager.get_field("Density")

    K = datamanager.get_stiffness_matrix()
    # create M^-1*K for linear run
    for (block, nodes) in pairs(block_nodes)
        for node in nodes
            for idof in 1:dof
                K[(node-1)*dof+dof, (node-1)*dof+dof]/=density[node]
            end
        end
    end
    datamanager.set_stiffness_matrix(K)
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
    dt::Float64 = solver_options["dt"]
    nsteps::Int64 = solver_options["Number of Steps"]
    time::Float64 = solver_options["Initial Time"]
    step_time::Float64 = 0
    max_cancel_damage::Float64 = solver_options["Maximum Damage"]
    numerical_damping::Float64 = solver_options["Numerical Damping"]
    max_damage::Float64 = 0
    damage_init::Bool = false
    uN = datamanager.get_field("Displacements", "N")
    uNP1 = datamanager.get_field("Displacements", "NP1")
    deformed_coorNP1 = datamanager.get_field("Deformed Coordinates", "NP1")
    forces = datamanager.get_field("Forces", "NP1")
    force_densities_N = datamanager.get_field("Force Densities", "N")
    force_densities_NP1 = datamanager.get_field("Force Densities", "NP1")
    volume = datamanager.get_field("Volume")
    forces = datamanager.get_field("Forces", "NP1")
    vN = datamanager.get_field("Velocity", "N")
    vNP1 = datamanager.get_field("Velocity", "NP1")
    external_force_densities = datamanager.get_field("External Force Densities")
    a = datamanager.get_field("Acceleration")
    perm = create_permutation(datamanager.get_nnodes(), datamanager.get_dof())

    rank = datamanager.get_rank()
    iter = progress_bar(rank, nsteps, silent)
    K = datamanager.get_stiffness_matrix()

    #nodes::Vector{Int64} = Vector{Int64}(1:datamanager.get_nnodes())
    @inbounds @fastmath for idt in iter
        datamanager = apply_bc_dirichlet(["Displacements", "Forces", "Force Densities"],
                                         bcs,
                                         datamanager, time,
                                         step_time)

        vNP1[active_nodes, :] = (1 - numerical_damping) .*
                                vN[active_nodes, :] .+
                                0.5 * dt .* a[active_nodes, :]
        uNP1[active_nodes, :] = uN[active_nodes, :] .+
                                dt .* vNP1[active_nodes, :]
        @views deformed_coorNP1 .= coor .+ uNP1

        # thermal model, usw.

        mul!(a[active_nodes], K[active_nodes, active_nodes], uNP1[active_nodes])

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
        #a+= + fexternal / density
    end
end

function te(a, K, u)
    a.=K*u
end

end
