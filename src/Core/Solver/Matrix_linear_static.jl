# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Matrix_linear_static
using SparseArrays
using LinearAlgebra

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
                K = datamanager.get_stiffness_matrix()
            end

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

            datamanager = apply_bc_dirichlet(["Displacements", "Temperature"],
                                             bcs,
                                             datamanager, time,
                                             step_time) #-> Dirichlet
            #needed because of optional deformation_gradient, Deformed bonds, etc.
            # all points to guarantee that the neighbors have coor as coordinates if they are not active

            #@timeit to "upload_to_cores" datamanager.synch_manager(synchronise_field,
            #                                                       "upload_to_cores")
            # synch
            F_external += K*U_extrnal
            # non zeros u weg von K; F_exteral an den stellen = u_external

            uNP1 = K/F
            if "Material" in solver_options["Models"]
                @views deformed_coorNP1[active_nodes,
                                        :] = coor[active_nodes, :] .+
                                             uNP1[active_nodes, :]
            end
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

end
