# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Matrix_Verlet
using LinearAlgebra
using TimerOutputs
using ProgressBars: set_multiline_postfix, set_postfix
using Printf
using LoopVectorization
using PrettyTables
using Logging

using ...Data_Manager
using ...Helpers: check_inf_or_nan, find_active_nodes, progress_bar, matrix_style,
                  create_permutation
using ...Parameter_Handling:
                             get_initial_time,
                             get_fixed_dt,
                             get_final_time,
                             get_numerical_damping,
                             get_safety_factor,
                             get_max_damage,
                             get_nsteps
using ...MPI_Communication: find_and_set_core_value_min, find_and_set_core_value_max,
                            barrier
using ..Model_Factory: compute_models
using ..Boundary_Conditions: apply_bc_dirichlet, apply_bc_neumann, find_bc_free_dof
using ...Logging_Module: print_table
include("../Model_reduction/Model_reduction.jl")
using .Model_reduction: guyan_reduction
include("../../Compute/compute_field_values.jl")
include("../../Models/Material/Material_Models/Correspondence/Correspondence_matrix_based.jl")
using .Correspondence_matrix_based
export init_solver, run_solver
function init_solver(solver_options::Dict{Any,Any},
                     params::Dict,
                     bcs::Dict{Any,Any},
                     block_nodes::Dict{Int64,Vector{Int64}})
    horizon = Data_Manager.get_field("Horizon")
    if Data_Manager.get_rank()>1
        @warn "implementation might not work for MPI. Especially for coupling. It has to be tested."
    end
    find_bc_free_dof(bcs)
    solver_options["Initial Time"] = get_initial_time(params)
    solver_options["Final Time"] = get_final_time(params)

    solver_options["Number of Steps"] = get_nsteps(params)
    solver_options["dt"] = (solver_options["Final Time"] - solver_options["Initial Time"]) /
                           solver_options["Number of Steps"]

    for (block, nodes) in pairs(block_nodes)
        model_param = Data_Manager.get_properties(block, "Material Model")
        Correspondence_matrix_based.init_model(nodes, model_param)
    end
    @timeit "init_matrix" Correspondence_matrix_based.init_matrix()
    dof = Data_Manager.get_dof()
    density = Data_Manager.get_field("Density")
    model_reduction = get(params["Verlet Matrix Based"], "Model Reduction", false)
    solver_options["Model Reduction"] = model_reduction
    solver_options["Numerical Damping"] = get_numerical_damping(params)
    K = Data_Manager.get_stiffness_matrix()
    reduction_blocks=[]
    if model_reduction!=false
        reduction_blocks = get(solver_options["Model Reduction"], "Reduction Blocks", [])

        if reduction_blocks==[]
            model_reduction==false
            @error "No reduction blocks defined for model reduction. If you want to use a reduced model please define 'Reduction Blocks' in the yaml input deck."
        else
            if reduction_blocks isa Float64
                @error "Type Float is not supported for Reduction Blocks"
            end
            if reduction_blocks isa Int64
                reduction_blocks = [reduction_blocks]
            else
                reduction_blocks = parse.(Int64,
                                          filter(!isempty,
                                                 split(solver_options["Model Reduction"]["Reduction Blocks"],
                                                       r"[,\s]")))
            end
        end
    end

    if reduction_blocks!=[]
        coor = Data_Manager.get_field("Coordinates")
        master_nodes = Int64[]
        slave_nodes = Int64[]
        pd_nodes = Int64[]
        nlist = Data_Manager.get_nlist()
        cn=Data_Manager.create_constant_node_scalar_field("Coupling Nodes", Int64)
        full_blocks = setdiff(collect(keys(block_nodes)), reduction_blocks)
        for block in full_blocks
            append!(pd_nodes, block_nodes[block])
        end

        # create coupling zone
        for node in pd_nodes
            append!(master_nodes, nlist[node])
        end

        append!(master_nodes, pd_nodes)
        master_nodes = sort(unique(master_nodes))
        slave_nodes = setdiff(collect(1:Data_Manager.get_nnodes()), master_nodes)

        cn[slave_nodes].=6
        cn[master_nodes].=2
        cn[pd_nodes].=1

        for block in eachindex(block_nodes)
            block_nodes[block] = intersect(block_nodes[block], pd_nodes)
        end

        if !(master_nodes==[])
            perm_master = create_permutation(master_nodes, Data_Manager.get_dof())
            perm_slave = create_permutation(slave_nodes, Data_Manager.get_dof())
            perm_pd_nodes = create_permutation(pd_nodes, Data_Manager.get_dof())
            # create reduced M^-1*K for linear run
            #TODO adapt for mass matrix
            density_mass = zeros(length(density)*Data_Manager.get_dof())
            density_mass .= density[1]
            perm_pd_reduced = findall(in(perm_pd_nodes), perm_master)

            K_reduced,
            mass_reduced = guyan_reduction(K, density_mass, perm_master, perm_slave, dof)
            K_reduced[perm_pd_reduced, :].=0
            #K_reduced[:, perm_pd_reduced].=0
            Data_Manager.set_stiffness_matrix(K_reduced)
            Data_Manager.set_mass_matrix(mass_reduced)
            Data_Manager.set_reduced_model_pd(pd_nodes)
            Data_Manager.set_reduced_model_master(master_nodes)
            @info "Model reduction is applied \n
             Slaves: $(length(slave_nodes)), Master: $(length(master_nodes)), PD: $(length(pd_nodes))."

            return
        end
        @warn "No master nodes defined for model reduction. Using full stiffness matrix."
    end

    # create M^-1*K for linear run
    #for (block, nodes) in pairs(block_nodes)
    #	for node in nodes
    #		for idof in 1:dof
    #			K[(node-1)*dof+idof, (node-1)*dof+idof]/=density[node]
    #		end
    #	end
    #end
    Data_Manager.set_stiffness_matrix(K)
    return

    # reduced matrix

end

function run_solver(solver_options::Dict{Any,Any},
                    block_nodes::Dict{Int64,Vector{Int64}},
                    bcs::Dict{Any,Any},
                    outputs::Dict{Int64,Dict{}},
                    result_files::Vector{Dict},
                    synchronise_field,
                    write_results,
                    compute_parabolic_problems_before_model_evaluation,
                    compute_parabolic_problems_after_model_evaluation,
                    silent::Bool)
    dt::Float64 = solver_options["dt"]
    nsteps::Int64 = solver_options["Number of Steps"]
    time::Float64 = solver_options["Initial Time"]
    step_time::Float64 = 0
    #max_cancel_damage::Float64 = solver_options["Maximum Damage"]
    numerical_damping::Float64 = solver_options["Numerical Damping"]
    max_damage::Float64 = 0
    damage_init::Bool = false

    density = Data_Manager.get_field("Density")

    coor = Data_Manager.get_field("Coordinates")

    active_list = Data_Manager.get_field("Active")
    volume = Data_Manager.get_field("Volume")

    if "Material" in solver_options["Models"]
        external_forces = Data_Manager.get_field("External Forces")
        external_force_densities = Data_Manager.get_field("External Force Densities")
        a = Data_Manager.get_field("Acceleration")
    end

    comm = Data_Manager.get_comm()
    rank = Data_Manager.get_rank()
    iter = progress_bar(rank, nsteps, silent)
    if solver_options["Model Reduction"]!=false
        master_nodes=Data_Manager.get_reduced_model_master()
        perm = collect(1:(length(master_nodes) * Data_Manager.get_dof()))
        mass = Data_Manager.get_mass_matrix()
    else
        perm = create_permutation(Data_Manager.get_nnodes(), Data_Manager.get_dof())
    end
    K = Data_Manager.get_stiffness_matrix()

    #nodes::Vector{Int64} = Vector{Int64}(1:Data_Manager.get_nnodes())
    @inbounds @fastmath for idt in iter
        uN = Data_Manager.get_field("Displacements", "N")
        uNP1 = Data_Manager.get_field("Displacements", "NP1")
        forces = Data_Manager.get_field("Forces", "NP1")
        deformed_coorNP1 = Data_Manager.get_field("Deformed Coordinates", "NP1")
        vN = Data_Manager.get_field("Velocity", "N")
        vNP1 = Data_Manager.get_field("Velocity", "NP1")
        force_densities_NP1 = Data_Manager.get_field("Force Densities", "NP1")
        active_nodes = Data_Manager.get_field("Active Nodes")
        active_nodes = find_active_nodes(active_list, active_nodes,
                                         1:Data_Manager.get_nnodes())
        apply_bc_dirichlet(["Displacements", "Forces", "Force Densities"],
                           bcs, time,
                           step_time)

        if solver_options["Model Reduction"]!=false
            active_nodes = master_nodes
            active_list.=false
            active_list[Data_Manager.get_reduced_model_pd()].=true
        end

        vNP1[active_nodes, :] = (1 - numerical_damping) .*
                                vN[active_nodes, :] .+
                                0.5 * dt .* a[active_nodes, :]

        uNP1[active_nodes, :] = uN[active_nodes, :] .+
                                dt .* vNP1[active_nodes, :]

        apply_bc_dirichlet(["Displacements", "Temperature"],
                           bcs,
                           time,
                           step_time) #-> Dirichlet
        @views deformed_coorNP1 .= coor .+ uNP1
        # thermal model, usw.

        #mul!(a[active_nodes], K[active_nodes, active_nodes], uNP1[active_nodes])

        sa = size(a[active_nodes, :])

        #mul!(a, K, vec(uNP1))
        #force_densities_NP1=reshape(K[perm, perm]*vec(uNP1), sa...)

        compute_models(block_nodes,
                       dt,
                       time,
                       solver_options["Models"],
                       synchronise_field)

        force_densities_NP1[active_nodes, :] -= reshape(K[perm,
                                                          perm]*vec(uNP1[active_nodes, :]),
                                                        sa...)
        .+ external_force_densities[active_nodes, :]
        .+ external_forces[active_nodes, :] ./ volume[active_nodes]

        # @timeit "download_from_cores" Data_Manager.synch_manager(synchronise_field,
        #                                                            "download_from_cores")
        #println(a)
        #@views a[active_nodes, :] = external_force_densities[active_nodes, :] .+
        #							external_forces[active_nodes, :] ./
        #							volume[active_nodes]
        if solver_options["Model Reduction"]!=false
            a[active_nodes, :] = reshape(vec(force_densities_NP1[active_nodes, :]) ./ mass,
                                         sa...)
        else
            a[active_nodes, :] = force_densities_NP1[active_nodes, :] ./
                                 (density[active_nodes])
        end

        @timeit "write_results" result_files=write_results(result_files, time,
                                                           max_damage, outputs)
        # for file in result_files
        #     flush(file)
        # end
        if rank == 0 && !silent && Data_Manager.get_cancel()
            set_multiline_postfix(iter, "Simulation canceled!")
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
        #a+= + fexternal / density
    end
end

function te(a, K, u)
    a.=K*u
end

end
