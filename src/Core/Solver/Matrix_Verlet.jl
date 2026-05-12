# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Matrix_Verlet
using LinearAlgebra
using SparseArrays
using TimerOutputs
using ProgressBars: set_multiline_postfix, set_postfix
using Printf
using LoopVectorization
using PrettyTables
using Logging
using LinearAlgebra: lu
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
using ...MPI_Communication: find_and_set_core_value_min, find_and_set_core_value_max
using ..Model_Factory: compute_models, compute_crititical_time_step
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
    if Data_Manager.get_rank() > 1
        @warn "Implementation might not work for MPI. Especially for coupling. It has to be tested."
    end
    find_bc_free_dof(bcs)
    initial_time = get_initial_time(params)
    final_time = get_final_time(params)
    nsteps = get_nsteps(params)

    mechanical = "Material" in solver_options["Models"]
    thermal = "Thermal" in solver_options["Models"]

    solver_options["dt"] = compute_crititical_time_step(block_nodes, mechanical, thermal)
    #solver_options["dt"] = 1e-8
    safety_factor = get_safety_factor(params)
    fixed_dt = get_fixed_dt(params)
    min_dt = compute_crititical_time_step(block_nodes, mechanical, thermal)
    if fixed_dt == -1.0
        nsteps, dt = get_integration_steps(initial_time, final_time, safety_factor * min_dt)
    else
        nsteps, dt = get_integration_steps(initial_time, final_time, fixed_dt)
    end
    K = Data_Manager.get_stiffness_matrix()
    if isnothing(K)
        for (block, nodes) in pairs(block_nodes)
            model_param = Data_Manager.get_properties(block, "Material Model")
            Correspondence_matrix_based.init_model(nodes, model_param, block)
        end
        @timeit "init_matrix" Correspondence_matrix_based.init_matrix()
        K = Data_Manager.get_stiffness_matrix()
    end
    dof = Data_Manager.get_dof()
    density = Data_Manager.get_field("Density")
    model_reduction = get(params["Verlet Matrix Based"], "Model Reduction", false)
    solver_options["Model Reduction"] = model_reduction
    solver_options["Numerical Damping"] = get_numerical_damping(params)
    solver_options["dt"] = dt
    solver_options["Initial Time"] = initial_time
    solver_options["Final Time"] = final_time
    solver_options["Number of Steps"] = nsteps
    @info "Solver parameter"
    @info "Initial Time: $(solver_options["Initial Time"])"
    @info "Final Time: $(solver_options["Final Time"])"
    @info "Number of Steps: $(solver_options["Number of Steps"])"
    @info "dt: $(solver_options["dt"])"
    @info "Numerical Damping: $(solver_options["Numerical Damping"])"

    reduction_blocks = []
    if model_reduction != false
        reduction_blocks = get(solver_options["Model Reduction"], "Reduction Blocks", [])

        if reduction_blocks == []
            model_reduction == false
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
    else
        # might lead to issues, because neighbors outside the block are included in
        for block in eachindex(block_nodes)
            block_nodes[block] = []
        end
        @warn "No other models work with full matrix Verlet right now."
    end

    if reduction_blocks != []
        master_nodes = Int64[]
        slave_nodes = Int64[]
        pd_nodes = Int64[]

        nlist = Data_Manager.get_nlist()
        cn = Data_Manager.create_constant_node_scalar_field("Coupling Nodes", Int64)
        full_blocks = setdiff(collect(keys(block_nodes)), reduction_blocks)
        for block in full_blocks
            append!(pd_nodes, block_nodes[block])
        end

        for node in pd_nodes
            append!(master_nodes, nlist[node])
        end

        append!(master_nodes, pd_nodes)

        master_nodes = sort(unique(master_nodes))

        slave_nodes = setdiff(1:Data_Manager.get_nnodes(), master_nodes)

        if !(get(solver_options["Model Reduction"], "Material Point Region", true))
            pd_nodes::Vector{Int64} = []
        end
        cn[slave_nodes] .= 6
        cn[master_nodes] .= 2
        cn[pd_nodes] .= 1

        if !(master_nodes == [])
            nnodes = Data_Manager.get_nnodes()
            master_nodes = Vector{Int64}(1:nnodes)
            # complete region of analysis with coupling
            perm_master = create_permutation(master_nodes, Data_Manager.get_dof(), nnodes)
            # complete region of reduction
            perm_slave = create_permutation(slave_nodes, Data_Manager.get_dof(), nnodes)

            # indices of pd nodes in global system of PD nodes for crack propagation and BCs application
            perm_pd_nodes = create_permutation(pd_nodes, Data_Manager.get_dof(), nnodes)
            # indices of deletes in global system
            # create reduced M^-1*K for linear run
            density_mass = zeros(length(density) * dof)
            for iID in eachindex(density_mass)
                density_mass[iID] = density[Int(ceil(iID / dof))]
            end

            # find indices of pd nodes and coupling nodes in reduced system; using master ids vector
            perm_pd_reduced = findall(in(perm_pd_nodes), perm_master)

            #  @timeit "Guyan Condensation" K_reduced,
            #  mass_reduced=guyan_reduction(K, density_mass,
            #                               perm_master, perm_slave)
            #
            #  # K_reduced[:, perm_pd_reduced] .= 0
            K_reduced = K
            mass_reduced = Diagonal(density_mass)
            K_reduced[perm_pd_reduced, :] .= 0
            #K_reduced[:, perm_pd_reduced] .= 0
            #   dropzeros!(mass_reduced)
            dropzeros!(K_reduced)

            Data_Manager.set_stiffness_matrix(sparse(K_reduced))
            Data_Manager.set_mass_matrix(lu(sparse(mass_reduced)))

            Data_Manager.set_reduced_model_pd(pd_nodes)
            Data_Manager.set_reduced_model_master(master_nodes)

            @info "Model reduction is applied"
            @info "Condensed: $(length(slave_nodes)), Master: $(length(master_nodes)), PD: $(length(pd_nodes))."
            return
        end
        @warn "No master nodes defined for model reduction. Using full stiffness matrix."
    end

    density_mass = zeros(length(density) * Data_Manager.get_dof())
    for iID in eachindex(density_mass)
        density_mass[iID] = density[Int(ceil(iID / dof))]
    end
    mass = Diagonal(density_mass)
    perm = collect(1:(Data_Manager.get_nnodes() * Data_Manager.get_dof()))
    Data_Manager.set_mass_matrix(lu(sparse(mass[perm, perm])))
    return

    # reduced matrix

end

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
    dt::Float64 = solver_options["dt"]
    nsteps::Int64 = solver_options["Number of Steps"]
    time::Float64 = solver_options["Initial Time"]
    step_time::Float64 = 0
    #max_cancel_damage::Float64 = solver_options["Maximum Damage"]
    numerical_damping::Float64 = solver_options["Numerical Damping"]
    max_damage::Float64 = 0
    damage_init::Bool = false

    coor::NodeVectorField{Float64} = Data_Manager.get_field("Coordinates")

    active_list::NodeScalarField{Bool} = Data_Manager.get_field("Active")
    volume::NodeScalarField{Float64} = Data_Manager.get_field("Volume")

    if "Material" in solver_options["Models"]
        external_forces::NodeVectorField{Float64} = Data_Manager.get_field("External Forces")
        external_force_densities::NodeVectorField{Float64} = Data_Manager.get_field("External Force Densities")
        a::NodeVectorField{Float64} = Data_Manager.get_field("Acceleration")
    end

    comm = Data_Manager.get_comm()
    rank = Data_Manager.get_rank()
    iter = progress_bar(rank, nsteps, silent)
    M_fact = Data_Manager.get_mass_matrix()
    if solver_options["Model Reduction"] != false
        master_nodes = Data_Manager.get_reduced_model_master()
        temp = zeros(length(master_nodes) * Data_Manager.get_dof())
        diff_nodes = setdiff(master_nodes, Data_Manager.get_reduced_model_pd())
    else
        temp = zeros(Data_Manager.get_dof() * Data_Manager.get_nnodes())
    end

    K::AbstractMatrix{Float64} = Data_Manager.get_stiffness_matrix()

    @timeit "Matrix Verlet" begin
        @inbounds @fastmath for idt in iter
            @timeit "solver start" begin
                uN::Matrix{Float64} = Data_Manager.get_field("Displacements", "N")
                uNP1::Matrix{Float64} = Data_Manager.get_field("Displacements", "NP1")
                forces::Matrix{Float64} = Data_Manager.get_field("Forces", "NP1")
                deformed_coorNP1::Matrix{Float64} = Data_Manager.get_field("Deformed Coordinates",
                                                                           "NP1")
                vN::Matrix{Float64} = Data_Manager.get_field("Velocity", "N")
                vNP1::Matrix{Float64} = Data_Manager.get_field("Velocity", "NP1")
                force_densities_NP1::Matrix{Float64} = Data_Manager.get_field("Force Densities",
                                                                              "NP1")
                active_nodes::Vector{Int64} = Data_Manager.get_field("Active Nodes")

                if solver_options["Model Reduction"] != false
                    active_nodes = master_nodes
                    active_list .= false
                    active_list[Data_Manager.get_reduced_model_pd()] .= true

                else
                    @timeit "active nodes" active_nodes=find_active_nodes(active_list,
                                                                          active_nodes,
                                                                          1:Data_Manager.get_nnodes())
                end
            end
            @timeit "compute Velocity" begin
                @. @views vNP1[active_nodes, :] = (1 - numerical_damping) *
                                                  vN[active_nodes, :] +
                                                  0.5 * dt * a[active_nodes, :]

                @. @views uNP1[active_nodes, :] = uN[active_nodes, :] +
                                                  dt * vNP1[active_nodes, :]
            end

            @timeit "apply BC" apply_bc_dirichlet(["Displacements", "Temperature"],
                                                  bcs,
                                                  time,
                                                  step_time) #-> Dirichlet
            @timeit "deformed coor" @views deformed_coorNP1 .= coor .+ uNP1

            @timeit "compute models" compute_models(block_nodes,
                                                    dt,
                                                    time,
                                                    solver_options["Models"],
                                                    synchronise_field)
            @timeit "Force computations" begin
                @timeit "sa" sa=size(a[active_nodes, :])
                @views fNP1 = force_densities_NP1[active_nodes, :]
                # masternode
                @info sum(force_densities_NP1[diff_nodes, :])
                @timeit "Force matrix computations" f_int_inplace!(fNP1, temp, K,
                                                                   vec(uNP1[active_nodes,
                                                                            :]), sa)
                @info sum(force_densities_NP1[diff_nodes, :])
                readline()
                @. @views fNP1 .+= external_force_densities[active_nodes, :] +
                                   external_forces[active_nodes, :] / volume[active_nodes]

                # forces = force_densities * volume (external bereits in force_densities drin)
                @. @views forces[active_nodes, :] .= fNP1 *
                                                     volume[active_nodes]
            end
            @timeit "Accelaration computation" begin
                @views a[active_nodes, :] = reshape(M_fact \
                                                    vec(fNP1),
                                                    sa...)
            end
            @timeit "write_results" result_files=write_results(result_files, time,
                                                               max_damage, outputs)

            if rank == 0 && !silent && Data_Manager.get_cancel()
                set_multiline_postfix(iter, "Simulation canceled!")
                break
            end
            @timeit "switch_NP1_to_N" Data_Manager.switch_NP1_to_N()
            @timeit "final solver steps" begin
                time += dt
                step_time += dt
                Data_Manager.set_current_time(time)

                if idt % ceil(nsteps / 100) == 0
                    @info "Step: $idt / $(nsteps+1) [$(@sprintf("%.3e", time)) s]"
                end
                if rank == 0 && !silent
                    set_postfix(iter, t = @sprintf("%.3e", time))
                end

                # barrier(comm)
            end
            #a+= + fexternal / density
        end
    end
    Data_Manager.set_current_time(time - dt)
    return result_files
end

function f_int_inplace!(F::AbstractMatrix{Float64},
                        temp::AbstractVector{Float64},
                        K,
                        u::AbstractVector{Float64}, sa::Tuple{Int64,Int64})
    mul!(temp, K, u)
    F .+= reshape(temp, sa)
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

end
