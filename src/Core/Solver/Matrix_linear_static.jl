# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Linear_static_matrix_based
using ProgressBars: set_multiline_postfix, set_postfix
using Printf
using SparseArrays
using LinearAlgebra
using TimerOutputs
include("../../MPI_communication/MPI_communication.jl")
using .MPI_communication: barrier
include("../../Support/Parameters/parameter_handling.jl")
using .Parameter_Handling: get_initial_time,
                           get_final_time
include("../../Support/Helpers.jl")
using .Helpers: check_inf_or_nan, find_active_nodes, progress_bar
include("../BC_manager.jl")
using .Boundary_conditions: apply_bc_dirichlet, apply_bc_neumann, find_bc_free_dof
include("../../Models/Material/Material_Models/Correspondence/Correspondence_matrix_based.jl")
using .Correspondence_matrix_based: assemble_stiffness_contributions_sparse,
                                    add_zero_energy_stiff

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
    # hier wird erstmal K bestimmt, bevor es woanders hinkommt.
    # indizes prüfen
    # vec und reshape prüfen
    #first test constant material
    find_bc_free_dof(datamanager, bcs)
    solver_options["Initial Time"] = get_initial_time(params, datamanager)
    solver_options["Final Time"] = get_final_time(params, datamanager)
    @warn "Number of steps is set to 1."
    solver_options["Number of Steps"] = 1

    for (block, nodes) in pairs(block_nodes)
        model_param = datamanager.get_properties(block, "Material Model")
        Correspondence_matrix_based.init_model(datamanager, nodes, model_param)
    end

    dof = datamanager.get_dof()
    C_voigt = datamanager.get_field("Hooke Matrix")
    inverse_shape_tensor = datamanager.get_field("Inverse Shape Tensor")
    nlist = datamanager.get_nlist()
    volume = datamanager.get_field("Volume")
    bond_geometry = datamanager.get_field("Bond Geometry")
    omega = datamanager.get_field("Influence Function")

    # verify K
    K_sparse = assemble_stiffness_contributions_sparse(datamanager.get_nnodes(),              # nnodes
                                                       dof,                # dof
                                                       C_voigt,            # C_Voigt
                                                       inverse_shape_tensor, # inverse_shape_tensor
                                                       nlist,              # nlist
                                                       volume,      # volume
                                                       bond_geometry,      # bond_geometry
                                                       omega)

    #add_zero_energy_stiff(
    #	K_sparse,
    #	datamanager.get_nnodes(),              # nnodes
    #	dof,                # dof
    #	C_voigt,            # C_Voigt
    #	inverse_shape_tensor, # inverse_shape_tensor
    #	nlist,              # nlist
    #	volume,      # volume
    #	bond_geometry,      # bond_geometry
    #	omega,
    #)
    # passt nicht zu dem Beispiel. matrix erstellung ist utnerschiedlich
    # - shape tensor prüfen
    # - funktionen verleichen

    datamanager.set_stiffness_matrix(K_sparse)
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

    @info "Run Linear Static Solver"
    max_damage::Float64 = 0
    volume = datamanager.get_field("Volume")
    coor = datamanager.get_field("Coordinates")
    comm = datamanager.get_comm()
    # needed to bring K_sparse into the consistent style of the transformation matrix -> vector
    perm = create_permutation(datamanager.get_nnodes(), datamanager.get_dof())
    if "Material" in solver_options["Models"]
        external_forces = datamanager.get_field("External Forces")
        external_force_densities = datamanager.get_field("External Force Densities")
        a = datamanager.get_field("Acceleration")
    end

    dt::Float64 = solver_options["Final Time"]
    nsteps::Int64 = solver_options["Number of Steps"]
    time::Float64 = solver_options["Initial Time"]
    step_time::Float64 = 0
    rank = datamanager.get_rank()
    iter = progress_bar(rank, nsteps, silent)
    #nodes::Vector{Int64} = Vector{Int64}(1:datamanager.get_nnodes())

    @inbounds @fastmath for idt in iter
        datamanager.set_iteration(idt)
        @timeit to "Linear Static" begin
            uNP1 = datamanager.get_field("Displacements", "NP1")
            deformed_coorNP1 = datamanager.get_field("Deformed Coordinates", "NP1")
            forces = datamanager.get_field("Forces", "NP1")
            force_densities = datamanager.get_field("Force Densities", "NP1")
            volume = datamanager.get_field("Volume")
            forces = datamanager.get_field("Forces", "NP1")
            K = datamanager.get_stiffness_matrix()
            external_force_densities = datamanager.get_field("External Force Densities")

            datamanager = apply_bc_dirichlet(["Displacements", "Forces", "Force Densities"],
                                             bcs,
                                             datamanager, time,
                                             step_time)
            external_force_densities.=external_forces ./ volume
            # reshape
            non_BCs = datamanager.get_bc_free_dof()
            compute_displacements!(K[perm, perm],
                                   non_BCs,
                                   uNP1,
                                   force_densities,
                                   external_force_densities)
            nodes = length(volume)
            for iID in 1:nodes
                @views forces[iID, :] .= force_densities[iID, :] .* volume[iID]
            end
            @views deformed_coorNP1 .= coor .+ uNP1

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

"""
	compute_displacements!(K, non_BCs, u, F, F_temp, K_reduced, lu_fact, temp)

Compute displacements with prescribed displacement boundary conditions.
No memory allocations.

Arguments:
- K: Global stiffness matrix
- non_BCs: Free DOF indices
- u: Displacement vector (contains prescribed values at fixed DOFs)
- F: External force vector
- F_temp: Temporary force vector (pre-allocated)

"""
function compute_displacements!(K::SparseMatrixCSC{Float64,Int64},
                                non_BCs::Vector{Int64},
                                u::Matrix{Float64},
                                F_int::Matrix{Float64},
                                F_ext::Matrix{Float64})

    # Compute modified force: F_modified = F - K * u_prescribed
    # (u contains prescribed displacements at fixed DOFs)

    #mul!(vec(F_int), K, vec(u))  # F_temp = K * u (in-place, no allocation)
    #F_int .+= F_ext   # F_temp = F - K*u (in-place)
    ## Update displacement vector at free DOFs
    #u[non_BCs] = K[non_BCs, non_BCs] \ F_int[non_BCs]

    F_modified = copy(vec(F_ext))
    BCs = setdiff(1:length(vec(u)), non_BCs)

    @views F_modified[non_BCs] .-= (K[non_BCs, BCs] * vec(u)[BCs])
    @views vec(u)[non_BCs] .= K[non_BCs, non_BCs] \ F_modified[non_BCs]

    #return reshape(F_int, 2, :)

end

function create_permutation(nnodes::Int, dof::Int)
    perm = Vector{Int}(undef, nnodes * dof)
    idx = 1
    for d in 1:dof           # Für jeden DOF
        for n in 1:nnodes     # Für jeden Knoten
            old_idx = (n-1)*dof + d  # Row-major: node, dann dof
            perm[idx] = old_idx
            idx += 1
        end
    end
    return perm
end
end
