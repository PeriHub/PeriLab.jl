# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

"""
start julia
]
activate .
run the code in the main structure in the REPL
if not in the REPL, make sure you that you adapt the pathes

"""
# for REPL
include("./test/developement_area/corr_stiff/corr_stiffnessmatrix.jl")
include("./test/developement_area/corr_stiff/init_fields.jl")
using .CorrespondenceStiffnessMatrix: assemble_stiffness_contributions_sparse,
                                      voigt_to_tensor4_2d,
                                      elasticity_matrix_2d_plane_stress, create_B_tensor,
                                      contraction, compute_linearized_operator
using .Data_manager
using SparseArrays
using LinearAlgebra
dm.create_constant_node_field("Hooke Matrix", Float64,
                              Int64((dof * (dof + 1)) / 2),
                              VectorOrMatrix = "Matrix")

material = Dict("Shear Modulus"=>20, "Bulk Modulus"=>40)

hm = elasticity_matrix_2d_plane_stress(material) # is in PeriLab
#    else
#        C_voigt = elasticity_matrix_2d_plane_strain(mat)
#    end
C_voigt = dm.get_field("Hooke Matrix")

for i in 1:nodes
    C_voigt[i, :, :] = hm
end

coor = dm.get_field("Coordinates")
bond_geometry = dm.create_constant_bond_field("Bond Geometry", Float64, dof)
inverse_shape_tensor = dm.create_constant_node_field("Inverse Shape Tensor", Float64, dof,
                                                     VectorOrMatrix = "Matrix")

omega = dm.create_constant_bond_field("Influence Function", Float64, 1, 1.0)

for iID in 1:nodes
    inverse_shape_tensor[iID, 1, 1]=1
    inverse_shape_tensor[iID, 2, 2]=1
    for (jID, nID) in enumerate(nlist[iID])
        bond_geometry[iID][jID] = coor[nID, :] - coor[iID, :]
    end
end

function assemble_stiffness_contributions_sparse_edit(nnodes::Int64,
                                                      dof::Int64,
                                                      C_Voigt::Array{Float64,3},
                                                      inverse_shape_tensor::Array{Float64,
                                                                                  3},
                                                      nlist::Vector{Vector{Int64}},
                                                      volume::Vector{Float64},
                                                      bond_geometry::Vector{Vector{Vector{Float64}}},
                                                      omega::Vector{Vector{Float64}})
    K = zeros(nnodes*dof, nnodes*dof)
    # Process each node i
    for i in 1:nnodes
        D_inv = inverse_shape_tensor[i, :, :]
        @views C_tensor = voigt_to_tensor4_2d(C_Voigt[i, :, :])

        K_ijk = compute_all_linearized_operators(i, C_tensor, D_inv,
                                                 volume, bond_geometry,
                                                 omega, nlist)

        # Process each neighbor j of node i
        for (j_idx, j) in enumerate(nlist[i])
            omega_ij = omega[i][j_idx]
            V_i = volume[i]
            V_j = volume[j]

            # Local K_ij matrix for this bond
            K_ij = zeros(dof, nnodes * dof)

            # Sum over all k in neighborhood of i
            for (k_idx, k) in enumerate(nlist[i])
                omega_ik = omega[i][k_idx]
                V_k = volume[k]

                # Get pre-computed K_ijk block
                # This already contains temp = CB_ik : D_i^(-1) ⊗ X_ij
                K_block = K_ijk[(i, j, k)]

                # Fill K_ij based on K_block
                for m in 1:dof
                    for o in 1:dof
                        K_ij[m, (i-1)*dof+o] -= omega_ij * V_k * omega_ik * K_block[m, o]
                        K_ij[m, (k-1)*dof+o] += omega_ij * V_k * omega_ik * K_block[m, o]
                    end
                end
            end

            # Add K_ij to global matrix K
            for m in 1:dof
                for n in 1:(nnodes * dof)
                    K[(j-1)*dof+m, n] -= K_ij[m, n] * V_i
                    K[(i-1)*dof+m, n] += K_ij[m, n] * V_j
                end
            end
        end
    end

    return sparse(K)
end

function compute_all_linearized_operators(iID::Int64, C_tensor::Array{Float64,4},
                                          D_inv::Matrix{Float64}, volume::Vector{Float64},
                                          bond_geometry::Vector{Vector{Vector{Float64}}},
                                          omega::Vector{Vector{Float64}},
                                          nlist::Vector{Vector{Int64}})
    K_operators = Dict{Tuple{Int64,Int64,Int64},Matrix{Float64}}()

    # Loop over all neighbors jID of node iID
    for (jID, njID) in enumerate(nlist[iID])
        X_ij = bond_geometry[iID][jID]
        omega_ij = omega[iID][jID]

        # Loop over all neighbors kID of node iID
        for (kID, nkID) in enumerate(nlist[iID])
            V_k = volume[kID]
            omega_ik = omega[iID][kID]
            X_ik = bond_geometry[iID][kID]

            # Create B-tensor with omega_ik * V_k / 2 factor
            B_ik = create_B_tensor(D_inv, X_ik, V_k, omega_ik)
            CB = contraction(C_tensor, B_ik)

            # Compute K_ijk operator
            K_ijk = compute_linearized_operator(CB, D_inv,
                                                X_ij, omega_ij, V_k, omega_ik)
            K_operators[(iID, njID, nkID)] = K_ijk
        end
    end

    return K_operators
end

# Assemble contributions from each node
K_sparse = assemble_stiffness_contributions_sparse_edit(nodes,              # nnodes
                                                        dof,                # dof
                                                        C_voigt,            # C_Voigt
                                                        inverse_shape_tensor, # inverse_shape_tensor
                                                        nlist,              # nlist
                                                        volume,      # volume
                                                        bond_geometry,      # bond_geometry
                                                        omega)

function compute_displacments!(K::SparseMatrixCSC{Float64,Int64}, non_BCs::Vector{Int64},
                               u::Vector{Float64}, F::Vector{Float64})
    @views F = K*u
    @views u[non_BCs] = K[non_BCs, non_BCs] / F[non_BCs]'
end

F = zeros(nodes*dof)
u = zeros(nodes*dof)
u[16:18] .= 0.1
non_BCs=collect(5:15)

compute_displacments!(K, non_BCs, u, F)

function run_solver(solver_options, datamanager)
    @info "Run Linear Static Solver"
    volume = datamanager.get_field("Volume")
    coor = datamanager.get_field("Coordinates")

    if "Material" in solver_options["Models"]
        external_forces = datamanager.get_field("External Forces")
        external_force_densities = datamanager.get_field("External Force Densities")
    end

    dt::Float64 = solver_options["dt"]
    nsteps::Int64 = solver_options["Number of Steps"]
    time::Float64 = solver_options["Initial Time"]
    step_time::Float64 = 0
    rank = datamanager.get_rank()
    iter = progress_bar(rank, nsteps, silent)

    #u[freenodes]=K[free_nodes, free_nodes]/F[free_nodes]

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
