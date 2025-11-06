# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Correspondence_matrix_based
using Base.Threads
using LinearAlgebra
using LoopVectorization: @turbo
using ...Material_Basis: get_Hooke_matrix
using ....Helpers: get_fourth_order, progress_bar
using ...Global_Zero_Energy_Control: create_zero_energy_mode_stiffness!
using SparseArrays
using StaticArrays: @MMatrix
export init_model
export add_zero_energy_stiff!
export compute_bond_force

"""
Peridynamic Correspondence Stiffness Matrix Assembly

Memory-optimized with pre-allocated mapping structure:
- mapping[i][1] → start index in V for block (i, i)
- mapping[i][k_idx+1] → start index in V for block (i, nlist[i][k_idx])
"""

kronecker_delta(x::Int64, y::Int64) = ==(x, y)

function contraction!(C::Array{Float64,4}, B::Array{Float64,3}, dof::Int64,
                      CB::AbstractArray{Float64,3})
    fill!(CB, 0.0)

    # Unroll for common cases
    if dof == 2
        @turbo for q in 1:2, p in 1:2, o in 1:2
            B_opq = B[o, p, q]
            for n in 1:2, m in 1:2
                CB[m, n, q] += C[m, n, o, p] * B_opq
            end
        end
    elseif dof == 3
        @turbo for q in 1:3, p in 1:3, o in 1:3
            B_opq = B[o, p, q]
            for n in 1:3, m in 1:3
                CB[m, n, q] += C[m, n, o, p] * B_opq
            end
        end
    else
        @inbounds for q in 1:dof, p in 1:dof, o in 1:dof
            B_opq = B[o, p, q]
            for n in 1:dof, m in 1:dof
                CB[m, n, q] += C[m, n, o, p] * B_opq
            end
        end
    end
end
function create_B_tensor!(D_inv::AbstractMatrix{Float64}, X_ik::Vector{Float64},
                          V_k::Float64, omega_ik::Float64, ::Val{2},
                          B_tensor::Array{Float64,3})
    factor = omega_ik * V_k * 0.5

    B1_1 = D_inv[1, 1] * X_ik[1] + D_inv[2, 1] * X_ik[2]
    B1_2 = D_inv[1, 2] * X_ik[1] + D_inv[2, 2] * X_ik[2]
    B2_1 = D_inv[1, 1] * X_ik[1] + D_inv[1, 2] * X_ik[2]
    B2_2 = D_inv[2, 1] * X_ik[1] + D_inv[2, 2] * X_ik[2]

    B_tensor[1, 1, 1] = factor * (B1_1 + B2_1)
    B_tensor[1, 1, 2] = 0.0
    B_tensor[1, 2, 1] = factor * B1_2
    B_tensor[1, 2, 2] = factor * B2_1
    B_tensor[2, 1, 1] = factor * B2_2
    B_tensor[2, 1, 2] = factor * B1_1
    B_tensor[2, 2, 1] = 0.0
    B_tensor[2, 2, 2] = factor * (B1_2 + B2_2)

    return nothing
end

function create_B_tensor!(D_inv::AbstractMatrix{Float64}, X_ik::Vector{Float64},
                          V_k::Float64, omega_ik::Float64, ::Val{3},
                          B_tensor::Array{Float64,3})
    factor = omega_ik * V_k * 0.5

    B1_1 = D_inv[1, 1] * X_ik[1] + D_inv[2, 1] * X_ik[2] + D_inv[3, 1] * X_ik[3]
    B1_2 = D_inv[1, 2] * X_ik[1] + D_inv[2, 2] * X_ik[2] + D_inv[3, 2] * X_ik[3]
    B1_3 = D_inv[1, 3] * X_ik[1] + D_inv[2, 3] * X_ik[2] + D_inv[3, 3] * X_ik[3]
    B2_1 = D_inv[1, 1] * X_ik[1] + D_inv[1, 2] * X_ik[2] + D_inv[1, 3] * X_ik[3]
    B2_2 = D_inv[2, 1] * X_ik[1] + D_inv[2, 2] * X_ik[2] + D_inv[2, 3] * X_ik[3]
    B2_3 = D_inv[3, 1] * X_ik[1] + D_inv[3, 2] * X_ik[2] + D_inv[3, 3] * X_ik[3]

    B_tensor[1, 1, 1] = factor * (B1_1 + B2_1)
    B_tensor[1, 1, 2] = 0.0
    B_tensor[1, 1, 3] = 0.0
    B_tensor[1, 2, 1] = factor * B1_2
    B_tensor[1, 2, 2] = factor * B2_1
    B_tensor[1, 2, 3] = 0.0
    B_tensor[1, 3, 1] = factor * B1_3
    B_tensor[1, 3, 2] = 0.0
    B_tensor[1, 3, 3] = factor * B2_1
    B_tensor[2, 1, 1] = factor * B2_2
    B_tensor[2, 1, 2] = factor * B1_1
    B_tensor[2, 1, 3] = 0.0
    B_tensor[2, 2, 1] = 0.0
    B_tensor[2, 2, 2] = factor * (B1_2 + B2_2)
    B_tensor[2, 2, 3] = 0.0
    B_tensor[2, 3, 1] = 0.0
    B_tensor[2, 3, 2] = factor * B1_3
    B_tensor[2, 3, 3] = factor * B2_2
    B_tensor[3, 1, 1] = factor * B2_3
    B_tensor[3, 1, 2] = 0.0
    B_tensor[3, 1, 3] = factor * B1_1
    B_tensor[3, 2, 1] = 0.0
    B_tensor[3, 2, 2] = factor * B2_3
    B_tensor[3, 2, 3] = factor * B1_2
    B_tensor[3, 3, 1] = 0.0
    B_tensor[3, 3, 2] = 0.0
    B_tensor[3, 3, 3] = factor * (B1_3 + B2_3)

    return nothing
end

function compute_linearized_operator(CB_ik::AbstractArray{Float64,3},
                                     D_inv_i::AbstractMatrix{Float64},
                                     X_ij::Vector{Float64}, omega_ij::Float64,
                                     V_k::Float64, omega_ik::Float64, dof::Int64,
                                     K_ijk::AbstractMatrix{Float64})
    if dof == 2
        DX_1 = D_inv_i[1, 1] * X_ij[1] + D_inv_i[1, 2] * X_ij[2]
        DX_2 = D_inv_i[2, 1] * X_ij[1] + D_inv_i[2, 2] * X_ij[2]
        factor = omega_ij * V_k * omega_ik

        @inbounds begin
            K_ijk[1, 1] = factor * (CB_ik[1, 1, 1]*DX_1 + CB_ik[1, 2, 1]*DX_2)
            K_ijk[1, 2] = factor * (CB_ik[1, 1, 2]*DX_1 + CB_ik[1, 2, 2]*DX_2)
            K_ijk[2, 1] = factor * (CB_ik[2, 1, 1]*DX_1 + CB_ik[2, 2, 1]*DX_2)
            K_ijk[2, 2] = factor * (CB_ik[2, 1, 2]*DX_1 + CB_ik[2, 2, 2]*DX_2)
        end
    elseif dof == 3
        DX_1 = D_inv_i[1, 1]*X_ij[1] + D_inv_i[1, 2]*X_ij[2] + D_inv_i[1, 3]*X_ij[3]
        DX_2 = D_inv_i[2, 1]*X_ij[1] + D_inv_i[2, 2]*X_ij[2] + D_inv_i[2, 3]*X_ij[3]
        DX_3 = D_inv_i[3, 1]*X_ij[1] + D_inv_i[3, 2]*X_ij[2] + D_inv_i[3, 3]*X_ij[3]
        factor = omega_ij * V_k * omega_ik

        @inbounds begin
            K_ijk[1, 1] = factor *
                          (CB_ik[1, 1, 1]*DX_1 + CB_ik[1, 2, 1]*DX_2 + CB_ik[1, 3, 1]*DX_3)
            K_ijk[1, 2] = factor *
                          (CB_ik[1, 1, 2]*DX_1 + CB_ik[1, 2, 2]*DX_2 + CB_ik[1, 3, 2]*DX_3)
            K_ijk[1, 3] = factor *
                          (CB_ik[1, 1, 3]*DX_1 + CB_ik[1, 2, 3]*DX_2 + CB_ik[1, 3, 3]*DX_3)
            K_ijk[2, 1] = factor *
                          (CB_ik[2, 1, 1]*DX_1 + CB_ik[2, 2, 1]*DX_2 + CB_ik[2, 3, 1]*DX_3)
            K_ijk[2, 2] = factor *
                          (CB_ik[2, 1, 2]*DX_1 + CB_ik[2, 2, 2]*DX_2 + CB_ik[2, 3, 2]*DX_3)
            K_ijk[2, 3] = factor *
                          (CB_ik[2, 1, 3]*DX_1 + CB_ik[2, 2, 3]*DX_2 + CB_ik[2, 3, 3]*DX_3)
            K_ijk[3, 1] = factor *
                          (CB_ik[3, 1, 1]*DX_1 + CB_ik[3, 2, 1]*DX_2 + CB_ik[3, 3, 1]*DX_3)
            K_ijk[3, 2] = factor *
                          (CB_ik[3, 1, 2]*DX_1 + CB_ik[3, 2, 2]*DX_2 + CB_ik[3, 3, 2]*DX_3)
            K_ijk[3, 3] = factor *
                          (CB_ik[3, 1, 3]*DX_1 + CB_ik[3, 2, 3]*DX_2 + CB_ik[3, 3, 3]*DX_3)
        end
    end
end

function precompute_CB_tensors!(CB_k::AbstractArray{Float64,4},
                                C_tensor::Array{Float64,4},
                                D_inv::AbstractMatrix{Float64},
                                neighbors::Vector{Int64},
                                volume::Vector{Float64},
                                bond_geometry_i::Vector{Vector{Float64}},
                                omega_i::Vector{Float64},
                                bond_damage_i::Vector{Float64},
                                dof::Int64,
                                B_ik::Array{Float64,3})
    @inbounds for (k_idx, k) in enumerate(neighbors)
        V_k = volume[k]
        omega_ik = omega_i[k_idx] * bond_damage_i[k_idx]
        X_ik = bond_geometry_i[k_idx]
        create_B_tensor!(D_inv, X_ik, V_k, omega_ik, Val(dof), B_ik)
        @views contraction!(C_tensor, B_ik, dof, CB_k[k_idx, :, :, :])
    end
    return nothing
end

"""
	create_mapping_structure(nodes, nlist, number_of_neighbors, dof)

Create pre-allocated mapping structure and I, J, V vectors.

Structure:
- mapping[i][1] → start index in V for block (i, i)
- mapping[i][k_idx+1] → start index in V for block (i, nlist[i][k_idx])

Size: (sum(number_of_neighbors) + length(nodes)) * dof * dof
"""
function create_mapping_structure(nodes::AbstractVector{Int64},
                                  nlist::Vector{Vector{Int64}},
                                  number_of_neighbors::Vector{Int64},
                                  dof::Int64)

    # Calculate total entries
    total_blocks = sum(number_of_neighbors) + length(nodes)
    total_entries = total_blocks * dof * dof

    # Pre-allocate I, J, V
    I = zeros(Int64, total_entries)
    J = zeros(Int64, total_entries)
    V = zeros(Float64, total_entries)

    # Create mapping structure
    mapping = Vector{Vector{Int64}}(undef, nodes[end])

    idx = 1
    for i in nodes
        ni = nlist[i]
        n_neighbors = number_of_neighbors[i]

        # Allocate mapping vector: [1] for i, [2..n+1] for nlist[i]
        mapping[i] = Vector{Int64}(undef, n_neighbors + 1)

        # Block (i, i)
        mapping[i][1] = idx
        for n in 1:dof
            for m in 1:dof
                I[idx] = (i-1)*dof + m
                J[idx] = (i-1)*dof + n
                idx += 1
            end
        end

        # Blocks (i, k) for k ∈ nlist[i]
        for (k_idx, k) in enumerate(ni)
            mapping[i][k_idx+1] = idx
            for n in 1:dof
                for m in 1:dof
                    I[idx] = (i-1)*dof + m
                    J[idx] = (k-1)*dof + n
                    idx += 1
                end
            end
        end
    end

    n_total = nodes[end] * dof

    return mapping, I, J, V, n_total
end

"""
	add_block_direct!(V, start_idx, block, dof)

Add dof×dof block directly to V at start_idx.
"""
@inline function add_block_direct!(V::Vector{Float64},
                                   start_idx::Int64,
                                   block::AbstractMatrix{Float64},
                                   dof::Int64)
    idx = start_idx
    @inbounds for n in 1:dof
        for m in 1:dof
            V[idx] += block[m, n]
            idx += 1
        end
    end
end

"""
	assemble_stiffness(...)

Optimized assembly with mapping vector structure.

Uses:
- mapping[i][1] for column i
- mapping[i][k_idx+1] for column nlist[i][k_idx]
- mapping[j][...] for writing row j when j ∈ nlist[i]
"""
function assemble_stiffness(nodes::AbstractVector{Int64},
                            active_nodes::AbstractVector{Int64},
                            dof::Int64,
                            C_Voigt::Array{Float64,3},
                            inverse_shape_tensor::Array{Float64,3},
                            number_of_neighbors::Vector{Int64},
                            nlist::Vector{Vector{Int64}},
                            volume::Vector{Float64},
                            bond_geometry::Vector{Vector{Vector{Float64}}},
                            omega::Vector{Vector{Float64}},
                            bond_damage::Vector{Vector{Float64}})

    # TODO: reuse mapping structure from initialization
    # right now nodes are needed to get the full mapping structure
    mapping, I, J, V,
    n_total = create_mapping_structure(nodes, nlist,
                                       number_of_neighbors, dof)

    # Neighbor lookup (array-based for speed)
    max_nodes = maximum(nodes)
    neighbor_positions = [zeros(Int64, max_nodes) for _ in 1:max_nodes]
    for node in nodes
        for (idx, neighbor) in enumerate(nlist[node])
            neighbor_positions[node][neighbor] = idx
        end
    end

    # Partition nodes by workload
    n_threads = Threads.nthreads()
    #workloads = [number_of_neighbors[i]^2 for i in nodes]
    workloads = [number_of_neighbors[i]^2 for i in active_nodes]
    #    node_partitions = partition_by_workload(nodes, workloads, n_threads)

    node_partitions = partition_by_workload(active_nodes, workloads, n_threads)

    # Thread-local V arrays
    thread_V = [zeros(Float64, length(V)) for _ in 1:n_threads]

    # Thread-local buffers
    max_neighbors = maximum(number_of_neighbors)
    thread_buffers = [(CB_k = zeros(max_neighbors, dof, dof, dof),
                       K_block = dof == 2 ? @MMatrix(zeros(2, 2)) : @MMatrix(zeros(3, 3)),
                       B_ik = dof == 2 ? zeros(2, 2, 2) : zeros(3, 3, 3))
                      for _ in 1:n_threads]

    # Process each partition in parallel
    Threads.@threads for tid in 1:n_threads
        buffers = thread_buffers[tid]
        V_thread = thread_V[tid]
        for i in node_partitions[tid]
            D_inv = @view inverse_shape_tensor[i, :, :]
            C_tensor = get_fourth_order(@view(C_Voigt[i, :, :]), dof)
            V_i = volume[i]
            ni = nlist[i]

            CB_k = @view buffers.CB_k[eachindex(ni), :, :, :]
            precompute_CB_tensors!(CB_k, C_tensor, D_inv, ni, volume,
                                   bond_geometry[i], omega[i], bond_damage[i], dof,
                                   buffers.B_ik)

            i_diag_idx = mapping[i][1]
            mapping_i = mapping[i]

            for (j_idx, j) in enumerate(ni)
                omega_ij = omega[i][j_idx]*bond_damage[i][j_idx]
                V_j = volume[j]
                X_ij = bond_geometry[i][j_idx]

                neighbor_pos_j = neighbor_positions[j]
                j_diag_idx = mapping[j][1]
                mapping_j = mapping[j]

                for (k_idx, k) in enumerate(ni)
                    omega_ik = omega[i][k_idx]*bond_damage[i][k_idx]
                    V_k = volume[k]

                    compute_linearized_operator(@view(CB_k[k_idx, :, :, :]),
                                                D_inv, X_ij, omega_ij, V_k, omega_ik,
                                                dof, buffers.K_block)

                    # Row j entries
                    i_pos_in_j = neighbor_pos_j[i]
                    if i_pos_in_j != 0
                        add_block_scaled!(V_thread, mapping_j[i_pos_in_j+1],
                                          buffers.K_block, V_i, dof)
                    end

                    if k == j
                        add_block_scaled!(V_thread, j_diag_idx, buffers.K_block, -V_i, dof)
                    else
                        k_pos_in_j = neighbor_pos_j[k]
                        if k_pos_in_j != 0
                            add_block_scaled!(V_thread, mapping_j[k_pos_in_j+1],
                                              buffers.K_block, -V_i, dof)
                        end
                    end

                    # Row i entries
                    add_block_scaled!(V_thread, i_diag_idx, buffers.K_block, -V_j, dof)
                    add_block_scaled!(V_thread, mapping_i[k_idx+1], buffers.K_block, V_j,
                                      dof)
                end
            end
        end
    end

    @inbounds Threads.@threads for i in eachindex(V)
        for tid in 1:n_threads
            V[i] += thread_V[tid][i]
        end
    end

    return I, J, V, n_total
end

# Helper function for load balancing
function partition_by_workload(nodes, workloads, n_threads)
    total_work = sum(workloads)
    target_work = total_work / n_threads

    partitions = [Int64[] for _ in 1:n_threads]
    current_thread = 1
    current_work = 0.0

    for (i, node) in enumerate(nodes)
        if current_thread < n_threads && current_work >= target_work
            current_thread += 1
            current_work = 0.0
        end
        push!(partitions[current_thread], node)
        current_work += workloads[i]
    end

    return partitions
end

@inline function add_block_scaled!(V, start_idx, block, scale, dof)
    @turbo for i in 1:dof, j in 1:dof
        V[start_idx+(j-1)*dof+i-1] += scale * block[i, j]
    end
end

# [Include all other functions from working version]

function add_zero_energy_stiff!(K::SparseMatrixCSC{Float64,Int64},
                                active_nodes::AbstractVector{Int64},
                                dof::Int64,
                                zStiff::Array{Float64,3},
                                inverse_shape_tensor::Array{Float64,3},
                                nlist::Vector{Vector{Int}},
                                volume::Vector{Float64},
                                bond_geometry::Vector{Vector{Vector{Float64}}},
                                bond_damage::Vector{Vector{Float64}},
                                omega::Vector{Vector{Float64}})
    I_indices = Int[]
    J_indices = Int[]
    values = Float64[]
    # TODO: not very memory efficient
    S_matrices = Vector{Matrix{Float64}}(undef, maximum(active_nodes))
    D_inv_X = Vector{Vector{Vector{Float64}}}(undef, maximum(active_nodes))
    1
    for i in active_nodes
        neighbors = nlist[i]
        n_neighbors = length(neighbors)
        @views D_inv_i = inverse_shape_tensor[i, :, :]
        D_inv_X[i] = [D_inv_i * bond_geometry[i][idx] for idx in 1:n_neighbors]

        S_matrices[i] = zeros(n_neighbors, n_neighbors)
        for p in 1:n_neighbors
            if bond_damage[i][p] == 0
                continue
            end
            for q in 1:n_neighbors
                if bond_damage[i][q] == 0
                    continue
                end
                S_matrices[i][p, q] = dot(bond_geometry[i][p], D_inv_X[i][q])
            end
        end
    end

    for i in active_nodes
        neighbors = nlist[i]
        @views Z_i = zStiff[i, :, :]
        V_i = volume[i]

        for (idx_j, j) in enumerate(neighbors)
            if i == j
                continue
            end

            scalar_sum_j = 0.0
            for idx_k in eachindex(neighbors)
                k = neighbors[idx_k]
                scalar_sum_j += omega[i][idx_k] * bond_damage[i][idx_k] * volume[k] *
                                S_matrices[i][idx_k, idx_j]
            end

            K_ji = -V_i * (1.0 - scalar_sum_j) * Z_i
            add_block_to_coo!(I_indices, J_indices, values, K_ji, j, i, dof)

            for (idx_k, k) in enumerate(neighbors)
                if k == i || k == j
                    continue
                end

                s_value = omega[i][idx_k] * bond_damage[i][idx_k] * volume[k] *
                          S_matrices[i][idx_k, idx_j]
                K_jk = -V_i * s_value * Z_i
                add_block_to_coo!(I_indices, J_indices, values, K_jk, j, k, dof)
            end
        end
    end

    for i in active_nodes
        neighbors = nlist[i]
        @views Z_i = zStiff[i, :, :]
        V_i = volume[i]

        for (idx_j, j) in enumerate(neighbors)
            if i == j
                continue
            end

            K_ij_total = zeros(dof, dof)

            for (idx_bond, bond_node) in enumerate(neighbors)
                if bond_node == i
                    continue
                end

                if idx_bond == idx_j
                    s_jj = S_matrices[i][idx_j, idx_j]
                    K_ij_total .-= V_i *
                                   (1.0 -
                                    omega[i][idx_j] * bond_damage[i][idx_j] * volume[j] *
                                    s_jj) * Z_i
                else
                    s_kj = S_matrices[i][idx_bond, idx_j]
                    K_ij_total .+= V_i * omega[i][idx_j] * bond_damage[i][idx_j] *
                                   volume[j] * s_kj * Z_i
                end
            end

            add_block_to_coo!(I_indices, J_indices, values, K_ij_total, i, j, dof)
        end
    end

    K_stab = sparse(I_indices, J_indices, values, size(K)...)

    for i in active_nodes
        K_ii_block = zeros(dof, dof)

        for d1 in 1:dof
            row = (i-1)*dof + d1
            for j in active_nodes
                if i != j
                    for d2 in 1:dof
                        col = (j-1)*dof + d2
                        K_ii_block[d1, d2] -= K_stab[row, col]
                    end
                end
            end
        end

        for d1 in 1:dof
            for d2 in 1:dof
                row = (i-1)*dof + d1
                col = (i-1)*dof + d2
                K[row, col] += K_ii_block[d1, d2]
            end
        end
    end
    K .+= K_stab

    return nothing
end

function add_block_to_coo!(I_indices::Vector{Int},
                           J_indices::Vector{Int},
                           values::Vector{Float64},
                           block::Matrix{Float64},
                           i::Int,
                           j::Int,
                           dof::Int)
    for d1 in 1:dof
        for d2 in 1:dof
            if abs(block[d1, d2]) > 1e-14
                push!(I_indices, (i-1)*dof + d1)
                push!(J_indices, (j-1)*dof + d2)
                push!(values, block[d1, d2])
            end
        end
    end
end

function compute_bond_force(bond_force::Vector{Vector{Vector{Float64}}},
                            K::SparseMatrixCSC{Float64,Int64}, u::Matrix{Float64},
                            nodes::AbstractVector{Int64}, nlist::Vector{Vector{Int64}},
                            dof::Int64)
    @inbounds for i in nodes
        ni = nlist[i]
        i_offset = (i-1)*dof

        for (j_idx, j) in enumerate(ni)
            if i == j
                continue
            end

            j_offset = (j-1)*dof
            f_ij = zeros(dof)

            for m in 1:dof
                for n in 1:dof
                    f_ij[m] += K[i_offset+m, j_offset+n] * (u[j, n] - u[i, n])
                end
            end

            bond_force[i][j_idx] = f_ij
        end
    end
end

function init_model(datamanager::Module,
                    nodes::AbstractVector{Int64},
                    material_parameter::Dict)
    if datamanager.get_max_rank()>1
        @error "Correspondence matrix based not implemented for parallel runs."
    end
    return datamanager
end

function init_matrix(datamanager::Module)
    nodes = collect(1:datamanager.get_nnodes())
    dof = datamanager.get_dof()
    zStiff = datamanager.create_constant_node_field("Zero Energy Stiffness",
                                                    Float64,
                                                    dof,
                                                    VectorOrMatrix = "Matrix")
    bond_geometry = datamanager.get_field("Bond Geometry")
    inverse_shape_tensor = datamanager.get_field("Inverse Shape Tensor")
    nlist = datamanager.get_nlist()
    number_of_neighbors=datamanager.get_field("Number of Neighbors")
    volume = datamanager.get_field("Volume")
    omega = datamanager.get_field("Influence Function")
    C_voigt = datamanager.get_field("Hooke Matrix")
    bond_geometry_N = datamanager.get_field("Deformed Bond Geometry", "N")
    bond_damage = datamanager.get_field("Bond Damage", "NP1")
    @info "Initializing stiffness matrix (mapping optimized)"
    index_x, index_y, vals,
    total_dof = assemble_stiffness(nodes,
                                   nodes,
                                   dof,
                                   C_voigt,
                                   inverse_shape_tensor,
                                   number_of_neighbors,
                                   nlist,
                                   volume,
                                   bond_geometry,
                                   omega,
                                   bond_damage)

    datamanager.init_stiffness_matrix(index_x, index_y, vals, total_dof)
    K_sparse = datamanager.get_stiffness_matrix()
    create_zero_energy_mode_stiffness!(nodes, dof, C_voigt,
                                       inverse_shape_tensor, zStiff)

    add_zero_energy_stiff!(K_sparse,
                           nodes,
                           dof,
                           zStiff,
                           inverse_shape_tensor,
                           nlist,
                           volume,
                           bond_geometry_N,
                           bond_damage,
                           omega)
end

function compute_model(datamanager::Module, nodes::AbstractVector{Int64})
    dof::Int64 = datamanager.get_dof()
    C_voigt::Array{Float64,3} = datamanager.get_field("Hooke Matrix")
    inverse_shape_tensor::Array{Float64,3} = datamanager.get_field("Inverse Shape Tensor")
    nlist::Vector{Vector{Int64}} = datamanager.get_nlist()
    volume::Vector{Float64} = datamanager.get_field("Volume")
    bond_geometry_N::Vector{Vector{Vector{Float64}}} = datamanager.get_field("Deformed Bond Geometry",
                                                                             "N")
    number_of_neighbors::Vector{Int64} = datamanager.get_field("Number of Neighbors")
    omega::Vector{Vector{Float64}} = datamanager.get_field("Influence Function")
    bond_damage::Vector{Vector{Float64}} = datamanager.get_field("Bond Damage", "NP1")

    zStiff = datamanager.get_field("Zero Energy Stiffness")
    # TODO: optimize update
    index_x, index_y, vals,
    total_dof = assemble_stiffness(collect(1:datamanager.get_nnodes()),
                                   nodes,
                                   dof,
                                   C_voigt,
                                   inverse_shape_tensor,
                                   number_of_neighbors,
                                   nlist,
                                   volume,
                                   bond_geometry_N,
                                   omega,
                                   bond_damage)

    datamanager.init_stiffness_matrix(index_x, index_y, vals, total_dof)

    create_zero_energy_mode_stiffness!(nodes, dof, C_voigt,
                                       inverse_shape_tensor, zStiff)

    K_sparse = datamanager.get_stiffness_matrix()
    add_zero_energy_stiff!(K_sparse,
                           nodes,
                           dof,
                           zStiff,
                           inverse_shape_tensor,
                           nlist,
                           volume,
                           bond_geometry_N,
                           bond_damage,
                           omega)

    return K_sparse
end

end # module
