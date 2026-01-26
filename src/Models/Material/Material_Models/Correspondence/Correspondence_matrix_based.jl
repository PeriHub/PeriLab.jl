# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Correspondence_matrix_based
using Base.Threads
using LinearAlgebra
using LoopVectorization: @turbo
using SparseArrays
using StaticArrays: @MMatrix
using Random, Statistics
using ...Data_Manager
using ...Material_Basis: get_Hooke_matrix
using ....Helpers: get_fourth_order, progress_bar
using ...Zero_Energy_Control

export init_model
export add_zero_energy_stiff!
export compute_bond_force
export init_matrix
export get_dof_index_block_style

function contraction!(C::Array{Float64,4}, B::Array{Float64,3}, dof::Int64,
                      CB::AbstractArray{Float64,3})
    fill!(CB, 0.0)

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
    end
end

function create_B_tensor!(D_inv::AbstractMatrix{Float64}, X_ik::Vector{Float64},
                          V_k::Float64, omega_ik::Float64, ::Val{2},
                          B_tensor::Array{Float64,3})
    factor = 0.5
    X_x = X_ik[1]
    X_y = X_ik[2]
    B1_1 = D_inv[1, 1] * X_x + D_inv[2, 1] * X_y
    B1_2 = D_inv[1, 2] * X_x + D_inv[2, 2] * X_y
    B2_1 = D_inv[1, 1] * X_x + D_inv[1, 2] * X_y
    B2_2 = D_inv[2, 1] * X_x + D_inv[2, 2] * X_y

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
    factor = 0.5

    X_x = X_ik[1]
    X_y = X_ik[2]
    X_z = X_ik[3]

    B1_1 = D_inv[1, 1] * X_x + D_inv[2, 1] * X_y + D_inv[3, 1] * X_z
    B1_2 = D_inv[1, 2] * X_x + D_inv[2, 2] * X_y + D_inv[3, 2] * X_z
    B1_3 = D_inv[1, 3] * X_x + D_inv[2, 3] * X_y + D_inv[3, 3] * X_z
    B2_1 = D_inv[1, 1] * X_x + D_inv[1, 2] * X_y + D_inv[1, 3] * X_z
    B2_2 = D_inv[2, 1] * X_x + D_inv[2, 2] * X_y + D_inv[2, 3] * X_z
    B2_3 = D_inv[3, 1] * X_x + D_inv[3, 2] * X_y + D_inv[3, 3] * X_z

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

function compute_stiffness_contribution(CB_ik::AbstractArray{Float64,3},
                                        D_inv_i::AbstractMatrix{Float64},
                                        X_ij::Vector{Float64},
                                        omega_ij::Float64,
                                        V_j::Float64,
                                        omega_ik::Float64,
                                        V_k::Float64,
                                        dof::Int64,
                                        K_block::AbstractMatrix{Float64})
    if dof == 2
        DX_1 = D_inv_i[1, 1] * X_ij[1] + D_inv_i[1, 2] * X_ij[2]
        DX_2 = D_inv_i[2, 1] * X_ij[1] + D_inv_i[2, 2] * X_ij[2]
        factor = omega_ij * V_j * omega_ik * V_k

        @inbounds begin
            K_block[1, 1] = factor * (CB_ik[1, 1, 1] * DX_1 + CB_ik[1, 2, 1] * DX_2)
            K_block[1, 2] = factor * (CB_ik[1, 1, 2] * DX_1 + CB_ik[1, 2, 2] * DX_2)
            K_block[2, 1] = factor * (CB_ik[2, 1, 1] * DX_1 + CB_ik[2, 2, 1] * DX_2)
            K_block[2, 2] = factor * (CB_ik[2, 1, 2] * DX_1 + CB_ik[2, 2, 2] * DX_2)
        end
    elseif dof == 3
        DX_1 = D_inv_i[1, 1] * X_ij[1] + D_inv_i[1, 2] * X_ij[2] + D_inv_i[1, 3] * X_ij[3]
        DX_2 = D_inv_i[2, 1] * X_ij[1] + D_inv_i[2, 2] * X_ij[2] + D_inv_i[2, 3] * X_ij[3]
        DX_3 = D_inv_i[3, 1] * X_ij[1] + D_inv_i[3, 2] * X_ij[2] + D_inv_i[3, 3] * X_ij[3]
        factor = omega_ij * V_j * omega_ik * V_k

        @inbounds begin
            K_block[1, 1] = factor * (CB_ik[1, 1, 1] * DX_1 + CB_ik[1, 2, 1] * DX_2 +
                             CB_ik[1, 3, 1] * DX_3)
            K_block[1, 2] = factor * (CB_ik[1, 1, 2] * DX_1 + CB_ik[1, 2, 2] * DX_2 +
                             CB_ik[1, 3, 2] * DX_3)
            K_block[1, 3] = factor * (CB_ik[1, 1, 3] * DX_1 + CB_ik[1, 2, 3] * DX_2 +
                             CB_ik[1, 3, 3] * DX_3)
            K_block[2, 1] = factor * (CB_ik[2, 1, 1] * DX_1 + CB_ik[2, 2, 1] * DX_2 +
                             CB_ik[2, 3, 1] * DX_3)
            K_block[2, 2] = factor * (CB_ik[2, 1, 2] * DX_1 + CB_ik[2, 2, 2] * DX_2 +
                             CB_ik[2, 3, 2] * DX_3)
            K_block[2, 3] = factor * (CB_ik[2, 1, 3] * DX_1 + CB_ik[2, 2, 3] * DX_2 +
                             CB_ik[2, 3, 3] * DX_3)
            K_block[3, 1] = factor * (CB_ik[3, 1, 1] * DX_1 + CB_ik[3, 2, 1] * DX_2 +
                             CB_ik[3, 3, 1] * DX_3)
            K_block[3, 2] = factor * (CB_ik[3, 1, 2] * DX_1 + CB_ik[3, 2, 2] * DX_2 +
                             CB_ik[3, 3, 2] * DX_3)
            K_block[3, 3] = factor * (CB_ik[3, 1, 3] * DX_1 + CB_ik[3, 2, 3] * DX_2 +
                             CB_ik[3, 3, 3] * DX_3)
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

@inline function get_dof_index_block_style(node_id::Int, dof_component::Int,
                                           nnodes::Int, dof::Int)
    return (dof_component - 1) * nnodes + node_id
end

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
                            bond_damage::Vector{Vector{Float64}},
                            use_block_style::Bool = true)
    nnodes = maximum(nodes)
    n_total = nnodes * dof

    I_indices = Int[]
    J_indices = Int[]
    values = Float64[]

    # Precompute CB tensors
    max_neighbors = maximum(number_of_neighbors)
    all_CB_tensors = [zeros(max_neighbors, dof, dof, dof) for _ in nodes]

    for i in nodes
        D_inv = @view inverse_shape_tensor[i, :, :]
        C_tensor = get_fourth_order(@view(C_Voigt[i, :, :]), dof)
        ni = nlist[i]

        B_temp = dof == 2 ? zeros(2, 2, 2) : zeros(3, 3, 3)
        CB_k = @view all_CB_tensors[i][eachindex(ni), :, :, :]
        @views precompute_CB_tensors!(CB_k, C_tensor, D_inv, ni, volume,
                                      bond_geometry[i], omega[i], bond_damage[i], dof,
                                      B_temp)
    end

    # FORWARD bonds assembly - exactly like development code
    for i in active_nodes
        D_inv_i = @view inverse_shape_tensor[i, :, :]
        V_i = volume[i]
        ni = nlist[i]
        CB_k_i = @view all_CB_tensors[i][eachindex(ni), :, :, :]

        # For each neighbor j, build K_ij matrix
        for (j_idx, j) in enumerate(ni)
            if bond_damage[i][j_idx] == 0.0
                continue
            end

            ω_ij = omega[i][j_idx] * bond_damage[i][j_idx]
            V_j = volume[j]
            X_ij = bond_geometry[i][j_idx]

            # Build local K_ij matrix for this bond
            K_ij = zeros(dof, nnodes * dof)
            K_block = zeros(dof, dof)

            for (k_idx, k) in enumerate(ni)
                if bond_damage[i][k_idx] == 0.0
                    continue
                end

                ω_ik = omega[i][k_idx] * bond_damage[i][k_idx]
                V_k = volume[k]

                compute_stiffness_contribution(@view(CB_k_i[k_idx, :, :, :]),
                                               D_inv_i, X_ij, ω_ij, 1.0,
                                               ω_ik, V_k, dof, K_block)

                # Fill K_ij with +/- pattern
                for m in 1:dof, o in 1:dof
                    if use_block_style
                        col_i = get_dof_index_block_style(i, o, nnodes, dof)
                        col_k = get_dof_index_block_style(k, o, nnodes, dof)
                    else
                        col_i = get_dof_index_interleaved(i, o, dof)
                        col_k = get_dof_index_interleaved(k, o, dof)
                    end

                    K_ij[m, col_i] += K_block[m, o]
                    K_ij[m, col_k] -= K_block[m, o]
                end
            end

            # Add to global matrix (scaled by V_j)
            for m in 1:dof
                if use_block_style
                    row = get_dof_index_block_style(i, m, nnodes, dof)
                else
                    row = get_dof_index_interleaved(i, m, dof)
                end

                for n_col in 1:n_total
                    val = K_ij[m, n_col] * V_j
                    if abs(val) > 1e-14
                        push!(I_indices, row)
                        push!(J_indices, n_col)
                        push!(values, val)
                    end
                end
            end
        end
    end

    # BACKWARD bonds assembly - exactly like development code
    for i in active_nodes
        V_i = volume[i]
        ni = nlist[i]

        for j in active_nodes
            if j == i
                continue
            end

            j_neighbors = nlist[j]
            i_pos_in_j = findfirst(==(i), j_neighbors)

            if i_pos_in_j === nothing ||
               omega[j][i_pos_in_j] * bond_damage[j][i_pos_in_j] == 0.0
                continue
            end

            ω_ji = omega[j][i_pos_in_j] * bond_damage[j][i_pos_in_j]
            X_ji = bond_geometry[j][i_pos_in_j]
            D_inv_j = @view inverse_shape_tensor[j, :, :]
            CB_k_j = @view all_CB_tensors[j][eachindex(j_neighbors), :, :, :]

            K_ji = zeros(dof, nnodes * dof)
            K_block = zeros(dof, dof)

            for (k_idx, k) in enumerate(j_neighbors)
                if omega[j][k_idx] * bond_damage[j][k_idx] == 0.0
                    continue
                end

                ω_jk = omega[j][k_idx] * bond_damage[j][k_idx]
                V_k = volume[k]

                compute_stiffness_contribution(@view(CB_k_j[k_idx, :, :, :]),
                                               D_inv_j, X_ji, ω_ji, 1.0,
                                               ω_jk, V_k, dof, K_block)

                for m in 1:dof, o in 1:dof
                    if use_block_style
                        col_j = get_dof_index_block_style(j, o, nnodes, dof)
                        col_k = get_dof_index_block_style(k, o, nnodes, dof)
                    else
                        col_j = get_dof_index_interleaved(j, o, dof)
                        col_k = get_dof_index_interleaved(k, o, dof)
                    end

                    K_ji[m, col_j] += K_block[m, o]
                    K_ji[m, col_k] -= K_block[m, o]
                end
            end

            # Subtract from global matrix (scaled by V_i)
            for m in 1:dof
                if use_block_style
                    row = get_dof_index_block_style(i, m, nnodes, dof)
                else
                    row = get_dof_index_interleaved(i, m, dof)
                end

                for n_col in 1:n_total
                    val = -K_ji[m, n_col] * V_i
                    if abs(val) > 1e-14
                        push!(I_indices, row)
                        push!(J_indices, n_col)
                        push!(values, val)
                    end
                end
            end
        end
    end

    return I_indices, J_indices, values, n_total
end

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

function contract_C_Dinv(C::Array{Float64,4}, D_inv::AbstractMatrix{Float64}, d::Int)
    Z = zeros(d, d)
    for m in 1:d, n in 1:d
        for o in 1:d, p in 1:d
            Z[m, n] += C[m, n, o, p] * D_inv[o, p]
        end
    end
    return Z
end

function add_zero_energy_stiff!(K::SparseMatrixCSC{Float64,Int64},
                                active_nodes::AbstractVector{Int64},
                                dof::Int64,
                                zStiff::Array{Float64,3},
                                inverse_shape_tensor::Array{Float64,3},
                                nlist::Vector{Vector{Int}},
                                volume::Vector{Float64},
                                bond_geometry::Vector{Vector{Vector{Float64}}},
                                bond_damage::Vector{Vector{Float64}},
                                omega::Vector{Vector{Float64}},
                                use_block_style::Bool = true)
    nnodes = Data_Manager.get_nnodes()
    n_total = nnodes * dof

    I_indices = Int[]
    J_indices = Int[]
    values = Float64[]

    # Z_i = C : D_i^(-1)
    Z = Vector{Matrix{Float64}}(undef, maximum(active_nodes))
    for i in active_nodes
        Z[i] = @view zStiff[i, :, :]
    end

    # =================================================================
    # FORWARD bonds assembly
    # =================================================================

    @inbounds for i in active_nodes
        neighbors_i = nlist[i]
        D_inv_i = @view inverse_shape_tensor[i, :, :]
        Z_i = Z[i]

        for (j_idx, j) in enumerate(neighbors_i)
            if bond_damage[i][j_idx] == 0.0
                continue
            end

            ω_ij = omega[i][j_idx] * bond_damage[i][j_idx]
            V_j = volume[j]
            X_ij = bond_geometry[i][j_idx]

            K_S_ij = zeros(dof, n_total)
            K_block = zeros(dof, dof)

            D_inv_X_ij = D_inv_i * X_ij

            # ================================================
            # Compute γ_ij = 1 - α_ijj
            # where α_ijj = ω_ij * V_j * (X_ij^T * D_inv * X_ij)
            # ================================================

            α_ijj = ω_ij * V_j * dot(X_ij, D_inv_X_ij)
            γ_ij = 1.0 - α_ijj

            # ================================================
            # Term 1: γ_ij * U_ij contribution
            # z_ij has: +γ_ij * U_ij = +γ_ij * (u_j - u_i)
            # ================================================

            K_block = ω_ij * γ_ij * Z_i

            for m in 1:dof, o in 1:dof
                if use_block_style
                    col_i = get_dof_index_block_style(i, o, nnodes, dof)
                    col_j = get_dof_index_block_style(j, o, nnodes, dof)
                else
                    col_i = get_dof_index_interleaved(i, o, dof)
                    col_j = get_dof_index_interleaved(j, o, dof)
                end

                K_S_ij[m, col_i] += K_block[m, o]  # -u_i
                K_S_ij[m, col_j] -= K_block[m, o]  # +u_j
            end

            # ================================================
            # Term 2: -Σ_{p≠j} [ω_ip * V_p * (X_ip^T D_inv X_ij) * U_ip]
            # ================================================

            for (p_idx, p) in enumerate(neighbors_i)
                if bond_damage[i][p_idx] == 0.0 || p == j
                    continue
                end

                X_ip = bond_geometry[i][p_idx]
                ω_ip = omega[i][p_idx] * bond_damage[i][p_idx]
                V_p = volume[p]

                coeff = dot(X_ip, D_inv_X_ij)
                K_block = -ω_ij * ω_ip * V_p * coeff * Z_i

                for m in 1:dof, o in 1:dof
                    if use_block_style
                        col_i = get_dof_index_block_style(i, o, nnodes, dof)
                        col_p = get_dof_index_block_style(p, o, nnodes, dof)
                    else
                        col_i = get_dof_index_interleaved(i, o, dof)
                        col_p = get_dof_index_interleaved(p, o, dof)
                    end

                    K_S_ij[m, col_i] += K_block[m, o]
                    K_S_ij[m, col_p] -= K_block[m, o]
                end
            end

            # Add to global matrix (scaled by V_j)
            for m in 1:dof
                if use_block_style
                    row = get_dof_index_block_style(i, m, nnodes, dof)
                else
                    row = get_dof_index_interleaved(i, m, dof)
                end

                for n_col in 1:n_total
                    val = K_S_ij[m, n_col] * V_j
                    if abs(val) > 1e-14
                        push!(I_indices, row)
                        push!(J_indices, n_col)
                        push!(values, val)
                    end
                end
            end
        end
    end

    # =================================================================
    # BACKWARD bonds assembly
    # =================================================================

    @inbounds for i in active_nodes
        V_i = volume[i]

        for j in active_nodes
            if j == i
                continue
            end

            neighbors_j = nlist[j]
            i_pos_in_j = findfirst(==(i), neighbors_j)

            if i_pos_in_j === nothing || bond_damage[j][i_pos_in_j] == 0.0
                continue
            end

            D_inv_j = @view inverse_shape_tensor[j, :, :]
            Z_j = Z[j]
            ω_ji = omega[j][i_pos_in_j] * bond_damage[j][i_pos_in_j]
            X_ji = bond_geometry[j][i_pos_in_j]

            K_S_ji = zeros(dof, n_total)
            K_block = zeros(dof, dof)

            D_inv_X_ji = D_inv_j * X_ji

            # ================================================
            # Compute γ_ji = 1 - α_jii
            # ================================================

            α_jii = ω_ji * V_i * dot(X_ji, D_inv_X_ji)
            γ_ji = 1.0 - α_jii

            # ================================================
            # Term 1: γ_ji * U_ji contribution
            # ================================================

            K_block = ω_ji * γ_ji * Z_j

            for m in 1:dof, o in 1:dof
                if use_block_style
                    col_j = get_dof_index_block_style(j, o, nnodes, dof)
                    col_i = get_dof_index_block_style(i, o, nnodes, dof)
                else
                    col_j = get_dof_index_interleaved(j, o, dof)
                    col_i = get_dof_index_interleaved(i, o, dof)
                end

                K_S_ji[m, col_j] += K_block[m, o]
                K_S_ji[m, col_i] -= K_block[m, o]
            end

            # ================================================
            # Term 2: -Σ_{p≠i} [ω_jp * V_p * (X_jp^T D_inv X_ji) * U_jp]
            # ================================================

            for (p_idx, p) in enumerate(neighbors_j)
                if bond_damage[j][p_idx] == 0.0 || p == i
                    continue
                end

                X_jp = bond_geometry[j][p_idx]
                ω_jp = omega[j][p_idx] * bond_damage[j][p_idx]
                V_p = volume[p]

                coeff = dot(X_jp, D_inv_X_ji)
                K_block = -ω_ji * ω_jp * V_p * coeff * Z_j

                for m in 1:dof, o in 1:dof
                    if use_block_style
                        col_j = get_dof_index_block_style(j, o, nnodes, dof)
                        col_p = get_dof_index_block_style(p, o, nnodes, dof)
                    else
                        col_j = get_dof_index_interleaved(j, o, dof)
                        col_p = get_dof_index_interleaved(p, o, dof)
                    end

                    K_S_ji[m, col_j] += K_block[m, o]
                    K_S_ji[m, col_p] -= K_block[m, o]
                end
            end

            # Subtract from global matrix
            for m in 1:dof
                if use_block_style
                    row = get_dof_index_block_style(i, m, nnodes, dof)
                else
                    row = get_dof_index_interleaved(i, m, dof)
                end

                for n_col in 1:n_total
                    val = -K_S_ji[m, n_col] * V_i
                    if abs(val) > 1e-14
                        push!(I_indices, row)
                        push!(J_indices, n_col)
                        push!(values, val)
                    end
                end
            end
        end
    end

    # =================================================================
    # Assemble and verify
    # =================================================================

    K_stab = sparse(I_indices, J_indices, values, size(K)...)
    K_combined = K + K_stab

    # Check eigenvalues
    λ_stab = eigvals(Matrix(K_stab))
    small_eigs_stab = sort(abs.(λ_stab))[1:5]

    @info "K_stab: 5 smallest |eigenvalues|:"
    for (idx, val) in enumerate(small_eigs_stab)
        @info "  λ[$idx] = $val"
    end

    expected_rank = nnodes * dof - Int(dof * (dof + 1) / 2)
    actual_rank = rank(K_combined)

    if actual_rank != expected_rank
        @warn "Stiffness matrix rank is $actual_rank, expected $expected_rank"
    else
        @info "✓ Stiffness matrix assembly successful: rank = $actual_rank"
    end

    Data_Manager.set_stiffness_matrix(-K_combined)
    return nothing
end

function compute_bond_force(bond_force::Vector{Vector{Vector{Float64}}},
                            K::SparseMatrixCSC{Float64,Int64},
                            u::Matrix{Float64},
                            nodes::AbstractVector{Int64},
                            nlist::Vector{Vector{Int64}},
                            dof::Int64,
                            use_block_style::Bool = true)
    nnodes = size(u, 1)

    @inbounds for i in nodes
        ni = nlist[i]

        for (j_idx, j) in enumerate(ni)
            f_ij = zeros(dof)

            for m in 1:dof
                if use_block_style
                    i_idx = get_dof_index_block_style(i, m, nnodes, dof)
                else
                    i_idx = get_dof_index_interleaved(i, m, dof)
                end

                for n in 1:dof
                    if use_block_style
                        j_idx_global = get_dof_index_block_style(j, n, nnodes, dof)
                    else
                        j_idx_global = get_dof_index_interleaved(j, n, dof)
                    end

                    f_ij[m] += K[i_idx, j_idx_global] * (u[j, n] - u[i, n])
                end
            end

            bond_force[i][j_idx] = f_ij
        end
    end
end

function init_model(nodes::AbstractVector{Int64},
                    material_parameter::Dict, block_id::Int64)
    Zero_Energy_Control.init_model(nodes, material_parameter, block_id)
    if Data_Manager.get_max_rank() > 1
        @error "Correspondence matrix based not implemented for parallel runs."
    end
end

function init_matrix(use_block_style::Bool = true)
    nodes = collect(1:Data_Manager.get_nnodes())
    dof = Data_Manager.get_dof()
    nnodes = length(nodes)
    zStiff = Data_Manager.create_constant_node_tensor_field("Zero Energy Stiffness",
                                                            Float64,
                                                            dof)
    bond_geometry = Data_Manager.get_field("Bond Geometry")
    inverse_shape_tensor = Data_Manager.get_field("Inverse Shape Tensor")
    nlist = Data_Manager.get_nlist()
    number_of_neighbors = Data_Manager.get_field("Number of Neighbors")
    volume = Data_Manager.get_field("Volume")
    omega = Data_Manager.get_field("Influence Function")
    C_voigt = Data_Manager.get_field("Elasticity Matrix")
    bond_damage = Data_Manager.get_field("Bond Damage", "NP1")

    @info "Initializing stiffness matrix with $(use_block_style ? "block-style" : "interleaved") layout"
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
                                   bond_damage,
                                   use_block_style)

    Data_Manager.init_stiffness_matrix(index_x, index_y, vals, total_dof)
    K_sparse = Data_Manager.get_stiffness_matrix()

    if haskey(Data_Manager.get_properties(1, "Material Model"), "Zero Energy Control")
        if Data_Manager.get_properties(1, "Material Model")["Zero Energy Control"] ==
           "Global"
            Zero_Energy_Control.create_zero_energy_mode_stiffness!(nodes, dof, C_voigt,
                                                                   inverse_shape_tensor,
                                                                   zStiff)
            add_zero_energy_stiff!(K_sparse,
                                   nodes,
                                   dof,
                                   zStiff,
                                   inverse_shape_tensor,
                                   nlist,
                                   volume,
                                   bond_geometry,
                                   bond_damage,
                                   omega,
                                   use_block_style)
        end
    else
        @warn "Global Energy Control Model is not active for block 1. Must be active for all blocks when used."
    end
end

function compute_model(nodes::AbstractVector{Int64}, use_block_style::Bool = true)
    dof::Int64 = Data_Manager.get_dof()
    nnodes = Data_Manager.get_nnodes()
    C_voigt::NodeTensorField{Float64} = Data_Manager.get_field("Elasticity Matrix")
    inverse_shape_tensor::NodeTensorField{Float64} = Data_Manager.get_field("Inverse Shape Tensor")
    nlist::BondScalarState{Int64} = Data_Manager.get_nlist()
    volume::NodeScalarField{Float64} = Data_Manager.get_field("Volume")
    bond_geometry_N = Data_Manager.get_field("Bond Geometry")
    number_of_neighbors::NodeScalarField{Int64} = Data_Manager.get_field("Number of Neighbors")
    omega::BondScalarState{Float64} = Data_Manager.get_field("Influence Function")
    bond_damage::BondScalarState{Float64} = Data_Manager.get_field("Bond Damage", "NP1")
    zStiff::NodeTensorField{Float64} = Data_Manager.get_field("Zero Energy Stiffness")

    index_x, index_y, vals,
    total_dof = assemble_stiffness(collect(1:Data_Manager.get_nnodes()),
                                   nodes,
                                   dof,
                                   C_voigt,
                                   inverse_shape_tensor,
                                   number_of_neighbors,
                                   nlist,
                                   volume,
                                   bond_geometry_N,
                                   omega,
                                   bond_damage,
                                   use_block_style)

    Data_Manager.init_stiffness_matrix(index_x, index_y, vals, total_dof)

    Zero_Energy_Control.create_zero_energy_mode_stiffness!(nodes, dof, C_voigt,
                                                           inverse_shape_tensor, zStiff)
    K_sparse = Data_Manager.get_stiffness_matrix()

    add_zero_energy_stiff!(K_sparse,
                           nodes,
                           dof,
                           zStiff,
                           inverse_shape_tensor,
                           nlist,
                           volume,
                           bond_geometry_N,
                           bond_damage,
                           omega,
                           use_block_style)

    Data_Manager.set_stiffness_matrix(K_sparse)
end

function voigt_to_tensor(sigma_voigt::Vector{Float64}, dof::Int64)
    if dof == 2
        return [sigma_voigt[1] sigma_voigt[3];
                sigma_voigt[3] sigma_voigt[2]]
    else
        return [sigma_voigt[1] sigma_voigt[6] sigma_voigt[5];
                sigma_voigt[6] sigma_voigt[2] sigma_voigt[4];
                sigma_voigt[5] sigma_voigt[4] sigma_voigt[3]]
    end
end

end # module
