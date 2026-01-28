# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Correspondence_matrix_based
using Base.Threads
using LinearAlgebra
using LoopVectorization: @turbo
using SparseArrays
using StaticArrays: @MMatrix
using TimerOutputs: @timeit
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
    @timeit "zero_energy_stiff" begin
        nnodes = Data_Manager.get_nnodes()
        n_total = nnodes * dof

        # Pre-allocate with estimated size (avg ~20 neighbors * 4 entries per bond * nnodes)
        @timeit "initialization" begin
            est_size = nnodes * 80 * dof * dof
            I_indices = sizehint!(Int[], est_size)
            J_indices = sizehint!(Int[], est_size)
            values = sizehint!(Float64[], est_size)

            # Z_i = C : D_i^(-1) - use views instead of copies
            Z = Vector{SubArray{Float64,2}}(undef, maximum(active_nodes))
            for i in active_nodes
                Z[i] = @view zStiff[i, :, :]
            end
        end

        # Reusable buffers to avoid allocations
        @timeit "buffer_allocation" begin
            K_block = zeros(dof, dof)
            D_inv_X = Vector{Float64}(undef, dof)
        end

        # =================================================================
        # FORWARD bonds assembly
        # =================================================================

        @timeit "forward_assembly" begin
            @inbounds for i in active_nodes
                neighbors_i = nlist[i]
                D_inv_i = @view inverse_shape_tensor[i, :, :]
                Z_i = Z[i]
                n_neighbors = length(neighbors_i)

                for (j_idx, j) in enumerate(neighbors_i)
                    if bond_damage[i][j_idx] == 0.0
                        continue
                    end

                    ω_ij = omega[i][j_idx] * bond_damage[i][j_idx]
                    V_j = volume[j]
                    X_ij = bond_geometry[i][j_idx]

                    # Precompute D_inv * X_ij once - reuse buffer
                    mul!(D_inv_X, D_inv_i, X_ij)

                    # Compute γ_ij = 1 - α_ijj
                    α_ijj = ω_ij * V_j * dot(X_ij, D_inv_X)
                    γ_ij = 1.0 - α_ijj

                    # ================================================
                    # Term 1: γ_ij * U_ij - direct assembly
                    # ================================================

                    factor1 = ω_ij * γ_ij * V_j

                    for m in 1:dof, o in 1:dof
                        val = factor1 * Z_i[m, o]
                        if abs(val) > 1e-14
                            row = get_dof_index_block_style(i, m, nnodes, dof)
                            col_i = get_dof_index_block_style(i, o, nnodes, dof)
                            col_j = get_dof_index_block_style(j, o, nnodes, dof)

                            push!(I_indices, row)
                            push!(J_indices, col_i)
                            push!(values, val)

                            push!(I_indices, row)
                            push!(J_indices, col_j)
                            push!(values, -val)
                        end
                    end

                    # ================================================
                    # Term 2: -Σ_{p≠j} - direct assembly
                    # ================================================

                    for (p_idx, p) in enumerate(neighbors_i)
                        if bond_damage[i][p_idx] == 0.0 || p == j
                            continue
                        end

                        X_ip = bond_geometry[i][p_idx]
                        ω_ip = omega[i][p_idx] * bond_damage[i][p_idx]
                        V_p = volume[p]

                        coeff = dot(X_ip, D_inv_X)
                        factor2 = -ω_ij * V_j * ω_ip * V_p * coeff

                        for m in 1:dof, o in 1:dof
                            val = factor2 * Z_i[m, o]
                            if abs(val) > 1e-14
                                row = get_dof_index_block_style(i, m, nnodes, dof)
                                col_i = get_dof_index_block_style(i, o, nnodes, dof)
                                col_p = get_dof_index_block_style(p, o, nnodes, dof)

                                push!(I_indices, row)
                                push!(J_indices, col_i)
                                push!(values, val)

                                push!(I_indices, row)
                                push!(J_indices, col_p)
                                push!(values, -val)
                            end
                        end
                    end
                end
            end
        end

        # =================================================================
        # BACKWARD bonds assembly
        # =================================================================

        @timeit "backward_assembly" begin
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

                    # Reuse buffer
                    mul!(D_inv_X, D_inv_j, X_ji)

                    # Compute γ_ji
                    α_jii = ω_ji * V_i * dot(X_ji, D_inv_X)
                    γ_ji = 1.0 - α_jii

                    # Term 1: Direct assembly
                    factor1 = -ω_ji * γ_ji * V_i

                    for m in 1:dof, o in 1:dof
                        val = factor1 * Z_j[m, o]
                        if abs(val) > 1e-14
                            row = get_dof_index_block_style(i, m, nnodes, dof)
                            col_j = get_dof_index_block_style(j, o, nnodes, dof)
                            col_i = get_dof_index_block_style(i, o, nnodes, dof)

                            push!(I_indices, row)
                            push!(J_indices, col_j)
                            push!(values, val)

                            push!(I_indices, row)
                            push!(J_indices, col_i)
                            push!(values, -val)
                        end
                    end

                    # Term 2: Direct assembly
                    for (p_idx, p) in enumerate(neighbors_j)
                        if bond_damage[j][p_idx] == 0.0 || p == i
                            continue
                        end

                        X_jp = bond_geometry[j][p_idx]
                        ω_jp = omega[j][p_idx] * bond_damage[j][p_idx]
                        V_p = volume[p]

                        coeff = dot(X_jp, D_inv_X)
                        factor2 = ω_ji * V_i * ω_jp * V_p * coeff

                        for m in 1:dof, o in 1:dof
                            val = factor2 * Z_j[m, o]
                            if abs(val) > 1e-14
                                row = get_dof_index_block_style(i, m, nnodes, dof)
                                col_j = get_dof_index_block_style(j, o, nnodes, dof)
                                col_p = get_dof_index_block_style(p, o, nnodes, dof)

                                push!(I_indices, row)
                                push!(J_indices, col_j)
                                push!(values, val)

                                push!(I_indices, row)
                                push!(J_indices, col_p)
                                push!(values, -val)
                            end
                        end
                    end
                end
            end
        end

        # =================================================================
        # Sparse assembly
        # =================================================================

        @timeit "sparse_assembly" begin
            K_stab = sparse(I_indices, J_indices, values, size(K)...)
            K_combined = K + K_stab
        end

        @timeit "verification" begin
            expected_rank = nnodes * dof - Int(dof * (dof + 1) / 2)
            actual_rank = rank(K_combined)

            if actual_rank != expected_rank
                @warn "Stiffness matrix rank is $actual_rank, expected $expected_rank"
            else
                @info "✓ Zero-energy stabilization: rank = $actual_rank"
            end
        end

        Data_Manager.set_stiffness_matrix(-K_combined)
    end

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

"""
    assemble_stiffness_with_zero_energy(nodes, active_nodes, dof, C_Voigt,
                                        inverse_shape_tensor, number_of_neighbors,
                                        nlist, volume, bond_geometry, omega,
                                        bond_damage, zStiff,
                                        include_zero_energy::Bool,
                                        use_block_style::Bool = true)

Assemble stiffness matrix with optional zero-energy stabilization in a single pass.
Memory-optimized by combining both contributions into a single triplet list assembly.

# Arguments
- `nodes`: All node IDs in the system
- `active_nodes`: Subset of nodes to assemble (typically active/damaged nodes)
- `dof`: Degrees of freedom per node (2 or 3)
- `C_Voigt`: Elasticity matrix in Voigt notation [nnodes × dof × dof]
- `inverse_shape_tensor`: Inverse shape tensor for each node [nnodes × dof × dof]
- `number_of_neighbors`: Number of neighbors per node
- `nlist`: Neighbor list for each node
- `volume`: Nodal volumes
- `bond_geometry`: Bond geometry vectors
- `omega`: Influence function values
- `bond_damage`: Bond damage values (0 = broken, 1 = intact)
- `zStiff`: Zero-energy stiffness tensor Z = C : D⁻¹ [nnodes × dof × dof]
- `include_zero_energy`: If true, adds zero-energy mode stabilization
- `use_block_style`: If true, uses block-style DOF indexing, otherwise interleaved

# Returns
- `I_indices, J_indices, values, n_total`: Triplet lists for sparse matrix assembly

# Memory optimization strategy
- Single triplet list for both contributions
- Reusable buffers for tensor operations
- Views instead of copies for submatrices
- Pre-allocated CB tensors shared between forward/backward passes
- Direct assembly without intermediate K_ij matrices for zero-energy terms
"""

function assemble_stiffness_with_zero_energy(nodes::AbstractVector{Int64},
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
                                             zStiff::Union{Array{Float64,3},Nothing} = nothing,
                                             include_zero_energy::Bool = false,
                                             use_block_style::Bool = true)
    nnodes = maximum(nodes)
    n_total = nnodes * dof

    # Validate zero-energy input
    if include_zero_energy && (zStiff === nothing)
        @error "Zero-energy stabilization requested but zStiff not provided"
        return Int[], Int[], Float64[], n_total
    end

    # Calculate accurate size estimate
    est_size = sum((length(nlist[i]) + 1)^2 * dof^2 for i in active_nodes)

    # Use Dictionary to accumulate entries - automatically handles duplicates
    stiffness_dict = Dict{Tuple{Int,Int},Float64}()
    sizehint!(stiffness_dict, est_size)

    # =================================================================
    # PRECOMPUTATION PHASE
    # =================================================================

    # Precompute CB tensors for normal stiffness
    max_neighbors = maximum(number_of_neighbors)
    all_CB_tensors = [zeros(max_neighbors, dof, dof, dof) for _ in 1:nnodes]

    @timeit "precompute_CB_tensors" begin
        @inbounds for i in nodes
            D_inv = @view inverse_shape_tensor[i, :, :]
            C_tensor = get_fourth_order(@view(C_Voigt[i, :, :]), dof)
            ni = nlist[i]

            B_temp = dof == 2 ? zeros(2, 2, 2) : zeros(3, 3, 3)
            CB_k = @view all_CB_tensors[i][eachindex(ni), :, :, :]
            @views precompute_CB_tensors!(CB_k, C_tensor, D_inv, ni, volume,
                                          bond_geometry[i], omega[i], bond_damage[i], dof,
                                          B_temp)
        end
    end

    @timeit "precompute_Z_views" begin
        # Precompute Z_i views for zero-energy (if needed)
        Z = if include_zero_energy
            Z_views = Vector{Union{SubArray{Float64,2},Nothing}}(nothing, nnodes)
            for i in active_nodes
                Z_views[i] = @view zStiff[i, :, :]
            end
            Z_views
        else
            nothing
        end
    end

    # Reusable buffers to minimize allocations
    K_block = zeros(dof, dof)
    D_inv_X = Vector{Float64}(undef, dof)

    # Helper function to add/accumulate entry
    @inline function add_entry!(row::Int, col::Int, val::Float64)
        if abs(val) > 1e-14
            key = (row, col)
            stiffness_dict[key] = get(stiffness_dict, key, 0.0) + val
        end
    end

    # =================================================================
    # FORWARD BONDS ASSEMBLY - Memory optimized
    # =================================================================
    @timeit "forward_assembly" begin
        @inbounds for i in active_nodes
            D_inv_i = @view inverse_shape_tensor[i, :, :]
            V_i = volume[i]
            ni = nlist[i]
            CB_k_i = @view all_CB_tensors[i][eachindex(ni), :, :, :]
            Z_i = include_zero_energy ? Z[i] : nothing

            # Loop over neighbors j
            for (j_idx, j) in enumerate(ni)
                if bond_damage[i][j_idx] == 0.0
                    continue
                end

                ω_ij = omega[i][j_idx] * bond_damage[i][j_idx]
                V_j = volume[j]
                X_ij = bond_geometry[i][j_idx]

                # ============================================================
                # NORMAL STIFFNESS CONTRIBUTION - Direct assembly
                # ============================================================

                for (k_idx, k) in enumerate(ni)
                    if bond_damage[i][k_idx] == 0.0
                        continue
                    end

                    ω_ik = omega[i][k_idx] * bond_damage[i][k_idx]
                    V_k = volume[k]

                    compute_stiffness_contribution(@view(CB_k_i[k_idx, :, :, :]),
                                                   D_inv_i, X_ij, ω_ij, 1.0, ω_ik, V_k, dof,
                                                   K_block)

                    # Direct assembly via dictionary
                    factor = V_j
                    for m in 1:dof
                        row = get_dof_index_block_style(i, m, nnodes, dof)
                        for o in 1:dof
                            val_block = K_block[m, o] * factor
                            if abs(val_block) > 1e-14
                                col_i = get_dof_index_block_style(i, o, nnodes, dof)
                                col_k = get_dof_index_block_style(k, o, nnodes, dof)

                                add_entry!(row, col_i, val_block)
                                add_entry!(row, col_k, -val_block)
                            end
                        end
                    end
                end

                # ============================================================
                # ZERO-ENERGY CONTRIBUTION (if enabled)
                # ============================================================

                if include_zero_energy
                    # Precompute D_inv * X_ij
                    mul!(D_inv_X, D_inv_i, X_ij)

                    # Compute γ_ij = 1 - α_ijj
                    α_ijj = ω_ij * V_j * dot(X_ij, D_inv_X)
                    γ_ij = 1.0 - α_ijj

                    # Term 1: γ_ij * U_ij - direct assembly
                    factor1 = ω_ij * γ_ij * V_j

                    for m in 1:dof
                        row = get_dof_index_block_style(i, m, nnodes, dof)
                        for o in 1:dof
                            val = factor1 * Z_i[m, o]
                            if abs(val) > 1e-14
                                col_i = get_dof_index_block_style(i, o, nnodes, dof)
                                col_j = get_dof_index_block_style(j, o, nnodes, dof)

                                add_entry!(row, col_i, val)
                                add_entry!(row, col_j, -val)
                            end
                        end
                    end

                    # Term 2: -Σ_{p≠j} coupling terms
                    for (p_idx, p) in enumerate(ni)
                        if bond_damage[i][p_idx] == 0.0 || p == j
                            continue
                        end

                        X_ip = bond_geometry[i][p_idx]
                        ω_ip = omega[i][p_idx] * bond_damage[i][p_idx]
                        V_p = volume[p]

                        coeff = dot(X_ip, D_inv_X)
                        factor2 = -ω_ij * V_j * ω_ip * V_p * coeff

                        for m in 1:dof
                            row = get_dof_index_block_style(i, m, nnodes, dof)
                            for o in 1:dof
                                val = factor2 * Z_i[m, o]
                                if abs(val) > 1e-14
                                    col_i = get_dof_index_block_style(i, o, nnodes, dof)
                                    col_p = get_dof_index_block_style(p, o, nnodes, dof)

                                    add_entry!(row, col_i, val)
                                    add_entry!(row, col_p, -val)
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    # =================================================================
    # BACKWARD BONDS ASSEMBLY - Memory optimized
    # =================================================================
    @timeit "backward_assembly" begin
        @inbounds for i in active_nodes
            V_i = volume[i]
            Z_i = include_zero_energy ? Z[i] : nothing

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
                Z_j = include_zero_energy ? Z[j] : nothing

                # ============================================================
                # NORMAL STIFFNESS CONTRIBUTION - Direct assembly
                # ============================================================

                for (k_idx, k) in enumerate(j_neighbors)
                    if omega[j][k_idx] * bond_damage[j][k_idx] == 0.0
                        continue
                    end

                    ω_jk = omega[j][k_idx] * bond_damage[j][k_idx]
                    V_k = volume[k]

                    compute_stiffness_contribution(@view(CB_k_j[k_idx, :, :, :]),
                                                   D_inv_j, X_ji, ω_ji, 1.0, ω_jk, V_k, dof,
                                                   K_block)

                    # Direct assembly via dictionary
                    factor = -V_i
                    for m in 1:dof
                        row = get_dof_index_block_style(i, m, nnodes, dof)
                        for o in 1:dof
                            val_block = K_block[m, o] * factor
                            if abs(val_block) > 1e-14
                                col_j = get_dof_index_block_style(j, o, nnodes, dof)
                                col_k = get_dof_index_block_style(k, o, nnodes, dof)

                                add_entry!(row, col_j, val_block)
                                add_entry!(row, col_k, -val_block)
                            end
                        end
                    end
                end

                # ============================================================
                # ZERO-ENERGY CONTRIBUTION (if enabled)
                # ============================================================

                if include_zero_energy
                    # Precompute D_inv * X_ji (reuse buffer)
                    mul!(D_inv_X, D_inv_j, X_ji)

                    # Compute γ_ji = 1 - α_jii
                    α_jii = ω_ji * V_i * dot(X_ji, D_inv_X)
                    γ_ji = 1.0 - α_jii

                    # Term 1: backward bond contribution
                    factor1 = -ω_ji * γ_ji * V_i

                    for m in 1:dof
                        row = get_dof_index_block_style(i, m, nnodes, dof)
                        for o in 1:dof
                            val = factor1 * Z_j[m, o]
                            if abs(val) > 1e-14
                                col_j = get_dof_index_block_style(j, o, nnodes, dof)
                                col_i = get_dof_index_block_style(i, o, nnodes, dof)

                                add_entry!(row, col_j, val)
                                add_entry!(row, col_i, -val)
                            end
                        end
                    end

                    # Term 2: coupling with j's other neighbors
                    for (p_idx, p) in enumerate(j_neighbors)
                        if bond_damage[j][p_idx] == 0.0 || p == i
                            continue
                        end

                        X_jp = bond_geometry[j][p_idx]
                        ω_jp = omega[j][p_idx] * bond_damage[j][p_idx]
                        V_p = volume[p]

                        coeff = dot(X_jp, D_inv_X)
                        factor2 = ω_ji * V_i * ω_jp * V_p * coeff

                        for m in 1:dof
                            row = get_dof_index_block_style(i, m, nnodes, dof)
                            for o in 1:dof
                                val = factor2 * Z_j[m, o]
                                if abs(val) > 1e-14
                                    col_j = get_dof_index_block_style(j, o, nnodes, dof)
                                    col_p = get_dof_index_block_style(p, o, nnodes, dof)

                                    add_entry!(row, col_j, val)
                                    add_entry!(row, col_p, -val)
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    # Convert dictionary to triplet format
    @timeit "dict_to_triplets" begin
        n_entries = length(stiffness_dict)
        I_indices = Vector{Int}(undef, n_entries)
        J_indices = Vector{Int}(undef, n_entries)
        values = Vector{Float64}(undef, n_entries)

        idx = 1
        for ((row, col), val) in stiffness_dict
            I_indices[idx] = row
            J_indices[idx] = col
            values[idx] = val
            idx += 1
        end
    end

    return I_indices, J_indices, values, n_total
end

function init_matrix(use_block_style::Bool = true, include_zero_energy::Bool = true)
    nodes = collect(1:Data_Manager.get_nnodes())
    dof = Data_Manager.get_dof()
    nnodes = length(nodes)

    # Get all required fields
    bond_geometry = Data_Manager.get_field("Bond Geometry")
    inverse_shape_tensor = Data_Manager.get_field("Inverse Shape Tensor")
    nlist = Data_Manager.get_nlist()
    number_of_neighbors = Data_Manager.get_field("Number of Neighbors")
    volume = Data_Manager.get_field("Volume")
    omega = Data_Manager.get_field("Influence Function")
    C_voigt = Data_Manager.get_field("Elasticity Matrix")
    bond_damage = Data_Manager.get_field("Bond Damage", "NP1")

    # Check if zero-energy control is requested
    use_zero_energy = false
    zStiff = nothing

    if haskey(Data_Manager.get_properties(1, "Material Model"), "Zero Energy Control")
        if Data_Manager.get_properties(1, "Material Model")["Zero Energy Control"] ==
           "Global"
            use_zero_energy = include_zero_energy
            if use_zero_energy
                zStiff = Data_Manager.create_constant_node_tensor_field("Zero Energy Stiffness",
                                                                        Float64, dof)
                Zero_Energy_Control.create_zero_energy_mode_stiffness!(nodes, dof, C_voigt,
                                                                       inverse_shape_tensor,
                                                                       zStiff)
            end
        end
    else
        if include_zero_energy
            @warn "Zero-energy control requested but not active in material model"
        end
    end

    @timeit "assemble_stiffness_with_zero_energy" begin
        index_x, index_y, vals,
        total_dof = assemble_stiffness_with_zero_energy(nodes,
                                                        nodes, dof,
                                                        C_voigt,
                                                        inverse_shape_tensor,
                                                        number_of_neighbors,
                                                        nlist,
                                                        volume,
                                                        bond_geometry,
                                                        omega,
                                                        bond_damage,
                                                        zStiff,
                                                        use_zero_energy,
                                                        use_block_style)
    end

    Data_Manager.init_stiffness_matrix(index_x, index_y, vals, total_dof)
    K_sparse = Data_Manager.get_stiffness_matrix()

    # Apply negative sign convention
    Data_Manager.set_stiffness_matrix(-K_sparse)

    # Verify rank if zero-energy stabilization was used
    if use_zero_energy
        expected_rank = nnodes * dof - Int(dof * (dof + 1) / 2)
        actual_rank = rank(-K_sparse)

        if actual_rank != expected_rank
            @warn "Stiffness matrix rank is $actual_rank, expected $expected_rank"
        else
            @info "✓ Zero-energy stabilization successful: rank = $actual_rank"
        end
    end
end

function compute_model(nodes::AbstractVector{Int64},
                       use_block_style::Bool = true,
                       include_zero_energy::Bool = true)
    dof::Int64 = Data_Manager.get_dof()
    nnodes = Data_Manager.get_nnodes()

    # Get fields
    C_voigt = Data_Manager.get_field("Elasticity Matrix")
    inverse_shape_tensor = Data_Manager.get_field("Inverse Shape Tensor")
    nlist = Data_Manager.get_nlist()
    volume = Data_Manager.get_field("Volume")
    bond_geometry_N = Data_Manager.get_field("Bond Geometry")
    number_of_neighbors = Data_Manager.get_field("Number of Neighbors")
    omega = Data_Manager.get_field("Influence Function")
    bond_damage = Data_Manager.get_field("Bond Damage", "NP1")

    # Setup zero-energy if needed
    zStiff = include_zero_energy ? Data_Manager.get_field("Zero Energy Stiffness") : nothing

    if include_zero_energy && zStiff !== nothing
        Zero_Energy_Control.create_zero_energy_mode_stiffness!(nodes, dof, C_voigt,
                                                               inverse_shape_tensor, zStiff)
    end

    # Assemble combined stiffness matrix
    @timeit "assemble_stiffness_with_zero_energy update" index_x, index_y, vals,
                                                         total_dof=assemble_stiffness_with_zero_energy(collect(1:nnodes),
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
                                                                                                       zStiff,
                                                                                                       include_zero_energy,
                                                                                                       use_block_style)

    Data_Manager.init_stiffness_matrix(index_x, index_y, vals, total_dof)
    K_sparse = Data_Manager.get_stiffness_matrix()
    Data_Manager.set_stiffness_matrix(-K_sparse)
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
