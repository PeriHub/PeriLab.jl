# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Correspondence_matrix_based
using LinearAlgebra
using LoopVectorization: @turbo
using SparseArrays
using ExtendableSparse
using TimerOutputs: @timeit
using ...Data_Manager
using ....Helpers: get_fourth_order
using ...Zero_Energy_Control

export init_model
export compute_model
export init_matrix

# =============================================================================
# Type-stable helper functions
# =============================================================================

@inline function get_dof_index_block_style(node_id::Int, dof_component::Int,
                                           nnodes::Int, ::Int)
    return (dof_component - 1) * nnodes + node_id
end

# =============================================================================
# Tensor computation kernels
# =============================================================================

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

function precompute_CB_all_nodes!(all_CB_tensors::Vector{Array{Float64,4}},
                                  nodes::AbstractVector{Int64},
                                  C_Voigt::Array{Float64,3},
                                  inverse_shape_tensor::Array{Float64,3},
                                  nlist::Vector{Vector{Int64}},
                                  volume::Vector{Float64},
                                  bond_geometry::Vector{Vector{Vector{Float64}}},
                                  omega::Vector{Vector{Float64}},
                                  bond_damage::Vector{Vector{Float64}},
                                  dof::Int64)
    B_temp = dof == 2 ? zeros(2, 2, 2) : zeros(3, 3, 3)

    @inbounds for i in nodes
        D_inv = @view inverse_shape_tensor[i, :, :]
        C_tensor = get_fourth_order(@view(C_Voigt[i, :, :]), dof)
        ni = nlist[i]
        CB_k = @view all_CB_tensors[i][eachindex(ni), :, :, :]
        @views precompute_CB_tensors!(CB_k, C_tensor, D_inv, ni, volume,
                                      bond_geometry[i], omega[i], bond_damage[i], dof,
                                      B_temp)
    end
    return nothing
end

function setup_zero_energy_views(include_zero_energy::Bool,
                                 zStiff::Union{Array{Float64,3},Nothing},
                                 active_nodes::AbstractVector{Int64},
                                 nnodes::Int)::Union{Vector{Union{SubArray{Float64,2},
                                                                  Nothing}},Nothing}
    if !include_zero_energy || zStiff === nothing
        return nothing
    end
    Z_views = Vector{Union{SubArray{Float64,2},Nothing}}(nothing, nnodes)
    @inbounds for i in active_nodes
        Z_views[i] = @view zStiff[i, :, :]
    end
    return Z_views
end

# =============================================================================
# Assembly subfunctions — now writing directly into ExtendableSparseMatrix
# =============================================================================

function assemble_normal_stiffness_forward!(K::ExtendableSparseMatrix{Float64,Int},
                                            i::Int,
                                            ni::Vector{Int64},
                                            CB_k_i::SubArray,
                                            D_inv_i::SubArray,
                                            volume::Vector{Float64},
                                            bond_geometry_i::Vector{Vector{Float64}},
                                            omega_i::Vector{Float64},
                                            bond_damage_i::Vector{Float64},
                                            nnodes::Int,
                                            dof::Int,
                                            K_block::Matrix{Float64})
    @inbounds for (j_idx, j) in enumerate(ni)
        bond_damage_i[j_idx] == 0.0 && continue
        ω_ij = omega_i[j_idx] * bond_damage_i[j_idx]
        V_j = volume[j]
        X_ij = bond_geometry_i[j_idx]

        for (k_idx, k) in enumerate(ni)
            bond_damage_i[k_idx] == 0.0 && continue
            ω_ik = omega_i[k_idx] * bond_damage_i[k_idx]
            V_k = volume[k]

            compute_stiffness_contribution(@view(CB_k_i[k_idx, :, :, :]),
                                           D_inv_i, X_ij, ω_ij, 1.0, ω_ik, V_k, dof,
                                           K_block)

            factor = V_j
            for m in 1:dof
                row = (m - 1) * nnodes + i
                for o in 1:dof
                    val = K_block[m, o] * factor
                    abs(val) <= 1e-14 && continue
                    col_i = (o - 1) * nnodes + i
                    col_k = (o - 1) * nnodes + k
                    updateindex!(K, +, val, row, col_i)
                    updateindex!(K, +, -val, row, col_k)
                end
            end
        end
    end
    return nothing
end

function assemble_zero_energy_forward!(K::ExtendableSparseMatrix{Float64,Int},
                                       i::Int,
                                       j::Int,
                                       ni::Vector{Int64},
                                       Z_i::SubArray{Float64,2},
                                       D_inv_i::SubArray{Float64,2},
                                       X_ij::Vector{Float64},
                                       ω_ij::Float64,
                                       V_j::Float64,
                                       bond_geometry_i::Vector{Vector{Float64}},
                                       omega_i::Vector{Float64},
                                       bond_damage_i::Vector{Float64},
                                       volume::Vector{Float64},
                                       nnodes::Int,
                                       dof::Int,
                                       D_inv_X::Vector{Float64})
    mul!(D_inv_X, D_inv_i, X_ij)
    α_ijj = ω_ij * V_j * dot(X_ij, D_inv_X)
    γ_ij = 1.0 - α_ijj
    factor1 = ω_ij * γ_ij * V_j

    @inbounds for m in 1:dof
        row = (m - 1) * nnodes + i
        for o in 1:dof
            val = factor1 * Z_i[m, o]
            abs(val) <= 1e-14 && continue
            col_i = (o - 1) * nnodes + i
            col_j = (o - 1) * nnodes + j
            updateindex!(K, +, val, row, col_i)
            updateindex!(K, +, -val, row, col_j)
        end
    end

    @inbounds for (p_idx, p) in enumerate(ni)
        (bond_damage_i[p_idx] == 0.0 || p == j) && continue
        X_ip = bond_geometry_i[p_idx]
        ω_ip = omega_i[p_idx] * bond_damage_i[p_idx]
        V_p = volume[p]
        coeff = dot(X_ip, D_inv_X)
        factor2 = -ω_ij * V_j * ω_ip * V_p * coeff

        for m in 1:dof
            row = (m - 1) * nnodes + i
            for o in 1:dof
                val = factor2 * Z_i[m, o]
                abs(val) <= 1e-14 && continue
                col_i = (o - 1) * nnodes + i
                col_p = (o - 1) * nnodes + p
                updateindex!(K, +, val, row, col_i)
                updateindex!(K, +, -val, row, col_p)
            end
        end
    end
    return nothing
end

function assemble_normal_stiffness_backward!(K::ExtendableSparseMatrix{Float64,Int},
                                             i::Int,
                                             j::Int,
                                             CB_k_j::SubArray,
                                             D_inv_j::SubArray,
                                             X_ji::Vector{Float64},
                                             ω_ji::Float64,
                                             V_i::Float64,
                                             j_neighbors::Vector{Int64},
                                             volume::Vector{Float64},
                                             omega_j::Vector{Float64},
                                             bond_damage_j::Vector{Float64},
                                             nnodes::Int,
                                             dof::Int,
                                             K_block::Matrix{Float64})
    @inbounds for (k_idx, k) in enumerate(j_neighbors)
        (omega_j[k_idx] * bond_damage_j[k_idx] == 0.0) && continue
        ω_jk = omega_j[k_idx] * bond_damage_j[k_idx]
        V_k = volume[k]

        compute_stiffness_contribution(@view(CB_k_j[k_idx, :, :, :]),
                                       D_inv_j, X_ji, ω_ji, 1.0, ω_jk, V_k, dof, K_block)

        factor = -V_i
        for m in 1:dof
            row = (m - 1) * nnodes + i
            for o in 1:dof
                val = K_block[m, o] * factor
                abs(val) <= 1e-14 && continue
                col_j = (o - 1) * nnodes + j
                col_k = (o - 1) * nnodes + k
                updateindex!(K, +, val, row, col_j)
                updateindex!(K, +, -val, row, col_k)
            end
        end
    end
    return nothing
end

function assemble_zero_energy_backward!(K::ExtendableSparseMatrix{Float64,Int},
                                        i::Int,
                                        j::Int,
                                        Z_j::SubArray{Float64,2},
                                        D_inv_j::SubArray,
                                        X_ji::Vector{Float64},
                                        ω_ji::Float64,
                                        V_i::Float64,
                                        j_neighbors::Vector{Int64},
                                        bond_geometry_j::Vector{Vector{Float64}},
                                        omega_j::Vector{Float64},
                                        bond_damage_j::Vector{Float64},
                                        volume::Vector{Float64},
                                        nnodes::Int,
                                        dof::Int,
                                        D_inv_X::Vector{Float64})
    mul!(D_inv_X, D_inv_j, X_ji)
    α_jii = ω_ji * V_i * dot(X_ji, D_inv_X)
    γ_ji = 1.0 - α_jii
    factor1 = -ω_ji * γ_ji * V_i

    @inbounds for m in 1:dof
        row = (m - 1) * nnodes + i
        for o in 1:dof
            val = factor1 * Z_j[m, o]
            abs(val) <= 1e-14 && continue
            col_j = (o - 1) * nnodes + j
            col_i = (o - 1) * nnodes + i
            updateindex!(K, +, val, row, col_j)
            updateindex!(K, +, -val, row, col_i)
        end
    end

    @inbounds for (p_idx, p) in enumerate(j_neighbors)
        (bond_damage_j[p_idx] == 0.0 || p == i) && continue
        X_jp = bond_geometry_j[p_idx]
        ω_jp = omega_j[p_idx] * bond_damage_j[p_idx]
        V_p = volume[p]
        coeff = dot(X_jp, D_inv_X)
        factor2 = ω_ji * V_i * ω_jp * V_p * coeff

        for m in 1:dof
            row = (m - 1) * nnodes + i
            for o in 1:dof
                val = factor2 * Z_j[m, o]
                abs(val) <= 1e-14 && continue
                col_j = (o - 1) * nnodes + j
                col_p = (o - 1) * nnodes + p
                updateindex!(K, +, val, row, col_j)
                updateindex!(K, +, -val, row, col_p)
            end
        end
    end
    return nothing
end

# =============================================================================
# Main assembly — returns SparseMatrixCSC directly via flush!
# =============================================================================

function assemble_stiffness_with_zero_energy(K::AbstractMatrix{Float64},
                                             nodes::AbstractVector{Int64},
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
                                             zStiff::Union{Array{Float64,3},Nothing},
                                             include_zero_energy::Bool,
                                             use_block_style::Bool)::SparseMatrixCSC{Float64,
                                                                                     Int}
    nnodes::Int = maximum(nodes)
    n_total::Int = nnodes * dof

    if include_zero_energy && zStiff === nothing
        @error "Zero-energy stabilization requested but zStiff not provided"
    end

    # --- ExtendableSparseMatrix statt Dict ---

    # Precompute CB tensors
    max_neighbors::Int = maximum(number_of_neighbors)
    all_CB_tensors = Vector{Array{Float64,4}}(undef, nnodes)
    @timeit "all_CB_tensors" begin
        @inbounds for i in 1:nnodes
            all_CB_tensors[i] = zeros(Float64, max_neighbors, dof, dof, dof)
        end
    end

    @timeit "precompute_CB" precompute_CB_all_nodes!(all_CB_tensors, nodes, C_Voigt,
                                                     inverse_shape_tensor, nlist, volume,
                                                     bond_geometry, omega, bond_damage, dof)

    @timeit "setup_Z" Z=setup_zero_energy_views(include_zero_energy, zStiff,
                                                active_nodes, nnodes)

    # Reusable buffers
    K_block = zeros(dof, dof)
    D_inv_X = Vector{Float64}(undef, dof)

    # Forward assembly
    @timeit "forward" begin
        @inbounds for i in active_nodes
            ni = nlist[i]
            D_inv_i = @view inverse_shape_tensor[i, :, :]
            CB_k_i = @view all_CB_tensors[i][eachindex(ni), :, :, :]

            assemble_normal_stiffness_forward!(K, i, ni, CB_k_i, D_inv_i,
                                               volume, bond_geometry[i],
                                               omega[i], bond_damage[i],
                                               nnodes, dof, K_block)

            if include_zero_energy
                Z_i = Z[i]
                for (j_idx, j) in enumerate(ni)
                    bond_damage[i][j_idx] == 0.0 && continue
                    ω_ij = omega[i][j_idx] * bond_damage[i][j_idx]
                    V_j = volume[j]
                    X_ij = bond_geometry[i][j_idx]
                    assemble_zero_energy_forward!(K, i, j, ni, Z_i, D_inv_i,
                                                  X_ij, ω_ij, V_j,
                                                  bond_geometry[i], omega[i],
                                                  bond_damage[i], volume,
                                                  nnodes, dof, D_inv_X)
                end
            end
        end
    end

    # Backward assembly
    @timeit "backward" begin
        @inbounds for i in active_nodes
            V_i = volume[i]
            for j in active_nodes
                j == i && continue
                j_neighbors = nlist[j]
                i_pos_in_j = findfirst(==(i), j_neighbors)
                (i_pos_in_j === nothing ||
                 omega[j][i_pos_in_j] * bond_damage[j][i_pos_in_j] == 0.0) && continue

                ω_ji = omega[j][i_pos_in_j] * bond_damage[j][i_pos_in_j]
                X_ji = bond_geometry[j][i_pos_in_j]
                D_inv_j = @view inverse_shape_tensor[j, :, :]
                CB_k_j = @view all_CB_tensors[j][eachindex(j_neighbors), :, :, :]

                assemble_normal_stiffness_backward!(K, i, j, CB_k_j, D_inv_j,
                                                    X_ji, ω_ji, V_i, j_neighbors,
                                                    volume, omega[j], bond_damage[j],
                                                    nnodes, dof, K_block)

                if include_zero_energy
                    Z_j = Z[j]
                    assemble_zero_energy_backward!(K, i, j, Z_j, D_inv_j,
                                                   X_ji, ω_ji, V_i, j_neighbors,
                                                   bond_geometry[j], omega[j],
                                                   bond_damage[j], volume,
                                                   nnodes, dof, D_inv_X)
                end
            end
        end
    end

    # flush! konvertiert LNK → CSC intern, sparse() gibt SparseMatrixCSC zurück
    @timeit "flush" flush!(K)
    return sparse(K)
end

# =============================================================================
# Voigt rotation helpers
# =============================================================================

function rotation_voigt_3D(R::AbstractMatrix{Float64}, C::AbstractMatrix{Float64})
    l1, m1, n1 = R[1, 1], R[1, 2], R[1, 3]
    l2, m2, n2 = R[2, 1], R[2, 2], R[2, 3]
    l3, m3, n3 = R[3, 1], R[3, 2], R[3, 3]
    M = [l1^2 m1^2 n1^2 2m1*n1 2l1*n1 2l1*m1;
         l2^2 m2^2 n2^2 2m2*n2 2l2*n2 2l2*m2;
         l3^2 m3^2 n3^2 2m3*n3 2l3*n3 2l3*m3;
         l2*l3 m2*m3 n2*n3 m2 * n3+m3 * n2 l2 * n3+l3 * n2 l2 * m3+l3 * m2;
         l1*l3 m1*m3 n1*n3 m1 * n3+m3 * n1 l1 * n3+l3 * n1 l1 * m3+l3 * m1;
         l1*l2 m1*m2 n1*n2 m1 * n2+m2 * n1 l1 * n2+l2 * n1 l1 * m2+l2 * m1]
    return M * C * M'
end

function rotation_voigt_2D(R::AbstractMatrix{Float64}, C::AbstractMatrix{Float64})
    l1, m1 = R[1, 1], R[1, 2]
    l2, m2 = R[2, 1], R[2, 2]
    M = [l1^2 m1^2 2l1*m1;
         l2^2 m2^2 2l2*m2;
         l1*l2 m1*m2 l1 * m2+l2 * m1]
    return M * C * M'
end

# =============================================================================
# Public interface
# =============================================================================

function init_model(nodes::AbstractVector{Int64},
                    material_parameter::Dict, block_id::Int64)
    Zero_Energy_Control.init_model(nodes, material_parameter, block_id)
    if Data_Manager.get_max_rank() > 1
        @error "Correspondence matrix based not implemented for parallel runs."
    end
end

function init_matrix(use_block_style::Bool = true, include_zero_energy::Bool = true)
    nodes = collect(1:Data_Manager.get_nnodes())
    dof = Data_Manager.get_dof()
    nnodes = length(nodes)

    bond_geometry = Data_Manager.get_field("Bond Geometry")
    inverse_shape_tensor = Data_Manager.get_field("Inverse Shape Tensor")
    nlist = Data_Manager.get_nlist()
    number_of_neighbors = Data_Manager.get_field("Number of Neighbors")
    volume = Data_Manager.get_field("Volume")
    omega = Data_Manager.get_field("Influence Function")
    C_voigt = Data_Manager.get_field("Elasticity Matrix")
    bond_damage = Data_Manager.get_field("Bond Damage", "NP1")

    if Data_Manager.get_rotation()
        C_voigt_trafo = similar(C_voigt)
        rotation_tensor = Data_Manager.get_field("Rotation Tensor")
        if dof == 2
            for i in 1:nnodes
                @views C_voigt_trafo[i, :, :] = rotation_voigt_2D(rotation_tensor[i, :, :],
                                                                  C_voigt[i, :, :])
            end
        elseif dof == 3
            for i in 1:nnodes
                @views C_voigt_trafo[i, :, :] = rotation_voigt_3D(rotation_tensor[i, :, :],
                                                                  C_voigt[i, :, :])
            end
        end
    else
        C_voigt_trafo = C_voigt
    end

    use_zero_energy = false
    zStiff = nothing
    if haskey(Data_Manager.get_properties(1, "Material Model"), "Zero Energy Control")
        if Data_Manager.get_properties(1, "Material Model")["Zero Energy Control"] ==
           "Global"
            use_zero_energy = include_zero_energy
            if use_zero_energy
                zStiff = Data_Manager.create_constant_node_tensor_field("Zero Energy Stiffness",
                                                                        Float64, dof)
                Zero_Energy_Control.create_zero_energy_mode_stiffness!(nodes, dof,
                                                                       C_voigt_trafo,
                                                                       inverse_shape_tensor,
                                                                       zStiff)
            end
        end
    end
    Data_Manager.init_stiffness_matrix(nnodes, dof)
    K_sparse = Data_Manager.get_stiffness_matrix()
    @timeit "assemble" K_sparse=assemble_stiffness_with_zero_energy(K_sparse, nodes, nodes,
                                                                    dof,
                                                                    C_voigt_trafo,
                                                                    inverse_shape_tensor,
                                                                    number_of_neighbors,
                                                                    nlist, volume,
                                                                    bond_geometry, omega,
                                                                    bond_damage, zStiff,
                                                                    use_zero_energy,
                                                                    use_block_style)

    Data_Manager.set_stiffness_matrix(-K_sparse)
end

function compute_model(nodes::AbstractVector{Int64},
                       use_block_style::Bool = true,
                       include_zero_energy::Bool = true)
    dof = Data_Manager.get_dof()
    nnodes = Data_Manager.get_nnodes()

    C_voigt = Data_Manager.get_field("Elasticity Matrix")
    inverse_shape_tensor = Data_Manager.get_field("Inverse Shape Tensor")
    nlist = Data_Manager.get_nlist()
    volume = Data_Manager.get_field("Volume")
    bond_geometry_N = Data_Manager.get_field("Bond Geometry")
    number_of_neighbors = Data_Manager.get_field("Number of Neighbors")
    omega = Data_Manager.get_field("Influence Function")
    bond_damage = Data_Manager.get_field("Bond Damage", "NP1")
    zStiff = include_zero_energy ? Data_Manager.get_field("Zero Energy Stiffness") : nothing

    if include_zero_energy && zStiff !== nothing
        Zero_Energy_Control.create_zero_energy_mode_stiffness!(nodes, dof, C_voigt,
                                                               inverse_shape_tensor, zStiff)
    end
    Data_Manager.init_stiffness_matrix(nnodes, dof)
    K_sparse = Data_Manager.get_stiffness_matrix()
    @timeit "assemble" K_sparse=assemble_stiffness_with_zero_energy(K_sparse,
                                                                    collect(1:nnodes),
                                                                    nodes, dof, C_voigt,
                                                                    inverse_shape_tensor,
                                                                    number_of_neighbors,
                                                                    nlist, volume,
                                                                    bond_geometry_N, omega,
                                                                    bond_damage, zStiff,
                                                                    include_zero_energy,
                                                                    use_block_style)

    Data_Manager.set_stiffness_matrix(-K_sparse)
end

end # module
