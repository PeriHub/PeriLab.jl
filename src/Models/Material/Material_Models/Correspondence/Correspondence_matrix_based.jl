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
# Index map structure
# =============================================================================

"""
    NZValIndexMap

Stores for each bond b = (iID → jID) with j_idx in nlist[iID] the direct
positions in nzval of the global stiffness matrix:

  idx_ii[di,dj,b]: position of K[row_i, col_i]  (iID diagonal)
  idx_ij[di,dj,b]: position of K[row_i, col_j]  (iID → jID)
  idx_ji[di,dj,b]: position of K[row_j, col_i]  (jID → iID)
  idx_jj[di,dj,b]: position of K[row_j, col_j]  (jID diagonal)

For the inner k-loop:
  bond_lookup[iID][k_idx] = bond_idx of bond (iID → nlist[iID][k_idx])

This allows O(1) lookup of nzval positions for every (j,k) contribution,
identical to the K_local logic in the fallback assembler.
"""
struct NZValIndexMap
    bond_node::Vector{Int32}          # bond_idx → iID
    bond_neighbor::Vector{Int32}      # bond_idx → jID
    bond_neighbor_idx::Vector{Int32}  # bond_idx → j_idx in nlist[iID]

    # [di, dj, bond_idx] → index in nzval
    idx_ii::Array{Int32,3}
    idx_ij::Array{Int32,3}
    idx_ji::Array{Int32,3}
    idx_jj::Array{Int32,3}

    # bond_lookup[iID][k_idx] = bond_idx of bond (iID → nlist[iID][k_idx])
    bond_lookup::Vector{Vector{Int32}}

    n_bonds::Int
    dof::Int
    nnodes::Int
end

"""
    find_nzval_idx(csc, row, col)

Find the index of K[row, col] in csc.nzval via binary search within the column.
Throws an assertion error if the entry does not exist in the sparsity pattern.
"""
@inline function find_nzval_idx(csc::SparseMatrixCSC, row::Int, col::Int)
    cs = csc.colptr[col]
    ce = csc.colptr[col + 1] - 1
    idx = searchsortedfirst(csc.rowval, row, cs, ce, Base.Order.Forward)
    @assert idx <= ce&&csc.rowval[idx] == row "Entry ($row,$col) not in sparsity pattern!"
    return idx
end

"""
    build_nzval_index_map(csc, active_nodes, nlist, dof, nnodes)

Build the complete index map once after the first flush!.
Must be called after the CSC sparsity pattern is known.
Uses Int32 for index storage (half the memory of Int64, sufficient for all
realistic problem sizes: max ~2.1 billion entries).
"""
function build_nzval_index_map(csc::SparseMatrixCSC{Float64,Int},
                               active_nodes::AbstractVector{Int},
                               nlist::Vector{Vector{Int}},
                               dof::Int,
                               nnodes::Int)::NZValIndexMap
    n_bonds = sum(length(nlist[i]) for i in active_nodes)

    bond_node = Vector{Int32}(undef, n_bonds)
    bond_neighbor = Vector{Int32}(undef, n_bonds)
    bond_neighbor_idx = Vector{Int32}(undef, n_bonds)

    idx_ii = Array{Int32}(undef, dof, dof, n_bonds)
    idx_ij = Array{Int32}(undef, dof, dof, n_bonds)
    idx_ji = Array{Int32}(undef, dof, dof, n_bonds)
    idx_jj = Array{Int32}(undef, dof, dof, n_bonds)

    # bond_lookup[iID] is empty for nodes not in active_nodes
    bond_lookup = [Vector{Int32}() for _ in 1:nnodes]

    bond_idx = 0
    @inbounds for iID in active_nodes
        nj = nlist[iID]
        mj = length(nj)
        bond_lookup[iID] = Vector{Int32}(undef, mj)

        for (j_idx, jID) in enumerate(nj)
            bond_idx += 1
            bond_node[bond_idx] = iID
            bond_neighbor[bond_idx] = jID
            bond_neighbor_idx[bond_idx] = j_idx
            bond_lookup[iID][j_idx] = bond_idx  # O(1) lookup for k-loop

            for di in 1:dof, dj in 1:dof
                # Global row/col in block-style DOF layout: (dof-1)*nnodes + node
                row_i = (di - 1) * nnodes + iID
                row_j = (di - 1) * nnodes + jID
                col_i = (dj - 1) * nnodes + iID
                col_j = (dj - 1) * nnodes + jID

                idx_ii[di, dj, bond_idx] = find_nzval_idx(csc, row_i, col_i)
                idx_ij[di, dj, bond_idx] = find_nzval_idx(csc, row_i, col_j)
                idx_ji[di, dj, bond_idx] = find_nzval_idx(csc, row_j, col_i)
                idx_jj[di, dj, bond_idx] = find_nzval_idx(csc, row_j, col_j)
            end
        end
    end

    return NZValIndexMap(bond_node, bond_neighbor, bond_neighbor_idx,
                         idx_ii, idx_ij, idx_ji, idx_jj,
                         bond_lookup, n_bonds, dof, nnodes)
end

# =============================================================================
# Tensor computation kernels
# =============================================================================

@inline function get_dof_index_block_style(node_id::Int, dof_component::Int,
                                           nnodes::Int, ::Int)
    return (dof_component - 1) * nnodes + node_id
end

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

function compute_stiffness_contribution(CB_k::Array{Float64,4}, k_idx::Int,
                                        D_inv_i::AbstractMatrix{Float64},
                                        X_ij::Vector{Float64},
                                        omega_ij::Float64, omega_ik::Float64,
                                        V_k::Float64, dof::Int,
                                        K_block::Matrix{Float64})
    if dof == 2
        DX_1 = D_inv_i[1, 1] * X_ij[1] + D_inv_i[1, 2] * X_ij[2]
        DX_2 = D_inv_i[2, 1] * X_ij[1] + D_inv_i[2, 2] * X_ij[2]
        f = omega_ij * omega_ik * V_k
        @inbounds begin
            K_block[1, 1] = f * (CB_k[k_idx, 1, 1, 1] * DX_1 + CB_k[k_idx, 1, 2, 1] * DX_2)
            K_block[1, 2] = f * (CB_k[k_idx, 1, 1, 2] * DX_1 + CB_k[k_idx, 1, 2, 2] * DX_2)
            K_block[2, 1] = f * (CB_k[k_idx, 2, 1, 1] * DX_1 + CB_k[k_idx, 2, 2, 1] * DX_2)
            K_block[2, 2] = f * (CB_k[k_idx, 2, 1, 2] * DX_1 + CB_k[k_idx, 2, 2, 2] * DX_2)
        end
    elseif dof == 3
        DX_1 = D_inv_i[1, 1] * X_ij[1] + D_inv_i[1, 2] * X_ij[2] + D_inv_i[1, 3] * X_ij[3]
        DX_2 = D_inv_i[2, 1] * X_ij[1] + D_inv_i[2, 2] * X_ij[2] + D_inv_i[2, 3] * X_ij[3]
        DX_3 = D_inv_i[3, 1] * X_ij[1] + D_inv_i[3, 2] * X_ij[2] + D_inv_i[3, 3] * X_ij[3]
        f = omega_ij * omega_ik * V_k
        @inbounds begin
            K_block[1, 1] = f *
                            (CB_k[k_idx, 1, 1, 1] * DX_1 + CB_k[k_idx, 1, 2, 1] * DX_2 +
                             CB_k[k_idx, 1, 3, 1] * DX_3)
            K_block[1, 2] = f *
                            (CB_k[k_idx, 1, 1, 2] * DX_1 + CB_k[k_idx, 1, 2, 2] * DX_2 +
                             CB_k[k_idx, 1, 3, 2] * DX_3)
            K_block[1, 3] = f *
                            (CB_k[k_idx, 1, 1, 3] * DX_1 + CB_k[k_idx, 1, 2, 3] * DX_2 +
                             CB_k[k_idx, 1, 3, 3] * DX_3)
            K_block[2, 1] = f *
                            (CB_k[k_idx, 2, 1, 1] * DX_1 + CB_k[k_idx, 2, 2, 1] * DX_2 +
                             CB_k[k_idx, 2, 3, 1] * DX_3)
            K_block[2, 2] = f *
                            (CB_k[k_idx, 2, 1, 2] * DX_1 + CB_k[k_idx, 2, 2, 2] * DX_2 +
                             CB_k[k_idx, 2, 3, 2] * DX_3)
            K_block[2, 3] = f *
                            (CB_k[k_idx, 2, 1, 3] * DX_1 + CB_k[k_idx, 2, 2, 3] * DX_2 +
                             CB_k[k_idx, 2, 3, 3] * DX_3)
            K_block[3, 1] = f *
                            (CB_k[k_idx, 3, 1, 1] * DX_1 + CB_k[k_idx, 3, 2, 1] * DX_2 +
                             CB_k[k_idx, 3, 3, 1] * DX_3)
            K_block[3, 2] = f *
                            (CB_k[k_idx, 3, 1, 2] * DX_1 + CB_k[k_idx, 3, 2, 2] * DX_2 +
                             CB_k[k_idx, 3, 3, 2] * DX_3)
            K_block[3, 3] = f *
                            (CB_k[k_idx, 3, 1, 3] * DX_1 + CB_k[k_idx, 3, 2, 3] * DX_2 +
                             CB_k[k_idx, 3, 3, 3] * DX_3)
        end
    end
end

function precompute_CB_tensors!(CB_k::AbstractArray{Float64,4},
                                C_tensor::Array{Float64,4},
                                D_inv::AbstractMatrix{Float64},
                                neighbors::Vector{Int64}, volume::Vector{Float64},
                                bond_geometry_i::Vector{Vector{Float64}},
                                omega_i::Vector{Float64}, bond_damage_i::Vector{Float64},
                                dof::Int64, B_ik::Array{Float64,3})
    @inbounds for (k_idx, k) in enumerate(neighbors)
        create_B_tensor!(D_inv, bond_geometry_i[k_idx], volume[k],
                         omega_i[k_idx] * bond_damage_i[k_idx], Val(dof), B_ik)
        @views contraction!(C_tensor, B_ik, dof, CB_k[k_idx, :, :, :])
    end
end

function precompute_CB_all_nodes!(all_CB_tensors::Vector{Array{Float64,4}},
                                  nodes::AbstractVector{Int64},
                                  C_Voigt::Array{Float64,3},
                                  inverse_shape_tensor::Array{Float64,3},
                                  nlist::Vector{Vector{Int64}}, volume::Vector{Float64},
                                  bond_geometry::Vector{Vector{Vector{Float64}}},
                                  omega::Vector{Vector{Float64}},
                                  bond_damage::Vector{Vector{Float64}}, dof::Int64)
    B_temp = dof == 2 ? zeros(2, 2, 2) : zeros(3, 3, 3)
    @inbounds for i in nodes
        D_inv = @view inverse_shape_tensor[i, :, :]
        C_tensor = get_fourth_order(@view(C_Voigt[i, :, :]), dof)
        ni = nlist[i]
        @views precompute_CB_tensors!(all_CB_tensors[i][eachindex(ni), :, :, :],
                                      C_tensor, D_inv, ni, volume,
                                      bond_geometry[i], omega[i], bond_damage[i],
                                      dof, B_temp)
    end
end

# =============================================================================
# Assembly variant 1: Classic K_local + rawupdateindex! (fallback / first call)
#
# Loop structure:
#   for iID in active_nodes           ← node i
#     for (j_idx, jID) in nj         ← bond (i→j)
#       for (k_idx, kID) in nj       ← k-loop: contributions to K[i,i], K[i,k],
#                                                                 K[j,i], K[j,k]
#
# K_local layout: rows/cols = [nj[1]..nj[mj], iID] × dof
#   row_i = (m-1)*(mj+1) + (mj+1)   ← iID in row m
#   row_j = (m-1)*(mj+1) + j_idx    ← jID in row m
#   col_i = (o-1)*(mj+1) + (mj+1)   ← iID in col o
#   col_k = (o-1)*(mj+1) + k_idx    ← kID in col o
# =============================================================================

function assemble_stiffness_with_zero_energy(K::AbstractMatrix{Float64},
                                             nodes::AbstractVector{Int},
                                             active_nodes::AbstractVector{Int},
                                             dof::Int,
                                             C_Voigt::Array{Float64,3},
                                             inverse_shape_tensor::Array{Float64,3},
                                             number_of_neighbors::Vector{Int},
                                             nlist::Vector{Vector{Int}},
                                             inverse_nlist::Vector{Dict{Int,Int}},
                                             volume::Vector{Float64},
                                             bond_geometry::Vector{Vector{Vector{Float64}}},
                                             omega::Vector{Vector{Float64}},
                                             bond_damage::Vector{Vector{Float64}},
                                             zStiff::Union{Array{Float64,3},Nothing},
                                             all_CB_tensors::Vector{Array{Float64,4}})::SparseMatrixCSC{Float64,
                                                                                                        Int}
    nnodes = maximum(nodes)
    max_neighbors = maximum(number_of_neighbors)

    @timeit "precompute_CB" precompute_CB_all_nodes!(all_CB_tensors, nodes, C_Voigt,
                                                     inverse_shape_tensor, nlist, volume,
                                                     bond_geometry, omega, bond_damage, dof)

    K_block = zeros(dof, dof)
    D_inv_X = Vector{Float64}(undef, dof)
    ldof = (max_neighbors + 1) * dof
    K_local = zeros(ldof, ldof)

    @timeit "forward+backward" @inbounds for iID in active_nodes
        nj = nlist[iID]
        mj = length(nj)
        D_inv_i = @view inverse_shape_tensor[iID, :, :]
        CB_k_i = all_CB_tensors[iID]
        Z_i = zStiff !== nothing ? (@view zStiff[iID, :, :]) : nothing
        V_i = volume[iID]

        fill!(@view(K_local[1:((mj + 1) * dof), 1:((mj + 1) * dof)]), 0.0)

        for (j_idx, jID) in enumerate(nj)
            bond_damage[iID][j_idx] == 0.0 && continue
            ω_ij = omega[iID][j_idx] * bond_damage[iID][j_idx]
            V_j = volume[jID]
            X_ij = bond_geometry[iID][j_idx]
            mul!(D_inv_X, D_inv_i, X_ij)

            # Forward k-loop: K[i,i] += val*V_j,  K[i,k] -= val*V_j
            for (k_idx, kID) in enumerate(nj)
                bond_damage[iID][k_idx] == 0.0 && continue
                compute_stiffness_contribution(CB_k_i, k_idx, D_inv_i, X_ij,
                                               ω_ij,
                                               omega[iID][k_idx] * bond_damage[iID][k_idx],
                                               volume[kID], dof, K_block)
                for m in 1:dof
                    row_i = (m - 1) * (mj + 1) + (mj + 1)
                    for o in 1:dof
                        val = K_block[m, o] * V_j
                        K_local[row_i, (o - 1) * (mj + 1) + (mj + 1)] += val  # col_i
                        K_local[row_i, (o - 1) * (mj + 1) + k_idx] -= val  # col_k
                    end
                end
            end

            # Zero-energy forward: K[i,i] += fac1*Z, K[i,j] -= fac1*Z
            #                      K[i,i] += fac2*Z, K[i,k] -= fac2*Z
            if Z_i !== nothing
                γ = 1.0 - ω_ij * V_j * dot(X_ij, D_inv_X)
                fac1 = ω_ij * γ * V_j
                for m in 1:dof
                    row_i = (m - 1) * (mj + 1) + (mj + 1)
                    for o in 1:dof
                        val = fac1 * Z_i[m, o]
                        K_local[row_i, (o - 1) * (mj + 1) + (mj + 1)] += val  # col_i
                        K_local[row_i, (o - 1) * (mj + 1) + j_idx] -= val  # col_j
                    end
                end
                for (k_idx, kID) in enumerate(nj)
                    (bond_damage[iID][k_idx] == 0.0 || kID == jID) && continue
                    fac2 = -ω_ij * V_j * omega[iID][k_idx] * bond_damage[iID][k_idx] *
                           volume[kID] * dot(bond_geometry[iID][k_idx], D_inv_X)
                    for m in 1:dof
                        row_i = (m - 1) * (mj + 1) + (mj + 1)
                        for o in 1:dof
                            val = fac2 * Z_i[m, o]
                            K_local[row_i, (o - 1) * (mj + 1) + (mj + 1)] += val  # col_i
                            K_local[row_i, (o - 1) * (mj + 1) + k_idx] -= val  # col_k
                        end
                    end
                end
            end

            # Backward k-loop: K[j,i] += val*(-V_i), K[j,k] -= val*(-V_i)
            for (k_idx, kID) in enumerate(nj)
                bond_damage[iID][k_idx] == 0.0 && continue
                compute_stiffness_contribution(CB_k_i, k_idx, D_inv_i, X_ij,
                                               ω_ij,
                                               omega[iID][k_idx] * bond_damage[iID][k_idx],
                                               volume[kID], dof, K_block)
                for m in 1:dof
                    row_j = (m - 1) * (mj + 1) + j_idx
                    for o in 1:dof
                        val = K_block[m, o] * (-V_i)
                        K_local[row_j, (o - 1) * (mj + 1) + (mj + 1)] += val  # col_i
                        K_local[row_j, (o - 1) * (mj + 1) + k_idx] -= val  # col_k
                    end
                end
            end

            # Zero-energy backward: K[j,i] += fac1*Z, K[j,j] -= fac1*Z
            #                       K[j,i] += fac2*Z, K[j,k] -= fac2*Z
            if Z_i !== nothing
                γ = 1.0 - ω_ij * V_i * dot(X_ij, D_inv_X)
                fac1 = -ω_ij * γ * V_i
                for m in 1:dof
                    row_j = (m - 1) * (mj + 1) + j_idx
                    for o in 1:dof
                        val = fac1 * Z_i[m, o]
                        K_local[row_j, (o - 1) * (mj + 1) + (mj + 1)] += val  # col_i
                        K_local[row_j, (o - 1) * (mj + 1) + j_idx] -= val  # col_j
                    end
                end
                for (k_idx, kID) in enumerate(nj)
                    (bond_damage[iID][k_idx] == 0.0 || kID == jID) && continue
                    fac2 = ω_ij * V_i * omega[iID][k_idx] * bond_damage[iID][k_idx] *
                           volume[kID] * dot(bond_geometry[iID][k_idx], D_inv_X)
                    for m in 1:dof
                        row_j = (m - 1) * (mj + 1) + j_idx
                        for o in 1:dof
                            val = fac2 * Z_i[m, o]
                            K_local[row_j, (o - 1) * (mj + 1) + (mj + 1)] += val  # col_i
                            K_local[row_j, (o - 1) * (mj + 1) + k_idx] -= val  # col_k
                        end
                    end
                end
            end
        end # j_idx loop

        # Scatter local → global stiffness matrix
        @timeit "rawupdateindex!" for row_local in 1:(mj + 1)
            row_node = row_local <= mj ? nj[row_local] : iID
            for m in 1:dof
                grow = (m - 1) * nnodes + row_node
                for col_local in 1:(mj + 1)
                    col_node = col_local <= mj ? nj[col_local] : iID
                    for o in 1:dof
                        val = K_local[(m - 1) * (mj + 1) + row_local,
                        (o - 1) * (mj + 1) + col_local]
                        abs(val) <= 1e-14 && continue
                        rawupdateindex!(K, +, val, grow, (o - 1) * nnodes + col_node)
                    end
                end
            end
        end
    end # iID loop

    @timeit "flush" flush!(K)
    return K
end

# =============================================================================
# Assembly variant 2: Fast direct nzval write
#
# Identical logic to variant 1. Mapping K_local → nzval:
#
#   K_local[row_i, col_i] ↔ nzval[idx_ii[m,o, b_ij]]
#   K_local[row_i, col_k] ↔ nzval[idx_ij[m,o, b_ik]]  b_ik = bond_lookup[iID][k_idx]
#   K_local[row_j, col_i] ↔ nzval[idx_ji[m,o, b_ij]]
#   K_local[row_j, col_k] ↔ nzval[idx_jj[m,o, b_ik]]  b_ik = bond_lookup[iID][k_idx]
#
# For col_j (zero-energy fac1): k_idx == j_idx, so b_ik = b_ij:
#   K_local[row_i, col_j] ↔ nzval[idx_ij[m,o, b_ij]]  ✓
#   K_local[row_j, col_j] ↔ nzval[idx_jj[m,o, b_ij]]  ✓
#
# No allocation, no LNK overhead, no flush! — only array index writes.
# nzval must be zeroed before calling this function.
# =============================================================================

function assemble_into_nzval!(nzval::Vector{Float64},
                              imap::NZValIndexMap,
                              all_CB_tensors::Vector{Array{Float64,4}},
                              nlist::Vector{Vector{Int64}},
                              inverse_shape_tensor::Array{Float64,3},
                              volume::Vector{Float64},
                              bond_geometry::Vector{Vector{Vector{Float64}}},
                              omega::Vector{Vector{Float64}},
                              bond_damage::Vector{Vector{Float64}},
                              zStiff::Union{Array{Float64,3},Nothing},
                              dof::Int, nnodes::Int)
    K_block = zeros(dof, dof)
    D_inv_X = Vector{Float64}(undef, dof)

    @inbounds for b_ij in 1:(imap.n_bonds)
        iID = Int(imap.bond_node[b_ij])
        jID = Int(imap.bond_neighbor[b_ij])
        j_idx = Int(imap.bond_neighbor_idx[b_ij])

        bond_damage[iID][j_idx] == 0.0 && continue

        nj = nlist[iID]
        D_inv_i = @view inverse_shape_tensor[iID, :, :]
        CB_k_i = all_CB_tensors[iID]
        Z_i = zStiff !== nothing ? (@view zStiff[iID, :, :]) : nothing
        V_i = volume[iID]
        V_j = volume[jID]
        X_ij = bond_geometry[iID][j_idx]
        ω_ij = omega[iID][j_idx] * bond_damage[iID][j_idx]

        mul!(D_inv_X, D_inv_i, X_ij)

        # Forward k-loop
        # K_local[row_i, col_i] += val*V_j → idx_ii[m,o, b_ij]
        # K_local[row_i, col_k] -= val*V_j → idx_ij[m,o, b_ik]
        for (k_idx, kID) in enumerate(nj)
            bond_damage[iID][k_idx] == 0.0 && continue
            b_ik = Int(imap.bond_lookup[iID][k_idx])
            compute_stiffness_contribution(CB_k_i, k_idx, D_inv_i, X_ij,
                                           ω_ij,
                                           omega[iID][k_idx] * bond_damage[iID][k_idx],
                                           volume[kID], dof, K_block)
            for m in 1:dof, o in 1:dof
                val = K_block[m, o] * V_j
                nzval[imap.idx_ii[m, o, b_ij]] += val
                nzval[imap.idx_ij[m, o, b_ik]] -= val
            end
        end

        # Zero-energy forward
        if Z_i !== nothing
            γ = 1.0 - ω_ij * V_j * dot(X_ij, D_inv_X)
            fac1 = ω_ij * γ * V_j
            # k == j case: b_ik = b_ij, idx_ij[b_ij] maps to K[i,j]
            for m in 1:dof, o in 1:dof
                val = fac1 * Z_i[m, o]
                nzval[imap.idx_ii[m, o, b_ij]] += val
                nzval[imap.idx_ij[m, o, b_ij]] -= val
            end
            for (k_idx, kID) in enumerate(nj)
                (bond_damage[iID][k_idx] == 0.0 || kID == jID) && continue
                b_ik = Int(imap.bond_lookup[iID][k_idx])
                fac2 = -ω_ij * V_j * omega[iID][k_idx] * bond_damage[iID][k_idx] *
                       volume[kID] * dot(bond_geometry[iID][k_idx], D_inv_X)
                for m in 1:dof, o in 1:dof
                    val = fac2 * Z_i[m, o]
                    nzval[imap.idx_ii[m, o, b_ij]] += val
                    nzval[imap.idx_ij[m, o, b_ik]] -= val
                end
            end
        end

        # Backward k-loop
        # K_local[row_j, col_i] += val*(-V_i) → idx_ji[m,o, b_ij]
        # K_local[row_j, col_k] -= val*(-V_i) → idx_jj[m,o, b_ik]
        for (k_idx, kID) in enumerate(nj)
            bond_damage[iID][k_idx] == 0.0 && continue
            b_ik = Int(imap.bond_lookup[iID][k_idx])
            compute_stiffness_contribution(CB_k_i, k_idx, D_inv_i, X_ij,
                                           ω_ij,
                                           omega[iID][k_idx] * bond_damage[iID][k_idx],
                                           volume[kID], dof, K_block)
            for m in 1:dof, o in 1:dof
                val = K_block[m, o] * (-V_i)
                nzval[imap.idx_ji[m, o, b_ij]] += val
                nzval[imap.idx_jj[m, o, b_ik]] -= val
            end
        end

        # Zero-energy backward
        if Z_i !== nothing
            γ = 1.0 - ω_ij * V_i * dot(X_ij, D_inv_X)
            fac1 = -ω_ij * γ * V_i
            # k == j case: idx_jj[b_ij] maps to K[j,j]
            for m in 1:dof, o in 1:dof
                val = fac1 * Z_i[m, o]
                nzval[imap.idx_ji[m, o, b_ij]] += val
                nzval[imap.idx_jj[m, o, b_ij]] -= val
            end
            for (k_idx, kID) in enumerate(nj)
                (bond_damage[iID][k_idx] == 0.0 || kID == jID) && continue
                b_ik = Int(imap.bond_lookup[iID][k_idx])
                fac2 = ω_ij * V_i * omega[iID][k_idx] * bond_damage[iID][k_idx] *
                       volume[kID] * dot(bond_geometry[iID][k_idx], D_inv_X)
                for m in 1:dof, o in 1:dof
                    val = fac2 * Z_i[m, o]
                    nzval[imap.idx_ji[m, o, b_ij]] += val
                    nzval[imap.idx_jj[m, o, b_ik]] -= val
                end
            end
        end
    end # b_ij loop
    return nothing
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
# Module-level persistent state
#
# These are allocated once in init_matrix and reused every compute_model call:
#   _K_sparse       : the ExtendableSparseMatrix — never recreated
#   _nzval          : direct reference to K_sparse.cscmatrix.nzval
#                     fill!(_nzval[], 0.0) zeros the matrix without allocation
#   _index_map      : precomputed nzval positions for all bonds
#   _all_CB_tensors : preallocated CB tensor storage, refilled each call
#   _pattern_nodes  : node set for which the index map was built
#
# Using Ref{Any} instead of Ref{ConcreteType} avoids world-age errors
# when the module is reloaded with Revise.jl.
# =============================================================================

const _K_sparse = Ref{Any}(nothing)
const _nzval = Ref{Union{Vector{Float64},Nothing}}(nothing)
const _index_map = Ref{Any}(nothing)
const _all_CB_tensors = Ref{Any}(nothing)
const _pattern_nodes = Ref{Vector{Int}}(Int[])

"""
Returns true if the stored index map is valid for the given active_nodes.
False if: not yet built, wrong type (world-age after Revise reload),
or node set has changed (additive manufacturing).
"""
function _index_map_valid(active_nodes::AbstractVector{Int})::Bool
    _index_map[] === nothing && return false
    !isa(_index_map[], NZValIndexMap) && return false  # world-age guard
    length(_pattern_nodes[]) != length(active_nodes) && return false
    return _pattern_nodes[] == collect(active_nodes)
end

"""
Build or rebuild all persistent state. Called from init_matrix and
from compute_model when the node set changes (additive manufacturing).
"""
function _build_persistent_state!(K::ExtendableSparseMatrix,
                                  nodes::AbstractVector{Int},
                                  nlist::Vector{Vector{Int}},
                                  dof::Int, nnodes::Int,
                                  max_neighbors::Int)
    # Store direct reference to nzval — fill! on _nzval[] zeros K automatically
    _K_sparse[] = K
    _nzval[] = K.cscmatrix.nzval
    _index_map[] = build_nzval_index_map(K.cscmatrix, nodes, nlist, dof, nnodes)
    _pattern_nodes[] = collect(nodes)

    # Preallocate CB tensors once — refilled each time, never reallocated
    cb = Vector{Array{Float64,4}}(undef, nnodes)
    @inbounds for i in 1:nnodes
        cb[i] = zeros(Float64, max_neighbors, dof, dof, dof)
    end
    _all_CB_tensors[] = cb
end

# =============================================================================
# Public interface
# =============================================================================

function init_model(nodes::AbstractVector{Int64},
                    material_parameter::Dict, block_id::Int64)
    Zero_Energy_Control.init_model(nodes, material_parameter, block_id)
    Data_Manager.get_max_rank() > 1 &&
        @error "Correspondence matrix based not implemented for parallel runs."
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
        trafo_fn = dof == 2 ? rotation_voigt_2D : rotation_voigt_3D
        for i in 1:nnodes
            @views C_voigt_trafo[i, :, :] = trafo_fn(rotation_tensor[i, :, :],
                                                     C_voigt[i, :, :])
        end
    else
        C_voigt_trafo = C_voigt
    end

    zStiff = nothing
    if haskey(Data_Manager.get_properties(1, "Material Model"), "Zero Energy Control") &&
       Data_Manager.get_properties(1, "Material Model")["Zero Energy Control"] ==
       "Global" &&
       include_zero_energy
        zStiff = Data_Manager.create_constant_node_tensor_field("Zero Energy Stiffness",
                                                                Float64, dof)
        Zero_Energy_Control.create_zero_energy_mode_stiffness!(nodes, dof, C_voigt_trafo,
                                                               inverse_shape_tensor, zStiff)
    end

    # Create the matrix once — this is the only call to init_stiffness_matrix
    Data_Manager.init_stiffness_matrix(nnodes, dof)
    K_sparse = Data_Manager.get_stiffness_matrix()
    inverse_nlist = Data_Manager.get_inverse_nlist()
    max_neighbors = maximum(number_of_neighbors)

    # Preallocate CB tensors for first assembly
    cb = Vector{Array{Float64,4}}(undef, nnodes)
    @inbounds for i in 1:nnodes
        cb[i] = zeros(Float64, max_neighbors, dof, dof, dof)
    end

    # First assembly: builds the sparsity pattern via variant 1
    @timeit "assemble_first" assemble_stiffness_with_zero_energy(K_sparse, nodes, nodes,
                                                                 dof, C_voigt_trafo,
                                                                 inverse_shape_tensor,
                                                                 number_of_neighbors,
                                                                 nlist, inverse_nlist,
                                                                 volume, bond_geometry,
                                                                 omega, bond_damage,
                                                                 zStiff, cb)

    # Build persistent state — index map requires a flushed CSC pattern
    @timeit "build_index_map" _build_persistent_state!(K_sparse, nodes, nlist,
                                                       dof, nnodes, max_neighbors)

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
    zStiff = include_zero_energy ?
             Data_Manager.get_field("Zero Energy Stiffness") : nothing

    if include_zero_energy && zStiff !== nothing
        Zero_Energy_Control.create_zero_energy_mode_stiffness!(nodes, dof, C_voigt,
                                                               inverse_shape_tensor, zStiff)
    end

    if _index_map_valid(nodes)
        # ── Fast path: write directly into persistent nzval ───────────────────
        # No init_stiffness_matrix, no flush!, no allocation.
        # _nzval[] is a direct reference to _K_sparse[].cscmatrix.nzval
        @timeit "precompute_CB" precompute_CB_all_nodes!(_all_CB_tensors[], nodes,
                                                         C_voigt, inverse_shape_tensor,
                                                         nlist, volume, bond_geometry_N,
                                                         omega, bond_damage, dof)
        @timeit "assemble_fast" begin
            fill!(_nzval[], 0.0)
            assemble_into_nzval!(_nzval[], _index_map[],
                                 _all_CB_tensors[], nlist, inverse_shape_tensor,
                                 volume, bond_geometry_N, omega, bond_damage,
                                 zStiff, dof, nnodes)
        end

        Data_Manager.set_stiffness_matrix(-_K_sparse[])

    else
        # ── Fallback: rebuild everything (e.g. new nodes in additive manufacturing)
        max_neighbors = maximum(number_of_neighbors)

        # Rebuild CB tensor storage for potentially larger node set
        cb = Vector{Array{Float64,4}}(undef, nnodes)
        @inbounds for i in 1:nnodes
            cb[i] = zeros(Float64, max_neighbors, dof, dof, dof)
        end

        # Recreate matrix with new sparsity pattern
        Data_Manager.init_stiffness_matrix(nnodes, dof)
        K_sparse = Data_Manager.get_stiffness_matrix()
        inverse_nlist = Data_Manager.get_inverse_nlist()

        @timeit "assemble_rebuild" assemble_stiffness_with_zero_energy(K_sparse,
                                                                       collect(1:nnodes),
                                                                       nodes, dof,
                                                                       C_voigt,
                                                                       inverse_shape_tensor,
                                                                       number_of_neighbors,
                                                                       nlist, inverse_nlist,
                                                                       volume,
                                                                       bond_geometry_N,
                                                                       omega, bond_damage,
                                                                       zStiff, cb)

        # Rebuild persistent state for subsequent fast-path calls
        @timeit "rebuild_index_map" _build_persistent_state!(K_sparse, nodes, nlist,
                                                             dof, nnodes, max_neighbors)

        Data_Manager.set_stiffness_matrix(-_K_sparse[])
    end
end

end # module
