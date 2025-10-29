# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Correspondence_matrix_based
using LinearAlgebra
using ProgressBars
using ...Material_Basis: get_Hooke_matrix
using ....Helpers: get_fourth_order, progress_bar
using ...Global_Zero_Energy_Control: create_zero_energy_mode_stiffness!
using SparseArrays
using StaticArrays: @MMatrix
export init_model
export add_zero_energy_stiff!
export compute_bond_force

"""
Peridynamic Correspondence Stiffness Matrix Assembly with Voigt Notation

Memory-optimized implementation using Dict-based sparse column accumulation.
Reduces memory from O(dof × n_total) to O(dof × n_neighbors) per iteration.
"""

kronecker_delta(x::Int64, y::Int64) = ==(x, y)

function contraction!(C::Array{Float64,4}, B::Array{Float64,3}, dof::Int64,
                      CB::AbstractArray{Float64,3})
    fill!(CB, 0.0)

    @inbounds for q in 1:dof, p in 1:dof, o in 1:dof
        B_opq = B[o, p, q]
        for n in 1:dof, m in 1:dof
            CB[m, n, q] += C[m, n, o, p] * B_opq
        end
    end
end

function create_B_tensor!(D_inv::AbstractMatrix{Float64}, X_ik::Vector{Float64},
                          V_k::Float64, omega_ik::Float64, dof::Int64,
                          B_tensor::Array{Float64,3})
    factor = omega_ik * V_k * 0.5

    @inbounds @fastmath if dof == 2
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

    else  # dof == 3
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
    end

    return nothing
end

function compute_linearized_operator(CB_ik::AbstractArray{Float64,3},
                                     D_inv_i::Matrix{Float64},
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
                                D_inv::Matrix{Float64},
                                neighbors::Vector{Int64},
                                volume::Vector{Float64},
                                bond_geometry_i::Vector{Vector{Float64}},
                                omega_i::Vector{Float64},
                                dof::Int64,
                                B_ik::Array{Float64,3})
    @inbounds for (k_idx, k) in enumerate(neighbors)
        V_k = volume[k]
        omega_ik = omega_i[k_idx]
        X_ik = bond_geometry_i[k_idx]
        create_B_tensor!(D_inv, X_ik, V_k, omega_ik, dof, B_ik)
        @views contraction!(C_tensor, B_ik, dof, CB_k[k_idx, :, :, :])
    end
    return nothing
end

"""
	create_sparse_matrix_mapping(nodes, nlist, dof)

Create mapping structure for efficient sparse matrix assembly.
"""
function create_sparse_matrix_mapping(nodes::AbstractVector{Int64},
                                      nlist::Vector{Vector{Int64}},
                                      dof::Int64)
    edge_pairs = Set{Tuple{Int64,Int64}}()

    @inbounds for i in nodes
        push!(edge_pairs, (i, i))

        for j in nlist[i]
            if i != j
                push!(edge_pairs, (i, j))
                push!(edge_pairs, (j, i))
            end

            for k in nlist[i]
                if k != i && k != j
                    push!(edge_pairs, (j, k))
                    push!(edge_pairs, (k, j))
                end
            end
        end
    end

    nnz = length(edge_pairs) * dof * dof

    I = Vector{Int64}(undef, nnz)
    J = Vector{Int64}(undef, nnz)
    V = Vector{Float64}(undef, nnz)
    fill!(V, 0.0)

    coo_map = Dict{Tuple{Int64,Int64},Int64}()

    idx = 1
    for (i_node, j_node) in sort(collect(edge_pairs))
        coo_map[(i_node, j_node)] = idx

        @inbounds for n in 1:dof
            for m in 1:dof
                I[idx] = (i_node - 1) * dof + m
                J[idx] = (j_node - 1) * dof + n
                idx += 1
            end
        end
    end

    return coo_map, I, J, V
end

"""
	add_block_mapped!(V, coo_map, block, i_node, j_node, dof)

Add dof×dof block to sparse matrix values using mapping structure.
"""
@inline function add_block_mapped!(V::Vector{Float64},
                                   coo_map::Dict{Tuple{Int64,Int64},Int64},
                                   block::AbstractMatrix{Float64},
                                   i_node::Int64,
                                   j_node::Int64,
                                   dof::Int64)
    start_idx = get(coo_map, (i_node, j_node), 0)
    if start_idx == 0
        return
    end

    idx = start_idx
    @inbounds for n in 1:dof
        for m in 1:dof
            V[idx] += block[m, n]
            idx += 1
        end
    end
end

"""
	init_assemble_stiffness_optimized(...)

Memory-optimized assembly using Dict for sparse column accumulation.

Key insight: Instead of K_ij_j = zeros(dof, n_total), use
K_ij_j = Dict{Int, Matrix{Float64}} to store only non-zero columns.

Memory reduction: O(dof × n_total) → O(dof² × n_neighbors)
"""
function init_assemble_stiffness_optimized(nodes::AbstractVector{Int64},
                                           dof::Int64,
                                           C_Voigt::Array{Float64,3},
                                           inverse_shape_tensor::Array{Float64,3},
                                           number_of_neighbors::Vector{Int64},
                                           nlist::Vector{Vector{Int64}},
                                           volume::Vector{Float64},
                                           bond_geometry::Vector{Vector{Vector{Float64}}},
                                           omega::Vector{Vector{Float64}})
    n_total = nodes[end] * dof

    @info "Building sparse matrix structure..."
    coo_map, I, J, V = create_sparse_matrix_mapping(nodes, nlist, dof)

    # Dict: column_node → dof×dof block
    K_ij_j = Dict{Int64,Matrix{Float64}}()
    K_ij_i = Dict{Int64,Matrix{Float64}}()

    if dof == 2
        K_block = @MMatrix zeros(2, 2)
        B_ik = zeros(2, 2, 2)
    else
        K_block = @MMatrix zeros(3, 3)
        B_ik = zeros(3, 3, 3)
    end

    max_neighbors = maximum(number_of_neighbors)
    CB_k_buffer = zeros(max_neighbors, dof, dof, dof)

    @info "Assembling stiffness contributions..."
    iter = progress_bar(0, length(nodes)-1, false)

    @inbounds for iter_dx in iter
        i = nodes[iter_dx]
        D_inv = inverse_shape_tensor[i, :, :]
        @views C_tensor = get_fourth_order(C_Voigt[i, :, :], dof)
        V_i = volume[i]
        ni = nlist[i]
        n_neighbors = length(ni)

        CB_k = @view CB_k_buffer[1:n_neighbors, :, :, :]

        precompute_CB_tensors!(CB_k, C_tensor, D_inv, ni, volume,
                               bond_geometry[i], omega[i], dof, B_ik)

        for (j_idx, j) in enumerate(ni)
            omega_ij = omega[i][j_idx]
            V_j = volume[j]
            X_ij = bond_geometry[i][j_idx]

            # Clear dictionaries for new (i,j) pair
            empty!(K_ij_j)
            empty!(K_ij_i)

            # Accumulate contributions from all k neighbors
            for (k_idx, k) in enumerate(ni)
                omega_ik = omega[i][k_idx]
                V_k = volume[k]

                @views compute_linearized_operator(CB_k[k_idx, :, :, :],
                                                   D_inv, X_ij,
                                                   omega_ij, V_k, omega_ik,
                                                   dof, K_block)

                # Initialize blocks if needed
                if !haskey(K_ij_j, i)
                    K_ij_j[i] = zeros(dof, dof)
                    K_ij_i[i] = zeros(dof, dof)
                end
                if !haskey(K_ij_j, k)
                    K_ij_j[k] = zeros(dof, dof)
                    K_ij_i[k] = zeros(dof, dof)
                end

                # Accumulate: subtract from column i, add to column k
                # This matches the original logic:
                # K_ij_j[m, i_offset+o] -= val
                # K_ij_j[m, k_offset+o] += val
                K_ij_j[i] .-= K_block
                K_ij_i[i] .-= K_block
                K_ij_j[k] .+= K_block
                K_ij_i[k] .+= K_block
            end

            # Write accumulated blocks to sparse matrix
            # Row j gets entries from K_ij_j (scaled by -V_i)
            for (col_node, block) in K_ij_j
                scaled_ji = -V_i .* block
                add_block_mapped!(V, coo_map, scaled_ji, j, col_node, dof)
            end

            # Row i gets entries from K_ij_i (scaled by V_j)
            for (col_node, block) in K_ij_i
                scaled_ij = V_j .* block
                add_block_mapped!(V, coo_map, scaled_ij, i, col_node, dof)
            end
        end
    end

    return I, J, V, n_total
end

function add_zero_energy_stiff!(K::SparseMatrixCSC{Float64,Int64},
                                nodes::AbstractVector{Int64},
                                dof::Int64,
                                zStiff::Array{Float64,3},
                                inverse_shape_tensor::Array{Float64,3},
                                nlist::Vector{Vector{Int}},
                                volume::Vector{Float64},
                                bond_geometry::Vector{Vector{Vector{Float64}}},
                                bond_damage::Vector{Vector{Float64}},
                                omega::Vector{Vector{Float64}})
    ndof_total = nodes[end] * dof

    I_indices = Int[]
    J_indices = Int[]
    values = Float64[]

    S_matrices = Vector{Matrix{Float64}}(undef, nodes[end])
    D_inv_X = Vector{Vector{Vector{Float64}}}(undef, nodes[end])

    for i in nodes
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

    for i in nodes
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
                scalar_sum_j += omega[i][idx_k] * volume[k] * S_matrices[i][idx_k, idx_j]
            end

            K_ji = -V_i * (1.0 - scalar_sum_j) * Z_i
            add_block_to_coo!(I_indices, J_indices, values, K_ji, j, i, dof)

            for (idx_k, k) in enumerate(neighbors)
                if k == i || k == j
                    continue
                end

                s_value = omega[i][idx_k] * volume[k] * S_matrices[i][idx_k, idx_j]
                K_jk = -V_i * s_value * Z_i
                add_block_to_coo!(I_indices, J_indices, values, K_jk, j, k, dof)
            end
        end
    end

    for i in nodes
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
                    K_ij_total .-= V_i * (1.0 - omega[i][idx_j] * volume[j] * s_jj) * Z_i
                else
                    s_kj = S_matrices[i][idx_bond, idx_j]
                    K_ij_total .+= V_i * omega[i][idx_j] * volume[j] * s_kj * Z_i
                end
            end

            add_block_to_coo!(I_indices, J_indices, values, K_ij_total, i, j, dof)
        end
    end

    K_stab = sparse(I_indices, J_indices, values, ndof_total, ndof_total)

    for i in nodes
        K_ii_block = zeros(dof, dof)

        for d1 in 1:dof
            row = (i-1)*dof + d1
            for j in nodes
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

    @info "Initializing stiffness matrix (optimized version)"
    index_x, index_y, vals,
    total_dof = init_assemble_stiffness_optimized(nodes,
                                                  dof,
                                                  C_voigt,
                                                  inverse_shape_tensor,
                                                  number_of_neighbors,
                                                  nlist,
                                                  volume,
                                                  bond_geometry,
                                                  omega)

    datamanager.init_stiffness_matrix(index_x, index_y, vals, total_dof)
end

function compute_model(datamanager::Module,
                       nodes::AbstractVector{Int64})
    dof = datamanager.get_dof()
    C_voigt = datamanager.get_field("Hooke Matrix")
    inverse_shape_tensor = datamanager.get_field("Inverse Shape Tensor")
    nlist = datamanager.get_nlist()
    volume = datamanager.get_field("Volume")
    bond_geometry_N = datamanager.get_field("Deformed Bond Geometry", "N")
    omega = datamanager.get_field("Influence Function")
    bond_damage = datamanager.get_field("Bond Damage", "NP1")

    zStiff = datamanager.get_field("Zero Energy Stiffness")
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

    return K_sparse
end

"""
	compute_bond_force(bond_force, K, u, nodes, nlist, dof)

Compute bond forces from stiffness matrix and displacements.

Calculates the linearized bond force component: f_ij = K_ij * (u_j - u_i)

# Arguments
- `bond_force`: Output vector for bond forces (modified in-place)
- `K`: Global stiffness matrix
- `u`: Displacement field (nnodes × dof)
- `nodes`: Vector of node indices
- `nlist`: Neighborhood list
- `dof`: Degrees of freedom per node

# Note
This computes only the linearized force component, not the complete
peridynamic bond force which includes nonlinear terms.
"""
function compute_bond_force(bond_force::Vector{Vector{Vector{Float64}}},
                            K::SparseMatrixCSC{Float64,Int64}, u::Matrix{Float64},
                            nodes::Vector{Int64}, nlist::Vector{Vector{Int64}}, dof::Int64)
    @inbounds for i in nodes
        ni = nlist[i]
        i_offset = (i-1)*dof

        for (j_idx, j) in enumerate(ni)
            if i == j
                continue
            end

            j_offset = (j-1)*dof
            f_ij = zeros(dof)

            # Extract K_ij block and compute force
            for m in 1:dof
                for n in 1:dof
                    f_ij[m] += K[i_offset+m, j_offset+n] * (u[j, n] - u[i, n])
                end
            end

            bond_force[i][j_idx] = f_ij
        end
    end
end

end # module
