# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module CorrespondenceStiffnessMatrix

include("correspondence_functions.jl")

"""
Peridynamic Correspondence Stiffness Matrix Assembly with Voigt Notation

This module implements corrected stiffness matrix assembly for peridynamics
correspondence theory based on the LaTeX formulation. It provides efficient
assembly of global stiffness matrices for implicit time integration and
static equilibrium solutions.

Key features:
- Linearized operator computation (K_ijk)
- Tensor contractions for elasticity
- Sparse matrix assembly
- 2D plane stress/strain support
- Memory-efficient implementations
"""

"""
Kronecker delta function
Returns 1 if x == y, 0 otherwise
"""
kronecker_delta(x::Int64, y::Int64) = ==(x, y)

"""
Compute tensor contraction C : B for elasticity and strain tensors

Performs the contraction: CB_{mnq} = C_{mnop} * B_{opq}
This operation combines the 4th order elasticity tensor C with the
3rd order strain-displacement tensor B.

Arguments:
- C: 4th order elasticity tensor (dim×dim×dim×dim)
- B: 3rd order strain-displacement tensor (dim×dim×dim)

Returns:
- CB: Contracted tensor (dim×dim×dim)

Mathematical background:
This contraction is fundamental in correspondence theory, relating
the material response (C) to the kinematic relationship (B).
"""
function contraction(C::Array{Float64,4}, B::Array{Float64,3})
    dim = size(C, 1)
    CB = zeros(dim, dim, dim)

    # CB_{mnq} = C_{mnop} * B_{opq}
    for m in 1:dim, n in 1:dim, q in 1:dim
        for o in 1:dim, p in 1:dim
            CB[m, n, q] += C[m, n, o, p] * B[o, p, q]
        end
    end

    return CB
end

"""
Create 3rd order B-tensor for strain-displacement relationship

The B-tensor relates displacement differences to strain components through:
B[m,n,p] = (ω_ik * V_k / 2) * (δ_mp * B1[n] + δ_np * B2[m])

Arguments:
- D_inv: Inverse shape tensor matrix (dim×dim)
- X_ik: Bond vector from node i to node k (dim×1)
- V_k: Volume of node k
- ω_ik: Influence function value for bond i-k

Returns:
- B_tensor: 3rd order strain-displacement tensor (dim×dim×dim)

Mathematical formulation:
B1 = X_ik' * D_inv (row vector)
B2 = D_inv' * X_ik (column vector)
The tensor incorporates both volumetric and influence function weights.
"""
function create_B_tensor(D_inv::Matrix{Float64}, X_ik::Vector{Float64},
                         V_k::Float64, ω_ik::Float64)
    dim = length(X_ik)
    factor = ω_ik * V_k / 2
    B_tensor = zeros(dim, dim, dim)
    B1 = X_ik' * D_inv  # Row vector (1×dim)
    B2 = D_inv' * X_ik  # Column vector (dim×1)

    for m in 1:dim, n in 1:dim, p in 1:dim
        B_tensor[m, n, p] = factor * (kronecker_delta(m, p) * B1[n] +
                             kronecker_delta(n, p) * B2[m])
    end

    return B_tensor
end

"""
Corrected multiplication: R[m,r] = sum_n CB[m,n,r] * v[n]

This function computes the matrix-vector product where v = D_inv * X.
It's used in the linearized operator computation.

Arguments:
- CB: Contracted elasticity-strain tensor (dim×dim×dim)
- D_inv: Inverse shape tensor matrix (dim×dim)
- X: Bond vector (dim×1)

Returns:
- R: Result matrix (dim×dim)

Mathematical operation:
v_n = (D_inv * X)_n
R[m,r] = Σ_n CB[m,n,r] * v[n]
"""
function multi(CB::Array{Float64,3}, D_inv::AbstractMatrix{Float64},
               X::AbstractVector{Float64})
    dim = length(X)
    v = D_inv * X                 # v_n = (D_inv * X)_n
    R = zeros(dim, dim)           # R[m,r]

    @inbounds for m in 1:dim, r in 1:dim
        s = 0.0
        for n in 1:dim
            s += CB[m, n, r] * v[n]
        end
        R[m, r] = s
    end
    return R
end

"""
Compute linearized operator K_ijk for correspondence theory

The linearized operator relates displacement at node k to bond force along i-j:
K_{ijk,mo} = ω_{ij} * V_k * ω_{ik} * Σ_{n,p} CB_{ik,mno} * (D_i^{-1})_{np} * X_{ij,p}

Arguments:
- CB_ik: Contracted elasticity-strain tensor for bond i-k
- D_inv_i: Inverse shape tensor at node i
- X_ij: Bond vector from node i to node j
- ω_ij: Influence function value for bond i-j
- V_k: Volume of node k
- ω_ik: Influence function value for bond i-k

Returns:
- K_ijk: Linearized operator matrix (dim×dim)

This operator is central to the linearized correspondence formulation,
enabling efficient implicit solution methods.
"""
function compute_linearized_operator(CB_ik::Array{Float64,3}, D_inv_i::Matrix{Float64},
                                     X_ij::Vector{Float64}, ω_ij::Float64,
                                     V_k::Float64, ω_ik::Float64)
    dim = size(CB_ik, 1)
    K_ijk = zeros(dim, dim)

    # K_{ijk,mo} = ω_{ij} * V_k * ω_{ik} * Σ_{n,p} CB_{ik,mno} * (D_i^{-1})_{np} * X_{ij,p}
    for m in 1:dim, o in 1:dim
        sum_val = 0.0
        for n in 1:dim, p in 1:dim
            sum_val += CB_ik[m, n, o] * D_inv_i[n, p] * X_ij[p]
        end
        K_ijk[m, o] = ω_ij * V_k * ω_ik * sum_val
    end

    return K_ijk
end

"""
Compute all linearized operators for node i

Calculates K_ijk operators for all neighbor combinations (j,k) of node i.
This function creates the complete set of linearized operators needed
for stiffness matrix assembly.

Arguments:
- system: Peridynamics system containing geometry and connectivity
- i: Central node index
- C_tensor: 4th order elasticity tensor
- D_inv: Inverse shape tensor for node i

Returns:
- K_operators: Dictionary mapping (i,j,k) tuples to K_ijk matrices

Note: This function computes B-tensors and CB contractions on-the-fly
for each bond combination.
"""
function compute_all_linearized_operators(system, i::Int, C_tensor, D_inv)
    K_operators = Dict{Tuple{Int,Int,Int},Matrix{Float64}}()

    # Loop over all neighbors j of node i
    for (j_idx, j) in enumerate(system.nlist[i])
        X_ij = system.positions[j, :] - system.positions[i, :]
        ω_ij = system.omega[i][j_idx]

        # Loop over all neighbors k of node i
        for (k_idx, k) in enumerate(system.nlist[i])
            V_k = system.volume[k]
            ω_ik = system.omega[i][k_idx]
            X_ik = system.positions[k, :] - system.positions[i, :]

            # Create B-tensor with ω_ik * V_k / 2 factor
            B_ik = create_B_tensor(D_inv, X_ik, V_k, ω_ik)
            CB = contraction(C_tensor, B_ik)

            # Compute K_ijk operator
            K_ijk = compute_linearized_operator(CB, D_inv,
                                                X_ij, ω_ij, V_k, ω_ik)
            K_operators[(i, j, k)] = K_ijk
        end
    end

    return K_operators
end

"""
Assemble stiffness contributions for specified nodes

Implements the stiffness matrix assembly based on linearized correspondence theory:
K_ik = -Σ_j ω_ij * V_j * X_ij^T * D_i^(-1) * [C : B_ik] * V_k

This function assembles the global stiffness matrix by computing contributions
from all bonds and their interactions through the linearized operators.

Arguments:
- K: Global stiffness matrix to be filled (modified in-place)
- nodes: Array of node indices to process
- model: Peridynamics model containing all system data
- C_tensor: 4th order elasticity tensor

Returns:
- Sparse representation of the assembled stiffness matrix

Mathematical background:
The stiffness matrix relates nodal forces to nodal displacements through
the linearized correspondence theory. Each entry K[i,j] represents the
force contribution at degree of freedom i due to unit displacement at
degree of freedom j.
"""
function assemble_stiffness_contributions!(K::Matrix{Float64}, nodes, model,
                                           C_tensor::Array{Float64,4})
    dim = model.dof
    n_nodes = model.n_nodes
    nnodes = model.n_nodes
    shape_tensor = zeros(nnodes, model.dof, model.dof)
    inverse_shape_tensor = zeros(nnodes, model.dof, model.dof)

    # Compute shape tensors for all nodes
    for iID in nodes
        compute_shape_tensor!(shape_tensor, model.volume, model.bond_damage[iID],
                              model.undeformed_bond[iID], iID)
        compute_shape_tensor_inverse!(inverse_shape_tensor, shape_tensor, iID)
    end

    # Process each node i
    for i in nodes
        D_inv = inverse_shape_tensor[i, :, :]
        K_ijk = compute_all_linearized_operators(model, i, C_tensor, D_inv)

        # Process each neighbor j of node i
        for (j_idx, j) in enumerate(model.nlist[i])
            X_ij = model.positions[j, :] - model.positions[i, :]
            ω_ij = model.omega[i][j_idx]
            V_i = model.volume[i]
            V_j = model.volume[j]

            # Local K_ij matrix for this bond
            K_ij = zeros(dim, n_nodes * dim)

            # Sum over all k in neighborhood of i
            for (k_idx, k) in enumerate(model.nlist[i])
                X_ik = model.positions[k, :] - model.positions[i, :]
                ω_ik = model.omega[i][k_idx]
                V_k = model.volume[k]

                # Get pre-computed K_ijk block
                # This already contains temp = CB_ik : D_i^(-1) ⊗ X_ij
                K_block = K_ijk[(i, j, k)]

                # Fill K_ij based on K_block
                for m in 1:dim
                    for o in 1:dim
                        K_ij[m, (i-1)*dim+o] -= ω_ij * V_k * ω_ik * K_block[m, o]
                        K_ij[m, (k-1)*dim+o] += ω_ij * V_k * ω_ik * K_block[m, o]
                    end
                end
            end

            # Add K_ij to global matrix K
            for m in 1:dim
                for n in 1:(n_nodes * dim)
                    K[(j-1)*dim+m, n] -= K_ij[m, n] * V_i
                    K[(i-1)*dim+m, n] += K_ij[m, n] * V_j
                end
            end
        end
    end

    return sparse(K)
end

"""
Assemble global stiffness matrix for peridynamics system

Main function for assembling the complete stiffness matrix using
correspondence theory. Supports different material models and
dimensionalities.

Arguments:
- model: Peridynamics model containing all system data
- nodes: Array of node indices to process
- material_type: Material model type (default: "plane_stress")
- is_3d: Boolean indicating 3D analysis (default: false)

Returns:
- K: Sparse global stiffness matrix

The resulting matrix can be used for:
- Implicit time integration schemes
- Static equilibrium analysis
- Modal analysis and eigenvalue problems
- Linearized stability analysis

Currently supports:
- 2D plane stress analysis
- Future extensions: plane strain, 3D analysis
"""
function assemble_global_stiffness_matrix(model, nodes,
                                          material_type::String = "plane_stress",
                                          is_3d::Bool = false)
    n_nodes = model.n_nodes
    dim = model.dof

    # Material-specific elasticity tensor setup
    # TODO: Add support for 3D and plane strain
    #if is_3d
    #    C_voigt = elasticity_matrix_3d(mat)
    #    C_tensor = voigt_to_tensor4_3d(C_voigt)
    #else
    #    if material_type == "plane_stress"
    C_voigt = elasticity_matrix_2d_plane_stress(model.material)
    #    else
    #        C_voigt = elasticity_matrix_2d_plane_strain(mat)
    #    end
    C_tensor = voigt_to_tensor4_2d(C_voigt)
    #end

    # Initialize global stiffness matrix
    K = zeros(n_nodes * dim, n_nodes * dim)

    # Assemble contributions from each node
    return assemble_stiffness_contributions!(K, nodes, model, C_tensor)
end

end # module
