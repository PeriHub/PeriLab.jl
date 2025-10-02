# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Correspondence_matrix_based

using SparseArrays

export init

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
- C: 4th order elasticity tensor (dof×dof×dof×dof)
- B: 3rd order strain-displacement tensor (dof×dof×dof)

Returns:
- CB: Contracted tensor (dof×dof×dof)

Mathematical background:
This contraction is fundamental in correspondence theory, relating
the material response (C) to the kinematic relationship (B).
"""
function contraction(C::Array{Float64,4}, B::Array{Float64,3})
    dof = size(C, 1)
    CB = zeros(dof, dof, dof)

    # CB_{mnq} = C_{mnop} * B_{opq}
    for m in 1:dof, n in 1:dof, q in 1:dof
        for o in 1:dof, p in 1:dof
            CB[m, n, q] += C[m, n, o, p] * B[o, p, q]
        end
    end

    return CB
end

"""
Create 3rd order B-tensor for strain-displacement relationship

The B-tensor relates displacement differences to strain components through:
B[m,n,p] = (omega_ik * V_k / 2) * (δ_mp * B1[n] + δ_np * B2[m])

Arguments:
- D_inv: Inverse shape tensor matrix (dof×dof)
- X_ik: Bond vector from node iID to node kID (dof×1)
- V_k: Volume of node kID
- omega_ik: Influence function value for bond iID-kID

Returns:
- B_tensor: 3rd order strain-displacement tensor (dof×dof×dof)

Mathematical formulation:
B1 = X_ik' * D_inv (row vector)
B2 = D_inv' * X_ik (column vector)
The tensor incorporates both volumetric and influence function weights.
"""
function create_B_tensor(D_inv::Matrix{Float64}, X_ik::Vector{Float64},
                         V_k::Float64, omega_ik::Float64)
    dof = length(X_ik)
    factor = omega_ik * V_k / 2
    B_tensor = zeros(dof, dof, dof)
    B1 = X_ik' * D_inv  # Row vector (1×dof)
    B2 = D_inv' * X_ik  # Column vector (dof×1)

    for m in 1:dof, n in 1:dof, p in 1:dof
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
- CB: Contracted elasticity-strain tensor (dof×dof×dof)
- D_inv: Inverse shape tensor matrix (dof×dof)
- X: Bond vector (dof×1)

Returns:
- R: Result matrix (dof×dof)

Mathematical operation:
v_n = (D_inv * X)_n
R[m,r] = Σ_n CB[m,n,r] * v[n]
"""
function multi(CB::Array{Float64,3}, D_inv::AbstractMatrix{Float64},
               X::AbstractVector{Float64})
    dof = length(X)
    v = D_inv * X                 # v_n = (D_inv * X)_n
    R = zeros(dof, dof)           # R[m,r]

    @inbounds for m in 1:dof, r in 1:dof
        s = 0.0
        for n in 1:dof
            s += CB[m, n, r] * v[n]
        end
        R[m, r] = s
    end
    return R
end

"""
Compute linearized operator K_ijk for correspondence theory

The linearized operator relates displacement at node kID to bond force along iID-jID:
K_{ijk,mo} = ω_{ij} * V_k * ω_{ik} * Σ_{n,p} CB_{ik,mno} * (D_i^{-1})_{np} * X_{ij,p}

Arguments:
- CB_ik: Contracted elasticity-strain tensor for bond iID-kID
- D_inv_i: Inverse shape tensor at node iID
- X_ij: Bond vector from node iID to node jID
- omega_ij: Influence function value for bond iID-jID
- V_k: Volume of node kID
- omega_ik: Influence function value for bond iID-kID

Returns:
- K_ijk: Linearized operator matrix (dof×dof)

This operator is central to the linearized correspondence formulation,
enabling efficient implicit solution methods.
"""
function compute_linearized_operator(CB_ik::Array{Float64,3}, D_inv_i::Matrix{Float64},
                                     X_ij::Vector{Float64}, omega_ij::Float64,
                                     V_k::Float64, omega_ik::Float64)
    dof = size(CB_ik, 1)
    K_ijk = zeros(dof, dof)

    # K_{ijk,mo} = ω_{ij} * V_k * ω_{ik} * Σ_{n,p} CB_{ik,mno} * (D_i^{-1})_{np} * X_{ij,p}
    for m in 1:dof, o in 1:dof
        sum_val = 0.0
        for n in 1:dof, p in 1:dof
            sum_val += CB_ik[m, n, o] * D_inv_i[n, p] * X_ij[p]
        end
        K_ijk[m, o] = omega_ij * V_k * omega_ik * sum_val
    end

    return K_ijk
end

"""
Compute all linearized operators for node iID

Calculates K_ijk operators for all neighbor combinations (jID,kID) of node iID.
This function creates the complete set of linearized operators needed
for stiffness matrix assembly.

Arguments:
- system: Peridynamics system containing geometry and connectivity
- iID: Central node index
- C_tensor: 4th order elasticity tensor
- D_inv: Inverse shape tensor for node iID

Returns:
- K_operators: Dictionary mapping (iID,jID,kID) tuples to K_ijk matrices

Note: This function computes B-tensors and CB contractions on-the-fly
for each bond combination.
"""
function compute_all_linearized_operators(iID::Int64, C_tensor::Array{Float64,4},
                                          D_inv::Matrix{Float64}, volume::Vector{Float64},
                                          bond_geometry::Vector{Vector{Vector{Float64}}},
                                          omega::Vector{Vector{Float64}},
                                          nlist::Vector{Vector{Int64}})
    K_operators = Dict{Tuple{Int64,Int64,Int64},Matrix{Float64}}()

    # Loop over all neighbors jID of node iID
    for (jID, jnID) in enumerate(nlist[iID])
        X_ij = bond_geometry[iID][jID]
        omega_ij = omega[iID][jID]

        # Loop over all neighbors kID of node iID
        for (kID, knID) in enumerate(nlist[iID])
            V_k = volume[knID]
            omega_ik = omega[iID][kID]
            X_ik = bond_geometry[iID][kID]

            # Create B-tensor with omega_ik * V_k / 2 factor
            B_ik = create_B_tensor(D_inv, X_ik, V_k, omega_ik)
            CB = contraction(C_tensor, B_ik)

            # Compute K_ijk operator
            K_ijk = compute_linearized_operator(CB, D_inv,
                                                X_ij, omega_ij, V_k, omega_ik)
            K_operators[(iID, jnID, knID)] = K_ijk
        end
    end

    return K_operators
end

"""
	assemble_stiffness_contributions_sparse(nnodes, dof, C_Voigt,
										   inverse_shape_tensor, nlist, volume,
										   bond_geometry, omega)

Assemble stiffness matrix directly as sparse matrix to save memory.
Uses I, J, V vectors to build sparse structure incrementally.

Arguments:
- nnodes: Number of nodes
- dof: Degrees of freedom per node
- C_Voigt: Voigt notation elasticity tensor per node (nnodes × voigt_dim × voigt_dim)
- inverse_shape_tensor: Inverse shape tensor per node (nnodes × dof × dof)
- nlist: Neighborhood list
- volume: Volume per node
- bond_geometry: Bond geometry vectors per node
- omega: Influence function values per node

Returns:
- K: Sparse global stiffness matrix (nnodes*dof × nnodes*dof)
"""

function assemble_stiffness_contributions_sparse(nnodes::Int64,
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

"""
Assemble global stiffness matrix for peridynamics system

Main function for assembling the complete stiffness matrix using
correspondence theory. Supports different material models and
dofensionalities.

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

"""
  init_model(datamanager::Module, nodes::AbstractVector{Int64}, material_parameter::Dict)

Initializes the material model.

# Arguments
  - `datamanager::Data_manager`: Datamanager.
  - `nodes::AbstractVector{Int64}`: List of block nodes.
  - `material_parameter::Dict(String, Any)`: Dictionary with material parameter.

# Returns
  - `datamanager::Data_manager`: Datamanager.
"""
function init_model(datamanager::Module,
                    nodes::AbstractVector{Int64},
                    material_parameter::Dict)
    dof = datamanager.get_dof()

    hooke_matrix = datamanager.create_constant_node_field("Hooke Matrix", Float64,
                                                          Int64((dof * (dof + 1)) / 2),
                                                          VectorOrMatrix = "Matrix")
    symmetry = get(material_parameter, "Symmetry", "default")::String
    for iID in nodes
        @views hooke_matrix[iID, :,
        :] = get_Hooke_matrix(datamanager,
                                                          material_parameter,
                                                          symmetry,
                                                          dof,
                                                          iID)
    end
    return datamanager
end

end # module
