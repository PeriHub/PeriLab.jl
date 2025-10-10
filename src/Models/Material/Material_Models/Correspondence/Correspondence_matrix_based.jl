# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Correspondence_matrix_based
include("../../Material_Basis.jl")
include("../../../../Support/Helpers.jl")
using LinearAlgebra
using .Helpers: get_fourth_order
include("../Zero_Energy_Control/global_control.jl")
using .Global_zero_energy_control: create_zero_energy_mode_stiffness!
using .Material_Basis: get_Hooke_matrix
using SparseArrays

export init_model
export add_zero_energy_stiff!
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
Adds Zero-Energy-Mode stabilization to the stiffness matrix (in-place).

# Arguments
- `K::SparseMatrixCSC{Float64, Int64}`: Stiffness matrix (will be modified)
- `nnodes::Int`: Number of nodes
- `dof::Int`: Degrees of freedom per node
- `C_voigt::Matrix{Float64}`: Elasticity tensor in Voigt notation (3×3 for 2D, 6×6 for 3D)
- `inverse_shape_tensor::Vector{Matrix{Float64}}`: Inverse shape tensors D^(-1) for each node
- `nlist::Vector{Vector{Int}}`: Neighborhood list [i] -> [j1, j2, ...]
- `volume::Vector{Float64}`: Node volumes
- `bond_geometry::Vector{Vector{Vector{Float64}}}`: Bond vectors X_ij for each node
- `omega::Vector{Vector{Float64}}`: Influence functions ω for each bond
"""

function add_zero_energy_stiff!(K::SparseMatrixCSC{Float64,Int64},
                                nodes::AbstractVector{Int64},
                                dof::Int64,
                                zStiff::Array{Float64,3},
                                inverse_shape_tensor::Array{Float64,3},
                                nlist::Vector{Vector{Int}},
                                volume::Vector{Float64},
                                bond_geometry::Vector{Vector{Vector{Float64}}},
                                omega::Vector{Vector{Float64}})
    ndof_total = nodes[end] * dof

    # COO format
    I_indices = Int[]
    J_indices = Int[]
    values = Float64[]

    # Pre-compute scalar product matrix S for each node
    # S[i][p,q] = X_ip^T * D_inv_i * X_iq
    S_matrices = Vector{Matrix{Float64}}(undef, nodes[end])
    D_inv_X = Vector{Vector{Vector{Float64}}}(undef, nodes[end])

    for i in nodes
        neighbors = nlist[i]
        n_neighbors = length(neighbors)

        @views D_inv_i = inverse_shape_tensor[i, :, :]

        # Pre-compute D_inv * X for all bonds
        D_inv_X[i] = [D_inv_i * bond_geometry[i][idx] for idx in 1:n_neighbors]

        # Compute S matrix
        S_matrices[i] = zeros(n_neighbors, n_neighbors)
        for p in 1:n_neighbors
            for q in 1:n_neighbors
                S_matrices[i][p, q] = dot(bond_geometry[i][p], D_inv_X[i][q])
            end
        end
    end

    # STEP 1: Fill rows for neighbor nodes (j, k, l, m, ...)
    # For each node i and each bond ij
    for i in nodes
        neighbors = nlist[i]
        n_neighbors = length(neighbors)

        if n_neighbors == 0
            continue
        end

        @views Z_i = zStiff[i, :, :]
        V_i = volume[i]

        for (idx_j, j) in enumerate(neighbors)
            if i == j
                continue
            end

            # Scalar sum for bond ij: Σ_k ω_ik V_k s_kj
            scalar_sum_j = 0.0
            for idx_k in 1:n_neighbors
                k = neighbors[idx_k]
                scalar_sum_j += omega[i][idx_k] * volume[k] * S_matrices[i][idx_k, idx_j]
            end

            # Row j: K_ji
            K_ji = -V_i * (1.0 - scalar_sum_j) * Z_i
            add_block_to_coo!(I_indices, J_indices, values, K_ji, j, i, dof)

            # Row j: K_jk for other neighbors k ≠ j
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

    # STEP 2: Fill row for central node i
    # Each row i gets contributions from ALL bonds of node i
    for i in nodes
        neighbors = nlist[i]
        n_neighbors = length(neighbors)

        if n_neighbors == 0
            continue
        end

        @views Z_i = zStiff[i, :, :]
        V_i = volume[i]

        # For each neighbor j, compute K_ij (sum over all bonds)
        for (idx_j, j) in enumerate(neighbors)
            if i == j
                continue
            end

            K_ij_total = zeros(dof, dof)

            # Each bond contributes to K_ij
            for (idx_bond, bond_node) in enumerate(neighbors)
                if bond_node == i
                    continue
                end

                if idx_bond == idx_j
                    # Contribution from bond ij to K_ij
                    s_jj = S_matrices[i][idx_j, idx_j]
                    K_ij_total .-= V_i * (1.0 - omega[i][idx_j] * volume[j] * s_jj) * Z_i
                else
                    # Contribution from bond ik (k≠j) to K_ij
                    s_kj = S_matrices[i][idx_bond, idx_j]
                    K_ij_total .+= V_i * omega[i][idx_j] * volume[j] * s_kj * Z_i
                end
            end

            add_block_to_coo!(I_indices, J_indices, values, K_ij_total, i, j, dof)
        end
    end

    # Create stabilization matrix from off-diagonal blocks
    K_stab = sparse(I_indices, J_indices, values, ndof_total, ndof_total)

    # STEP 3: Diagonal blocks K_ii = -Σ_(ℓ≠i) K_iℓ
    for i in nodes
        K_ii_block = zeros(dof, dof)

        # Sum all off-diagonal entries in row i
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

        # Add diagonal block to K
        for d1 in 1:dof
            for d2 in 1:dof
                row = (i-1)*dof + d1
                col = (i-1)*dof + d2
                K[row, col] += K_ii_block[d1, d2]
            end
        end
    end

    # Add off-diagonal stabilization to K
    K .+= K_stab

    return nothing
end

"""
Computes Z = C : D^(-1) for 2D or 3D.
"""
function compute_Z_tensor(C_voigt::AbstractMatrix{Float64}, D_inv::AbstractMatrix{Float64},
                          dof::Int)
    if dof == 2
        # 2D: Voigt notation [ε11, ε22, 2ε12] -> [σ11, σ22, σ12]
        # Z_ij = C_ijkl D^(-1)_kl
        Z = zeros(2, 2)
        #@warn "must be implementated for general C"
        # Simplified for isotropic materials in 2D:
        Z[1, 1] = C_voigt[1, 1]*D_inv[1, 1] + C_voigt[1, 2]*D_inv[2, 2] +
                  C_voigt[1, 3]*(D_inv[1, 2] + D_inv[2, 1])
        Z[1, 2] = C_voigt[1, 1]*D_inv[1, 2] + C_voigt[1, 2]*D_inv[2, 1] +
                  C_voigt[1, 3]*(D_inv[1, 1] + D_inv[2, 2])
        Z[2, 1] = C_voigt[2, 1]*D_inv[1, 1] + C_voigt[2, 2]*D_inv[2, 2] +
                  C_voigt[2, 3]*(D_inv[1, 2] + D_inv[2, 1])
        Z[2, 2] = C_voigt[2, 1]*D_inv[1, 2] + C_voigt[2, 2]*D_inv[2, 1] +
                  C_voigt[2, 3]*(D_inv[1, 1] + D_inv[2, 2])

        return Z

    elseif dof == 3
        # 3D: Voigt notation [ε11, ε22, ε33, 2ε23, 2ε13, 2ε12]
        Z = zeros(3, 3)

        # Full tensor contraction
        for i in 1:3
            for j in 1:3
                Z[i, j] = C_voigt[1, 1]*D_inv[i, j] + C_voigt[1, 2]*D_inv[i, j] +
                          C_voigt[1, 3]*D_inv[i, j] +
                          C_voigt[2, 1]*D_inv[i, j] + C_voigt[2, 2]*D_inv[i, j] +
                          C_voigt[2, 3]*D_inv[i, j] +
                          C_voigt[3, 1]*D_inv[i, j] + C_voigt[3, 2]*D_inv[i, j] +
                          C_voigt[3, 3]*D_inv[i, j]
            end
        end

        return Z
    else
        error("Only dof=2 or dof=3 supported")
    end
end
"""
Computes Γ_ij = I - Σ_k ω_ik V_k X_ik (D^(-1) X_ij)^T.

Note: This is a MATRIX equation, not scalar! The term X_ik (D^(-1) X_ij)^T
is an outer product creating a dof×dof matrix.
"""
function compute_gamma(bond_vectors::Vector{Vector{Float64}},
                       X_ij::Vector{Float64},
                       omega_i::Vector{Float64},
                       volume::Vector{Float64},
                       neighbors::Vector{Int},
                       D_inv::AbstractMatrix{Float64},
                       dof::Int)
    Γ = Matrix{Float64}(I, dof, dof)  # Identity matrix

    # Precompute D^(-1) * X_ij (column vector)
    D_inv_X_ij = D_inv * X_ij

    # Σ_k ω_ik V_k X_ik (D^(-1) X_ij)^T
    # This is an outer product: X_ik * (D_inv_X_ij)^T creates a dof×dof matrix
    for (idx_k, k) in enumerate(neighbors)
        X_ik = bond_vectors[idx_k]

        # Outer product: X_ik (column) times (D_inv_X_ij)^T (row)
        # Results in a dof×dof matrix
        outer_prod = X_ik * D_inv_X_ij'

        # Subtract from Γ
        Γ .-= omega_i[idx_k] * volume[k] * outer_prod
    end

    return Γ
end

"""
Adds a dof×dof block to COO sparse matrix format.
"""
function add_block_to_coo!(I_indices::Vector{Int},
                           J_indices::Vector{Int},
                           values::Vector{Float64},
                           block::Matrix{Float64},
                           i::Int,
                           j::Int,
                           dof::Int)
    for d1 in 1:dof
        for d2 in 1:dof
            if abs(block[d1, d2]) > 1e-14  # Only add non-zero entries
                push!(I_indices, (i-1)*dof + d1)
                push!(J_indices, (j-1)*dof + d2)
                push!(values, block[d1, d2])
            end
        end
    end
end

"""
	assemble_stiffness_contributions_sparse(nodes, dof, C_Voigt,
										   inverse_shape_tensor, nlist, volume,
										   bond_geometry, omega)

Assemble stiffness matrix directly as sparse matrix to save memory.
Uses I, J, V vectors to build sparse structure incrementally.

Arguments:
- nodes: list of nodes
- dof: Degrees of freedom per node
- C_Voigt: Voigt notation elasticity tensor per node (nnodes × voigt_dim × voigt_dim)
- inverse_shape_tensor: Inverse shape tensor per node (nnodes × dof × dof)
- nlist: Neighborhood list
- volume: Volume per node
- bond_geometry: Bond geometry vectors per node
- omega: Influence function values per node

Returns:
- K: Sparse global stiffness matrix (nodes[end]*dof × nodes[end]*dof)
"""

function assemble_stiffness_contributions_sparse(nodes::AbstractVector{Int64},
                                                 dof::Int64,
                                                 C_Voigt::Array{Float64,3},
                                                 inverse_shape_tensor::Array{Float64,
                                                                             3},
                                                 nlist::Vector{Vector{Int64}},
                                                 volume::Vector{Float64},
                                                 bond_geometry::Vector{Vector{Vector{Float64}}},
                                                 omega::Vector{Vector{Float64}})
    K = zeros(nodes[end]*dof, nodes[end]*dof)
    # Process each node i
    for i in nodes
        D_inv = inverse_shape_tensor[i, :, :]
        @views C_tensor = get_fourth_order(C_Voigt[i, :, :], dof)

        K_ijk = compute_all_linearized_operators(i, C_tensor, D_inv,
                                                 volume, bond_geometry,
                                                 omega, nlist)

        # Process each neighbor j of node i
        for (j_idx, j) in enumerate(nlist[i])
            omega_ij = omega[i][j_idx]
            V_i = volume[i]
            V_j = volume[j]

            # Local K_ij matrix for this bond
            K_ij = zeros(dof, nodes[end] * dof)

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
                for n in 1:(nodes[end] * dof)
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
    if datamanager.get_max_rank()>1
        @error "Correspondence matrix based not implemented for parallel runs."
    end

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

    zStiff = datamanager.create_constant_node_field("Zero Energy Stiffness",
                                                    Float64,
                                                    dof,
                                                    VectorOrMatrix = "Matrix")
    @warn "Memory inefficient version for testing, K_sparse is not initialized and sparsed after full initial"
    @warn "damages are not included yet in stiffness matrix"
    return datamanager
end
function compute_model(datamanager::Module,
                       nodes::AbstractVector{Int64})
    dof = datamanager.get_dof()
    C_voigt = datamanager.get_field("Hooke Matrix")
    inverse_shape_tensor = datamanager.get_field("Inverse Shape Tensor")
    nlist = datamanager.get_nlist()
    volume = datamanager.get_field("Volume")
    #bond_geometry = datamanager.get_field("Bond Geometry")
    bond_geometry_N = datamanager.get_field("Deformed Bond Geometry", "N")
    omega = datamanager.get_field("Influence Function")

    # verify K
    K_sparse = assemble_stiffness_contributions_sparse(nodes,              # nnodes
                                                       dof,                # dof
                                                       C_voigt,            # C_Voigt
                                                       inverse_shape_tensor, # inverse_shape_tensor
                                                       nlist,              # nlist
                                                       volume,      # volume
                                                       bond_geometry_N,      # bond_geometry
                                                       omega)

    zStiff = datamanager.get_field("Zero Energy Stiffness")
    create_zero_energy_mode_stiffness!(nodes, dof, C_voigt,
                                       inverse_shape_tensor, zStiff)
    add_zero_energy_stiff!(K_sparse,
                           nodes,              # nnodes
                           dof,                # dof
                           zStiff,            # C_Voigt
                           inverse_shape_tensor, # inverse_shape_tensor
                           nlist,              # nlist
                           volume,      # volume
                           bond_geometry_N,      # bond_geometry
                           omega)
    datamanager.set_stiffness_matrix(K_sparse)
    return K_sparse
end
end # module
