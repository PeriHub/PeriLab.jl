# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Guyan
using LinearAlgebra
using SparseArrays
export model_reduction_name
export reduce_matrices

"""
    model_reduction_name()

Name of the reduction scheme, used to identify it in the input deck and in log output.

# Returns
- `::String`: The name of the scheme
"""
function model_reduction_name()
    return "Static Condensation"
end

"""
    reduce_matrices(K, M_diag, m, s, n_modes = 1)

Guyan reduction (static condensation) [GuyanRJ1965](@cite) of a stiffness and a lumped mass matrix.

The slave degrees of freedom are expressed through the master ones by the static
relation `x_s = T x_m` with `T = -K_ss \\ K_sm`, which gives

    K_r = K_mm + K_ms T
    M_r = M_mm + T' M_ss T

The inertia of the slave degrees of freedom is only carried over through `T`, so the
reduction is exact for the static case and approximate for dynamics — the error grows
with frequency. Use Craig-Bampton where the dynamic behaviour matters.

The mass term omits the coupling blocks `M_ms T` and `T' M_sm`, which is correct here
because the mass matrix is lumped and therefore has no off-diagonal entries. With a
consistent mass matrix those blocks would have to be included.

`n_modes` is not used; it is part of the signature so that all reduction schemes share
one interface.

# Arguments
- `K::AbstractMatrix{Float64}`: Stiffness matrix
- `M_diag::Vector{Float64}`: Lumped mass matrix as a vector, one entry per degree of freedom
- `m::Vector{Int64}`: Indices of the master degrees of freedom
- `s::Vector{Int64}`: Indices of the slave degrees of freedom
- `n_modes::Int64`: Unused, kept for interface compatibility
# Returns
- `K_reduced::SparseMatrixCSC`: Reduced stiffness matrix, size `length(m)`
- `M_reduced::SparseMatrixCSC`: Reduced mass matrix, same size
"""
function reduce_matrices(K::AbstractMatrix{Float64}, M_diag::Vector{Float64},
                         m::Vector{Int64}, s::Vector{Int64}, n_modes::Int64 = 1)
    nm = length(m)
    ns = length(s)

    if !isempty(intersect(m, s))
        throw(ArgumentError("Master and slave index sets overlap."))
    end

    # Extract submatrices
    K_mm = K[m, m]
    K_ms = K[m, s]
    K_ss = K[s, s]
    K_sm = K[s, m]

    # Compute transformation matrix: T = -K_ss^(-1) * K_sm
    K_ss_fact = lu(K_ss)

    K_sm_dense = Matrix(K_sm)
    T = similar(K_sm_dense)
    ldiv!(T, K_ss_fact, K_sm_dense)
    T .*= -1.0

    # Stiffness reduction: K_reduced = K_mm + K_ms * T
    temp_K = zeros(nm, nm)
    mul!(temp_K, K_ms, T)
    K_reduced = K_mm + temp_K

    # Mass reduction: M_reduced = M_mm + T^T * M_ss * T
    # Step 1: temp_M = Diagonal(M_ss) * T
    M_ss_diag = Diagonal(M_diag[s])
    temp_M = zeros(ns, nm)
    mul!(temp_M, M_ss_diag, T)  # Diagonal * Matrix multiplication

    # Step 2: M_reduced = T^T * temp_M
    M_reduced = zeros(nm, nm)
    mul!(M_reduced, T', temp_M)

    # Step 3: Add M_mm diagonal
    @inbounds for i in 1:nm
        M_reduced[i, i] += M_diag[m[i]]
    end

    # Both are symmetric by construction; rounding leaves a small asymmetry that some
    # solvers reject.
    K_reduced = (K_reduced + K_reduced') / 2
    M_reduced = (M_reduced + M_reduced') / 2

    return sparse(K_reduced), sparse(M_reduced)
end

end
