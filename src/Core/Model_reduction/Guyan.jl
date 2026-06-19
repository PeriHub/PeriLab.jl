# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Guyan
using LinearAlgebra
using SparseArrays
export model_reduction_name
export reduce_matrices

function model_reduction_name()
    return "Static Condensation"
end
function reduce_matrices(K::AbstractMatrix{Float64}, M_diag::Vector{Float64},
                         m::Vector{Int64}, s::Vector{Int64}, n_modes::Int64 = 1)
    nm = length(m)
    ns = length(s)

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

    return sparse(K_reduced), sparse(M_reduced)
end

end
