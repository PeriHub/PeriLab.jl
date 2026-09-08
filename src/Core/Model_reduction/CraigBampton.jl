# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module CraigBampton
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
    return "Craig Bampton"
end

"""
    reduce_matrices(K, M_diag, m, s, n_modes = 1)

Craig-Bampton reduction of a stiffness and a lumped mass matrix.

The master degrees of freedom are kept as physical coordinates, the condensed ones are
represented by `n_modes` fixed-interface normal modes. The transformation is

    | x_m |   | I      0     | | x_m |
    | x_s | = | Phi_c  Phi_n | | eta |

with the constraint modes `Phi_c = -K_ss \\ K_sm`, the same static relation Guyan uses,
and the normal modes `Phi_n` from the fixed-interface eigenvalue problem
`K_ss Phi = M_ss Phi Lambda`. Guyan is the special case `n_modes = 0`: the modal part
is what carries the dynamics of the condensed degrees of freedom, which static condensation
drops.

The normal modes are mass normalised, so the lower right block of `K_r` holds the
eigenvalues and the one of `M_r` is the identity. The stiffness has no coupling between
the master and the modal part — the two mode sets are K-orthogonal — while the mass
does.

The reduced system has `length(m) + n_modes` degrees of freedom. The first `length(m)`
of them are the master degrees of freedom in the order given by `m`, the remaining ones
are modal coordinates without a physical meaning.

# Arguments
- `K::AbstractMatrix{Float64}`: Stiffness matrix
- `M_diag::Vector{Float64}`: Lumped mass matrix as a vector, one entry per degree of freedom
- `m::Vector{Int64}`: Indices of the master degrees of freedom
- `s::Vector{Int64}`: Indices of the condensed degrees of freedom
- `n_modes::Int64`: Number of fixed-interface normal modes to keep
# Returns
- `K_reduced::SparseMatrixCSC`: Reduced stiffness matrix, size `length(m) + n_modes`
- `M_reduced::SparseMatrixCSC`: Reduced mass matrix, same size
"""
function reduce_matrices(K::AbstractMatrix{Float64}, M_diag::Vector{Float64},
                         m::Vector{Int64}, s::Vector{Int64}, n_modes::Int64 = 1)
    nm = length(m)
    ns = length(s)

    if !isempty(intersect(m, s))
        throw(ArgumentError("Master and condensed index sets overlap."))
    end
    if n_modes < 0 || n_modes > ns
        throw(ArgumentError("n_modes = $n_modes, but only $ns condensed degrees of " *
                            "freedom are available."))
    end

    # Extract submatrices
    K_mm = K[m, m]
    K_ms = K[m, s]
    K_ss = K[s, s]
    K_sm = K[s, m]

    # 1. Constraint modes: Phi_c = -K_ss^(-1) * K_sm, identical to Guyan
    K_ss_fact = lu(K_ss)

    K_sm_dense = Matrix(K_sm)
    Phi_c = similar(K_sm_dense)
    ldiv!(Phi_c, K_ss_fact, K_sm_dense)
    Phi_c .*= -1.0

    M_ss_diag = Diagonal(M_diag[s])

    # 2. Fixed-interface normal modes. Symmetric() selects the symmetric solver, which
    #    returns real eigenvalues in ascending order, so the lowest modes come first.
    #    Dense matrices are needed here, eigen has no method for sparse input.
    factorization = eigen(Symmetric(Matrix(K_ss)), Symmetric(Matrix(M_ss_diag)))
    Phi_n = factorization.vectors[:, 1:n_modes]
    lambda_n = factorization.values[1:n_modes]

    # Mass normalisation: Phi_n' * M_ss * Phi_n = I. eigen does not guarantee this, and
    # without it the lower right blocks below lose their meaning.
    for k in 1:n_modes
        Phi_n[:, k] ./= sqrt(Phi_n[:, k]' * M_ss_diag * Phi_n[:, k])
    end

    # 3. Reduced matrices, block by block. Assembling the full system and multiplying it
    #    by the transformation matrix would allocate a matrix the size of the unreduced
    #    problem, which is what the reduction is meant to avoid.
    n_total = nm + n_modes

    K_reduced = zeros(n_total, n_total)
    K_reduced[1:nm, 1:nm] = K_mm + K_ms * Phi_c

    # The coupling blocks of the stiffness are zero, they are not left out by accident.
    # T' K T gives K_ms * Phi_n + Phi_c' * K_ss * Phi_n for the upper right block, and
    # since K_ss * Phi_c = -K_sm by the definition of the constraint modes, the second
    # term cancels the first exactly. Constraint modes and fixed-interface modes are
    # K-orthogonal; writing K_ms * Phi_n here would put twice the coupling into the
    # matrix. The mass has no such cancellation, its coupling block below is real.
    K_reduced[(nm + 1):end, (nm + 1):end] = diagm(lambda_n)

    # The coupling blocks M_ms and M_sm are zero because the mass matrix is lumped. With
    # a consistent mass matrix they would have to be included here.
    M_reduced = zeros(n_total, n_total)
    M_reduced[1:nm, 1:nm] = Phi_c' * M_ss_diag * Phi_c
    @inbounds for i in 1:nm
        M_reduced[i, i] += M_diag[m[i]]
    end
    M_reduced[1:nm, (nm + 1):end] = Phi_c' * M_ss_diag * Phi_n
    M_reduced[(nm + 1):end, 1:nm] = M_reduced[1:nm, (nm + 1):end]'
    M_reduced[(nm + 1):end, (nm + 1):end] = Matrix{Float64}(I, n_modes, n_modes)

    return sparse(K_reduced), sparse(M_reduced)
end

end
