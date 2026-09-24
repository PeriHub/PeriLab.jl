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
    reduce_matrices(K, M_diag, r, l, n_modes = 1)

Guyan reduction (static condensation) [GuyanRJ1965](@cite) of a stiffness and a lumped mass matrix.

The condensed degrees of freedom are expressed through the retained ones by the static
relation `x_l = T x_r` with `T = -K_ll \\ K_lr`, which gives

    K_R = K_rr + K_rl T
    M_R = M_rr + T' M_ll T

The inertia of the condensed degrees of freedom is only carried over through `T`, so the
reduction is exact for the static case and approximate for dynamics — the error grows
with frequency. Use Craig-Bampton where the dynamic behaviour matters.

The mass term omits the coupling blocks `M_rl T` and `T' M_lr`, which is correct here
because the mass matrix is lumped and therefore has no off-diagonal entries. With a
consistent mass matrix those blocks would have to be included.

A retained degree of freedom with no bond into the condensed region has an all-zero
column in `K_lr`, so its column of `T` is `K_ll^-1 * 0 = 0`, and — independently — an
all-zero row in `K_rl` makes its row of `K_rl T` zero regardless of `T`. `coupling` is the
union of both, so `K_rl T` and `T' M_ll T` are computed only on that `nc x nc` block
instead of densely over all `nr` retained degrees of freedom; every other entry of
`K_R`/`M_R` is exactly `K_rr`/`Diagonal(M_diag[r])`. This holds without assuming `K`
symmetric.

`n_modes` is not used; it is part of the signature so that all reduction schemes share
one interface.

# Arguments
- `K::AbstractMatrix{Float64}`: Stiffness matrix
- `M_diag::Vector{Float64}`: Lumped mass matrix as a vector, one entry per degree of freedom
- `r::Vector{Int64}`: Indices of the retained degrees of freedom
- `l::Vector{Int64}`: Indices of the condensed degrees of freedom
- `n_modes::Int64`: Unused, kept for interface compatibility
# Returns
- `K_reduced::SparseMatrixCSC`: Reduced stiffness matrix, size `length(r)`
- `M_reduced::SparseMatrixCSC`: Reduced mass matrix, same size
"""
function reduce_matrices(K::AbstractMatrix{Float64}, M_diag::Vector{Float64},
                         r::Vector{Int64}, l::Vector{Int64}, n_modes::Int64 = 1)
    nr = length(r)
    nl = length(l)

    if !isempty(intersect(r, l))
        throw(ArgumentError("Retained and condensed index sets overlap."))
    end

    K_rr = K[r, r]
    K_rl = K[r, l]
    K_ll = K[l, l]
    K_lr = K[l, r]

    coupling = union(findall(!iszero, vec(sum(abs, K_rl; dims = 2))),
                     findall(!iszero, vec(sum(abs, K_lr; dims = 1))))
    nc = length(coupling)

    K_ll_fact = lu(K_ll)
    T = Matrix(K_lr[:, coupling])
    ldiv!(K_ll_fact, T)
    T .*= -1.0

    Kbb_fill = Matrix{Float64}(undef, nc, nc)
    mul!(Kbb_fill, K_rl[coupling, :], T)

    rows, columns, values = findnz(K_rr)
    @inbounds for (j, cj) in enumerate(coupling), (i, ci) in enumerate(coupling)
        push!(rows, ci)
        push!(columns, cj)
        push!(values, Kbb_fill[i, j])
    end
    K_reduced = sparse(rows, columns, values, nr, nr)

    # M_ll is diagonal, so scaling T by sqrt(M_ll) turns T' M_ll T into a plain product.
    root_mass = sqrt.(M_diag[l])
    @inbounds for j in 1:nc, i in 1:nl
        T[i, j] *= root_mass[i]
    end
    Mbb_fill = Matrix{Float64}(undef, nc, nc)
    mul!(Mbb_fill, T', T)

    mrows = collect(1:nr)
    mcolumns = collect(1:nr)
    mvalues = M_diag[r]
    @inbounds for (j, cj) in enumerate(coupling), (i, ci) in enumerate(coupling)
        push!(mrows, ci)
        push!(mcolumns, cj)
        push!(mvalues, Mbb_fill[i, j])
    end
    M_reduced = sparse(mrows, mcolumns, mvalues, nr, nr)

    return K_reduced, M_reduced
end

end
