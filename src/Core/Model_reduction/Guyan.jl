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

The condensed degrees of freedom are expressed through the master ones by the static
relation `x_s = T x_m` with `T = -K_ss \\ K_sm`, which gives

    K_r = K_mm + K_ms T
    M_r = M_mm + T' M_ss T

The inertia of the condensed degrees of freedom is only carried over through `T`, so the
reduction is exact for the static case and approximate for dynamics — the error grows
with frequency. Use Craig-Bampton where the dynamic behaviour matters.

The mass term omits the coupling blocks `M_ms T` and `T' M_sm`, which is correct here
because the mass matrix is lumped and therefore has no off-diagonal entries. With a
consistent mass matrix those blocks would have to be included.

A master degree of freedom with no bond into the condensed region has an all-zero column
in `K_sm`, so its column of `T` is `K_ss^-1 * 0 = 0`, and — independently — an all-zero
row in `K_ms` makes its row of `K_ms T` zero regardless of `T`. `coupling` is the union
of both, so `K_ms T` and `T' M_ss T` are computed only on that `nc x nc` block instead of
densely over all `nm` master degrees of freedom; every other entry of `K_r`/`M_r` is
exactly `K_mm`/`Diagonal(M_diag[m])`. This holds without assuming `K` symmetric.

`n_modes` is not used; it is part of the signature so that all reduction schemes share
one interface.

# Arguments
- `K::AbstractMatrix{Float64}`: Stiffness matrix
- `M_diag::Vector{Float64}`: Lumped mass matrix as a vector, one entry per degree of freedom
- `m::Vector{Int64}`: Indices of the master degrees of freedom
- `s::Vector{Int64}`: Indices of the condensed degrees of freedom
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
        throw(ArgumentError("Master and condensed index sets overlap."))
    end

    K_mm = K[m, m]
    K_ms = K[m, s]
    K_ss = K[s, s]
    K_sm = K[s, m]

    coupling = union(findall(!iszero, vec(sum(abs, K_ms; dims = 2))),
                     findall(!iszero, vec(sum(abs, K_sm; dims = 1))))
    nc = length(coupling)

    K_ss_fact = lu(K_ss)
    T = Matrix(K_sm[:, coupling])
    ldiv!(K_ss_fact, T)
    T .*= -1.0

    Kbb_fill = Matrix{Float64}(undef, nc, nc)
    mul!(Kbb_fill, K_ms[coupling, :], T)

    rows, columns, values = findnz(K_mm)
    @inbounds for (j, cj) in enumerate(coupling), (i, ci) in enumerate(coupling)
        push!(rows, ci)
        push!(columns, cj)
        push!(values, Kbb_fill[i, j])
    end
    K_reduced = sparse(rows, columns, values, nm, nm)

    # M_ss is diagonal, so scaling T by sqrt(M_ss) turns T' M_ss T into a plain product.
    root_mass = sqrt.(M_diag[s])
    @inbounds for j in 1:nc, i in 1:ns
        T[i, j] *= root_mass[i]
    end
    Mbb_fill = Matrix{Float64}(undef, nc, nc)
    mul!(Mbb_fill, T', T)

    m_rows = collect(1:nm)
    m_columns = collect(1:nm)
    m_values = M_diag[m]
    @inbounds for (j, cj) in enumerate(coupling), (i, ci) in enumerate(coupling)
        push!(m_rows, ci)
        push!(m_columns, cj)
        push!(m_values, Mbb_fill[i, j])
    end
    M_reduced = sparse(m_rows, m_columns, m_values, nm, nm)

    return K_reduced, M_reduced
end

end
