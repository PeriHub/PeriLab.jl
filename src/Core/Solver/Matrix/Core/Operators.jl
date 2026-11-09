# SPDX-FileCopyrightText: 2026 Christian Willberg
#
# SPDX-License-Identifier: BSD-3-Clause

module Operators

using KernelAbstractions: @kernel, @index, @Const
using LinearAlgebra: LinearAlgebra, mul!
using TimerOutputs: @timeit
using ..Sparse_CSR: CSRMatrix, launch!
using ..Dof_Sets: DofSet

export ConstrainedOperator, ShiftedOperator, constrained_rhs!, apply_prescribed!

"""
    ConstrainedOperator(K, prescribed)

Operator of the system with prescribed degrees of freedom `p` and free ones `f`, applied
directly on the full matrix without forming sub-matrices (README A6.1, A6.4):

    y = M_f K M_f x + M_p x

i.e. `K_ff` on the free part and the identity on the prescribed part. The matrix stays
unchanged, so the prescribed set can change between load steps (README A6.3).

The same operator restricted to the complement of a master set is `K_ss` of a model
reduction (README A4a.3). On an MPI rank `n_own` is the number of own rows; ghost rows
(partial sums sent back to the owners) get no identity part.
"""
struct ConstrainedOperator{M<:CSRMatrix,S<:DofSet}
    K::M
    prescribed::S
    n_own::Int64
end
function ConstrainedOperator(K::CSRMatrix, prescribed::DofSet; n_own::Int64 = K.n_rows)
    return ConstrainedOperator(K, prescribed, n_own)
end

Base.size(A::ConstrainedOperator) = size(A.K)
Base.size(A::ConstrainedOperator, d::Integer) = size(A.K, d)
Base.eltype(A::ConstrainedOperator) = eltype(A.K)

@kernel function constrained_spmv_kernel!(y, @Const(rowptr), @Const(colval), @Const(nzval),
                                          @Const(prescribed), @Const(x), n_own::Int64)
    i = @index(Global)
    if prescribed[i]
        # identity on own prescribed rows; ghost rows contribute nothing to them
        @inbounds y[i] = i <= n_own ? x[i] : zero(eltype(y))
    else
        s = zero(eltype(y))
        @inbounds for k in rowptr[i]:(rowptr[i + 1] - 1)
            c = colval[k]
            s += ifelse(prescribed[c], zero(s), nzval[k] * x[c])
        end
        @inbounds y[i] = s
    end
end

function LinearAlgebra.mul!(y::AbstractVector, A::ConstrainedOperator, x::AbstractVector)
    K = A.K
    @timeit "SpMV constrained" launch!(constrained_spmv_kernel!, y, K.n_rows, y, K.rowptr,
                                       K.colval, K.nzval, A.prescribed.mask, x,
                                       A.n_own)
    return y
end

"""
    ShiftedOperator(A::ConstrainedOperator, shift)

`A` plus a diagonal `shift` on the free own degrees of freedom (e.g. `K + a₀ M` of an
implicit time integration with a lumped mass):

    y = A x + M_f · shift .* x   (own rows only)
"""
struct ShiftedOperator{O<:ConstrainedOperator,V<:AbstractVector}
    A::O
    shift::V
end
Base.size(S::ShiftedOperator) = size(S.A)
Base.size(S::ShiftedOperator, d::Integer) = size(S.A, d)
Base.eltype(S::ShiftedOperator) = eltype(S.A)

@kernel function shift_kernel!(y, @Const(prescribed), @Const(shift), @Const(x))
    i = @index(Global)
    @inbounds y[i] += ifelse(prescribed[i], zero(eltype(y)), shift[i] * x[i])
end

function LinearAlgebra.mul!(y::AbstractVector, S::ShiftedOperator, x::AbstractVector)
    mul!(y, S.A, x)
    @timeit "shift" launch!(shift_kernel!, y, S.A.n_own, y, S.A.prescribed.mask, S.shift,
                            x)
    return y
end

@kernel function rhs_kernel!(b, @Const(rowptr), @Const(colval), @Const(nzval),
                             @Const(prescribed), @Const(u), @Const(f))
    i = @index(Global)
    if prescribed[i]
        @inbounds b[i] = zero(eltype(b))
    else
        s = f[i]
        @inbounds for k in rowptr[i]:(rowptr[i + 1] - 1)
            c = colval[k]
            s -= ifelse(prescribed[c], nzval[k] * u[c], zero(s))
        end
        @inbounds b[i] = s
    end
end

"""
    constrained_rhs!(b, K, prescribed, u, f)

Right-hand side for `ConstrainedOperator` (README A6.4), homogeneous form:

    b = M_f (f − K M_p u)

`u` holds the prescribed values at the prescribed degrees of freedom (other entries are
ignored). One fused kernel, no temporary vector. Because `b_p = 0` and the operator is the
identity on the prescribed part, every Krylov vector is exactly zero there; the prescribed
values are inserted exactly afterwards with `apply_prescribed!`.
"""
function constrained_rhs!(b::AbstractVector, K::CSRMatrix, prescribed::DofSet,
                          u::AbstractVector, f::AbstractVector)
    @timeit "rhs" launch!(rhs_kernel!, b, K.n_rows, b, K.rowptr, K.colval, K.nzval,
                          prescribed.mask, u, f)
    return b
end

@kernel function apply_prescribed_kernel!(x, @Const(prescribed), @Const(u))
    i = @index(Global)
    @inbounds x[i] = prescribed[i] ? u[i] : x[i]
end

"""
    apply_prescribed!(x, prescribed, u)

Inserts the prescribed values `u` into the solution `x` (exact, no rounding).
"""
function apply_prescribed!(x::AbstractVector, prescribed::DofSet, u::AbstractVector)
    @timeit "apply prescribed" launch!(apply_prescribed_kernel!, x, length(x), x,
                                       prescribed.mask, u)
    return x
end

end
