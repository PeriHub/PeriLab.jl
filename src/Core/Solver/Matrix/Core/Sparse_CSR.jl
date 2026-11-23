# SPDX-FileCopyrightText: 2026 Christian Willberg
#
# SPDX-License-Identifier: BSD-3-Clause

module Sparse_CSR

using Adapt: Adapt, @adapt_structure
using KernelAbstractions: @kernel, @index, @Const, get_backend, synchronize
using LinearAlgebra: LinearAlgebra
using TimerOutputs: @timeit
using SparseArrays: SparseMatrixCSC, sparse

export CSRMatrix, n_rows, n_cols, nnz_csr, to_sparse, launch!

"""
    CSRMatrix{Tv,Ti,Vv,Vi}

Sparse matrix in compressed sparse row format. The arrays may live on the host (`Vector`) or
on a device (e.g. `CuVector`); move it with `Adapt.adapt(CuArray, K)`.

The pattern (`rowptr`, `colval`) is fixed after the symbolic phase, only `nzval` changes
during a numerical assembly (README A2.5).

# Fields
- `n_rows::Int64`: Number of rows
- `n_cols::Int64`: Number of columns
- `rowptr::Vi`: Row pointers, length `n_rows + 1`, 1-based
- `colval::Vi`: Column index of every stored entry
- `nzval::Vv`: Value of every stored entry
"""
struct CSRMatrix{Tv,Ti<:Integer,Vv<:AbstractVector{Tv},Vi<:AbstractVector{Ti}}
    n_rows::Int64
    n_cols::Int64
    rowptr::Vi
    colval::Vi
    nzval::Vv
end
@adapt_structure CSRMatrix

n_rows(K::CSRMatrix) = K.n_rows
n_cols(K::CSRMatrix) = K.n_cols
nnz_csr(K::CSRMatrix) = length(K.nzval)
Base.size(K::CSRMatrix) = (K.n_rows, K.n_cols)
Base.size(K::CSRMatrix, d::Integer) = d == 1 ? K.n_rows : (d == 2 ? K.n_cols : 1)
Base.eltype(::CSRMatrix{Tv}) where {Tv} = Tv

"""
    launch!(kernel, y, n, args...)

Launches a KernelAbstractions kernel on the backend of `y` over `n` work items and waits
for it, so that the result can be used directly by the host (e.g. inside Krylov.jl).
"""
function launch!(kernel, y::AbstractArray, n::Integer, args...)
    n == 0 && return nothing
    backend = get_backend(y)
    kernel(backend)(args...; ndrange = n)
    synchronize(backend)
    return nothing
end

@kernel function spmv_kernel!(y, @Const(rowptr), @Const(colval), @Const(nzval), @Const(x))
    i = @index(Global)
    s = zero(eltype(y))
    @inbounds for k in rowptr[i]:(rowptr[i + 1] - 1)
        s += nzval[k] * x[colval[k]]
    end
    @inbounds y[i] = s
end

"""
    LinearAlgebra.mul!(y, K::CSRMatrix, x)

Computes `y = K x` with one work item per row (no write conflicts).
"""
function LinearAlgebra.mul!(y::AbstractVector, K::CSRMatrix, x::AbstractVector)
    @timeit "SpMV" launch!(spmv_kernel!, y, K.n_rows, y, K.rowptr, K.colval, K.nzval, x)
    return y
end

"""
    to_sparse(K::CSRMatrix)

Copies the matrix into a host `SparseMatrixCSC`. For reference solutions and tests only.
"""
function to_sparse(K::CSRMatrix)
    rowptr = Array(K.rowptr)
    colval = Array(K.colval)
    nzval = Array(K.nzval)
    rows = Vector{Int64}(undef, length(nzval))
    for i in 1:(K.n_rows)
        for k in rowptr[i]:(rowptr[i + 1] - 1)
            rows[k] = i
        end
    end
    return sparse(rows, colval, nzval, K.n_rows, K.n_cols)
end

end
