# SPDX-FileCopyrightText: 2026 Christian Willberg
#
# SPDX-License-Identifier: BSD-3-Clause

module Matrix_Assembly

using Adapt: Adapt, @adapt_structure
using ArraysOfArrays: VectorOfVectors, flatview
using KernelAbstractions: @kernel, @index, @Const
using TimerOutputs: @timeit
using ..Sparse_CSR: CSRMatrix, launch!

export AssemblyPattern, create_assembly_pattern, create_matrix, create_element_buffer,
       assemble!, n_element_values

"""
    AssemblyPattern{Vi}

Result of the symbolic phase (README A4.0a). It knows nothing about the material model,
only where every value of every element matrix ends up.

Element `e` writes its dense matrix `K_e` (size `|row_dofs[e]| × |col_dofs[e]|`, column
major) into the flat element buffer starting at `ke_ptr[e]`. Every stored entry `p` of the
matrix is the sum of the buffer values `src[src_ptr[p]:(src_ptr[p+1]-1)]`.

Assembly is therefore a *gather*: one work item per matrix entry, no atomics, no write
conflicts, and the summation order is fixed (bitwise reproducible), also when elements
write into the same rows (correspondence, ordinary state-based, IGA-PD; README A3.6).

On an MPI rank the local matrix has the own rows followed by the ghost rows (partial sums
of elements that write into rows of neighbours, README A4.0b); the columns are the own
plus the ghost degrees of freedom in local numbering (README A4.2). Columns are sorted
within each row.

# Fields
- `n_rows::Int64`, `n_cols::Int64`: Matrix size
- `ke_ptr::Vi`: Start of each element matrix in the element buffer, length `n_elements + 1`
- `src_ptr::Vi`: Start of the contributions of each matrix entry in `src`, length `nnz + 1`
- `src::Vi`: Element buffer positions, sorted by matrix entry
- `rowptr::Vi`, `colval::Vi`: CSR pattern
"""
struct AssemblyPattern{Vi<:AbstractVector{Int64}}
    n_rows::Int64
    n_cols::Int64
    ke_ptr::Vi
    src_ptr::Vi
    src::Vi
    rowptr::Vi
    colval::Vi
end
@adapt_structure AssemblyPattern

n_elements(pattern::AssemblyPattern) = length(pattern.ke_ptr) - 1
n_element_values(pattern::AssemblyPattern) = Int64(Array(pattern.ke_ptr[end:end])[1] - 1)

"""
    create_assembly_pattern(row_dofs, col_dofs, n_rows, n_cols = n_rows)

Symbolic phase on the host. Builds the CSR pattern and the gather map from the degrees of
freedom each element writes into (`row_dofs`) and depends on (`col_dofs`).

# Arguments
- `row_dofs::VectorOfVectors{Int64}`: Row degrees of freedom per element
- `col_dofs::VectorOfVectors{Int64}`: Column degrees of freedom per element
- `n_rows::Int64`: Number of rows (own degrees of freedom)
- `n_cols::Int64`: Number of columns (own + ghost degrees of freedom)
# Returns
- `pattern::AssemblyPattern`: Host pattern; move it to a device with `Adapt.adapt`
"""
function create_assembly_pattern(row_dofs::AbstractVector{<:AbstractVector{Int64}},
                                 col_dofs::AbstractVector{<:AbstractVector{Int64}},
                                 n_rows::Int64, n_cols::Int64 = n_rows)
    n_el = length(row_dofs)
    if length(col_dofs) != n_el
        throw(ArgumentError("row_dofs and col_dofs must have one entry per element."))
    end
    @timeit "symbolic phase" begin
        ke_ptr = Vector{Int64}(undef, n_el + 1)
        ke_ptr[1] = 1
        for e in 1:n_el
            ke_ptr[e + 1] = ke_ptr[e] + length(row_dofs[e]) * length(col_dofs[e])
        end
        n_values = ke_ptr[end] - 1

        # one key per element value; sorting groups the values of each matrix entry
        keys = Vector{Int64}(undef, n_values)
        for e in 1:n_el
            fill_keys!(keys, ke_ptr[e], row_dofs[e], col_dofs[e], n_rows, n_cols)
        end
        src = sortperm(keys)
        rowptr, colval, src_ptr = compress_keys(keys, src, n_rows, n_cols)
    end
    return AssemblyPattern(n_rows, n_cols, ke_ptr, src_ptr, src, rowptr, colval)
end

function fill_keys!(keys::Vector{Int64}, start::Int64, rows::AbstractVector{Int64},
                    cols::AbstractVector{Int64}, n_rows::Int64, n_cols::Int64)
    n_r = length(rows)
    for (jc, c) in enumerate(cols), (ir, r) in enumerate(rows)
        if !(1 <= r <= n_rows)
            throw(ArgumentError("Row $r outside 1:$n_rows."))
        end
        if !(1 <= c <= n_cols)
            throw(ArgumentError("Column $c outside 1:$n_cols."))
        end
        keys[start + (jc - 1) * n_r + ir - 1] = (r - 1) * n_cols + c
    end
end

function compress_keys(keys::Vector{Int64}, src::Vector{Int64}, n_rows::Int64,
                       n_cols::Int64)
    rowptr = zeros(Int64, n_rows + 1)
    colval = Int64[]
    src_ptr = Int64[]
    last_key = -1
    for (q, k) in enumerate(src)
        key = keys[k]
        key == last_key && continue
        last_key = key
        row = div(key - 1, n_cols) + 1
        push!(colval, key - (row - 1) * n_cols)
        push!(src_ptr, q)
        rowptr[row + 1] += 1
    end
    push!(src_ptr, length(src) + 1)
    rowptr[1] = 1
    cumsum!(rowptr, rowptr)
    return rowptr, colval, src_ptr
end

"""
    create_matrix(pattern)

Allocates the matrix for a pattern on the pattern's backend (values zero).
"""
function create_matrix(pattern::AssemblyPattern; T::Type = Float64)
    nzval = similar(pattern.colval, T)
    fill!(nzval, zero(T))
    return CSRMatrix(pattern.n_rows, pattern.n_cols, pattern.rowptr, pattern.colval, nzval)
end

"""
    create_element_buffer(pattern)

Allocates the flat element buffer that the element kernels fill (README A3.4).
"""
function create_element_buffer(pattern::AssemblyPattern; T::Type = Float64)
    ke = similar(pattern.colval, T, n_element_values(pattern))
    fill!(ke, zero(T))
    return ke
end

@kernel function gather_kernel!(nzval, @Const(src_ptr), @Const(src), @Const(ke))
    p = @index(Global)
    s = zero(eltype(nzval))
    @inbounds for q in src_ptr[p]:(src_ptr[p + 1] - 1)
        s += ke[src[q]]
    end
    @inbounds nzval[p] = s
end

@kernel function mask_kernel!(nzval, @Const(rowptr), @Const(colval), @Const(active),
                              n_own::Int64)
    i = @index(Global)
    @inbounds for k in rowptr[i]:(rowptr[i + 1] - 1)
        c = colval[k]
        if !active[i] || !active[c]
            # identity on own inactive rows; ghost rows are partial sums of the owner's row
            nzval[k] = ifelse(c == i && i <= n_own, one(eltype(nzval)),
                              zero(eltype(nzval)))
        end
    end
end

"""
    assemble!(K, pattern, ke; active = nothing, n_own = K.n_rows)

Numerical phase: writes the sums of the element values into `K.nzval`. No allocation,
pattern unchanged (README A2.5, A4.4).

With the status mask `active` (one `Bool` per local degree of freedom, on the backend of
`K`, README A6a.2) rows and columns of inactive degrees of freedom are zero and the
diagonal of the own inactive rows (`≤ n_own`) is one, so the matrix stays regular.
"""
function assemble!(K::CSRMatrix, pattern::AssemblyPattern, ke::AbstractVector;
                   active = nothing, n_own::Int64 = K.n_rows)
    @timeit "assemble (gather)" launch!(gather_kernel!, K.nzval, length(K.nzval), K.nzval,
                                        pattern.src_ptr, pattern.src, ke)
    active === nothing ||
        @timeit "status mask" launch!(mask_kernel!, K.nzval, K.n_rows, K.nzval, K.rowptr,
                                      K.colval, active, n_own)
    return K
end

end
