# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module CraigBampton_cascade
using LinearAlgebra
using SparseArrays
using TimerOutputs: @timeit
using Arpack: eigs
export model_reduction_name
export reduce_matrices
export reduce_stage
export group_condensed_by_level

"""
    model_reduction_name()

Name of the reduction scheme, used to identify it in the input deck and in log output.

# Returns
- `::String`: The name of the scheme
"""
function model_reduction_name()
    return "Craig Bampton Cascade"
end

"""
    nonzero_rows(A)

Indices of the rows of a matrix that hold at least one entry.

A condensed degree of freedom whose row is empty carries no stiffness at all — an
isolated node, or one whose bonds were all cut by a bond filter. Such a row makes `Kcc`
singular, and the eigenvalue problem cannot separate its rigid body mode from the modes
being looked for. They are removed before the factorization and the eigenvalue problem,
and the corresponding entries of the modes stay zero.

# Arguments
- `A::AbstractMatrix`: The matrix
# Returns
- `::Vector{Int64}`: Indices of the non-empty rows
"""
function nonzero_rows(A::SparseMatrixCSC)
    occupied = falses(size(A, 1))
    rows = rowvals(A)
    values = nonzeros(A)
    @inbounds for column in axes(A, 2)
        for index in nzrange(A, column)
            if values[index] != 0.0
                occupied[rows[index]] = true
            end
        end
    end
    return findall(occupied)
end

nonzero_rows(A::AbstractMatrix) = findall(!iszero, vec(sum(abs, A; dims = 2)))

"""
    nonzero_columns(A)

Indices of the columns of a matrix that hold at least one entry.

For `Kcr` these are the retained degrees of freedom with a bond into the condensed
region — the coupling layer. Every other column is empty, so its recovery mode is zero
and never has to be solved for. That is the difference between one triangular solve per
retained degree of freedom and one per coupling degree of freedom, and it also bounds
where the reduced matrices can fill in: `Krc * B` has entries only in the rows and
columns of the coupling layer, everything else keeps the sparsity of `Krr`.

# Arguments
- `A::AbstractMatrix`: The matrix
# Returns
- `::Vector{Int64}`: Indices of the non-empty columns
"""
function nonzero_columns(A::SparseMatrixCSC)
    values = nonzeros(A)
    columns = Int64[]
    for column in axes(A, 2)
        for index in nzrange(A, column)
            if values[index] != 0.0
                push!(columns, column)
                break
            end
        end
    end
    return columns
end

nonzero_columns(A::AbstractMatrix) = findall(!iszero, vec(sum(abs, A; dims = 1)))

"""
    _warn_not_pd(err = nothing)

Warns that the condensed stiffness failed to factorize as positive definite, shared
between the dense in-place path and the sparse fallback below.
"""
function _warn_not_pd(err = nothing)
    @warn "Cholesky of the condensed stiffness failed, falling back to LU. Kcc is not " *
          "positive definite: check whether every condensed node is tied to the " *
          "retained region -- an isolated node or a detached cluster has a rigid body " *
          "mode." exception=err
end

"""
    factorize_condensed!(Kcc)

Factorization of the condensed stiffness, Cholesky where possible.

With the interface held fixed `Kcc` is symmetric positive definite. A sparse Cholesky is
faster to build than an LU and considerably faster to apply to a dense block of right
hand sides, because CHOLMOD goes through supernodal BLAS3 kernels while UMFPACK works
column by column. With one right hand side per coupling degree of freedom that
difference dominates the whole reduction.

A failure is informative in itself: `Kcc` is then not positive definite, which for a
fixed interface points at a condensed region falling apart into pieces not tied to any
retained degree of freedom — or, in peridynamics, simply at a stiffness that is not
symmetric.

Dense `Kcc` (a `StridedMatrix` of a BLAS float type -- what a single, small condensed
region typically comes down to) factorizes in place: LAPACK's `potrf` touches only the
upper triangle it is asked for, so on failure the lower triangle plus the saved diagonal
is enough to rebuild the original matrix for the LU fallback, without a second copy of
`Kcc`. Anything else -- in particular a `SparseMatrixCSC`, which `cholesky!`/`lu!` do not
support in place -- goes through the ordinary, allocating `cholesky`/`lu`.

# Arguments
- `Kcc::AbstractMatrix`: Stiffness of the condensed part
# Returns
- A factorization object supporting `ldiv!`
"""
function factorize_condensed!(Kcc::StridedMatrix{<:LinearAlgebra.BlasFloat})
    n = LinearAlgebra.checksquare(Kcc)
    d = diag(Kcc)                                   # nur die Diagonale sichern
    F = cholesky!(Symmetric(Kcc, :U); check = false)
    issuccess(F) && return F

    _warn_not_pd()
    # potrf hat nur das obere Dreieck inkl. Diagonale verändert; das untere ist intakt
    @inbounds for j in 1:n
        Kcc[j, j] = d[j]
        for i in 1:(j - 1)
            Kcc[i, j] = Kcc[j, i]
        end
    end
    return lu!(Kcc)
end

function factorize_condensed!(Kcc::AbstractMatrix)
    try
        return cholesky(Symmetric(Kcc))
    catch err
        _warn_not_pd(err)
        return lu(Kcc)
    end
end

"""
    ShiftInvertOperator

Applies `D Kcc^-1 D` with `D = sqrt(Mcc)`, as a matrix free operator for Arpack.

The generalised problem `Kcc X = Mcc X W` becomes the standard symmetric one
`A u = w u` with `A = D^-1 Kcc D^-1` and `u = D x`. Its smallest eigenvalues are the
largest of `A^-1 = D Kcc^-1 D`, which is what this operator applies. Two things follow:
shift-invert without asking Arpack for it, so an existing factorization is reused
instead of a second one being built, and eigenvectors orthonormal in the standard sense,
so `x = D^-1 u` is mass normalised without a further pass.

# Fields
- `factorization`: Factorization of `Kcc`
- `scale::Vector{Float64}`: Diagonal of `D`, the square root of the lumped mass
- `buffer::Vector{Float64}`: Work space
"""
struct ShiftInvertOperator{F}
    factorization::F
    scale::Vector{Float64}
    buffer::Vector{Float64}
end

Base.size(op::ShiftInvertOperator) = (length(op.scale), length(op.scale))
Base.size(op::ShiftInvertOperator, ::Integer) = length(op.scale)
Base.eltype(::ShiftInvertOperator) = Float64
LinearAlgebra.issymmetric(::ShiftInvertOperator) = true
LinearAlgebra.ishermitian(::ShiftInvertOperator) = true

function LinearAlgebra.mul!(y::AbstractVector{Float64}, op::ShiftInvertOperator,
                            x::AbstractVector{Float64})
    @inbounds @simd for i in eachindex(x)
        op.buffer[i] = op.scale[i] * x[i]
    end
    ldiv!(y, op.factorization, op.buffer)
    @inbounds @simd for i in eachindex(y)
        y[i] *= op.scale[i]
    end
    return y
end

Base.:*(op::ShiftInvertOperator, x::AbstractVector{Float64}) = mul!(similar(x), op, x)

"""
    GeneralizedShiftInvertOperator

Applies `Kcc^-1 Mcc`, as a matrix free operator for Arpack, for a `Mcc` that need not be
diagonal.

Unlike [`ShiftInvertOperator`](@ref), there is no `D = sqrt(Mcc)` scaling that turns the
generalised problem `Kcc X = Mcc X W` into a standard symmetric one -- that scaling only
works for a diagonal `Mcc`, true only for a cascade's very first level (see
[`reduce_matrices`](@ref)). The generalised problem still reduces to shift-invert: its
smallest eigenvalues `w` are the largest of `A = Kcc^-1 Mcc`, `A x = (1/w) x`, reusing
whatever factorization of `Kcc` the caller already holds instead of a second one being
built. `A` is not symmetric even though `Kcc` and `Mcc` both are, so Arpack runs its
general (Arnoldi) iteration here, not the symmetric Lanczos one -- slower per iteration
than [`ShiftInvertOperator`](@ref), and the eigenvectors come out `Mcc`-orthogonal only up
to Arnoldi's own tolerance, not exactly, so mass normalisation is a separate pass
afterwards rather than automatic.

# Fields
- `factorization`: Factorization of `Kcc`
- `Mcc::AbstractMatrix`: Mass of the condensed part, possibly non-diagonal
- `buffer::Vector{Float64}`: Work space, `Mcc * x` before the solve
"""
struct GeneralizedShiftInvertOperator{F,M<:AbstractMatrix}
    factorization::F
    Mcc::M
    buffer::Vector{Float64}
end

Base.size(op::GeneralizedShiftInvertOperator) = (length(op.buffer), length(op.buffer))
Base.size(op::GeneralizedShiftInvertOperator, ::Integer) = length(op.buffer)
Base.eltype(::GeneralizedShiftInvertOperator) = Float64
LinearAlgebra.issymmetric(::GeneralizedShiftInvertOperator) = false
LinearAlgebra.ishermitian(::GeneralizedShiftInvertOperator) = false

function LinearAlgebra.mul!(y::AbstractVector{Float64}, op::GeneralizedShiftInvertOperator,
                            x::AbstractVector{Float64})
    mul!(op.buffer, op.Mcc, x)
    ldiv!(y, op.factorization, op.buffer)
    return y
end

function Base.:*(op::GeneralizedShiftInvertOperator, x::AbstractVector{Float64})
    mul!(similar(x),
         op, x)
end

"""
    fixed_interface_modes_general(factorization, Kcc, Mcc, n_modes; tol, maxiter)

The `n_modes` lowest eigenpairs of `Kcc X = Mcc X W`, mass normalised, for a general (not
necessarily diagonal) `Mcc`, reusing a factorization of `Kcc` the caller already computed.

Below a few hundred degrees of freedom, or for a dense `Kcc`, the plain dense
`eigen(Symmetric(Kcc), Symmetric(Mcc))` from [`fixed_interface_modes_dense`](@ref) is
cheaper and simpler and is used instead. Above that,
[`GeneralizedShiftInvertOperator`](@ref) avoids ever densifying `Kcc`/`Mcc` or solving for
more than `n_modes` eigenpairs -- what made a cascade level with an outlier-sized shell
(see `dense_shell_limit` in [`reduce_matrices`](@ref)) and `n_modes > 0` cost `O(nc^3)`
time and `O(nc^2)` memory before, dense or not.

# Arguments
- `factorization`: Factorization of `Kcc`, as already computed by
  [`factorize_condensed!`](@ref)
- `Kcc::AbstractMatrix`: Stiffness of the condensed part
- `Mcc::AbstractMatrix`: Mass of the condensed part, possibly non-diagonal
- `n_modes::Int64`: Number of modes
# Keywords
- `tol::Float64`: Relative tolerance handed to Arpack
- `maxiter::Int64`: Iteration limit per attempt
# Returns
- `X::Matrix{Float64}`: Modes, `nc x n_modes`, mass normalised (`X' Mcc X = I`)
- `w::Vector{Float64}`: Eigenvalues, ascending
"""
function fixed_interface_modes_general(factorization, Kcc::AbstractMatrix,
                                       Mcc::AbstractMatrix, n_modes::Int64;
                                       tol::Float64 = 1e-10, maxiter::Int64 = 1000)
    nc = size(Kcc, 1)
    n_modes == 0 && return zeros(nc, 0), Float64[]
    if !issparse(Kcc) || nc <= 500
        return fixed_interface_modes_dense(Kcc, Mcc, n_modes)
    end

    operator = GeneralizedShiftInvertOperator(factorization, Mcc,
                                              Vector{Float64}(undef, nc))
    last_error = nothing
    for ncv in (min(nc, max(20, 4 * n_modes + 1)),
        min(nc, max(60, 10 * n_modes + 1)),
        min(nc, max(150, 20 * n_modes + 1)))
        try
            theta,
            U = eigs(operator; nev = n_modes, which = :LM, ncv = ncv,
                     tol = tol, maxiter = maxiter)
            w = 1.0 ./ real.(theta)
            order = sortperm(w)
            X = real.(U)[:, order]
            # Arnoldi's Mcc-orthogonality only holds up to its own tolerance, not
            # exactly like the symmetric Lanczos path in fixed_interface_modes -- each
            # mode is normalised on its own rather than jointly orthogonalised, which is
            # enough since distinct eigenvalues already make X Mcc-orthogonal in exact
            # arithmetic.
            for k in 1:n_modes
                @views nrm = sqrt(dot(X[:, k], Mcc * X[:, k]))
                X[:, k] ./= nrm
            end
            return X, w[order]
        catch err
            last_error = err
            @debug "Arpack did not converge, retrying with a larger subspace" ncv
        end
    end

    @warn "Arpack did not converge for $n_modes modes." exception=last_error
    throw(ErrorException("Could not compute $n_modes fixed-interface modes."))
end

"""
    fixed_interface_modes(Kcc, Mcc, n_modes; tol, maxiter)

The `n_modes` lowest eigenpairs of `Kcc X = Mcc X W`, mass normalised.

Rows of `Kcc` without any entry are excluded from the eigenvalue problem and their
entries in the modes stay zero, so that an isolated condensed degree of freedom does not
make the problem singular.

The eigenpairs come from shift-invert Lanczos through [`ShiftInvertOperator`](@ref), so
`Kcc` is never densified. A dense `eigen` would cost `O(nc^2)` memory and `O(nc^3)` time
for all `nc` eigenpairs, of which `n_modes` are wanted. Below a few hundred degrees of
freedom the dense route is the cheaper one and is taken instead.

The eigenvectors are mass normalised, `X' Mcc X = I`. Without that the modal blocks of
the reduced matrices have no defined scaling.

Both routes treat `Kcc` as symmetric. For an unsymmetric `Kcc` the modes are those of
its symmetric part, which is an approximation — the coupling blocks formed in
[`reduce_matrices`](@ref) no longer cancel then, which is why they are formed
explicitly.

# Arguments
- `Kcc::AbstractMatrix`: Stiffness of the condensed part
- `Mcc::Diagonal`: Lumped mass of the condensed part
- `n_modes::Int64`: Number of modes
# Keywords
- `tol::Float64`: Relative tolerance handed to Arpack
- `maxiter::Int64`: Iteration limit per attempt
# Returns
- `X::Matrix{Float64}`: Modes, `nc x n_modes`
- `w::Vector{Float64}`: Eigenvalues, ascending
"""
function fixed_interface_modes(Kcc::AbstractMatrix, Mcc::Diagonal, n_modes::Int64;
                               tol::Float64 = 1e-10, maxiter::Int64 = 1000)
    nc = size(Kcc, 1)
    X = zeros(nc, n_modes)
    n_modes == 0 && return X, Float64[]

    active = nonzero_rows(Kcc)
    if length(active) < nc
        @warn "$(nc - length(active)) of $nc condensed degrees of freedom carry no " *
              "stiffness and are excluded from the eigenvalue problem. Check whether " *
              "every condensed node is still tied to the retained region."
    end
    if n_modes > length(active)
        throw(ArgumentError("n_modes = $n_modes, but only $(length(active)) condensed " *
                            "degrees of freedom carry stiffness."))
    end

    K_active = Kcc[active, active]
    mass_active = diag(Mcc)[active]
    scale = sqrt.(mass_active)

    if !issparse(Kcc) || length(active) <= 500
        # The scaling by sqrt(M) turns the generalised problem into a standard symmetric
        # one, whose eigenvectors come out orthonormal — which is the mass
        # normalisation, for free.
        inverse_scale = Diagonal(inv.(scale))
        factorization = eigen(Symmetric(inverse_scale * Matrix(K_active) *
                                        inverse_scale))
        X[active, :] = inverse_scale * factorization.vectors[:, 1:n_modes]
        return X, factorization.values[1:n_modes]
    end

    operator = ShiftInvertOperator(factorize_condensed!(K_active), scale,
                                   Vector{Float64}(undef, length(active)))

    last_error = nothing
    for ncv in (min(length(active), max(20, 4 * n_modes + 1)),
        min(length(active), max(60, 10 * n_modes + 1)),
        min(length(active), max(150, 20 * n_modes + 1)))
        try
            theta,
            U = eigs(operator; nev = n_modes, which = :LM, ncv = ncv,
                     tol = tol, maxiter = maxiter)
            w = 1.0 ./ real.(theta)
            order = sortperm(w)
            V = real.(U)[:, order]
            @inbounds for k in 1:n_modes, i in eachindex(active)
                V[i, k] /= scale[i]
            end
            X[active, :] = V
            return X, w[order]
        catch err
            last_error = err
            @debug "Arpack did not converge, retrying with a larger subspace" ncv
        end
    end

    @warn "Arpack did not converge for $n_modes modes." exception=last_error
    throw(ErrorException("Could not compute $n_modes fixed-interface modes."))
end

"""
    recovery_modes(factorization, Kcr, coupling; block_size = 512)

Recovery modes `B = Kcc^-1 Kcr` for the coupling degrees of freedom only.

Columns of `Kcr` outside `coupling` are empty, so their recovery modes are zero and are
neither solved for nor stored. The result is `nc x length(coupling)` instead of
`nc x nr`.

Solved in column blocks so that the dense right hand side buffer stays at
`nc x block_size` rather than being a second matrix of the full size.

# Arguments
- `factorization`: Factorization of `Kcc`
- `Kcr::AbstractMatrix`: Coupling block, condensed to retained
- `coupling::Vector{Int64}`: Columns to solve for
# Keywords
- `block_size::Int64`: Right hand sides solved at once
# Returns
- `::Matrix{Float64}`: The recovery modes, `nc x length(coupling)`
"""
function recovery_modes(factorization, Kcr::AbstractMatrix, coupling::Vector{Int64};
                        block_size::Int64 = 512)
    nc = size(Kcr, 1)
    nrc = length(coupling)
    B = Matrix{Float64}(undef, nc, nrc)
    nrc == 0 && return B

    width = min(block_size, nrc)
    rhs = Matrix{Float64}(undef, nc, width)

    for first in 1:width:nrc
        last = min(first + width - 1, nrc)
        local_columns = first:last
        buffer = view(rhs, :, 1:length(local_columns))
        fill!(buffer, 0.0)
        copy_sparse_columns!(buffer, Kcr, view(coupling, local_columns))
        ldiv!(view(B, :, local_columns), factorization, buffer)
    end

    return B
end

"""
    copy_sparse_columns!(target, A, columns)

Copies selected columns of a sparse matrix into a dense buffer, without the temporary
`A[:, columns]` would allocate.
"""
function copy_sparse_columns!(target::AbstractMatrix{Float64},
                              A::SparseMatrixCSC{Float64,<:Integer},
                              columns::AbstractVector{Int64})
    rows = rowvals(A)
    values = nonzeros(A)
    @inbounds for (local_column, column) in enumerate(columns)
        for index in nzrange(A, column)
            target[rows[index], local_column] = values[index]
        end
    end
    return target
end

function copy_sparse_columns!(target::AbstractMatrix{Float64}, A::AbstractMatrix,
                              columns::AbstractVector{Int64})
    copyto!(target, A[:, columns])
end

"""
    check_symmetry(K; samples = 2000, tolerance = 1e-8)

Measures how far the stiffness is from symmetric, sampling stored entries.

Craig-Bampton is built on symmetry: the fixed-interface modes use the symmetric solver,
and the classical derivation drops the stiffness coupling blocks because `Kcc B = Kcr`
makes them cancel — which needs `Krc = Kcr'`. In peridynamics the stiffness is not
symmetric in general, two points with different horizons do not contribute equally to
each other, so the deviation is measured rather than assumed away. The coupling blocks
are formed explicitly either way; this only reports how much they will carry.

Sampling runs over the stored entries of `K`, not over random index pairs. For a sparse
matrix a random pair is almost surely a structural zero on both sides, which would make
any deviation invisible and report perfect symmetry for every input.

# Arguments
- `K::AbstractMatrix`: The full stiffness matrix
# Keywords
- `samples::Int64`: Number of stored entries drawn
- `tolerance::Float64`: Relative deviation accepted as rounding
# Returns
- `::Float64`: The largest relative deviation found
"""
function check_symmetry(K::SparseMatrixCSC; samples::Int64 = 2000,
                        tolerance::Float64 = 1e-8)
    rows = rowvals(K)
    values = nonzeros(K)
    total = length(values)
    total == 0 && return 0.0

    stride = max(1, total ÷ samples)
    scale = 0.0
    deviation = 0.0
    drawn = 0

    @inbounds for column in axes(K, 2)
        for index in nzrange(K, column)
            index % stride == 0 || continue
            i = rows[index]
            upper = values[index]
            lower = K[column, i]
            scale = max(scale, abs(upper), abs(lower))
            deviation = max(deviation, abs(upper - lower))
            drawn += 1
        end
    end

    drawn == 0 && return 0.0
    relative = deviation / max(scale, eps())

    if relative > tolerance
        @warn "Stiffness matrix is not symmetric (largest relative deviation " *
              "$(round(relative; sigdigits = 3)) over $drawn sampled entries). The " *
              "modal coupling blocks are kept, so the reduction stays consistent, but " *
              "the fixed-interface modes are those of the symmetric part."
    end
    return relative
end

check_symmetry(K::AbstractMatrix; kwargs...) = check_symmetry(sparse(K); kwargs...)

"""
    add_dense_block!(rows, columns, values, block, positions)

Appends a dense square block to a coordinate list, mapped to the given global positions.

Entries below the rounding level of the block are skipped: they are fill-in that carries
no information and would only make the result denser.
"""
function add_dense_block!(rows::Vector{Int64}, columns::Vector{Int64},
                          values::Vector{Float64}, block::AbstractMatrix{Float64},
                          positions::AbstractVector{Int64})
    threshold = 1.0e-12 * max(maximum(abs, block; init = 0.0), eps())
    @inbounds for j in axes(block, 2), i in axes(block, 1)
        value = block[i, j]
        abs(value) <= threshold && continue
        push!(rows, positions[i])
        push!(columns, positions[j])
        push!(values, value)
    end
    return nothing
end

"""
    add_dense_block!(rows, columns, values, block, row_positions, column_positions)

Appends a rectangular block, rows and columns mapped independently. Used for the
coupling blocks between the retained and the modal part, which are `nc x n_modes`.
"""
function add_dense_block!(rows::Vector{Int64}, columns::Vector{Int64},
                          values::Vector{Float64}, block::AbstractMatrix{Float64},
                          row_positions::AbstractVector{Int64},
                          column_positions::AbstractVector{Int64};
                          threshold::Float64 = -1.0)
    if threshold < 0
        threshold = 1.0e-12 * max(maximum(abs, block; init = 0.0), eps())
    end
    @inbounds for j in axes(block, 2), i in axes(block, 1)
        value = block[i, j]
        abs(value) <= threshold && continue
        push!(rows, row_positions[i])
        push!(columns, column_positions[j])
        push!(values, value)
    end
    return nothing
end

"""
    report_modal_coupling(K_reduced, M_reduced, nr, n_modes)

Reports the coupling between the retained and the modal part, and warns if it is empty.

A mode that no block couples to stays at rest for the whole simulation. It still carries
mass and shifts the eigenfrequencies, but never responds, which makes the result worse
than the static condensation it was meant to improve on — and, characteristically,
independent of the number of modes.

There are two paths. `Mrcm` drives the modes through the acceleration of the retained
degrees of freedom and is always present. `Krcm` drives them through the displacement and
vanishes exactly for a symmetric stiffness, so an empty one is only expected there.

# Arguments
- `K_reduced::SparseMatrixCSC`, `M_reduced::SparseMatrixCSC`: The reduced matrices
- `nr::Int64`: Number of retained degrees of freedom
- `n_modes::Int64`: Number of modes
# Returns
- `::Bool`: Whether any coupling holds entries
"""
function report_modal_coupling(K_reduced::SparseMatrixCSC, M_reduced::SparseMatrixCSC,
                               nr::Int64, n_modes::Int64)
    n_modes == 0 && return true

    modal = (nr + 1):(nr + n_modes)
    mass_coupling = M_reduced[1:nr, modal]
    stiffness_coupling = K_reduced[1:nr, modal]

    mass_entries = nnz(mass_coupling)
    stiffness_entries = nnz(stiffness_coupling)

    if mass_entries == 0 && stiffness_entries == 0
        @warn "Neither the mass nor the stiffness couples the retained and the modal " *
              "part. The modes will stay at rest and the reduction behaves worse than " *
              "a static condensation."
        return false
    end

    retained_mass = maximum(abs, nonzeros(M_reduced[1:nr, 1:nr]); init = 0.0)
    retained_stiffness = maximum(abs, nonzeros(K_reduced[1:nr, 1:nr]); init = 0.0)

    @info "Modal coupling: mass $mass_entries entries, largest " *
          "$(round(maximum(abs, nonzeros(mass_coupling); init = 0.0); sigdigits = 3)) " *
          "against $(round(retained_mass; sigdigits = 3)); stiffness " *
          "$stiffness_entries entries, largest " *
          "$(round(maximum(abs, nonzeros(stiffness_coupling); init = 0.0); sigdigits = 3)) " *
          "against $(round(retained_stiffness; sigdigits = 3))"

    return true
end

"""
    reduce_stage(K, M_diag, r, c, n_modes = 1; block_size = 512,
                check_symmetry_sample = 2000)

Craig-Bampton reduction of a stiffness and a lumped mass matrix -- one stage.

The single-region reduction, unchanged from `CraigBampton.reduce_matrices`; the cascade's
`reduce_matrices` below calls this once per stage on a local slice of `K`, carrying the
previous stage's reduced boundary forward into the next slice instead of ever factorizing
the full condensed region at once.

The degrees of freedom are split into the retained ones `r` and the condensed ones `c`.
The retained ones stay physical coordinates, the condensed ones are represented by the
`n_modes` lowest fixed-interface modes `X` of `Kcc X = Mcc X W`:

    | u_r |   | I   0 | | u_r |
    | u_c | = | -B  X | | eta |

with the recovery modes `B = Kcc^-1 Kcr`, the same static relation Guyan uses. The
reduced blocks are

    Krcrc = Krr - Krc B         Kmm = X' Kcc X
    Krcm = Krc X - B' Kcc X     Kmrc = X' Kcr - X' Kcc B
    Mrcrc = Mrr + B' Mcc B      Mmm = I
    Mrcm = -B' Mcc X

The textbook derivation drops `Krcm` and `Kmrc`: with `Kcc B = Kcr` the two contributions
to `Krcm` become `Krc X - Kcr' X`, which is zero for `Krc = Kcr'`. That holds for a
symmetric stiffness only. Peridynamic stiffness matrices are not symmetric in general —
two points with different horizons do not contribute equally to each other — and
dropping the blocks then removes the path through which a displacement of the retained
region drives the modes. What remains is the mass coupling alone, the modes stay
underexcited, and adding modes no longer changes the answer. The blocks are therefore
formed explicitly; for a symmetric `K` they come out at rounding level and are dropped
by the threshold, leaving the classical result unchanged.

`Kmm` is formed as `X' Kcc X` for the same reason, rather than being set to the
eigenvalues, which are the right answer only for a symmetric `Kcc`.

The mass coupling terms `Mrc X` vanish for a genuine reason and are not formed: the mass
matrix is diagonal and the two index sets are disjoint.

The reduced matrices stay sparse. Only retained degrees of freedom with a bond into the
condensed region appear in `Kcr`; the reduction fills in exactly their rows and columns,
and everything else keeps the sparsity of `Krr`. Both the solves for `B` and the dense
blocks therefore scale with the size of the coupling layer, not with the number of
retained degrees of freedom.

`n_modes = 0` drops the modal part and leaves Guyan condensation.

# Arguments
- `K::AbstractMatrix`: Stiffness matrix
- `M_diag::AbstractVector`: Lumped mass matrix as a vector, one entry per degree of freedom.
  `K*u` delivers force densities (see `Guyan.reduce_matrices`), so this is a mass
  density, not a physical mass; it is used as-is, with no volume weighting.
- `r::AbstractVector{<:Integer}`: Indices of the retained degrees of freedom
- `c::AbstractVector{<:Integer}`: Indices of the condensed degrees of freedom
- `n_modes::Integer`: Number of fixed-interface modes to keep
# Keywords
- `block_size::Int64`: Right hand sides solved at once for the recovery modes
- `check_symmetry_sample::Int64`: Entries drawn for the symmetry check, 0 disables it
# Returns
- `K_reduced::SparseMatrixCSC`: Reduced stiffness, size `length(r) + n_modes`
- `M_reduced::SparseMatrixCSC`: Reduced mass, same size
"""
function reduce_stage(K::AbstractMatrix,
                      M_diag::AbstractVector,
                      r::AbstractVector{<:Integer},
                      c::AbstractVector{<:Integer},
                      n_modes::Integer = 1;
                      block_size::Int64 = 512,
                      check_symmetry_sample::Int64 = 2000)
    r = collect(Int64, r)
    c = collect(Int64, c)
    n_modes = Int64(n_modes)

    nr = length(r)
    nc = length(c)

    if !isempty(intersect(r, c))
        throw(ArgumentError("Retained and condensed index sets overlap."))
    end
    if n_modes < 0 || n_modes > nc
        throw(ArgumentError("n_modes = $n_modes, but only $nc condensed degrees of " *
                            "freedom are available."))
    end

    if check_symmetry_sample > 0
        @timeit "CB symmetry check" check_symmetry(K; samples = check_symmetry_sample)
    end

    @timeit "CB extract submatrices" begin
        Krr = K[r, r]
        Krc = K[r, c]
        Kcr = K[c, r]
        Kcc = K[c, c]
        # M_diag is a mass density, matching the force densities K*u delivers; see
        # Guyan.reduce_matrices for why no volume weighting belongs here.
        mass_c = M_diag[c]
        mass_r = M_diag[r]
        Mcc = Diagonal(mass_c)
    end

    # A retained dof with a nonzero Krc row but (for an unsymmetric K) a zero Kcr
    # column would otherwise be silently dropped from the coupling correction; see
    # Guyan.reduce_matrices for the same fix on the same asymmetry.
    coupling = union(nonzero_rows(Krc), nonzero_columns(Kcr))
    nrc = length(coupling)

    @info "Craig-Bampton: $nr retained, $nc condensed, $nrc coupling, $n_modes modes; " *
          "reduced size $(nr + n_modes), dense blocks of " *
          "$(round(nrc^2 * 8 / 2^20; digits = 1)) MiB"

    @timeit "CB factorize Kcc" factorization=factorize_condensed!(Kcc)

    @timeit "CB recovery modes" B=recovery_modes(factorization, Kcr, coupling;
                                                 block_size = block_size)

    # The modes are normalised against Mcc, so X' Mcc X = I.
    @timeit "CB fixed interface modes" X, w=fixed_interface_modes(Kcc, Mcc, n_modes)

    if n_modes > 0
        f_lo = sqrt(w[1]) / (2 * pi)
        f_hi = sqrt(w[end]) / (2 * pi)
        @info "Craig-Bampton fixed-interface frequencies: $(round(f_lo; sigdigits = 4)) Hz " *
              "(mode 1) to $(round(f_hi; sigdigits = 4)) Hz (mode $n_modes)"
    end

    n_total = nr + n_modes
    modal_positions = collect((nr + 1):(nr + n_modes))

    @timeit "CB reduced stiffness" begin
        rows, columns, values = findnz(sparse(Krr))

        if nrc > 0
            # Only the coupling rows of Krc carry entries, so the product is nrc x nrc.
            Krcrc_fill = Matrix{Float64}(undef, nrc, nrc)
            mul!(Krcrc_fill, Krc[coupling, :], B, -1.0, 0.0)
            add_dense_block!(rows, columns, values, Krcrc_fill, coupling)
        end

        if n_modes > 0
            Kcc_X = Kcc * X                                  # nc x n_modes

            if nrc > 0
                # Krcm = Krc X - B' Kcc X
                Krcm = Matrix{Float64}(undef, nrc, n_modes)
                mul!(Krcm, Krc[coupling, :], X)
                mul!(Krcm, transpose(B), Kcc_X, -1.0, 1.0)

                # For an unsymmetric K the transposed block is not Krcm', so it is
                # formed separately: Kmrc = X' Kcr - (Kcc' X)' B
                Kmrc = Matrix{Float64}(undef, n_modes, nrc)
                mul!(Kmrc, transpose(X), Kcr[:, coupling])
                mul!(Kmrc, transpose(transpose(Kcc) * X), B, -1.0, 1.0)

                threshold = 1.0e-12 *
                            max(maximum(abs, Krcm; init = 0.0),
                                maximum(abs, Kmrc; init = 0.0), eps())
                add_dense_block!(rows, columns, values, Krcm, coupling, modal_positions;
                                 threshold = threshold)
                add_dense_block!(rows, columns, values, Kmrc, modal_positions, coupling;
                                 threshold = threshold)

                reference = maximum(abs, nonzeros(Krr); init = eps())
                @info "Craig-Bampton modal stiffness coupling: largest |Krcm| = " *
                      "$(round(maximum(abs, Krcm; init = 0.0); sigdigits = 3)), " *
                      "|Kmrc| = $(round(maximum(abs, Kmrc; init = 0.0); sigdigits = 3)), " *
                      "against $(round(reference; sigdigits = 3)) in Krr"
            end

            # Kmm = X' Kcc X. Equal to diagm(w) for a symmetric Kcc; formed explicitly
            # so that an unsymmetric one does not silently fall back to the eigenvalues.
            Kmm = transpose(X) * Kcc_X
            add_dense_block!(rows, columns, values, Kmm, modal_positions,
                             modal_positions)
        end

        K_reduced = sparse(rows, columns, values, n_total, n_total)
    end

    @timeit "CB reduced mass" begin
        # Both mass blocks have the form Y' Mcc Z with a diagonal Mcc, so scaling the
        # factors by sqrt(Mcc) turns them into plain products. B and X are scaled in
        # place, which is why the stiffness block above had to be formed first.
        root_mass = sqrt.(mass_c)
        @inbounds for j in 1:nrc, i in 1:nc
            B[i, j] *= root_mass[i]
        end

        rows = Int64[]
        columns = Int64[]
        values = Float64[]

        if nrc > 0
            Mrcrc_fill = Matrix{Float64}(undef, nrc, nrc)
            mul!(Mrcrc_fill, transpose(B), B)
            add_dense_block!(rows, columns, values, Mrcrc_fill, coupling)
        end

        for local_index in eachindex(r)
            push!(rows, local_index)
            push!(columns, local_index)
            push!(values, mass_r[local_index])
        end

        if n_modes > 0 && nrc > 0
            @inbounds for k in 1:n_modes, i in 1:nc
                X[i, k] *= root_mass[i]
            end
            # Mrcm = -B' Mcc X, the minus coming from u_c = -B u_r + X eta.
            Mrcm = Matrix{Float64}(undef, nrc, n_modes)
            mul!(Mrcm, transpose(B), X, -1.0, 0.0)

            threshold = 1.0e-12 * max(maximum(abs, Mrcm; init = 0.0), eps())
            add_dense_block!(rows, columns, values, Mrcm, coupling, modal_positions;
                             threshold = threshold)
            add_dense_block!(rows, columns, values, transpose(Mrcm), modal_positions,
                             coupling; threshold = threshold)
        end

        for k in 1:n_modes
            push!(rows, nr + k)
            push!(columns, nr + k)
            push!(values, 1.0)
        end

        M_reduced = sparse(rows, columns, values, n_total, n_total)
    end

    report_modal_coupling(K_reduced, M_reduced, nr, n_modes)

    @debug "Reduced matrices: K with $(nnz(K_reduced)) entries, M with " *
           "$(nnz(M_reduced)) entries, out of $(n_total^2) possible"

    return K_reduced, M_reduced
end

"""
    fixed_interface_modes_dense(Kcc, Mcc, n_modes)

The `n_modes` lowest eigenpairs of `Kcc X = Mcc X W`, mass normalised, for a general
(not necessarily diagonal) symmetric positive definite `Mcc`.

Used instead of [`fixed_interface_modes`](@ref) for a cascade level's own local modes:
`Mcc` there is a cascade level's own shell plus whatever it absorbs from the previous
level, small by construction (bounded by the shell width transverse to `r`, not by the
condensed region's total size -- see [`reduce_matrices`](@ref)), so the dense
generalised eigenproblem -- `LinearAlgebra.eigen` on the
matrix pencil, `B`-orthonormal eigenvectors by construction (`X' Mcc X = I`) -- is both
simpler and cheap enough, unlike [`fixed_interface_modes`](@ref)'s sparse Arpack path
built for a single, potentially large, whole condensed region. `Mcc` needs the general
form here specifically because a cascade stage's own group can itself contain a degree
of freedom carried forward from an earlier stage (see [`reduce_matrices`](@ref)), whose
mass there is a full `Mrcrc` block, not part of the original lumped diagonal.

# Arguments
- `Kcc::AbstractMatrix`: Stiffness of the condensed part
- `Mcc::AbstractMatrix`: Mass of the condensed part, possibly non-diagonal
- `n_modes::Int64`: Number of modes
# Returns
- `X::Matrix{Float64}`: Modes, `nc x n_modes`
- `w::Vector{Float64}`: Eigenvalues, ascending
"""
function fixed_interface_modes_dense(Kcc::AbstractMatrix, Mcc::AbstractMatrix,
                                     n_modes::Int64)
    nc = size(Kcc, 1)
    n_modes == 0 && return zeros(nc, 0), Float64[]
    if n_modes > nc
        throw(ArgumentError("n_modes = $n_modes, but this stage's group only has $nc " *
                            "degrees of freedom."))
    end
    factorization = eigen(Symmetric(Matrix(Kcc)), Symmetric(Matrix(Mcc)))
    return factorization.vectors[:, 1:n_modes], factorization.values[1:n_modes]
end

"""
    reduce_stage_dense(K, Mcc_source, r, c, n_modes = 0; block_size = 512)

Condensation of one cascade stage where the condensed block's mass, `Mcc`, is not
necessarily diagonal, optionally retaining that stage's own fixed-interface normal
modes.

Used instead of [`reduce_stage`](@ref) exactly when a cascade stage's own group contains
a degree of freedom that was carried forward, retained, from an earlier stage: its mass
there is `Mrcrc = Mrr + B' Mcc B` from that earlier stage, generally full, not the lumped
diagonal `reduce_stage` assumes. Retaining modes here is what makes a cascade a
Craig-Bampton one and not a purely static (Guyan) one: since a stage's own condensed
group is discarded once eliminated, its internal dynamics can only be captured *during*
that elimination, as this stage's own local fixed-interface modes -- there is no
retrieving them afterwards from the boundary-only matrices further stages see. Those
modes are carried forward exactly like a retained physical degree of freedom, as extra
rows/columns of the boundary block handed to the next stage, and are never eliminated
again.

The same congruence transform as [`reduce_stage`](@ref) applies, generalised for a
possibly non-diagonal `Mrc`/`Mcr`: `Mrcm = Mrc X - B' Mcc X` (`reduce_stage`'s `-B' Mcc X`
alone assumes `Mrc = 0`, true only for the original lumped mass, not for a carried-in
block). `M` stays symmetric throughout the cascade, so `Mmrc = Mrcm'` is used directly
rather than formed separately, unlike the stiffness coupling blocks, which need both
`Krcm` and `Kmrc` explicitly since `K` need not be symmetric.

# Arguments
- `K::AbstractMatrix`: Stiffness matrix, this stage's local slice
- `Mcc_source::AbstractMatrix{Float64}`: Mass matrix, this stage's local slice, `c` block
  possibly non-diagonal
- `r::AbstractVector{<:Integer}`: Indices of the retained degrees of freedom
- `c::AbstractVector{<:Integer}`: Indices of the condensed degrees of freedom
- `n_modes::Integer`: Number of this stage's own fixed-interface modes to retain
# Keywords
- `block_size::Int64`: Right hand sides solved at once for the recovery modes
# Returns
- `K_reduced::SparseMatrixCSC`, `M_reduced::SparseMatrixCSC`: Reduced stiffness and mass,
  size `length(r) + n_modes`
"""
function reduce_stage_dense(K::AbstractMatrix, Mcc_source::AbstractMatrix{Float64},
                            r::AbstractVector{<:Integer}, c::AbstractVector{<:Integer},
                            n_modes::Integer = 0; block_size::Int64 = 512)
    r = collect(Int64, r)
    c = collect(Int64, c)
    n_modes = Int64(n_modes)
    nr = length(r)
    nc = length(c)

    Krr = K[r, r]
    Krc = K[r, c]
    Kcr = K[c, r]
    Kcc = K[c, c]
    Mrr = Mcc_source[r, r]
    Mrc = Mcc_source[r, c]
    Mcr = Mcc_source[c, r]
    Mcc = Mcc_source[c, c]

    coupling = union(nonzero_rows(Krc), nonzero_columns(Kcr))
    nrc = length(coupling)
    # Printed before any of the expensive dense work below (the Cholesky of Kcc, the
    # eigenvalue solve for modes, the Krcrc/Mrcrc fill-in) so the sizes that actually
    # drove a crash mid-stage are on record even if nothing after this line ever gets to
    # run.
    @info "Craig-Bampton Cascade stage: nr=$nr (this stage's boundary), nc=$nc " *
          "(eliminated now), nrc=$nrc (of nr actually coupled to nc)"

    # Kcc is read again below (Kcc_X, fixed_interface_modes_dense), so factorizing it
    # needs a defensive copy only on the dense path: factorize_condensed! mutates a
    # dense StridedMatrix in place (cholesky!/lu!), but its sparse fallback
    # (cholesky/lu, no `!`) never touches Kcc itself -- forcing Matrix(Kcc)
    # unconditionally would both discard a genuinely sparse Kcc's sparsity (see
    # dense_shell_limit in reduce_matrices) and pay for a copy that a sparse Kcc does
    # not need in the first place.
    factorization = Kcc isa AbstractSparseMatrix ? factorize_condensed!(Kcc) :
                    factorize_condensed!(Matrix(Kcc))
    B = recovery_modes(factorization, Kcr, coupling; block_size = block_size)

    # Reuses `factorization` (already built above for B) rather than densifying Kcc a
    # second time here -- see fixed_interface_modes_general for why that matters once
    # nc is large (an outlier-sized shell, see dense_shell_limit in reduce_matrices).
    X, w = fixed_interface_modes_general(factorization, Kcc, Mcc, n_modes)
    if n_modes > 0
        f_lo = sqrt(abs(w[1])) / (2 * pi)
        f_hi = sqrt(abs(w[end])) / (2 * pi)
        @info "Craig-Bampton Cascade stage: fixed-interface frequencies " *
              "$(round(f_lo; sigdigits = 4)) Hz (mode 1) to $(round(f_hi; sigdigits = 4)) " *
              "Hz (mode $n_modes)"
    end

    n_total = nr + n_modes
    modal_positions = collect((nr + 1):(nr + n_modes))

    rows, columns, values = findnz(sparse(Krr))
    if nrc > 0
        Krcrc_fill = Matrix{Float64}(undef, nrc, nrc)
        mul!(Krcrc_fill, Krc[coupling, :], B, -1.0, 0.0)
        add_dense_block!(rows, columns, values, Krcrc_fill, coupling)
    end
    if n_modes > 0
        Kcc_X = Kcc * X
        if nrc > 0
            Krcm = Matrix{Float64}(undef, nrc, n_modes)
            mul!(Krcm, Krc[coupling, :], X)
            mul!(Krcm, transpose(B), Kcc_X, -1.0, 1.0)

            Kmrc = Matrix{Float64}(undef, n_modes, nrc)
            mul!(Kmrc, transpose(X), Kcr[:, coupling])
            mul!(Kmrc, transpose(transpose(Kcc) * X), B, -1.0, 1.0)

            threshold = 1.0e-12 *
                        max(maximum(abs, Krcm; init = 0.0),
                            maximum(abs, Kmrc; init = 0.0), eps())
            add_dense_block!(rows, columns, values, Krcm, coupling, modal_positions;
                             threshold = threshold)
            add_dense_block!(rows, columns, values, Kmrc, modal_positions, coupling;
                             threshold = threshold)
        end
        Kmm = transpose(X) * Kcc_X
        add_dense_block!(rows, columns, values, Kmm, modal_positions, modal_positions)
    end
    K_reduced = sparse(rows, columns, values, n_total, n_total)

    mrows, mcolumns, mvalues = findnz(sparse(Mrr))
    if nrc > 0
        # Mcc_source can carry a non-diagonal Mrc/Mcr coupling here (unlike the raw
        # lumped mass): this block comes from a previous cascade stage's Mrcrc, which is
        # Mrr + B'MccB in general position, not block-diagonal. The full congruence
        # transform is therefore needed, not just the B'MccB term -- and, unlike Krc
        # and Kcr, Mrc and Mcr are *not* guaranteed to vanish outside `coupling`: that
        # exclusion is only valid for a quantity B itself annihilates there (Krc/Kcr,
        # zero outside `coupling` by definition of `coupling`), and B being zero outside
        # `coupling` says nothing about Mrc/Mcr's own support, which comes from whatever
        # a previous stage's Mrcrc carried, unrelated to this stage's own K sparsity. The
        # `Mrc B` / `B' Mcr` terms therefore use every retained row/column; only the
        # `B' Mcc B` term, where both factors are B, stays restricted to `coupling`.
        Mrcrc_rc = Matrix{Float64}(undef, nr, nrc)
        mul!(Mrcrc_rc, Mrc, B, -1.0, 0.0)
        add_dense_block!(mrows, mcolumns, mvalues, Mrcrc_rc, collect(1:nr), coupling)

        Mrcrc_cr = Matrix{Float64}(undef, nrc, nr)
        mul!(Mrcrc_cr, transpose(B), Mcr, -1.0, 0.0)
        add_dense_block!(mrows, mcolumns, mvalues, Mrcrc_cr, coupling, collect(1:nr))

        Mrcrc_cc = Matrix{Float64}(undef, nrc, nrc)
        mul!(Mrcrc_cc, transpose(B), Mcc * B)
        add_dense_block!(mrows, mcolumns, mvalues, Mrcrc_cc, coupling)
    end
    if n_modes > 0
        # General form: Mrcm = Mrc X - B' Mcc X. Same reasoning as Mrcrc above: Mrc X
        # uses every retained row, the B' Mcc X correction only the coupling ones.
        Mrcm = Matrix{Float64}(undef, nr, n_modes)
        mul!(Mrcm, Mrc, X)
        if nrc > 0
            Mrcm_corr = Matrix{Float64}(undef, nrc, n_modes)
            mul!(Mrcm_corr, transpose(B), Mcc * X, -1.0, 0.0)
            @views Mrcm[coupling, :] .+= Mrcm_corr
        end

        threshold = 1.0e-12 * max(maximum(abs, Mrcm; init = 0.0), eps())
        add_dense_block!(mrows, mcolumns, mvalues, Mrcm, collect(1:nr), modal_positions;
                         threshold = threshold)
        add_dense_block!(mrows, mcolumns, mvalues, transpose(Mrcm), modal_positions,
                         collect(1:nr); threshold = threshold)
    end
    if n_modes > 0
        for k in 1:n_modes
            push!(mrows, nr + k)
            push!(mcolumns, nr + k)
            push!(mvalues, 1.0)
        end
    end
    M_reduced = sparse(mrows, mcolumns, mvalues, n_total, n_total)

    dropzeros!(K_reduced)
    dropzeros!(M_reduced)
    return K_reduced, M_reduced
end

"""
    group_condensed_by_level(K, r, c)

Partitions the condensed degrees of freedom `c` into shells of constant graph distance
(in bonds) from the retained set `r`, farthest first -- the levels the multi-level scheme
[`reduce_matrices`](@ref) implements processes one at a time.

Two condensed degrees of freedom at different distances from `r` can only be directly
coupled if those distances differ by exactly one -- a basic property of breadth-first
distance layers, not an assumption about the geometry: a node's own distance is the
length of its *shortest* path to `r`, so a neighbour one hop closer or farther is
possible, but a neighbour two or more hops away in either direction never is (that would
mean a shorter path than the one that set that neighbour's own distance). Consequently a
shell can only be coupled to its two immediate neighbouring shells and, for the shell
adjacent to `r` itself, to `r`. Eliminating a whole shell together with the *entire*
previous interface is therefore always exact: nothing outside the current shell and the
interface it replaces can be reachable from what is being eliminated, however wide or
narrow the region is transverse to `r` -- unlike slicing by a fixed size, which has no
such guarantee and so cannot safely eliminate more than the one group it was handed.

# Arguments
- `K::AbstractMatrix`: Stiffness matrix, used only for its sparsity pattern
- `r::AbstractVector{<:Integer}`: Indices of the retained degrees of freedom
- `c::AbstractVector{<:Integer}`: Indices of the condensed degrees of freedom
# Returns
- `::Vector{Vector{Int64}}`: `c`, partitioned into shells, farthest first. Degrees of
  freedom with no path to `r` through `c` at all (an isolated cluster) form their own
  leading group, processed before the wavefront proper -- safe for the same reason any
  order is (see [`reduce_matrices`](@ref)), since they have no coupling to anything this
  ever eliminates alongside them either.
"""
function group_condensed_by_level(K::AbstractMatrix, r::AbstractVector{<:Integer},
                                  c::AbstractVector{<:Integer})
    c_set = Set(c)
    level = Dict{Int64,Int64}()

    # Neighbours of a single node come from one CSC column via nzrange -- O(its own
    # bond count), not from slicing K[frontier, :] / K[:, frontier] into a fresh
    # submatrix every BFS round, which for a wide frontier (thousands of degrees of
    # freedom at the same distance from r, exactly the domains this function exists
    # for) copies a large chunk of K's sparsity on every single level.
    Ksp = K isa SparseMatrixCSC ? K : sparse(K)
    rows = rowvals(Ksp)
    vals = nonzeros(Ksp)

    r_neighbors = Set{Int64}()
    for node in r
        for idx in nzrange(Ksp, node)
            vals[idx] != 0.0 && push!(r_neighbors, rows[idx])
        end
    end
    frontier = collect(Int64, intersect(r_neighbors, c_set))
    depth = 1
    for node in frontier
        level[node] = depth
    end
    while !isempty(frontier)
        depth += 1
        next_frontier = Int64[]
        for parent in frontier
            for idx in nzrange(Ksp, parent)
                vals[idx] == 0.0 && continue
                node = rows[idx]
                if (node in c_set) && !haskey(level, node)
                    level[node] = depth
                    push!(next_frontier, node)
                end
            end
        end
        frontier = next_frontier
    end
    @info "Craig-Bampton Cascade: condensed region spans $depth bond-hops from the " *
          "retained set, $(length(c)) condensed degrees of freedom to process"

    groups = Vector{Int64}[]
    unreached = [node for node in c if !haskey(level, node)]
    isempty(unreached) || push!(groups, unreached)
    for d in depth:-1:1
        shell = [node for node in c if get(level, node, -1) == d]
        isempty(shell) || push!(groups, shell)
    end
    return groups
end

"""
    transpose_structure(A)

The sparsity *pattern* of `transpose(A)`, as `(colptr, rowval)`, without ever allocating
a values array.

[`reduce_matrices`](@ref)'s level bookkeeping needs, for a node `A` stores by column,
which other nodes have an entry in its *row* -- the one direction a `SparseMatrixCSC`
cannot answer without either scanning every column or transposing. It never needs the
transposed values themselves: whether such an entry is actually a nonzero coupling (as
opposed to an explicitly stored zero) is instead checked with a direct lookup back into
`A`, once a candidate column is known (see the call site). Building the transpose the
usual way, `SparseMatrixCSC(transpose(A))`, allocates a second `nnz(A)`-length value
array purely to be discarded again right after -- half of what that step needs to
allocate, for nothing. Built from scratch here rather than reusing `A`'s own index type,
the result is `Int32` whenever `A`'s size and `nnz` fit -- ample range for any cascade
level in practice -- halving this step's own allocation again, and only ever `A`'s own
(usually `Int64`) index type for a problem actually too large for that.

# Arguments
- `A::SparseMatrixCSC`: The matrix
# Returns
- `colptr::Vector`, `rowval::Vector`: `transpose(A)`'s column pointers and row indices,
  in the same layout `SparseMatrixCSC` itself uses
"""
function transpose_structure(A::SparseMatrixCSC{Tv,Ti}) where {Tv,Ti}
    m, n = size(A)
    rows = rowvals(A)
    nz = length(rows)
    # Built from scratch here, so the index type is ours to choose, independent of A's
    # own Ti (usually Int64): Int32 halves this step's own allocation again on top of
    # dropping the values array, and is plenty of range for any cascade level this ever
    # sees in practice -- falls back to Ti only if the problem is actually that large.
    IdxT = (nz <= typemax(Int32) && m <= typemax(Int32) && n <= typemax(Int32)) ?
           Int32 : Ti

    row_counts = zeros(IdxT, m)
    for row in rows
        row_counts[row] += 1
    end
    colptr = Vector{IdxT}(undef, m + 1)
    colptr[1] = 1
    for row in 1:m
        colptr[row+1] = colptr[row] + row_counts[row]
    end
    rowval = Vector{IdxT}(undef, nz)
    next = copy(colptr)
    for col in 1:n
        for idx in nzrange(A, col)
            row = rows[idx]
            rowval[next[row]] = col
            next[row] += 1
        end
    end
    return colptr, rowval
end

"""
    reduce_matrices(K, M_diag, r, c, n_modes = 0; block_size = 512,
                    check_symmetry_sample = 2000)

Craig-Bampton reduction done as a cascade of local eliminations instead of one
factorization of the whole condensed region -- the same math as [`reduce_stage`](@ref),
called once per level on a local slice of `K`. This implements the accompanying paper's
multi-level scheme (its multi-level section): the condensed region is grown shell by
shell from the side farthest from `r` inward, and at every level the *entire*
superelement built so far is re-condensed together with the next shell, rather than
merely carried alongside it.

`c` is first partitioned by [`group_condensed_by_level`](@ref) into shells of constant
bond-distance from `r`, farthest first -- see its docstring for why a whole shell, not an
arbitrary chunk of it, is the unit this cascade can safely eliminate all at once. At each
level, the shell being eliminated is joined by the *entire* interface and modal state
carried in from the previous level (see below), a small local matrix is built for exactly
that union, and [`reduce_stage_dense`](@ref) is called on it. `Kcc` for every level is
therefore only ever the size of one shell plus whatever interface preceded it, never the
whole condensed region -- the memory and factorization cost this avoids scales with the
shell width transverse to `r`, not with the number of condensed degrees of freedom, and
does not accumulate from one level to the next.

Both the physical interface and the modal coordinates are *replaced*, not accumulated, at
every level -- this is what bounds the local system size regardless of how many levels
run. [`group_condensed_by_level`](@ref)'s docstring gives the reason this is exact for
the physical interface: a shell can only be coupled to its two neighbouring shells, so
the entire previous interface has no remaining coupling to anything outside the current
shell and can be eliminated alongside it without loss. For the modal coordinates, this is
the same tradeoff the paper states explicitly for its own levels: fixed-interface modes
of the *whole* condensed region would require its whole `Kcc`/`Mcc` at once, exactly what
this cascade avoids factorizing, so `n_modes` is instead the number of *local* modes
computed at every level, re-absorbed into the *next* level's condensation together with
the shell it eliminates, and replaced there by a fresh set of that level's own modes.
Truncation error therefore still accumulates across levels (the same cutoff/mode count
should be used throughout, well above whatever frequency band actually matters), even
though the state carried forward does not. `n_modes = 0` remains the exact, static
(Guyan) limit, unaffected by any of this.

# Arguments
- `K::AbstractMatrix`: Stiffness matrix
- `M_diag::AbstractVector`: Lumped mass matrix as a vector, one entry per degree of freedom
- `r::AbstractVector{<:Integer}`: Indices of the retained degrees of freedom
- `c::AbstractVector{<:Integer}`: Indices of the condensed degrees of freedom
- `n_modes::Integer`: Fixed-interface modes computed at *every* level; 0 is pure Guyan
# Keywords
- `block_size::Int64`: Right hand sides solved at once for the recovery modes, per level
- `check_symmetry_sample::Int64`: Entries drawn for the symmetry check on the full `K`,
  once, before the cascade starts; 0 disables it
- `max_rss_mib::Union{Nothing,Float64}`: Abort with an informative error once this
  process' peak resident memory exceeds this many MiB, checked once per level. A SIGKILL
  from the OS (or a job scheduler's own memory limit) cannot be caught by any Julia code
  -- the process is gone before anything, including a `try`/`catch`, gets to run -- so a
  run that is actually going to be killed only ever produces a level number and an
  interface size to debug from if something inside the process chooses to stop first, on
  its own terms. Defaults to 8192.0 (8 GiB); `nothing` disables this.
- `dense_shell_limit::Int64`: Above this many physical degrees of freedom in a single
  shell, that level's local system is built and factorized sparse instead of dense (see
  the loop body below) -- `Kcc`'s own physical part is exactly as sparse as `K` itself
  (peridynamic bonds are local), so a dense `zeros(n_full, n_full)` is `O(n_full^2)`
  memory spent on what is overwhelmingly zero once a shell gets wide. Below the limit,
  dense stays faster for the same reason [`fixed_interface_modes_dense`](@ref) prefers it
  for a typical, small shell. Defaults to 1500.
# Returns
- `K_reduced::SparseMatrixCSC`: Reduced stiffness, size `length(r)` plus `n_modes` -- the
  last level's own modes; earlier levels' modes were re-absorbed along the way, not
  accumulated (see above)
- `M_reduced::SparseMatrixCSC`: Reduced mass, same size
"""
function reduce_matrices(K::AbstractMatrix,
                         M_diag::AbstractVector,
                         r::AbstractVector{<:Integer},
                         c::AbstractVector{<:Integer},
                         n_modes::Integer = 0;
                         block_size::Int64 = 512,
                         check_symmetry_sample::Int64 = 2000,
                         max_rss_mib::Union{Nothing,Float64} = 8192.0,
                         dense_shell_limit::Int64 = 1500)
    @timeit "cascade setup" begin
        n_modes = Int64(n_modes)
        if n_modes < 0
            throw(ArgumentError("n_modes = $n_modes, must be >= 0."))
        end

        r = collect(Int64, r)
        c = collect(Int64, c)
        k_nnz = K isa AbstractSparseMatrix ? nnz(K) : length(K)
        @info "Craig-Bampton Cascade: entering reduce_matrices with $(length(r)) retained, " *
              "$(length(c)) condensed degrees of freedom, n_modes=$n_modes, nnz(K)=$k_nnz"
    end

    if check_symmetry_sample > 0
        @timeit "cascade symmetry check" check_symmetry(K; samples = check_symmetry_sample)
    end

    @timeit "cascade level grouping" begin
        groups = group_condensed_by_level(K, r, c)
        group_sizes = length.(groups)
        @info "Craig-Bampton Cascade: $(length(groups)) levels, shell size " *
              "min=$(minimum(group_sizes)) max=$(maximum(group_sizes)) " *
              "mean=$(round(sum(group_sizes) / length(group_sizes); digits = 1))"
    end

    # A shell many times the median, or simply large in absolute terms, is exactly what
    # drives the local system built for it in the loop below past dense_shell_limit --
    # checked here, before any of that memory is actually allocated, so an outlier shows
    # up in the log even on a run that goes on to complete.
    @timeit "cascade outlier check" begin
        sorted_sizes = sort(group_sizes)
        median_size = sorted_sizes[(length(sorted_sizes)+1)÷2]

        outlier_levels = findall(n -> n > max(4 * median_size, dense_shell_limit),
                                 group_sizes)
        if !isempty(outlier_levels)
            @warn "Craig-Bampton Cascade: $(length(outlier_levels)) of $(length(groups)) " *
                  "levels have an outlier-sized shell (median $median_size): levels " *
                  "$outlier_levels, sizes $(group_sizes[outlier_levels]). Their local " *
                  "system is built and factorized sparse instead of dense (see " *
                  "dense_shell_limit) to bound memory."
        end
    end

    # Every level needs, for its own shell, exactly which other degrees of freedom it
    # couples to in *either* direction (K need not be numerically symmetric, see
    # check_symmetry above) -- unlike group_condensed_by_level's ordering, this
    # determines the next interface itself and so must stay exact, not merely a good
    # heuristic. K[shell, :] is the expensive direction on a column-major
    # SparseMatrixCSC (it has to scan every stored entry of the whole matrix, not just
    # the shell's own), and doing that once per level made it effectively
    # O(levels * nnz(K)).
    #
    # K itself is the *whole system's* stiffness matrix, not just this reduction's part
    # of it -- Model_reduction.jl passes it in straight from Data_Manager, before it
    # knows anything about which blocks this particular reduction even touches (other
    # PD/FEM regions the reduction has nothing to do with, a contact pair, whatever else
    # the deck has, all still show up in K's sparsity). Restricting to r union c -- the
    # only degrees of freedom `touched` can ever keep, everything else is discarded by
    # the intersect with r_set/remaining_c right below regardless -- before transposing
    # bounds that one-time cost by this reduction's own size instead of the whole
    # system's.
    restricted = vcat(r, c)
    # Slicing to r union c only helps when it actually shrinks anything: once other
    # code (update_material_point_part rebuilding K without the material point nodes,
    # say) has already zeroed out everything K would otherwise have outside this
    # reduction's own degrees of freedom, restricted's own entry count comes out equal
    # to K's, and slicing it out is then a second full copy paid for nothing -- on top
    # of the transpose, which is needed regardless. Below half of K's own dimension is
    # a cheap, correctness-irrelevant proxy for "worth doing": it only ever affects how
    # this one-time cost is paid, never the result (both branches build exactly the same
    # sparsity information, just addressed differently).
    n_total = size(K, 1)
    restrict_worthwhile = length(restricted) < n_total ÷ 2
    local to_local, to_global
    if restrict_worthwhile
        id_map = Dict{Int64,Int64}(g => i for (i, g) in enumerate(restricted))
        to_local = node -> id_map[node]
        to_global = idx -> restricted[idx]
        @timeit "cascade restrict" Ksp=sparse(K[restricted, restricted])
    else
        to_local = identity
        to_global = identity
        @timeit "cascade restrict" Ksp=(K isa SparseMatrixCSC ? K : sparse(K))
    end
    @timeit "cascade transpose" begin
        (Kt_colptr, Kt_rowval) = transpose_structure(Ksp)
        Krows = rowvals(Ksp)
        Kvals = nonzeros(Ksp)
    end

    @timeit "cascade loop setup" begin
        r_set = Set(r)
        r_index = Dict(node => k for (k, node) in enumerate(r))
        nr = length(r)
        remaining_c = Set(c)
        incoming_pos = Int64[]        # physical dof indices carried forward (⊆ current shell)
        incoming_modal = 0            # count of modal coordinates carried forward
        incoming_K = zeros(0, 0)
        incoming_M = zeros(0, 0)
    end

    level_num = 0
    for shell in groups
        @timeit "cascade level bookkeeping" begin
            level_num += 1
            setdiff!(remaining_c, shell)
            shell_set = Set(shell)

            # What this shell actually couples to, outside itself: into the retained
            # set, and into whatever of c has not been eliminated yet. incoming_pos --
            # the entire previous interface -- is *not* unioned in here: by
            # group_condensed_by_level's BFS-layering guarantee it is already a subset
            # of this shell (the previous, deeper shell could only reach this one or
            # its own), so it is eliminated below together with the shell, not carried
            # past it.
            touched = Set{Int64}()
            for node in shell
                lnode = to_local(node)
                for idx in nzrange(Ksp, lnode)
                    Kvals[idx] != 0.0 && push!(touched, to_global(Krows[idx]))
                end
                # Kt_rowval holds, for row lnode, the columns Ksp stores an entry at --
                # transpose_structure never kept their values, so whether one is an
                # actual nonzero coupling (as opposed to an explicitly stored zero) is
                # checked directly against Ksp itself instead.
                for idx in Kt_colptr[lnode]:(Kt_colptr[lnode + 1] - 1)
                    col = Kt_rowval[idx]
                    Ksp[lnode, col] != 0.0 && push!(touched, to_global(col))
                end
            end
            setdiff!(touched, shell_set)
            new_interface = sort(collect(Int64,
                                         union(intersect(touched, r_set),
                                               intersect(touched, remaining_c))))
        end

        @timeit "cascade level local assembly" begin
            full_nodes = vcat(shell, new_interface)
            n_full_physical = length(full_nodes)
            n_full = n_full_physical + incoming_modal

            # incoming_pos/incoming_modal from the previous level are substituted as a
            # block, replacing (not adding to) whatever raw value sits there: that block
            # already is the complete effective stiffness/mass (own value plus
            # everything eliminated so far), so adding it would count the interface's
            # own stiffness twice. incoming_pos sits within `shell` itself (see above);
            # modal coordinates have no dof index of their own and so no position in
            # full_nodes, always occupying the last incoming_modal local positions.
            has_incoming = !isempty(incoming_pos) || incoming_modal > 0
            combined_idx = Int64[]
            if has_incoming
                local_idx_physical = [findfirst(==(p), full_nodes) for p in incoming_pos]
                local_idx_modal = (n_full_physical + 1):(n_full_physical + incoming_modal)
                combined_idx = vcat(local_idx_physical, local_idx_modal)
            end

            # Dense here costs O(n_full^2); fine for a typical shell, wasteful for an
            # outlier one (see the level-grouping check above), since the shell's own
            # physical part is exactly as sparse as K itself. Above dense_shell_limit,
            # both matrices are instead assembled once from triplets, exactly like the
            # final assembly at the end of this function: setindex! with a range or
            # vector index into an *existing* SparseMatrixCSC is one of the slowest
            # operations SparseArrays has -- every inserted structural nonzero can shift
            # every following column's storage -- so writing the (near-fully dense, see
            # incoming_K/incoming_M above) incoming block in that way, one entry at a
            # time, costs orders of magnitude more than building the whole sparse
            # matrix once from its final triplets.
            @timeit "cascade matrix alloc" begin
                if n_full_physical > dense_shell_limit
                    Kf_rows, Kf_cols, Kf_vals = findnz(sparse(K[full_nodes, full_nodes]))
                    Mf_diag = M_diag[full_nodes]
                    if has_incoming
                        covered = Set(combined_idx)
                        keep = [!(Kf_rows[k] in covered && Kf_cols[k] in covered)
                                for k in eachindex(Kf_rows)]
                        Ki_rows, Ki_cols, Ki_vals = findnz(sparse(incoming_K))
                        k_rows = vcat(Kf_rows[keep], combined_idx[Ki_rows])
                        k_cols = vcat(Kf_cols[keep], combined_idx[Ki_cols])
                        k_vals = vcat(Kf_vals[keep], Ki_vals)

                        m_uncovered = [i for i in 1:n_full_physical if !(i in covered)]
                        Mi_rows, Mi_cols, Mi_vals = findnz(sparse(incoming_M))
                        m_rows = vcat(m_uncovered, combined_idx[Mi_rows])
                        m_cols = vcat(m_uncovered, combined_idx[Mi_cols])
                        m_vals = vcat(Mf_diag[m_uncovered], Mi_vals)
                    else
                        k_rows, k_cols, k_vals = Kf_rows, Kf_cols, Kf_vals
                        m_rows = m_cols = collect(1:n_full_physical)
                        m_vals = Mf_diag
                    end
                    K_local = sparse(k_rows, k_cols, k_vals, n_full, n_full)
                    M_local = sparse(m_rows, m_cols, m_vals, n_full, n_full)
                else
                    K_local = zeros(n_full, n_full)
                    M_local = zeros(n_full, n_full)
                    K_local[1:n_full_physical, 1:n_full_physical] = K[full_nodes,
                                                                      full_nodes]
                    M_local[1:n_full_physical,
                            1:n_full_physical] = Diagonal(M_diag[full_nodes])
                    if has_incoming
                        K_local[combined_idx, combined_idx] .= incoming_K
                        M_local[combined_idx, combined_idx] .= incoming_M
                    end
                end
            end
        end
        @timeit "cascade level partition" begin
            # The entire shell -- including wherever incoming_pos sits within it -- plus
            # the previous level's modal coordinates are eliminated together; only the
            # fresh new_interface stays retained into the next level. See
            # group_condensed_by_level and this function's own docstring for why this is
            # exact.
            c_local = vcat(collect(1:length(shell)), (n_full_physical + 1):n_full)
            r_local = collect((length(shell) + 1):n_full_physical)
        end

        K_stage,
        M_stage=@timeit "cascade level reduction" if n_modes > 0
            reduce_stage_dense(K_local, M_local, r_local, c_local, n_modes;
                               block_size = block_size)
        elseif has_incoming
            reduce_stage_dense(K_local, M_local, r_local, c_local; block_size = block_size)
        else
            reduce_stage(K_local, diag(M_local), r_local, c_local, 0;
                         block_size = block_size, check_symmetry_sample = 0)
        end

        @timeit "cascade level carry update" begin
            incoming_pos = new_interface
            incoming_modal = n_modes
            incoming_K = Matrix(K_stage)
            incoming_M = Matrix(M_stage)
        end

        @timeit "cascade level memory check" begin
            # Peak resident memory, not just what this process itself allocated: an
            # external OOM kill depends on the whole machine's memory state at that
            # instant, so which level it happens to land on is not reproducible run to run
            # even for the exact same deck -- Sys.maxrss() is this process' own high-water
            # mark and gives a deterministic figure to correlate against interface size
            # instead.
            rss_mib = Sys.maxrss() / 2^20
            @info "Craig-Bampton Cascade: level $level_num, $(length(shell)) eliminated, " *
                  "$(length(new_interface)) carried forward, $(length(remaining_c)) left, " *
                  "$incoming_modal modes, peak RSS $(round(rss_mib; digits = 1)) MiB"

            # An OS (or job scheduler) OOM kill is a SIGKILL: no Julia code, including a
            # try/catch around this whole call, ever gets to run once it happens, so it
            # cannot be caught after the fact -- only avoided by stopping deliberately,
            # before the real limit is hit, with an ordinary catchable Julia exception that
            # at least says which level and how large the interface was.
            readline()
            #if !isnothing(max_rss_mib) && rss_mib > max_rss_mib
            #    throw(ErrorException("Craig-Bampton Cascade: peak resident memory " *
            #        "$(round(rss_mib; digits = 1)) MiB exceeded the " *
            #        "$(max_rss_mib) MiB limit at level $level_num " *
            #        "($(length(new_interface)) carried forward, " *
            #        "$(length(remaining_c)) left). Raise max_rss_mib."))
            #end
        end
    end
    @timeit "assemble final reduced matrices" begin
        # A retained degree of freedom with no bond into c at all (an isolated PD node, say)
        # never appears in any level's coupling and so never enters incoming_pos; its row and
        # column are simply the raw, untouched K[r,r]/M[r,r] -- exactly like reduce_stage,
        # where such a row of Krc is all zero and Krcrc reduces to Krr there. Starting from
        # the raw values and overwriting only what the cascade actually touched covers both
        # cases without a lookup that could miss an untouched degree of freedom.
        readline()
        # The result is sparse, and for a large retained set nr is itself large -- a
        # dense (nr+modes)x(nr+modes) intermediate (as a naive "start from raw, overwrite
        # what changed" would build) costs O(nr^2) memory for a matrix that is
        # overwhelmingly zero. Every piece below is already available as triplets or as
        # a block small enough to densify on its own (incoming_K/incoming_M, bounded
        # by the very last level's own interface, not by nr), so the result is
        # assembled directly as one.
        touched_order = [r_index[node] for node in incoming_pos]
        modal_order = (nr + 1):(nr + incoming_modal)
        combined_order = vcat(touched_order, modal_order)

        # A raw K[r,r]/M[r,r] entry is stale wherever the last level's own carried-forward
        # block replaced it -- that substitution is a correction, never an addition (see
        # above), so any raw entry with both its row and column covered by it must be
        # dropped, not summed alongside its replacement.
        covered = Set(touched_order)

        Kr_rows, Kr_cols, Kr_vals = findnz(sparse(K[r, r]))
        keep = [!(Kr_rows[k] in covered && Kr_cols[k] in covered)
                for k in eachindex(Kr_rows)]
        rows = Kr_rows[keep]
        cols = Kr_cols[keep]
        vals = Float64.(Kr_vals[keep])

        m_uncovered = [i for i in 1:nr if !(i in covered)]
        mrows = copy(m_uncovered)
        mcols = copy(m_uncovered)
        mvals = M_diag[r][m_uncovered]

        if !isempty(combined_order)
            Ki_rows, Ki_cols, Ki_vals = findnz(sparse(incoming_K))
            append!(rows, combined_order[Ki_rows])
            append!(cols, combined_order[Ki_cols])
            append!(vals, Ki_vals)

            Mi_rows, Mi_cols, Mi_vals = findnz(sparse(incoming_M))
            append!(mrows, combined_order[Mi_rows])
            append!(mcols, combined_order[Mi_cols])
            append!(mvals, Mi_vals)
        end

        n_total_reduced = nr + incoming_modal
        K_reduced = sparse(rows, cols, vals, n_total_reduced, n_total_reduced)
        M_reduced = sparse(mrows, mcols, mvals, n_total_reduced, n_total_reduced)
        dropzeros!(K_reduced)
        dropzeros!(M_reduced)
    end
    return K_reduced, M_reduced
end

end
