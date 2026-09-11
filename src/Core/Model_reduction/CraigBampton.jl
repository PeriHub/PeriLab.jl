# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module CraigBampton
using LinearAlgebra
using SparseArrays
using TimerOutputs: @timeit
using Arpack: eigs
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
    nonzero_rows(A)

Indices of the rows of a matrix that hold at least one entry.

A condensed degree of freedom whose row is empty carries no stiffness at all — an
isolated node, or one whose bonds were all cut by a bond filter. Such a row makes `Kll`
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

For `Klr` these are the retained degrees of freedom with a bond into the condensed
region — the coupling layer. Every other column is empty, so its recovery mode is zero
and never has to be solved for. That is the difference between one triangular solve per
retained degree of freedom and one per coupling degree of freedom, and it also bounds
where the reduced matrices can fill in: `Krl * B` has entries only in the rows and
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
    factorize_condensed(Kll)

Factorization of the condensed stiffness, Cholesky where possible.

With the interface held fixed `Kll` is symmetric positive definite. A sparse Cholesky is
faster to build than an LU and considerably faster to apply to a dense block of right
hand sides, because CHOLMOD goes through supernodal BLAS3 kernels while UMFPACK works
column by column. With one right hand side per coupling degree of freedom that
difference dominates the whole reduction.

A failure is informative in itself: `Kll` is then not positive definite, which for a
fixed interface points at a condensed region falling apart into pieces not tied to any
retained degree of freedom.

# Arguments
- `Kll::AbstractMatrix`: Stiffness of the condensed part
# Returns
- A factorization object supporting `ldiv!`
"""
function factorize_condensed(Kll::AbstractMatrix)
    try
        return cholesky(Symmetric(Kll))
    catch err
        @warn "Cholesky of the condensed stiffness failed, falling back to LU. Kll is " *
              "not positive definite: check whether every condensed node is tied to " *
              "the retained region, an isolated node or a detached cluster has a " *
              "rigid body mode." exception=err
        return lu(Kll)
    end
end

"""
    ShiftInvertOperator

Applies `D Kll^-1 D` with `D = sqrt(Mll)`, as a matrix free operator for Arpack.

The generalised problem `Kll X = Mll X W` becomes the standard symmetric one
`A u = w u` with `A = D^-1 Kll D^-1` and `u = D x`. Its smallest eigenvalues are the
largest of `A^-1 = D Kll^-1 D`, which is what this operator applies. Two things follow:
shift-invert without asking Arpack for it, so an existing factorization is reused
instead of a second one being built, and eigenvectors orthonormal in the standard sense,
so `x = D^-1 u` is mass normalised without a further pass.

# Fields
- `factorization`: Factorization of `Kll`
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
    fixed_interface_modes(Kll, Mll, n_modes; tol, maxiter)

The `n_modes` lowest eigenpairs of `Kll X = Mll X W`, mass normalised.

Rows of `Kll` without any entry are excluded from the eigenvalue problem and their
entries in the modes stay zero, so that an isolated condensed degree of freedom does not
make the problem singular.

The eigenpairs come from shift-invert Lanczos through [`ShiftInvertOperator`](@ref), so
`Kll` is never densified. A dense `eigen` would cost `O(nl^2)` memory and `O(nl^3)` time
for all `nl` eigenpairs, of which `n_modes` are wanted. Below a few hundred degrees of
freedom the dense route is the cheaper one and is taken instead.

The eigenvectors are mass normalised, `X' Mll X = I`. Without that the modal blocks of
the reduced matrices have no defined scaling.

# Arguments
- `Kll::AbstractMatrix`: Stiffness of the condensed part
- `Mll::Diagonal`: Lumped mass of the condensed part
- `n_modes::Int64`: Number of modes
# Keywords
- `tol::Float64`: Relative tolerance handed to Arpack
- `maxiter::Int64`: Iteration limit per attempt
# Returns
- `X::Matrix{Float64}`: Modes, `nl x n_modes`
- `w::Vector{Float64}`: Eigenvalues, ascending
"""
function fixed_interface_modes(Kll::AbstractMatrix, Mll::Diagonal, n_modes::Int64;
                               tol::Float64 = 1e-10, maxiter::Int64 = 1000)
    nl = size(Kll, 1)
    X = zeros(nl, n_modes)
    n_modes == 0 && return X, Float64[]

    active = nonzero_rows(Kll)
    if length(active) < nl
        @warn "$(nl - length(active)) of $nl condensed degrees of freedom carry no " *
              "stiffness and are excluded from the eigenvalue problem. Check whether " *
              "every condensed node is still tied to the retained region."
    end
    if n_modes > length(active)
        throw(ArgumentError("n_modes = $n_modes, but only $(length(active)) condensed " *
                            "degrees of freedom carry stiffness."))
    end

    K_active = Kll[active, active]
    mass_active = diag(Mll)[active]
    scale = sqrt.(mass_active)

    if !issparse(Kll) || length(active) <= 500
        # The scaling by sqrt(M) turns the generalised problem into a standard symmetric
        # one, whose eigenvectors come out orthonormal — which is the mass
        # normalisation, for free.
        inverse_scale = Diagonal(inv.(scale))
        factorization = eigen(Symmetric(inverse_scale * Matrix(K_active) *
                                        inverse_scale))
        X[active, :] = inverse_scale * factorization.vectors[:, 1:n_modes]
        return X, factorization.values[1:n_modes]
    end

    operator = ShiftInvertOperator(factorize_condensed(K_active), scale,
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
    recovery_modes(factorization, Klr, coupling; block_size = 512)

Recovery modes `B = Kll^-1 Klr` for the coupling degrees of freedom only.

Columns of `Klr` outside `coupling` are empty, so their recovery modes are zero and are
neither solved for nor stored. The result is `nl x length(coupling)` instead of
`nl x nr`.

Solved in column blocks so that the dense right hand side buffer stays at
`nl x block_size` rather than being a second matrix of the full size.

# Arguments
- `factorization`: Factorization of `Kll`
- `Klr::AbstractMatrix`: Coupling block, condensed to retained
- `coupling::Vector{Int64}`: Columns to solve for
# Keywords
- `block_size::Int64`: Right hand sides solved at once
# Returns
- `::Matrix{Float64}`: The recovery modes, `nl x length(coupling)`
"""
function recovery_modes(factorization, Klr::AbstractMatrix, coupling::Vector{Int64};
                        block_size::Int64 = 512)
    nl = size(Klr, 1)
    nc = length(coupling)
    B = Matrix{Float64}(undef, nl, nc)
    nc == 0 && return B

    width = min(block_size, nc)
    rhs = Matrix{Float64}(undef, nl, width)

    for first in 1:width:nc
        last = min(first + width - 1, nc)
        local_columns = first:last
        buffer = view(rhs, :, 1:length(local_columns))
        fill!(buffer, 0.0)
        copy_sparse_columns!(buffer, Klr, view(coupling, local_columns))
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
    check_symmetry(K, r, l; samples = 2000, tolerance = 1e-8)

Warns if the stiffness is not symmetric, judged from a random sample of entries.

Craig-Bampton assumes symmetry: the fixed-interface modes use the symmetric solver, and
the coupling blocks of the reduced stiffness are dropped because `Kll B = Klr` makes
them cancel, which needs `Kll' = Kll`. In peridynamics the stiffness is not symmetric in
general, two points with different horizons do not contribute equally to each other. A
sample is compared rather than the full matrix, which would cost as much memory as the
reduction saves.

# Arguments
- `K::AbstractMatrix`: The full stiffness matrix
- `r::Vector{Int64}`, `l::Vector{Int64}`: Retained and condensed index sets
# Keywords
- `samples::Int64`: Number of entries drawn
- `tolerance::Float64`: Relative deviation accepted as rounding
# Returns
- `::Bool`: Whether the sample looked symmetric
"""
function check_symmetry(K::AbstractMatrix, r::Vector{Int64}, l::Vector{Int64};
                        samples::Int64 = 2000, tolerance::Float64 = 1e-8)
    indices = vcat(r, l)
    n = length(indices)
    n < 2 && return true

    scale = 0.0
    deviation = 0.0
    for _ in 1:samples
        i = indices[rand(1:n)]
        j = indices[rand(1:n)]
        upper = K[i, j]
        lower = K[j, i]
        scale = max(scale, abs(upper), abs(lower))
        deviation = max(deviation, abs(upper - lower))
    end

    scale = max(scale, eps())
    if deviation / scale > tolerance
        @warn "Stiffness matrix appears not to be symmetric (relative deviation " *
              "$(round(deviation / scale; sigdigits = 3))). Craig-Bampton assumes " *
              "symmetry; unequal horizons in the condensed region are the usual cause."
        return false
    end
    return true
end

"""
    add_dense_block!(rows, columns, values, block, positions)

Appends a dense block to a coordinate list, mapped to the given global positions.

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
    report_modal_coupling(M_reduced, nr, n_modes)

Reports the mass coupling between the retained and the modal part, and warns if it is
empty.

The modes are driven through the mass, not through the stiffness: `Kbm` is zero because
recovery modes and fixed-interface modes are K-orthogonal, so the only path from a
moving retained degree of freedom to a modal amplitude is `Mbm`. An empty `Mbm` leaves
the modes at rest for the whole simulation — they still carry mass and shift the
eigenfrequencies, but never respond, which makes the result worse than the static
condensation it was meant to improve on.

An empty block usually means the mass coupling was dropped by the threshold in
[`add_dense_block!`](@ref), or that `Mbm` came out with the wrong sign convention and
cancelled itself.

# Arguments
- `M_reduced::SparseMatrixCSC`: The reduced mass matrix
- `nr::Int64`: Number of retained degrees of freedom
- `n_modes::Int64`: Number of modes
# Returns
- `::Bool`: Whether the coupling holds entries
"""
function report_modal_coupling(M_reduced::SparseMatrixCSC, nr::Int64, n_modes::Int64)
    n_modes == 0 && return true

    coupling_block = M_reduced[1:nr, (nr + 1):(nr + n_modes)]
    entries = nnz(coupling_block)

    if entries == 0
        @warn "The mass coupling between the retained and the modal part is empty. " *
              "The modes are driven through the mass only, so they will stay at rest " *
              "and the reduction behaves worse than a static condensation."
        return false
    end

    modal_mass = maximum(abs, nonzeros(coupling_block))
    retained_mass = maximum(abs, nonzeros(M_reduced[1:nr, 1:nr]); init = 0.0)
    @info "Modal coupling: $entries entries, largest $(round(modal_mass; sigdigits = 3)) " *
          "against $(round(retained_mass; sigdigits = 3)) in the retained block"

    return true
end

"""
    reduce_matrices(K, M_diag, r, l, n_modes = 1; block_size = 512,
                    check_symmetry_sample = 2000)

Craig-Bampton reduction of a stiffness and a lumped mass matrix.

The degrees of freedom are split into the retained ones `r` and the condensed ones `l`.
The retained ones stay physical coordinates, the condensed ones are represented by the
`n_modes` lowest fixed-interface modes `X` of `Kll X = Mll X W`:

    | u_r |   | I   0 | | u_r |
    | u_l | = | -B  X | | eta |

with the recovery modes `B = Kll^-1 Klr`, the same static relation Guyan uses. The
reduced blocks are

    Kbb = Krr - Krl B          Kmm = W
    Mbb = Mrr + B' Mll B       Mmm = I
    Mbm = -B' Mll X            Kbm = 0

The stiffness coupling blocks vanish: `Kll B = Klr` makes the two contributions
`Krl X` and `-B' Kll X` cancel, so recovery modes and fixed-interface modes are
K-orthogonal. The mass has no such cancellation. The mass coupling terms `Mrl X` vanish
as well, but for a different reason: the mass matrix is diagonal and the two index sets
are disjoint.

The reduced matrices stay sparse. Only retained degrees of freedom with a bond into the
condensed region appear in `Klr`; the reduction fills in exactly their rows and columns,
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
- `l::AbstractVector{<:Integer}`: Indices of the condensed degrees of freedom
- `n_modes::Integer`: Number of fixed-interface modes to keep
# Keywords
- `block_size::Int64`: Right hand sides solved at once for the recovery modes
- `check_symmetry_sample::Int64`: Entries drawn for the symmetry check, 0 disables it
# Returns
- `K_reduced::SparseMatrixCSC`: Reduced stiffness, size `length(r) + n_modes`
- `M_reduced::SparseMatrixCSC`: Reduced mass, same size
"""
function reduce_matrices(K::AbstractMatrix,
                         M_diag::AbstractVector,
                         r::AbstractVector{<:Integer},
                         l::AbstractVector{<:Integer},
                         n_modes::Integer = 1;
                         block_size::Int64 = 512,
                         check_symmetry_sample::Int64 = 2000)
    r = collect(Int64, r)
    l = collect(Int64, l)
    n_modes = Int64(n_modes)

    nr = length(r)
    nl = length(l)

    if !isempty(intersect(r, l))
        throw(ArgumentError("Retained and condensed index sets overlap."))
    end
    if n_modes < 0 || n_modes > nl
        throw(ArgumentError("n_modes = $n_modes, but only $nl condensed degrees of " *
                            "freedom are available."))
    end

    if check_symmetry_sample > 0
        @timeit "CB symmetry check" check_symmetry(K, r, l;
                                                   samples = check_symmetry_sample)
    end

    @timeit "CB extract submatrices" begin
        Krr = K[r, r]
        Krl = K[r, l]
        Klr = K[l, r]
        Kll = K[l, l]
        # M_diag is a mass density, matching the force densities K*u delivers; see
        # Guyan.reduce_matrices for why no volume weighting belongs here.
        mass_l = M_diag[l]
        mass_r = M_diag[r]
        Mll = Diagonal(mass_l)
    end

    coupling = nonzero_columns(Klr)
    nc = length(coupling)

    @info "Craig-Bampton: $nr retained, $nl condensed, $nc coupling, $n_modes modes; " *
          "reduced size $(nr + n_modes), dense blocks of " *
          "$(round(nc^2 * 8 / 2^20; digits = 1)) MiB"

    @timeit "CB factorize Kll" factorization=factorize_condensed(Kll)

    @timeit "CB recovery modes" B=recovery_modes(factorization, Klr, coupling;
                                                 block_size = block_size)

    # The modes are normalised against Mll, so X' Mll X = I.
    @timeit "CB fixed interface modes" X, w=fixed_interface_modes(Kll, Mll, n_modes)

    n_total = nr + n_modes

    @timeit "CB reduced stiffness" begin
        rows, columns, values = findnz(Krr)

        if nc > 0
            # Only the coupling rows of Krl carry entries, so the product is nc x nc.
            Kbb_fill = Matrix{Float64}(undef, nc, nc)
            mul!(Kbb_fill, Krl[coupling, :], B, -1.0, 0.0)
            add_dense_block!(rows, columns, values, Kbb_fill, coupling)
        end
        for k in 1:n_modes
            push!(rows, nr + k)
            push!(columns, nr + k)
            push!(values, w[k])
        end
        K_reduced = sparse(rows, columns, values, n_total, n_total)
    end

    @timeit "CB reduced mass" begin
        # Both mass blocks have the form Y' Mll Z with a diagonal Mll, so scaling the
        # factors by sqrt(Mll) turns them into plain products. B is scaled in place,
        # which is why the stiffness block above had to be formed first.
        root_mass = sqrt.(mass_l)
        @inbounds for j in 1:nc, i in 1:nl
            B[i, j] *= root_mass[i]
        end

        rows = Int64[]
        columns = Int64[]
        values = Float64[]

        if nc > 0
            Mbb_fill = Matrix{Float64}(undef, nc, nc)
            mul!(Mbb_fill, transpose(B), B)
            add_dense_block!(rows, columns, values, Mbb_fill, coupling)
        end

        for local_index in eachindex(r)
            push!(rows, local_index)
            push!(columns, local_index)
            push!(values, mass_r[local_index])
        end

        if n_modes > 0 && nc > 0
            # X is not needed unscaled again, so it is scaled in place like B was above.
            @inbounds for k in 1:n_modes, i in 1:nl
                X[i, k] *= root_mass[i]
            end
            # Mbm = -B' Mll X, the minus coming from u_l = -B u_r + X eta.
            Mbm = Matrix{Float64}(undef, nc, n_modes)
            mul!(Mbm, transpose(B), X, -1.0, 0.0)

            threshold = 1.0e-12 * max(maximum(abs, Mbm; init = 0.0), eps())
            @inbounds for k in 1:n_modes, i in 1:nc
                value = Mbm[i, k]
                abs(value) <= threshold && continue
                push!(rows, coupling[i])
                push!(columns, nr + k)
                push!(values, value)
                push!(rows, nr + k)
                push!(columns, coupling[i])
                push!(values, value)
            end
        end

        for k in 1:n_modes
            push!(rows, nr + k)
            push!(columns, nr + k)
            push!(values, 1.0)
        end

        M_reduced = sparse(rows, columns, values, n_total, n_total)
    end

    report_modal_coupling(M_reduced, nr, n_modes)

    @debug "Reduced matrices: K with $(nnz(K_reduced)) entries, M with " *
           "$(nnz(M_reduced)) entries, out of $(n_total^2) possible"

    return K_reduced, M_reduced
end

end
