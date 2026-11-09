# SPDX-FileCopyrightText: 2026 Christian Willberg
#
# SPDX-License-Identifier: BSD-3-Clause

module BiCGStab

using KernelAbstractions: @kernel, @index, @Const
using LinearAlgebra: mul!
using TimerOutputs: @timeit
using ..Sparse_CSR: CSRMatrix, launch!
using ..Dof_Sets: DofSet
using ..Distributed_Operators: AbstractSpace, inner, n_own
using ..Halo_Exchange: accumulate_ghosts!

export BiCGStabWorkspace, jacobi!, solve_bicgstab!

"""
    BiCGStabWorkspace(template)

Preallocated vectors for BiCGStab (README A5.1) with right Jacobi preconditioning. All
vectors have the local length of `template` (own + ghost degrees of freedom) and its
backend. 10 vectors, independent of the iteration count (GMRES(m) needs m + 4).
"""
struct BiCGStabWorkspace{V<:AbstractVector{Float64}}
    x::V
    r::V
    r0::V
    p::V
    v::V
    s::V
    t::V
    phat::V
    shat::V
    diag_inv::V
end
function BiCGStabWorkspace(template::AbstractVector{Float64})
    vectors = [fill!(similar(template), 0.0) for _ in 1:9]
    diag_inv = fill!(similar(template), 1.0)
    return BiCGStabWorkspace(vectors..., diag_inv)
end

@kernel function diagonal_kernel!(d, @Const(rowptr), @Const(colval), @Const(nzval))
    i = @index(Global)
    value = zero(eltype(d))
    @inbounds for k in rowptr[i]:(rowptr[i + 1] - 1)
        if colval[k] == i
            value = nzval[k]
        end
    end
    @inbounds d[i] = value
end

@kernel function add_shift_kernel!(d, @Const(shift))
    i = @index(Global)
    @inbounds d[i] += shift[i]
end

@kernel function invert_kernel!(d, @Const(prescribed))
    i = @index(Global)
    @inbounds d[i] = (prescribed[i] || d[i] == zero(eltype(d))) ? one(eltype(d)) :
                     one(eltype(d)) / d[i]
end

"""
    jacobi!(workspace, K, prescribed; halo = nothing, accumulate = false, n_own = K.n_rows,
            shift = nothing)

Jacobi preconditioner from the diagonal (identity on prescribed degrees of freedom).
Recompute after a new assembly. Locally the own degrees of freedom are numbered first, so
the diagonal entry of row `i` is in column `i`; ghost rows hold the diagonal entries of
the ghost degrees of freedom. With `accumulate` the partial diagonals of the ghost rows
are added to their owners first (README A4.0b) and then inverted. `shift` (own degrees of
freedom) is added to the diagonal before the inversion (`ShiftedOperator`).
"""
function jacobi!(ws::BiCGStabWorkspace, K::CSRMatrix, prescribed::DofSet;
                 halo = nothing, accumulate::Bool = false, n_own::Int64 = K.n_rows,
                 shift = nothing)
    @timeit "jacobi" begin
        launch!(diagonal_kernel!, ws.diag_inv, K.n_rows, ws.diag_inv, K.rowptr, K.colval,
                K.nzval)
        accumulate && accumulate_ghosts!(halo, ws.diag_inv)
        shift === nothing ||
            launch!(add_shift_kernel!, ws.diag_inv, n_own, ws.diag_inv, shift)
        launch!(invert_kernel!, ws.diag_inv, n_own, ws.diag_inv, prescribed.mask)
    end
    return ws
end

"""
    solve_bicgstab!(workspace, A, b, space; rtol, itmax, x0 = nothing)

Solves `A x = b` with right-preconditioned BiCGStab, starting from `x0` (local vector, zero
at prescribed degrees of freedom) or zero. `A` is any operator with `mul!`
on local vectors (e.g. `DistributedOperator`); `space` defines the own part and the inner
product (serial or MPI). No allocation of vectors inside the iteration.

# Returns
- `x`: Solution (owned by the workspace)
- `stats::NamedTuple`: `solved`, `niter`, `residual` (relative)
"""
function solve_bicgstab!(ws::BiCGStabWorkspace, A, b::AbstractVector, space::AbstractSpace;
                         rtol::Float64 = 1e-10, itmax::Int64 = 10_000, x0 = nothing)
    @timeit "solver (BiCGStab)" begin
        n = n_own(space)
        x = view(ws.x, 1:n)
        r = view(ws.r, 1:n)
        r0 = view(ws.r0, 1:n)
        p = view(ws.p, 1:n)
        v = view(ws.v, 1:n)
        s = view(ws.s, 1:n)
        t = view(ws.t, 1:n)
        phat = view(ws.phat, 1:n)
        shat = view(ws.shat, 1:n)
        d = view(ws.diag_inv, 1:n)
        bo = view(b, 1:n)

        fill!(ws.p, 0.0)
        fill!(ws.v, 0.0)
        b_norm = sqrt(inner(space, b, b))
        if b_norm == 0.0
            fill!(ws.x, 0.0)
            return ws.x, (solved = true, niter = 0, residual = 0.0)
        end
        if x0 === nothing
            fill!(ws.x, 0.0)
            r .= bo
        else
            copyto!(ws.x, x0)
            mul!(ws.v, A, ws.x)
            @. r = bo - v
            fill!(ws.v, 0.0)
        end
        r0 .= r
        tol = rtol * b_norm
        residual = sqrt(inner(space, ws.r, ws.r))
        residual <= tol &&
            return ws.x, (solved = true, niter = 0, residual = residual / b_norm)
        rho = alpha = omega = 1.0
        for k in 1:itmax
            rho_new = inner(space, ws.r0, ws.r)
            beta = (rho_new / rho) * (alpha / omega)
            @. p = r + beta * (p - omega * v)
            @. phat = d * p
            mul!(ws.v, A, ws.phat)
            alpha = rho_new / inner(space, ws.r0, ws.v)
            @. s = r - alpha * v
            residual = sqrt(inner(space, ws.s, ws.s))
            if residual <= tol
                @. x += alpha * phat
                return ws.x, (solved = true, niter = k, residual = residual / b_norm)
            end
            @. shat = d * s
            mul!(ws.t, A, ws.shat)
            omega = inner(space, ws.t, ws.s) / inner(space, ws.t, ws.t)
            @. x += alpha * phat + omega * shat
            @. r = s - omega * t
            residual = sqrt(inner(space, ws.r, ws.r))
            if residual <= tol
                return ws.x, (solved = true, niter = k, residual = residual / b_norm)
            end
            rho = rho_new
        end
        return ws.x, (solved = false, niter = itmax, residual = residual / b_norm)
    end
end

end
