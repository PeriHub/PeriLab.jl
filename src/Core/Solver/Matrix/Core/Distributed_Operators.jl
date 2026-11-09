# SPDX-FileCopyrightText: 2026 Christian Willberg
#
# SPDX-License-Identifier: BSD-3-Clause

module Distributed_Operators

using MPI
using LinearAlgebra: LinearAlgebra, dot, mul!
using ..Halo_Exchange: exchange!, accumulate_ghosts!

export AbstractSpace, SerialSpace, MPISpace, inner, n_own, DistributedOperator,
       needs_accumulate

"""
    AbstractSpace

Vector space the solver works in. Vectors always have the local length (own + ghost
degrees of freedom); the solver only touches the own part `1:n_own(space)`, the inner
product sums over all ranks.
"""
abstract type AbstractSpace end

"Single process (no MPI)."
struct SerialSpace <: AbstractSpace
    n::Int64
end

"Distributed over the ranks of `comm`; `n` own degrees of freedom on this rank."
struct MPISpace <: AbstractSpace
    comm::MPI.Comm
    n::Int64
end

n_own(space::AbstractSpace) = space.n

"""
    inner(space, x, y)

Inner product over the own degrees of freedom of all ranks.
"""
function inner(space::SerialSpace, x::AbstractVector, y::AbstractVector)
    return dot(view(x, 1:(space.n)), view(y, 1:(space.n)))
end
function inner(space::MPISpace, x::AbstractVector, y::AbstractVector)
    local_value = dot(view(x, 1:(space.n)), view(y, 1:(space.n)))
    return MPI.Allreduce(local_value, +, space.comm)
end

"""
    DistributedOperator(op, halo; accumulate = false)

Distributed application of a rank's local operator (README A4.0b, one-horizon ghosts):

1. `exchange!`: ghost entries of `x` from their owners (one horizon),
2. local `mul!` over all local rows (own rows and ghost rows),
3. `accumulate_ghosts!` (only if `accumulate`): partial sums of the ghost rows are sent to the
   owners and added to their own rows.

Step 3 is needed when elements write into rows of their neighbours (correspondence,
ordinary state-based, IGA-PD). For bond-based elements the ghost rows are empty and step 3
is skipped (`needs_accumulate`). With `halo = nothing` it is the serial operator.
"""
struct DistributedOperator{O,H}
    op::O
    halo::H
    accumulate::Bool
end
function DistributedOperator(op, halo; accumulate::Bool = false)
    DistributedOperator(op, halo,
                        accumulate)
end

function LinearAlgebra.mul!(y::AbstractVector, A::DistributedOperator, x::AbstractVector)
    exchange!(A.halo, x)
    mul!(y, A.op, x)
    A.accumulate && accumulate_ghosts!(A.halo, y)
    return y
end

"""
    needs_accumulate(comm, rowptr, n_own)

`true` if any rank has entries in ghost rows (rows `n_own+1:end` of its local matrix), i.e.
if partial row sums have to be sent back to the owners.
"""
function needs_accumulate(comm::MPI.Comm, rowptr::AbstractVector, n_own::Int64)
    host = Array(rowptr)
    has_ghost_rows = host[end] > host[n_own + 1]
    return MPI.Allreduce(has_ghost_rows, |, comm)
end

end
