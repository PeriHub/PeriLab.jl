# SPDX-FileCopyrightText: 2026 Christian Willberg
#
# SPDX-License-Identifier: BSD-3-Clause

module Halo_Exchange

using MPI
using KernelAbstractions: @kernel, @index, @Const
using TimerOutputs: @timeit
using ..Sparse_CSR: launch!

export HaloExchange, create_halo_exchange, exchange!, accumulate_ghosts!

"""
    HaloExchange

Communication pattern for the ghost degrees of freedom of one rank (README A4.2, A0.3a).

Local numbering of a rank: own degrees of freedom `1:n_own`, then the ghost degrees of
freedom `n_own+1:n_own+n_ghost` in the caller's order (PeriLab: neighbour-based, not
grouped by owner). Ghost values arrive per owner rank in a receive buffer and are
scattered to their local positions with `recv_idx`. `exchange!(halo, x)` fills the ghost
part of `x` with the current values of the owners. All buffers are allocated once.

Without CUDA-aware MPI the device buffers are copied through preallocated host buffers
(README F6); with `cuda_aware = true` MPI reads and writes the device buffers directly.
"""
struct HaloExchange{Vi<:AbstractVector{Int64},Vf<:AbstractVector{Float64}}
    comm::MPI.Comm
    n_own::Int64
    n_ghost::Int64
    send_ranks::Vector{Int64}
    send_ptr::Vector{Int64}
    recv_ranks::Vector{Int64}
    recv_ptr::Vector{Int64}
    send_idx::Vi
    send_buf::Vf
    recv_idx::Vi
    recv_buf::Vf
    send_host::Vector{Float64}
    recv_host::Vector{Float64}
    staged::Bool
    requests::MPI.MultiRequest
end

"""
    create_halo_exchange(comm, own_points, ghost_points, ghost_owner, dof; array_type, cuda_aware)

Builds the exchange pattern from the point distribution delivered by the caller (as in
PeriLab: own points, ghost points and their owner ranks). The ranks agree on the send lists
with one `Alltoall` and one `Alltoallv` (setup only).

# Arguments
- `comm::MPI.Comm`: Communicator
- `own_points::AbstractVector{Int64}`: Global ids of the own points (local order)
- `ghost_points::AbstractVector{Int64}`: Global ids of the ghost points (local order, any
  owner order)
- `ghost_owner::AbstractVector{Int64}`: Owner rank of each ghost point
- `dof::Int64`: Degrees of freedom per point
- `array_type`: `Array` or `CuArray`, backend of the vectors that are exchanged
- `cuda_aware::Bool`: MPI can read device memory directly
"""
function create_halo_exchange(comm::MPI.Comm, own_points::AbstractVector{Int64},
                              ghost_points::AbstractVector{Int64},
                              ghost_owner::AbstractVector{Int64}, dof::Int64;
                              array_type = Array, cuda_aware::Bool = false)
    n_ranks = MPI.Comm_size(comm)
    if length(ghost_owner) != length(ghost_points)
        throw(ArgumentError("ghost_owner needs one entry per ghost point."))
    end
    # receive order: grouped by owner, local order within an owner
    order = sortperm(ghost_owner; alg = Base.Sort.DEFAULT_STABLE)
    recv_count = zeros(Cint, n_ranks)
    for o in ghost_owner
        recv_count[o + 1] += 1
    end
    send_count = MPI.Alltoall(MPI.UBuffer(recv_count, 1), comm)
    requested = Vector{Int64}(undef, sum(send_count))
    MPI.Alltoallv!(MPI.VBuffer(collect(Int64, ghost_points[order]), recv_count),
                   MPI.VBuffer(requested, send_count), comm)

    local_of = Dict(g => i for (i, g) in enumerate(own_points))
    send_points = [local_of[g] for g in requested]
    send_idx = [(i - 1) * dof + d for i in send_points for d in 1:dof]
    n_own_dofs = dof * length(own_points)
    recv_idx = [n_own_dofs + (i - 1) * dof + d for i in order for d in 1:dof]

    send_ranks = [r - 1 for r in 1:n_ranks if send_count[r] > 0]
    recv_ranks = [r - 1 for r in 1:n_ranks if recv_count[r] > 0]
    send_ptr = cumsum([1; dof .* send_count[send_ranks .+ 1]])
    recv_ptr = cumsum([1; dof .* recv_count[recv_ranks .+ 1]])

    n_send = length(send_idx)
    n_ghost = dof * length(ghost_points)
    staged = !(array_type <: Array) && !cuda_aware
    return HaloExchange(comm, n_own_dofs, n_ghost, send_ranks, send_ptr, recv_ranks,
                        recv_ptr, array_type(send_idx), array_type(zeros(n_send)),
                        array_type(recv_idx), array_type(zeros(n_ghost)),
                        zeros(staged ? n_send : 0), zeros(staged ? n_ghost : 0), staged,
                        MPI.MultiRequest(length(send_ranks) + length(recv_ranks)))
end

@kernel function pack_kernel!(buf, @Const(idx), @Const(x))
    k = @index(Global)
    @inbounds buf[k] = x[idx[k]]
end

@kernel function scatter_kernel!(x, @Const(idx), @Const(buf))
    k = @index(Global)
    @inbounds x[idx[k]] = buf[k]
end

"""
    exchange!(halo, x)

Writes the owners' values into the ghost part of `x` (`x[n_own+1:end]`). No allocation
of vectors; works for host and device vectors.
"""
function exchange!(halo::HaloExchange, x::AbstractVector)
    @timeit "halo exchange" begin
        n_send = length(halo.send_idx)
        launch!(pack_kernel!, halo.send_buf, n_send, halo.send_buf, halo.send_idx, x)
        halo.staged && n_send > 0 && copyto!(halo.send_host, 1, halo.send_buf, 1, n_send)
        send_base = halo.staged ? halo.send_host : halo.send_buf
        recv_base = halo.staged ? halo.recv_host : halo.recv_buf
        n_recv = length(halo.recv_ranks)
        for (k, rank) in enumerate(halo.recv_ranks)
            range = halo.recv_ptr[k]:(halo.recv_ptr[k + 1] - 1)
            MPI.Irecv!(view(recv_base, range), halo.comm, halo.requests[k]; source = rank,
                       tag = 0)
        end
        for (k, rank) in enumerate(halo.send_ranks)
            range = halo.send_ptr[k]:(halo.send_ptr[k + 1] - 1)
            MPI.Isend(view(send_base, range), halo.comm, halo.requests[n_recv + k];
                      dest = rank, tag = 0)
        end
        MPI.Waitall(halo.requests)
        halo.staged && halo.n_ghost > 0 &&
            copyto!(halo.recv_buf, 1, halo.recv_host, 1, halo.n_ghost)
        launch!(scatter_kernel!, x, halo.n_ghost, x, halo.recv_idx, halo.recv_buf)
    end
    return x
end

exchange!(::Nothing, x::AbstractVector) = x

@kernel function add_kernel!(y, @Const(idx), @Const(buf), offset::Int64)
    k = @index(Global)
    @inbounds y[idx[offset + k]] += buf[offset + k]
end

"""
    accumulate_ghosts!(halo, y)

Reverse of `exchange!`: sends the ghost part of `y` (partial sums of rows that belong to
other ranks) to the owners and adds it to their own rows. Same pattern as PeriLab
`synch_responder_to_controller` with `add = true`. Afterwards the own part of `y` is
complete; the ghost part is not used.

The additions are done per neighbour rank in a fixed order, so the result is
reproducible and free of write conflicts (a rank requests every point only once).
"""
function accumulate_ghosts!(halo::HaloExchange, y::AbstractVector)
    @timeit "halo accumulate" begin
        launch!(pack_kernel!, halo.recv_buf, halo.n_ghost, halo.recv_buf, halo.recv_idx, y)
        halo.staged && halo.n_ghost > 0 &&
            copyto!(halo.recv_host, 1, halo.recv_buf, 1, halo.n_ghost)
        recv_base = halo.staged ? halo.recv_host : halo.recv_buf
        n_back = length(halo.send_ranks)
        for (k, rank) in enumerate(halo.send_ranks)
            range = halo.send_ptr[k]:(halo.send_ptr[k + 1] - 1)
            buf = halo.staged ? view(halo.send_host, range) : view(halo.send_buf, range)
            MPI.Irecv!(buf, halo.comm, halo.requests[k]; source = rank, tag = 1)
        end
        for (k, rank) in enumerate(halo.recv_ranks)
            range = halo.recv_ptr[k]:(halo.recv_ptr[k + 1] - 1)
            MPI.Isend(view(recv_base, range), halo.comm, halo.requests[n_back + k];
                      dest = rank, tag = 1)
        end
        MPI.Waitall(halo.requests)
        n_send = length(halo.send_idx)
        halo.staged && n_send > 0 && copyto!(halo.send_buf, 1, halo.send_host, 1, n_send)
        for k in eachindex(halo.send_ranks)
            offset = halo.send_ptr[k] - 1
            launch!(add_kernel!, y, halo.send_ptr[k + 1] - halo.send_ptr[k], y,
                    halo.send_idx, halo.send_buf, offset)
        end
    end
    return y
end

accumulate_ghosts!(::Nothing, y::AbstractVector) = y

end
