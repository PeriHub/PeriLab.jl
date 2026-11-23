# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
module MPI_Communication
import MPI
using ArraysOfArrays: flatview
using StaticArrays: StaticArray
using ...PeriLabExceptions: @abort
export send_single_value_from_vector
export synch_responder_to_controller
export synch_controller_to_responder
export synch_controller_bonds_to_responder_flattened
export send_vector_from_root_to_core_i
export send_rows_from_root_to_core_i
export broadcast_value
export find_and_set_core_value_min
export find_and_set_core_value_sum
export find_and_set_core_value_avg
export gather_values
# export barrier

"""
TODO
Contact
send all information to first core and synch to all otherwise
optimization is possible by reducing it to slave and master. Therefore its only the surface.

Master are known and their core. local to global is know and the sending can occur.

"""

"""
    send_single_value_from_vector(comm::MPI.Comm, controller::Int64, values::Union{Int64,Vector{Float64},Vector{Int64},Vector{Bool}}, type::Type)

Sends a single value from a vector to a controller

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `controller::Int64`: The controller
- `values::Union{Int64,Vector{Float64},Vector{Int64},Vector{Bool}}`: The values
- `type::Type`: The type
# Returns
- `recv_msg::Union{Int64,Vector{Float64},Vector{Int64},Vector{Bool}}`: The received message
"""
function send_single_value_from_vector(comm::MPI.Comm,
                                       controller::Int64,
                                       values::Union{Int64,Vector{Float64},Vector{Int64},
                                                     Vector{Bool}},
                                       type::Type)
    ncores = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    requests = Vector{MPI.Request}()
    if type == String
        @abort "Wrong type - String in function send_single_value_from_vector"
        return nothing
    end
    recv_msg = zeros(type, 1, 1)
    if rank == controller
        send_msg = zeros(type, 1, 1)
        for i in 0:(ncores - 1)
            # +1 because the index of cores is zero based and julia matrices are one based
            send_msg[1] = values[i + 1]
            if i != controller
                push!(requests, MPI.Isend(send_msg, comm; dest = i, tag = 0))
                # @debug "Sending   $rank -> $i"
            else
                recv_msg[1] = send_msg[1]
            end
        end

    else
        MPI.Recv!(recv_msg, comm; source = controller, tag = 0)
        # @debug "Receiving $controller -> $rank"
    end
    if !isempty(requests)
        MPI.Waitall(requests)
    end
    MPI.Barrier(comm)
    return recv_msg[1]
end

"""
function synch_overlapnodes(comm::MPI.Comm, topo, vector)
    currentRank = MPI.Comm_rank(comm)
    ncores = MPI.Comm_size(comm)
    overlapCurrentRank = topo[currentRank+1]
    for icore in 0:ncores-1
        if icore == currentRank
            continue
        end
        if overlapCurrentRank[icore+1]["Responder"] > 0
            send_msg = vector[overlapCurrentRank[icore+1]["Responder"]]
            MPI.Send(send_msg, comm; dest=icore+i, tag=0)
        end
        if overlapCurrentRank[icore+1]["Controller"] > 0
            recv_msg = vector[overlapCurrentRank[icore+1]["Controller"]]
            MPI.Recv!(recv_msg, comm; source=0, tag=0)
            vector[overlapCurrentRank[icore+1]["Controller"]]
        end
    end
    return vector

end
"""

"""
    SynchBuffers{T}

Send and receive buffers of the node and bond synchronisation, one per neighbour rank, and
their requests. Created once per element type `T` and reused by every synchronisation:
`resize!` within the capacity of a buffer does not allocate, so after the first steps no
buffer is created any more.
"""
struct SynchBuffers{T}
    send::Vector{Vector{T}}
    recv::Vector{Vector{T}}
    send_requests::MPI.MultiRequest
    recv_requests::MPI.MultiRequest
end

const SYNCH_BUFFERS = Dict{DataType,Any}()

"""
    synch_buffers(::Type{T}, ncores::Int64)

The `SynchBuffers{T}` of this process for `ncores` ranks (created on first use).
"""
function synch_buffers(::Type{T}, ncores::Int64)::SynchBuffers{T} where {T}
    buffers = get(SYNCH_BUFFERS, T, nothing)
    buffers isa SynchBuffers{T} && length(buffers.send) == ncores && return buffers
    buffers = SynchBuffers{T}([T[] for _ in 1:ncores], [T[] for _ in 1:ncores],
                              MPI.MultiRequest(ncores), MPI.MultiRequest(ncores))
    SYNCH_BUFFERS[T] = buffers
    return buffers
end

# values per row: all trailing dimensions (vector: 1, matrix: columns, tensor field: d·d)
columns(vector::AbstractArray) = length(vector) ÷ max(size(vector, 1), 1)

# rows `index` of `vector` (all values of a row, column major) into `buffer`
function pack!(buffer::Vector{T}, vector::AbstractArray{T},
               index::Vector{Int64}) where {T}
    n = length(index)
    rows = size(vector, 1)
    resize!(buffer, n * columns(vector))
    @inbounds for c in 1:columns(vector), (k, i) in enumerate(index)
        buffer[(c - 1) * n + k] = vector[(c - 1) * rows + i]
    end
    return buffer
end

# `buffer` (layout of `pack!`) into the rows `index` of `vector`; `add` sums instead
function unpack!(vector::AbstractArray{T}, index::Vector{Int64}, buffer::Vector{T},
                 add::Bool) where {T}
    n = length(index)
    rows = size(vector, 1)
    @inbounds for c in 1:columns(vector), (k, i) in enumerate(index)
        j = (c - 1) * rows + i
        vector[j] = add ? vector[j] + buffer[(c - 1) * n + k] : buffer[(c - 1) * n + k]
    end
    return vector
end

"""
    synch_responder_to_controller(comm::MPI.Comm, overlapnodes, vector, dof)

Sends the values of the responder rows to the ranks that control these nodes and adds them
there (Bool fields are only exchanged). Uses the reused `SynchBuffers`; no buffer is
allocated per call.

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `overlapnodes::Dict`: The overlap nodes
- `vector::Vector`: The vector
- `dof::Int`: The degree of freedom
# Returns
- `vector::Vector`: The vector
"""
function synch_responder_to_controller(comm::MPI.Comm, overlapnodes, vector, dof)
    ncores = MPI.Comm_size(comm)
    ncores == 1 && return vector
    rank = MPI.Comm_rank(comm)
    maps = overlapnodes[rank + 1]
    buffers = synch_buffers(eltype(vector), ncores)
    # all receives are posted before the sends, so no rank waits for the loop order of
    # another; the responder (sent) and controller (received) rows are disjoint
    for jcore in 1:ncores
        jcore == rank + 1 && continue
        recv_index = maps[jcore]["Controller"]
        isempty(recv_index) && continue
        buffer = resize!(buffers.recv[jcore], length(recv_index) * columns(vector))
        MPI.Irecv!(buffer, comm, buffers.recv_requests[jcore]; source = jcore - 1, tag = 0)
    end
    for jcore in 1:ncores
        jcore == rank + 1 && continue
        send_index = maps[jcore]["Responder"]
        isempty(send_index) && continue
        buffer = pack!(buffers.send[jcore], vector, send_index)
        MPI.Isend(buffer, comm, buffers.send_requests[jcore]; dest = jcore - 1, tag = 0)
    end
    MPI.Waitall(buffers.recv_requests)
    # summed in the rank order, as before; Bool fields are only exchanged
    if !(eltype(vector) <: Bool)
        for jcore in 1:ncores
            jcore == rank + 1 && continue
            recv_index = maps[jcore]["Controller"]
            isempty(recv_index) && continue
            unpack!(vector, recv_index, buffers.recv[jcore], true)
        end
    end
    MPI.Waitall(buffers.send_requests)
    return vector
end

"""
    synch_controller_to_responder(comm::MPI.Comm, overlapnodes, vector, dof)

Sends the values of the controller rows to the ranks where these nodes are responders.
Uses the reused `SynchBuffers`; no buffer is allocated per call.

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `overlapnodes::Dict`: The overlap nodes
- `vector::Vector`: The vector
- `dof::Int`: The degree of freedom
# Returns
- `vector::Vector`: The vector
"""
function synch_controller_to_responder(comm::MPI.Comm, overlapnodes, vector, dof)
    ncores = MPI.Comm_size(comm)
    ncores == 1 && return vector
    rank = MPI.Comm_rank(comm)
    maps = overlapnodes[rank + 1]
    buffers = synch_buffers(eltype(vector), ncores)
    # receives first, then sends (see `synch_responder_to_controller`)
    for jcore in 1:ncores
        jcore == rank + 1 && continue
        recv_index = maps[jcore]["Responder"]
        isempty(recv_index) && continue
        buffer = resize!(buffers.recv[jcore], length(recv_index) * columns(vector))
        MPI.Irecv!(buffer, comm, buffers.recv_requests[jcore]; source = jcore - 1, tag = 0)
    end
    for jcore in 1:ncores
        jcore == rank + 1 && continue
        send_index = maps[jcore]["Controller"]
        isempty(send_index) && continue
        buffer = pack!(buffers.send[jcore], vector, send_index)
        MPI.Isend(buffer, comm, buffers.send_requests[jcore]; dest = jcore - 1, tag = 0)
    end
    MPI.Waitall(buffers.recv_requests)
    for jcore in 1:ncores
        jcore == rank + 1 && continue
        recv_index = maps[jcore]["Responder"]
        isempty(recv_index) && continue
        unpack!(vector, recv_index, buffers.recv[jcore], false)
    end
    MPI.Waitall(buffers.send_requests)
    return vector
end

# scalar values of bond data (also of SVector / SMatrix bonds) in bond order, as MPI buffer
bond_values_view(bonds::AbstractVector{<:Number}) = bonds
function bond_values_view(bonds::AbstractVector{<:StaticArray{S,T}}) where {S,T}
    reinterpret(T, bonds)
end

"""
    synch_controller_bonds_to_responder_flattened(comm::MPI.Comm, overlapnodes, array)

Sends the bond values of the controller nodes to the ranks where they are responders. A
responder row only holds the bonds known on its rank; the flat bond indices of both sides
are the `"Bond Controller"` and `"Bond Responder"` entries of the overlap map (see
`create_bond_overlap_map!`). Uses the reused `SynchBuffers`; no buffer is allocated per
call.

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `overlapnodes::Dict`: The overlap nodes, including the bond indices
- `array`: The bond field
# Returns
- `array`: The bond field
"""
function synch_controller_bonds_to_responder_flattened(comm::MPI.Comm,
                                                       overlapnodes,
                                                       array)
    ncores = MPI.Comm_size(comm)
    ncores == 1 && return array
    rank = MPI.Comm_rank(comm)
    maps = overlapnodes[rank + 1]
    flat = flatview(array)
    isempty(flat) && return array
    # all scalar values of the bonds; bond b holds values (b - 1) * per_bond .+ (1:per_bond)
    values = bond_values_view(flat)
    per_bond = length(values) ÷ length(flat)
    buffers = synch_buffers(eltype(values), ncores)
    # receives first, then sends (see `synch_responder_to_controller`)
    for jcore in 1:ncores
        jcore == rank + 1 && continue
        recv_bonds = maps[jcore]["Bond Responder"]
        isempty(recv_bonds) && continue
        buffer = resize!(buffers.recv[jcore], length(recv_bonds) * per_bond)
        MPI.Irecv!(buffer, comm, buffers.recv_requests[jcore]; source = jcore - 1, tag = 0)
    end
    for jcore in 1:ncores
        jcore == rank + 1 && continue
        send_bonds = maps[jcore]["Bond Controller"]
        isempty(send_bonds) && continue
        buffer = resize!(buffers.send[jcore], length(send_bonds) * per_bond)
        @inbounds for (k, b) in enumerate(send_bonds), l in 1:per_bond
            buffer[(k - 1) * per_bond + l] = values[(b - 1) * per_bond + l]
        end
        MPI.Isend(buffer, comm, buffers.send_requests[jcore]; dest = jcore - 1, tag = 0)
    end
    MPI.Waitall(buffers.recv_requests)
    for jcore in 1:ncores
        jcore == rank + 1 && continue
        recv_bonds = maps[jcore]["Bond Responder"]
        isempty(recv_bonds) && continue
        buffer = buffers.recv[jcore]
        @inbounds for (k, b) in enumerate(recv_bonds), l in 1:per_bond
            values[(b - 1) * per_bond + l] = buffer[(k - 1) * per_bond + l]
        end
    end
    MPI.Waitall(buffers.send_requests)
    return array
end

"""
    send_vector_from_root_to_core_i(comm::MPI.Comm, send_msg, recv_msg, distribution)

Sends a vector from the root to the core i

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `send_msg::Union{Int64,Vector{Float64},Vector{Int64},Vector{Bool}}`: The send message
- `recv_msg::Union{Int64,Vector{Float64},Vector{Int64},Vector{Bool}}`: The receive message
- `distribution::Vector{Int64}`: The distribution
# Returns
- `recv_msg::Union{Int64,Vector{Float64},Vector{Int64},Vector{Bool}}`: The received message
"""
function send_vector_from_root_to_core_i(comm::MPI.Comm, send_msg, recv_msg, distribution)
    currentRank = MPI.Comm_rank(comm)
    # MPI.Barrier(comm)
    requests = Vector{MPI.Request}()
    if currentRank == 0
        for rank in 1:(MPI.Comm_size(comm) - 1)
            push!(requests,
                  MPI.Isend(send_msg[distribution[rank + 1]], comm; dest = rank, tag = 0))
            # @debug "Sending $currentRank -> $rank"
        end
        recv_msg .= send_msg[distribution[1]]
    else
        MPI.Recv!(recv_msg, comm; source = 0, tag = 0)
        # @debug "Receiving 0 -> $currentRank"
    end
    if !isempty(requests)
        MPI.Waitall(requests)
    end
    MPI.Barrier(comm)
    return recv_msg
end

"""
    send_rows_from_root_to_core_i(comm::MPI.Comm, rows, distribution)

Sends every rank only its rows `rows[distribution[rank]]` of a nested list held by the root,
flattened into two messages (row lengths and values) instead of a serialized broadcast.

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `rows::Vector{Vector{Int64}}`: All rows (root), ignored on the other ranks
- `distribution::Vector{Vector{Int64}}`: global nodes of all ranks
# Returns
- `Vector{Vector{Int64}}`: the rows of this rank
"""
function send_rows_from_root_to_core_i(comm::MPI.Comm, rows,
                                       distribution::Vector{Vector{Int64}})
    rank = MPI.Comm_rank(comm)
    if rank == 0
        requests = MPI.Request[]
        for jrank in 1:(MPI.Comm_size(comm) - 1)
            local_rows = view(rows, distribution[jrank + 1])
            lengths = length.(local_rows)
            values = Vector{Int64}(undef, sum(lengths; init = 0))
            pos = 1
            for row in local_rows
                copyto!(values, pos, row, 1, length(row))
                pos += length(row)
            end
            push!(requests, MPI.Isend(lengths, comm; dest = jrank, tag = 1))
            push!(requests, MPI.Isend(values, comm; dest = jrank, tag = 2))
        end
        MPI.Waitall(requests)
        return rows[distribution[1]]
    end
    lengths = Vector{Int64}(undef, length(distribution[rank + 1]))
    MPI.Recv!(lengths, comm; source = 0, tag = 1)
    values = Vector{Int64}(undef, sum(lengths; init = 0))
    MPI.Recv!(values, comm; source = 0, tag = 2)
    local_rows = Vector{Vector{Int64}}(undef, length(lengths))
    pos = 1
    for (i, n) in enumerate(lengths)
        local_rows[i] = values[pos:(pos + n - 1)]
        pos += n
    end
    return local_rows
end

"""
    broadcast_value(comm::MPI.Comm, send_msg)

Broadcast a value to all ranks

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `send_msg::Union{Int64,Vector{Float64},Vector{Int64},Vector{Bool}}`: The send message
# Returns
- `recv_msg::Union{Int64,Vector{Float64},Vector{Int64},Vector{Bool}}`: The received message
"""
function broadcast_value(comm::MPI.Comm,
                         send_msg::T) where {T<:Union{Float64,Int64,
                                                      Vector{Float64},
                                                      Vector{Vector{Float64}},
                                                      Vector{Int64},Vector{Vector{Int64}},
                                                      Matrix{Float64},
                                                      Matrix{Int64},
                                                      Dict,
                                                      Bool,
                                                      Nothing,
                                                      Any}}
    # recv_msg = MPI.Comm_rank(comm) == controller ? send_msg : nothing
    # recv_msg = MPI.bcast(send_msg, controller, comm)
    # return recv_msg
    return MPI.bcast(send_msg, 0, comm)
end

"""
    find_and_set_core_value_sum(comm::MPI.Comm, value::Union{Float64,Int64})

Find and set core value sum

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `value::Union{Float64,Int64}`: The value
# Returns
- `recv_msg::Union{Int64,Vector{Float64},Vector{Int64},Vector{Bool}}`: The received message
"""
function find_and_set_core_value_sum(comm::MPI.Comm,
                                     value::T) where {T<:Union{Float64,Int64}}
    return MPI.Allreduce(value, MPI.SUM, comm)
end
# arrays are summed in place (no new array per call)
function find_and_set_core_value_sum(comm::MPI.Comm,
                                     value::T) where {T<:Union{Vector{Float64},
                                                               Vector{Int64},
                                                               Matrix{Float64},
                                                               Matrix{Int64}}}
    return MPI.Allreduce!(value, MPI.SUM, comm)
end

"""
    find_and_set_core_value_max(comm::MPI.Comm, value::Union{Float64,Int64})

Find and set core value max

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `value::Union{Float64,Int64}`: The value
# Returns
- `recv_msg::Union{Int64,Vector{Float64},Vector{Int64},Vector{Bool}}`: The received message
"""
function find_and_set_core_value_max(comm::MPI.Comm,
                                     value::T) where {T<:Union{Float64,Int64}}
    return MPI.Allreduce(value, MPI.MAX, comm)
end

"""
    find_and_set_core_value_min(comm::MPI.Comm, value::Union{Float64,Int64})

Find and set core value sum

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `value::Union{Float64,Int64}`: The value
# Returns
- `recv_msg::Union{Int64,Vector{Float64},Vector{Int64},Vector{Bool}}`: The received message
"""
function find_and_set_core_value_min(comm::MPI.Comm,
                                     value::T) where {T<:Union{Float64,Int64}}
    return MPI.Allreduce(value, MPI.MIN, comm)
end
# arrays are reduced in place (no new array per call)
function find_and_set_core_value_min(comm::MPI.Comm,
                                     value::T) where {T<:Union{Vector{Float64},
                                                               Vector{Int64},
                                                               Matrix{Float64},
                                                               Matrix{Int64}}}
    return MPI.Allreduce!(value, MPI.MIN, comm)
end

"""
    find_and_set_core_value_avg(comm::MPI.Comm,
                                     value::T,
                                     nnodes::Int64) where {T<:Union{Float64,Int64}}

Find and set core value avg

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `value::Union{Float64,Int64}`: The value
# Returns
- `recv_msg::Float64`: The received a Float64 message
"""
function find_and_set_core_value_avg(comm::MPI.Comm,
                                     value::T,
                                     nnodes::Int64) where {T<:Union{Float64,Int64}}
    return MPI.Allreduce(value, MPI.SUM, comm) / MPI.Allreduce(nnodes, MPI.SUM, comm)
end

"""
    gather_values(comm::MPI.Comm, value::Any)

Gather values

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `value::Any`: The value
# Returns
- `recv_msg::Any`: The received message
"""
function gather_values(comm::MPI.Comm,
                       value::T) where {T<:Union{Float64,Int64,
                                                 Vector{Float64},
                                                 Vector{Int64},
                                                 Matrix{Float64},
                                                 Matrix{Int64},
                                                 Any}}
    return MPI.gather(value, comm; root = 0)
end

# function barrier(comm::MPI.Comm)
#     MPI.Barrier(comm)
# end

end
