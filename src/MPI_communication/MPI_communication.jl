# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
module MPI_Communication
import MPI
import ..Data_Manager: SynchBuffer, fill_send!, flush_recv!
export send_single_value_from_vector
export synch_responder_to_controller!
export synch_controller_to_responder!
export synch_controller_bonds_to_responder
export split_vector
export synch_controller_bonds_to_responder_flattened!
export send_vector_from_root_to_core_i
export broadcast_value
export find_and_set_core_value_min
export find_and_set_core_value_max
export find_and_set_core_value_sum
export find_and_set_core_value_avg
export gather_values
export barrier

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
    if type == String
        @error "Wrong type - String in function send_single_value_from_vector"
        return nothing
    end
    recv_msg = zeros(type, 1, 1)
    if rank == controller
        send_msg = zeros(type, 1, 1)
        for i in 0:(ncores - 1)
            send_msg[1] = values[i + 1]
            if i != controller
                MPI.Isend(send_msg, comm; dest = i, tag = 0)
            else
                recv_msg[1] = send_msg[1]
            end
        end
    else
        MPI.Recv!(recv_msg, comm; source = controller, tag = 0)
    end
    return recv_msg[1]
end

"""
    synch_responder_to_controller!(comm::MPI.Comm,  vector, buf::SynchBuffer{T,N}) where {T,N}

Synchronise responder nodes to controller nodes using pre-allocated buffers.
Responder values are summed into the controller nodes.

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `vector`: Field array to synchronise (modified in-place)
- `buf::SynchBuffer`: Pre-allocated send/receive buffers matching `vector`
"""
function synch_responder_to_controller!(comm::MPI.Comm, vector::Array{T,N},
                                        buf::SynchBuffer{T,N}) where {T,N}
    ncores = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    ncores == 1 && return

    for jcore in 1:ncores
        rank + 1 == jcore && continue
        key = (rank + 1, jcore)

        # Responder sends its values to the controller
        if haskey(buf.send_bufs, key)
            fill_send!(buf, key, vector)
            MPI.Send(buf.send_bufs[key], comm; dest = jcore - 1, tag = 0)
        end

        # Controller receives and accumulates responder values
        if haskey(buf.recv_bufs, key)
            MPI.Recv!(buf.recv_bufs[key], comm; source = jcore - 1, tag = 0)
            if N == 1
                vector[buf.recv_indices[key]] .+= buf.recv_bufs[key]
            else
                vector[buf.recv_indices[key], :] .+= buf.recv_bufs[key]
            end
        end
    end
end

"""
    synch_controller_to_responder!(comm::MPI.Comm, vector, buf::SynchBuffer{T,N}) where {T,N}

Synchronise controller nodes to responder nodes using pre-allocated buffers.

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `overlapnodes`: Overlap map defining controller/responder node indices per core pair
- `vector`: Field array to synchronise (modified in-place)
- `buf::SynchBuffer`: Pre-allocated send/receive buffers matching `vector`
"""
function synch_controller_to_responder!(comm::MPI.Comm, vector::Array{T,N},
                                        buf::SynchBuffer{T,N}) where {T,N}
    ncores = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    ncores == 1 && return

    for jcore in 1:ncores
        rank + 1 == jcore && continue
        key = (rank + 1, jcore)

        if haskey(buf.send_bufs, key)
            fill_send!(buf, key, vector)
            MPI.Send(buf.send_bufs[key], comm; dest = jcore - 1, tag = 0)
        end

        if haskey(buf.recv_bufs, key)
            MPI.Recv!(buf.recv_bufs[key], comm; source = jcore - 1, tag = 0)
            flush_recv!(buf, key, vector)
        end
    end
end

"""
    synch_controller_bonds_to_responder(comm::MPI.Comm, overlapnodes, array, dof)

Synch the controller bonds to the responder

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `overlapnodes::Dict`: The overlap nodes
- `array::Array`: The array
- `dof::Int`: The degree of freedom
# Returns
- `array::Array`: The array
"""
function synch_controller_bonds_to_responder(comm::MPI.Comm, overlapnodes, array, dof)
    ncores = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    ncores == 1 && return array

    for jcore in 1:ncores
        rank + 1 == jcore && continue

        if !isempty(overlapnodes[rank + 1][jcore]["Controller"])
            for iID in overlapnodes[rank + 1][jcore]["Controller"]
                if dof == 1
                    @views send_msg = array[iID]
                else
                    @views send_msg = mapreduce(permutedims, vcat, array[iID])
                end
                MPI.Isend(send_msg, comm; dest = jcore - 1, tag = 0)
            end
        end

        if !isempty(overlapnodes[rank + 1][jcore]["Responder"])
            for iID in overlapnodes[rank + 1][jcore]["Responder"]
                if dof == 1
                    @views recv_msg = similar(array[iID])
                else
                    @views recv_msg = similar(mapreduce(permutedims, vcat, array[iID]))
                end
                MPI.Recv!(recv_msg, comm; source = jcore - 1, tag = 0)
                recv_msg = reshape(recv_msg, :, 2)
                if dof == 1
                    array[iID] = recv_msg
                else
                    array[iID] = Vector{eltype(recv_msg)}[eachcol(recv_msg)...]
                end
            end
        end
    end
    return array
end

"""
    split_vector(input, row_nums, dof)

Split a vector into a vector of matrices

# Arguments
- `input::Vector`: The input vector
- `row_nums::Vector`: The row numbers
- `dof::Int`: The degree of freedom
# Returns
- `result::Vector`: The result vector
"""
function split_vector(input, row_nums, dof)
    result = Vector{Vector{Vector{eltype(input)}}}()
    start = firstindex(input)
    for (i, len) in enumerate(row_nums)
        push!(result, [input[(start + (j - 1) * dof):(start - 1 + j * dof)] for j in 1:len])
        start += dof * len
    end
    return result
end

"""
    synch_controller_bonds_to_responder_flattened!(comm::MPI.Comm, overlapnodes, array, dof)

Synchronise controller bond data to responder nodes in flattened form.

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `overlapnodes`: Overlap map defining controller/responder node indices per core pair
- `array`: Bond-level field array to synchronise (modified in-place)
- `dof::Int`: The degree of freedom
# Returns
- `array`: The updated array
"""
function synch_controller_bonds_to_responder_flattened!(comm::MPI.Comm, overlapnodes,
                                                        array, dof)
    ncores = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    ncores == 1 && return array

    for jcore in 1:ncores
        rank + 1 == jcore && continue

        if !isempty(overlapnodes[rank + 1][jcore]["Controller"])
            @views send_indices = overlapnodes[rank + 1][jcore]["Controller"]
            MPI.Send(vcat(vcat(array...)[send_indices]...), comm; dest = jcore - 1, tag = 0)
        end

        if !isempty(overlapnodes[rank + 1][jcore]["Responder"])
            @views recv_indices = overlapnodes[rank + 1][jcore]["Responder"]
            row_nums = [length(subarr) for subarr in array[recv_indices]]
            recv_msg = zeros(sum(row_nums * dof))
            MPI.Recv!(recv_msg, comm; source = jcore - 1, tag = 0)
            recv_msg = split_vector(recv_msg, row_nums, dof)
            array[recv_indices] .= recv_msg
        end
    end
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
    MPI.Barrier(comm)
    if currentRank == 0
        for rank in 1:(MPI.Comm_size(comm) - 1)
            MPI.Isend(send_msg[distribution[rank + 1]], comm; dest = rank, tag = 0)
        end
        recv_msg .= send_msg[distribution[1]]
    else
        MPI.Recv!(recv_msg, comm; source = 0, tag = 0)
    end
    return recv_msg
end

"""
    broadcast_value(comm::MPI.Comm, send_msg)

Broadcast a value to all ranks.

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `send_msg`: The value to broadcast
# Returns
- The broadcasted value on all ranks
"""
function broadcast_value(comm::MPI.Comm, send_msg::T) where {T}
    return MPI.bcast(send_msg, 0, comm)
end

"""
    find_and_set_core_value_sum(comm::MPI.Comm, value)

Find and set core value sum.

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `value`: The value
# Returns
- The global sum across all ranks
"""
function find_and_set_core_value_sum(comm::MPI.Comm,
                                     value::T) where {T<:Union{Float64,Int64,
                                                               Vector{Float64},
                                                               Vector{Int64},
                                                               Matrix{Float64},
                                                               Matrix{Int64}}}
    return MPI.Allreduce(value, MPI.SUM, comm)
end

"""
    find_and_set_core_value_max(comm::MPI.Comm, value::Union{Float64,Int64})

Find and set core value max.

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `value::Union{Float64,Int64}`: The value
# Returns
- The global maximum across all ranks
"""
function find_and_set_core_value_max(comm::MPI.Comm,
                                     value::T) where {T<:Union{Float64,Int64}}
    return MPI.Allreduce(value, MPI.MAX, comm)
end

"""
    find_and_set_core_value_min(comm::MPI.Comm, value)

Find and set core value min.

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `value`: The value
# Returns
- The global minimum across all ranks
"""
function find_and_set_core_value_min(comm::MPI.Comm,
                                     value::T) where {T<:Union{Float64,Int64,
                                                               Vector{Float64},
                                                               Vector{Int64},
                                                               Matrix{Float64},
                                                               Matrix{Int64}}}
    return MPI.Allreduce(value, MPI.MIN, comm)
end

"""
    find_and_set_core_value_avg(comm::MPI.Comm, value, nnodes::Int64)

Find and set core value avg.

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `value::Union{Float64,Int64}`: The value
- `nnodes::Int64`: The total number of nodes
# Returns
- `Float64`: The global average across all ranks
"""
function find_and_set_core_value_avg(comm::MPI.Comm,
                                     value::T,
                                     nnodes::Int64) where {T<:Union{Float64,Int64}}
    return MPI.Allreduce(value, MPI.SUM, comm) / MPI.Allreduce(nnodes, MPI.SUM, comm)
end

"""
    gather_values(comm::MPI.Comm, value)

Gather values from all ranks to root.

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `value`: The value to gather
# Returns
- Gathered values on root rank, `nothing` on other ranks
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

"""
    barrier(comm::MPI.Comm)

Wait until all MPI ranks have reached this point.

# Arguments
- `comm::MPI.Comm`: The MPI communicator
"""
function barrier(comm::MPI.Comm)
    MPI.Barrier(comm)
end

end
