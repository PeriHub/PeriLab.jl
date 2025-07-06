# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
module MPI_communication
import MPI
export send_single_value_from_vector
export synch_responder_to_controller
export synch_controller_to_responder
export synch_controller_bonds_to_responder
export split_vector
export synch_controller_bonds_to_responder_flattened
export send_vector_from_root_to_core_i
export send_value
export find_and_set_core_value_min
export find_and_set_core_value_sum
export find_and_set_core_value_avg
export gather_values
export barrier

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
    if type == String
        @error "Wrong type - String in function send_single_value_from_vector"
        return nothing
    end
    recv_msg = zeros(type, 1, 1)
    if rank == controller
        send_msg = zeros(type, 1, 1)
        for i in 0:(ncores - 1)
            # +1 because the index of cores is zero based and julia matrices are one based
            send_msg[1] = values[i + 1]
            if i != controller
                MPI.Isend(send_msg, comm; dest = i, tag = 0)
                # @debug "Sending   $rank -> $i"
            else
                recv_msg[1] = send_msg[1]
            end
        end

    else
        MPI.Recv!(recv_msg, comm; source = controller, tag = 0)
        # @debug "Receiving $controller -> $rank"
    end
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
    synch_responder_to_controller(comm::MPI.Comm, overlapnodes, vector, dof)

Synch the responder to the controller

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `overlapnodes::Dict`: The overlap nodes
- `vector::Vector`: The vector
- `dof::Int`: The degree of freedom
# Returns
- `vector::Vector`: The vector
"""
function synch_responder_to_controller(comm::MPI.Comm, overlapnodes, vector, dof)
    rank = MPI.Comm_rank(comm)
    ncores = MPI.Comm_size(comm)

    if ncores == 1
        return vector
    end

    # Create buffers for send and receive operations
    recv_buffers = Vector{Union{Nothing,Matrix,Vector}}(undef, ncores)
    send_buffers = Vector{Union{Nothing,Matrix,Vector}}(undef, ncores)

    # Prepare send and receive operations
    for jcore in 1:ncores
        if (rank + 1 == jcore)
            continue
        end
        if !isempty(overlapnodes[rank + 1][jcore]["Responder"])
            send_index = overlapnodes[rank + 1][jcore]["Responder"]
            if dof == 1
                send_buffers[jcore] = vector[send_index]
            else
                send_buffers[jcore] = vector[send_index, :]
            end
            # @debug "Sending $rank -> $(jcore-1)"
            # @debug size(send_buffers[jcore])
            # MPI.Isend(send_buffers[jcore], comm; dest = jcore - 1, tag = 0)
            MPI.Send(send_buffers[jcore], comm; dest = jcore - 1, tag = 0)
        end

        if !isempty(overlapnodes[rank + 1][jcore]["Controller"])
            recv_index = overlapnodes[rank + 1][jcore]["Controller"]
            if dof == 1
                recv_buffers[jcore] = similar(vector[recv_index])
            else
                recv_buffers[jcore] = similar(vector[recv_index, :])
            end
            # @debug "Receiving $(jcore-1) -> $rank"
            # @debug size(recv_buffers[jcore])
            MPI.Recv!(recv_buffers[jcore], comm; source = jcore - 1, tag = 0)
            if recv_buffers[jcore][1, 1] isa Bool
                continue
            end
            if dof == 1
                vector[recv_index] .+= recv_buffers[jcore]
            else
                vector[recv_index, :] .+= recv_buffers[jcore]
            end
        end
    end

    return vector
end

"""
    synch_controller_to_responder(comm::MPI.Comm, overlapnodes, vector, dof)

Synch the controller to the responder

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
    rank = MPI.Comm_rank(comm)

    if ncores == 1
        return vector
    end
    for jcore in 1:ncores
        if (rank + 1 == jcore)
            continue
        end
        if !isempty(overlapnodes[rank + 1][jcore]["Controller"])
            send_index = overlapnodes[rank + 1][jcore]["Controller"]
            if dof == 1
                MPI.Send(vector[send_index], comm; dest = jcore - 1, tag = 0)
            else
                MPI.Send(vector[send_index, :], comm; dest = jcore - 1, tag = 0)
            end
            # @debug "Sending $rank -> $(jcore-1)"
        end
        if !isempty(overlapnodes[rank + 1][jcore]["Responder"])
            recv_index = overlapnodes[rank + 1][jcore]["Responder"]
            if dof == 1
                vector[recv_index] = MPI.Recv!(vector[recv_index], comm; source = jcore - 1,
                                               tag = 0)
            else
                vector[recv_index,
                       :] = MPI.Recv!(vector[recv_index, :], comm;
                                      source = jcore - 1, tag = 0)
            end
            # @debug "Receiving $(jcore-1) -> $rank"

            # if dof == 1
            #     recv_msg = similar(vector[recv_index])
            # else
            #     recv_msg = similar(vector[recv_index, :])
            # end
            # MPI.Recv!(recv_msg, comm; source=jcore - 1, tag=0)
            # if dof == 1
            #     vector[recv_index] .= recv_msg
            # else
            #     vector[recv_index, :] .= recv_msg
            # end
            # @debug "Received $(jcore-1) -> $rank = $recv_msg"
        end
    end
    return vector
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

    if ncores == 1
        return array
    end
    for jcore in 1:ncores
        if (rank + 1 == jcore)
            continue
        end
        if !isempty(overlapnodes[rank + 1][jcore]["Controller"])
            for iID in overlapnodes[rank + 1][jcore]["Controller"]
                if dof == 1
                    @views send_msg = array[iID]
                else
                    #TODO: Check if we can remove the [:,:]
                    @views send_msg = mapreduce(permutedims, vcat, array[iID])
                end
                MPI.Isend(send_msg, comm; dest = jcore - 1, tag = 0)
                # @debug "Sending   $rank -> $(jcore-1)"
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
                # @debug "Receiving $(jcore-1) -> $rank"
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
    result
end

"""
    synch_controller_bonds_to_responder_flattened(comm::MPI.Comm, overlapnodes, array, dof)

Synch the controller bonds to the responder

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `overlapnodes::Dict`: The overlap nodes
- `array::Array`: The array
- `dof::Int`: The degree of freedom
# Returns
- `array::Array`: The array
"""
function synch_controller_bonds_to_responder_flattened(comm::MPI.Comm,
                                                       overlapnodes,
                                                       array,
                                                       dof)
    ncores = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)

    if ncores == 1
        return array
    end
    for jcore in 1:ncores
        if (rank + 1 == jcore)
            continue
        end
        if !isempty(overlapnodes[rank + 1][jcore]["Controller"])
            @views send_indices = overlapnodes[rank + 1][jcore]["Controller"]
            # @debug "Sending $rank -> $(jcore-1)"
            MPI.Send(vcat(vcat(array...)[send_indices]...), comm; dest = jcore - 1, tag = 0)
        end
        if !isempty(overlapnodes[rank + 1][jcore]["Responder"])
            @views recv_indices = overlapnodes[rank + 1][jcore]["Responder"]
            row_nums = [length(subarr) for subarr in array[recv_indices]]
            @views recv_msg = zeros(sum(row_nums * dof))
            # @debug "Receiving $(jcore-1) -> $rank"
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
    if currentRank == 0
        for rank in 1:(MPI.Comm_size(comm) - 1)
            MPI.Isend(send_msg[distribution[rank + 1]], comm; dest = rank, tag = 0)
            # @debug "Sending $currentRank -> $rank"
        end
        recv_msg .= send_msg[distribution[1]]
    else
        MPI.Recv!(recv_msg, comm; source = 0, tag = 0)
        # @debug "Receiving 0 -> $currentRank"
    end
    return recv_msg
end

"""
    send_value(comm::MPI.Comm, controller, send_msg)

Sends a value to a controller

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `controller::Int64`: The controller
- `send_msg::Union{Int64,Vector{Float64},Vector{Int64},Vector{Bool}}`: The send message
# Returns
- `recv_msg::Union{Int64,Vector{Float64},Vector{Int64},Vector{Bool}}`: The received message
"""
function send_value(comm::MPI.Comm, controller::Int64,
                    send_msg::T) where {T<:Union{Float64,Int64,
                                                 Vector{Float64},Vector{Vector{Float64}},
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
    return MPI.bcast(send_msg, controller, comm)
end

"""
    find_and_set_core_value_min(comm::MPI.Comm, value::Union{Float64,Int64})

Find and set core value min

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
                                     value::T) where {T<:Union{Float64,Int64,
                                                               Vector{Float64},
                                                               Vector{Int64},
                                                               Matrix{Float64},
                                                               Matrix{Int64}}}
    return MPI.Allreduce(value, MPI.SUM, comm)
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
                                     value::T) where {T<:Union{Float64,Int64,
                                                               Vector{Float64},
                                                               Vector{Int64},
                                                               Matrix{Float64},
                                                               Matrix{Int64}}}
    return MPI.Allreduce(value, MPI.MIN, comm)
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
    average = zero(Float64)
    # must be Float to avoid that at some cores it is Int and at some it is a Float
    if nnodes != 0
        average = value / nnodes
    end
    return MPI.Allreduce(average, MPI.SUM, comm) / MPI.Comm_size(comm)
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

function barrier(comm::MPI.Comm)
    MPI.Barrier(comm)
end

end
