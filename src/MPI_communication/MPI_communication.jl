# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

import MPI

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
function send_single_value_from_vector(comm::MPI.Comm, controller::Int64, values::Union{Int64,Vector{Float64},Vector{Int64},Vector{Bool}}, type::Type)
    ncores = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    if type == String
        @error "Wrong type - String in function send_single_value_from_vector"
        return nothing
    end
    recv_msg = zeros(type, 1, 1)
    if rank == controller
        send_msg = zeros(type, 1, 1)
        for i = 0:ncores-1
            # +1 because the index of cores is zero based and julia matrices are one based
            send_msg[1] = values[i+1]
            if i != controller
                MPI.Send(send_msg, comm; dest=i, tag=0)
            else
                recv_msg[1] = send_msg[1]
            end
        end

    else
        MPI.Recv!(recv_msg, comm; source=controller, tag=0)
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
    # does not work for bool fields
    ncores = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)

    if ncores == 1
        return vector
    end
    for jcore in 1:ncores
        if (rank + 1 == jcore)
            continue
        end
        # Send
        if overlapnodes[rank+1][jcore]["Responder"] != []
            # check ut_create_overlap_map test
            # the function is clear there
            if dof == 1
                send_msg = vector[overlapnodes[rank+1][jcore]["Responder"]]
            else
                send_msg = vector[overlapnodes[rank+1][jcore]["Responder"], :]
            end
            MPI.Send(send_msg, comm; dest=jcore - 1, tag=0)
        end
        # Receive
        if overlapnodes[rank+1][jcore]["Controller"] != []
            if dof == 1
                recv_msg = similar(vector[overlapnodes[rank+1][jcore]["Controller"]])
            else
                recv_msg = similar(vector[overlapnodes[rank+1][jcore]["Controller"], :])
            end
            MPI.Recv!(recv_msg, comm; source=jcore - 1, tag=0)
            if recv_msg[1, 1] isa Bool
                continue
            end
            #TODO add dot operator, memory leak?
            if dof == 1
                vector[overlapnodes[rank+1][jcore]["Controller"]] .+= recv_msg
            else
                vector[overlapnodes[rank+1][jcore]["Controller"], :] .+= recv_msg
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
        if overlapnodes[rank+1][jcore]["Controller"] != []
            # check ut_create_overlap_map test
            # the function is clear there
            if dof == 1
                send_msg = vector[overlapnodes[rank+1][jcore]["Controller"]]
            else
                send_msg = vector[overlapnodes[rank+1][jcore]["Controller"], :]
            end
            MPI.Send(send_msg, comm; dest=jcore - 1, tag=0)
        end
        if overlapnodes[rank+1][jcore]["Responder"] != []
            if dof == 1
                recv_msg = similar(vector[overlapnodes[rank+1][jcore]["Responder"]])
            else
                recv_msg = similar(vector[overlapnodes[rank+1][jcore]["Responder"], :])
            end
            MPI.Recv!(recv_msg, comm; source=jcore - 1, tag=0)
            if dof == 1
                vector[overlapnodes[rank+1][jcore]["Responder"]] = recv_msg
            else
                vector[overlapnodes[rank+1][jcore]["Responder"], :] = recv_msg
            end
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
        if overlapnodes[rank+1][jcore]["Controller"] != []
            for iID in overlapnodes[rank+1][jcore]["Controller"]
                if dof == 1
                    send_msg = array[iID][:]
                else
                    send_msg = array[iID][:, :]
                end
                MPI.Send(send_msg, comm; dest=jcore - 1, tag=0)
            end
        end
        if overlapnodes[rank+1][jcore]["Responder"] != []
            for iID in overlapnodes[rank+1][jcore]["Responder"]
                if dof == 1
                    recv_msg = similar(array[iID][:])
                else
                    recv_msg = similar(array[iID][:, :])
                end
                MPI.Recv!(recv_msg, comm; source=jcore - 1, tag=0)
                recv_msg = reshape(recv_msg, :, 2)
                if dof == 1
                    array[iID] = recv_msg
                else
                    array[iID][:, :] = recv_msg
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
    result = Vector{Matrix{eltype(input)}}()
    start = firstindex(input)
    for len in row_nums
        push!(result, reshape(input[start:(start+len-1)], (Int64(len / dof), dof)))
        start += len
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
function synch_controller_bonds_to_responder_flattened(comm::MPI.Comm, overlapnodes, array, dof)

    ncores = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)

    if ncores == 1
        return array
    end
    for jcore in 1:ncores
        if (rank + 1 == jcore)
            continue
        end
        if overlapnodes[rank+1][jcore]["Controller"] != []
            send_msg = array[overlapnodes[rank+1][jcore]["Controller"]]
            send_msg = vcat(send_msg...)
            MPI.Send(send_msg, comm; dest=jcore - 1, tag=0)
        end
        if overlapnodes[rank+1][jcore]["Responder"] != []
            row_nums = [length(subarr) for subarr in array[overlapnodes[rank+1][jcore]["Responder"]]]
            recv_msg = zeros(sum(row_nums))
            MPI.Recv!(recv_msg, comm; source=jcore - 1, tag=0)
            recv_msg = split_vector(recv_msg, row_nums, dof)
            if dof == 1
                array[overlapnodes[rank+1][jcore]["Responder"]] = recv_msg
            else
                array[overlapnodes[rank+1][jcore]["Responder"]] = recv_msg
            end
        end
    end
    return array
end

"""
    send_vectors_to_cores(comm::MPI.Comm, controller, values, type)

Sends a vector to a controller

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `controller::Int64`: The controller
- `values::Union{Int64,Vector{Float64},Vector{Int64},Vector{Bool}}`: The values
- `type::Type`: The type
# Returns
- `recv_msg::Union{Int64,Vector{Float64},Vector{Int64},Vector{Bool}}`: The received message
"""
function send_vectors_to_cores(comm::MPI.Comm, controller, values, type)
    #tbd
    ncores = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    if type == String
        @error "Wrong type - String in function send_vectors_to_cores"
    end
    recv_msg = zeros(type, length(values), 1)
    if rank == controller
        send_msg = zeros(type, length(values), 1)
        for i = 0:ncores-1
            # +1 because the index of cores is zero based and julia matrices are one based
            send_msg = values
            if i != controller
                MPI.Send(send_msg, comm; dest=i, tag=0)
            else
                recv_msg = send_msg
            end
        end

    else
        MPI.Recv!(recv_msg, comm; source=controller, tag=0)
    end
    return recv_msg
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
        recv_msg = send_msg[distribution[1]]
        for rank in 1:MPI.Comm_size(comm)-1
            MPI.Send(send_msg[distribution[rank+1]], comm; dest=rank, tag=0)
        end
    else
        MPI.Recv!(recv_msg, comm; source=0, tag=0)
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
function send_value(comm::MPI.Comm, controller, send_msg)

    if MPI.Comm_rank(comm) == controller
        recv_msg = send_msg
    else
        recv_msg = nothing
    end
    recv_msg = MPI.bcast(send_msg, controller, comm)
    return recv_msg
end

"""
    get_vector(comm::MPI.Comm, vector, topo)

Gets a vector from all cores

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `vector::Union{Int64,Vector{Float64},Vector{Int64},Vector{Bool}}`: The vector
- `topo::Vector{Vector{Tuple{Int64,Int64}}}`: The topology
# Returns
- `recv_msg::Union{Int64,Vector{Float64},Vector{Int64},Vector{Bool}}`: The received message
"""
function get_vector(comm::MPI.Comm, vector, topo)
    rank = MPI.Comm_rank(comm)
    ncores = MPI.Comm_size(comm)
    synchTopo = topo[rank]
    for i in 0:ncores-1
        if i != rank
            recv_msg = zeros(length(synchTopo[i][2]), 1)
            MPI.Send(vector[synchTopo[i][1]], comm; source=i, tag=0)
            MPI.Recv!(recv_msg, comm; source=i, tag=0)
            vector[synchTopo[i][2]] += recv_msg
        end
    end
    return vector
end

"""
    recv_vector_from_root(comm::MPI.Comm, recv_msg)

Receives a vector from the root

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `recv_msg::Union{Int64,Vector{Float64},Vector{Int64},Vector{Bool}}`: The receive message
# Returns
- `recv_msg::Union{Int64,Vector{Float64},Vector{Int64},Vector{Bool}}`: The received message
"""
function recv_vector_from_root(comm::MPI.Comm, recv_msg)

    #if rank == MPI.Comm_rank(comm)
    MPI.Recv!(recv_msg, comm; source=0, tag=0)
    # end
    return recv_msg
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
function find_and_set_core_value_min(comm::MPI.Comm, value::Union{Float64,Int64})
    return MPI.Allreduce!([value], MPI.MIN, comm)[1]
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
function find_and_set_core_value_max(comm::MPI.Comm, value::Union{Float64,Int64})
    return MPI.Allreduce!([value], MPI.MAX, comm)[1]
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
function find_and_set_core_value_sum(comm::MPI.Comm, value::Union{Float64,Int64,Bool})
    return MPI.Allreduce!([value], MPI.SUM, comm)[1]
end

"""
    find_and_set_core_value_avg(comm::MPI.Comm, value::Union{Float64,Int64})

Find and set core value avg

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `value::Union{Float64,Int64}`: The value
# Returns
- `recv_msg::Union{Int64,Vector{Float64},Vector{Int64},Vector{Bool}}`: The received message
"""
function find_and_set_core_value_avg(comm::MPI.Comm, value::Union{Float64,Int64}, nnodes::Int64)
    nnodes = MPI.Allreduce!([nnodes], MPI.SUM, comm)[1]
    return MPI.Allreduce!([value], MPI.SUM, comm)[1] / nnodes
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
function gather_values(comm::MPI.Comm, value::Any)
    return MPI.gather(value, comm; root=0)
end