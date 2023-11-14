# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

import MPI

function send_single_value_from_vector(comm::MPI.Comm, controller::Int64, values::Union{Int64,Vector{Float64},Vector{Int64},Vector{Bool}}, type::Type)
    ncores = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    if type == String
        @error "Wrong type - String in function send_single_value_from_vector"
    end
    recv_msg = zeros(type, 1, 1)
    if rank == controller
        send_msg = zeros(type, 1, 1)
        for i = 0:ncores-1
            # +1 because the index of cores is zero based and julia matrices are one based
            send_msg[1] = values[i+1]
            if i != controller
                MPI.Send(send_msg, i, 0, comm)
            else
                recv_msg[1] = send_msg[1]
            end
        end

    else
        MPI.Recv!(recv_msg, controller, 0, comm)
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
            MPI.Send(send_msg, icore+1, 0, comm)
        end
        if overlapCurrentRank[icore+1]["Controller"] > 0
            recv_msg = vector[overlapCurrentRank[icore+1]["Controller"]]
            MPI.Recv!(recv_msg, 0, 0, comm)
            vector[overlapCurrentRank[icore+1]["Controller"]]
        end
    end
    return vector
  
end
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
            MPI.Send(send_msg, jcore - 1, 0, comm)
        end
        # Receive
        if overlapnodes[rank+1][jcore]["Controller"] != []
            if dof == 1
                recv_msg = similar(vector[overlapnodes[rank+1][jcore]["Controller"]])
            else
                recv_msg = similar(vector[overlapnodes[rank+1][jcore]["Controller"], :])
            end
            MPI.Recv!(recv_msg, jcore - 1, 0, comm)
            if typeof(recv_msg[1, 1]) == Bool
                continue
            end
            if dof == 1
                vector[overlapnodes[rank+1][jcore]["Controller"]] += recv_msg
            else
                vector[overlapnodes[rank+1][jcore]["Controller"], :] += recv_msg
            end
        end
    end
    return vector
end


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
            MPI.Send(send_msg, jcore - 1, 0, comm)
        end
        if overlapnodes[rank+1][jcore]["Responder"] != []
            if dof == 1
                recv_msg = similar(vector[overlapnodes[rank+1][jcore]["Responder"]])
            else
                recv_msg = similar(vector[overlapnodes[rank+1][jcore]["Responder"], :])
            end
            MPI.Recv!(recv_msg, jcore - 1, 0, comm)
            if dof == 1
                vector[overlapnodes[rank+1][jcore]["Responder"]] = recv_msg
            else
                vector[overlapnodes[rank+1][jcore]["Responder"], :] = recv_msg
            end
        end
    end
    return vector
end

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
                MPI.Send(send_msg, i, 0, comm)
            else
                recv_msg = send_msg
            end
        end

    else
        MPI.Recv!(recv_msg, controller, 0, comm)
    end
    return recv_msg
end

function send_vector_from_root_to_core_i(comm::MPI.Comm, send_msg, recv_msg, distribution)
    currentRank = MPI.Comm_rank(comm)
    if currentRank == 0
        recv_msg = send_msg[distribution[1]]
        for rank in 1:MPI.Comm_size(comm)-1
            MPI.Send(send_msg[distribution[rank+1]], rank, 0, comm)
        end
    else
        MPI.Recv!(recv_msg, 0, 0, comm)
    end
    return recv_msg
end
function send_value(comm::MPI.Comm, controller, send_msg)

    if MPI.Comm_rank(comm) == controller
        recv_msg = send_msg
    else
        recv_msg = nothing
    end
    recv_msg = MPI.bcast(send_msg, controller, comm)
    return recv_msg
end

function get_vector(comm::MPI.Comm, vector, topo)
    rank = MPI.Comm_rank(comm)
    ncores = MPI.Comm_size(comm)
    synchTopo = topo[rank]
    for i in 0:ncores-1
        if i != rank
            recv_msg = zeros(length(synchTopo[i][2]), 1)
            MPI.Send(vector[synchTopo[i][1]], i, 0, comm)
            MPI.Recv!(recv_msg, i, 0, comm)
            vector[synchTopo[i][2]] += recv_msg
        end
    end
    return vector
end

function recv_vector_from_root(comm::MPI.Comm, recv_msg)

    #if rank == MPI.Comm_rank(comm)
    MPI.Recv!(recv_msg, 0, 0, comm)
    # end
    return recv_msg
end


function find_and_set_core_value_min(comm::MPI.Comm, value::Union{Float64,Int64})
    return MPI.Allreduce!([value], MPI.MIN, comm)[1]
end

function find_and_set_core_value_max(comm::MPI.Comm, value::Union{Float64,Int64})
    return MPI.Allreduce!([value], MPI.MAX, comm)[1]
end

function find_and_set_core_value_sum(comm::MPI.Comm, value::Union{Float64,Int64,Bool})
    return MPI.Allreduce!([value], MPI.SUM, comm)[1]
end

function find_and_set_core_value_avg(comm::MPI.Comm, value::Union{Float64,Int64}, nnodes::Int64)
    nnodes = MPI.Allreduce!([nnodes], MPI.SUM, comm)[1]
    return MPI.Allreduce!([value], MPI.SUM, comm)[1] / nnodes
end