import MPI

function send_single_value_from_vector(comm, master, values, type)
    size = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    if type == String
        @error "Wrong type - String in function send_single_value_from_vector"
    end
    recv_msg = zeros(type, 1, 1)
    if rank == master
        send_msg = zeros(type, 1, 1)
        for i = 0:size-1
            # +1 because the index of cores is zero based and julia matrices are one based
            send_msg[1] = values[i+1]
            if i != master
                MPI.Send(send_msg, i, 0, comm)
            else
                recv_msg[1] = send_msg[1]
            end
        end

    else
        MPI.Recv!(recv_msg, master, 0, comm)
    end
    return recv_msg[1]
end
"""
function synch_overlapnodes(comm, topo, vector)
    currentRank = MPI.Comm_rank(comm)
    size = MPI.Comm_size(comm)
    overlapCurrentRank = topo[currentRank+1]
    for icore in 0:size-1
        if icore == currentRank
            continue
        end
        if overlapCurrentRank[icore+1]["Send"] > 0
            send_msg = vector[overlapCurrentRank[icore+1]["Send"]]
            MPI.Send(send_msg, icore+1, 0, comm)
        end
        if overlapCurrentRank[icore+1]["Receive"] > 0
            recv_msg = vector[overlapCurrentRank[icore+1]["Receive"]]
            MPI.Recv!(recv_msg, 0, 0, comm)
            vector[overlapCurrentRank[icore+1]["Receive"]]
        end
    end
    return vector
  
end
"""

function synch_overlapnodes(comm, overlapnodes, vector)

    size = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)

    if size == 1
        return vector
    end
    for jcore in 1:size
        if (rank + 1 == jcore)
            continue
        end
        if length(overlapnodes[rank+1][jcore]["Send"]) > 0
            # check ut_create_overlap_map test
            # the function is clear there
            send_msg = vector[overlapnodes[rank+1][jcore]["Send"]]
            MPI.Send(send_msg, jcore - 1, 0, comm)
        end
        if length(overlapnodes[rank+1][jcore]["Receive"]) > 0
            recv_msg = vector[overlapnodes[rank+1][jcore]["Receive"]]
            MPI.Recv!(recv_msg, jcore - 1, 0, comm)
            vector[overlapnodes[rank+1][jcore]["Receive"]] = recv_msg
        end


    end
    return vector
end
function send_vectors_to_cores(comm, master, values, type)
    #tbd
    size = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    if type == String
        @error "Wrong type - String in function send_vectors_to_cores"
    end
    recv_msg = zeros(type, length(values), 1)
    if rank == master
        send_msg = zeros(type, length(values), 1)
        for i = 0:size-1
            # +1 because the index of cores is zero based and julia matrices are one based
            send_msg = values
            if i != master
                MPI.Send(send_msg, i, 0, comm)
            else
                recv_msg = send_msg
            end
        end

    else
        MPI.Recv!(recv_msg, master, 0, comm)
    end
    return recv_msg
end

function send_vector_from_root_to_core_i(comm, send_msg, recv_msg, distribution)
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
function send_value(comm, master, send_msg)

    if MPI.Comm_rank(comm) == master
        recv_msg = send_msg
    else
        recv_msg = nothing
    end
    recv_msg = MPI.bcast(send_msg, master, comm)
    return recv_msg
end

function get_vector(comm, vector, topo)
    rank = MPI.Comm_rank(comm)
    size = MPI.Comm_size(comm)
    synchTopo = topo[rank]
    for i in 0:size-1
        if i != rank
            recv_msg = zeros(length(synchTopo[i][2]), 1)
            MPI.Send(vector[synchTopo[i][1]], i, 0, comm)
            MPI.Recv!(recv_msg, i, 0, comm)
            vector[synchTopo[i][2]] += recv_msg
        end
    end
    return vector
end


function synch_vector(comm, vector, topo)
    vector = get_vector(comm, vector, topo)
    vector = set_vector(comm, vector, topo)
    return vector
end

function set_overlap_information(comm, topo, vector)
    rank = MPI.Comm_rank(comm)
    size = MPI.Comm_size(comm)
    recvTopo = topo[rank+1]
    for i in 0:size-1
        if i != rank
            if length(recvTopo[i+1]) == 0
                continue
            end
            # temp = vector[recvTopo[i+1][1]]
            #index = recvTopo[i+1]["Receive"]
            #println("Receive $rank - $i index $index value $temp")
            MPI.Recv!(vector[recvTopo[i+1]["Receive"]], i, 0, comm)

        else
            for j in 0:size-1
                if j == rank
                    continue
                end
                sendTopo = topo[j+1]
                if length(sendTopo[rank+1]) > 0
                    # temp = vector[sendTopo[rank+1][2]]
                    # index = sendTopo[rank+1]["Send"]
                    # println("Send $rank - $i index $index value $temp")
                    MPI.Send(vector[sendTopo[rank+1]["Send"]], j, 0, comm)
                end
            end
        end
    end
    return vector
end
function get_overlap_information(comm, vector, topo)
    rank = MPI.Comm_rank(comm)
    size = MPI.Comm_size(comm)
    recvTopo = topo[rank+1]
    vec = vector[rank+1]
    for i in 0:size-1
        if i != rank
            if length(recvTopo[i+1]) > 0
                temp = vec[recvTopo[i+1][1]]
                index = recvTopo[i+1][1]
                println("Receive $rank - $i index $index value $temp")
                recv_msg = zeros(length(recvTopo[i+1][1]), 1)
                MPI.Recv!(recv_msg, i, 0, comm)
                vec[recvTopo[i+1][1]] += recv_msg
            end
        else
            for j in 0:size-1
                sendTopo = topo[j+1]
                if length(sendTopo[rank+1]) > 0
                    temp = vec[sendTopo[rank+1][2]]
                    index = sendTopo[rank+1][2]
                    println("Send $rank - $i index $index value $temp")
                    MPI.Send(vec[sendTopo[rank+1][2]], j, 0, comm)
                end
            end
        end
    end
    return vector
end

function recv_vector_from_root(comm, recv_msg)

    #if rank == MPI.Comm_rank(comm)
    MPI.Recv!(recv_msg, 0, 0, comm)
    # end
    return recv_msg
end


function find_and_set_core_value_min(comm, value)
    return MPI.Allreduce!([value], MPI.MIN, comm)[1]
end

function find_and_set_core_value_max(comm, value)
    return MPI.Allreduce!([value], MPI.MAX, comm)[1]
end