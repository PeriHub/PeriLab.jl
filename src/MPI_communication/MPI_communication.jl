import MPI

function send_single_value_from_vector(comm, master, values, type)
    ncores = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    if type == String
        @error "Wrong type - String in function send_single_value_from_vector"
    end
    recv_msg = zeros(type, 1, 1)
    if rank == master
        send_msg = zeros(type, 1, 1)
        for i = 0:ncores-1
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
    ncores = MPI.Comm_size(comm)
    overlapCurrentRank = topo[currentRank+1]
    for icore in 0:ncores-1
        if icore == currentRank
            continue
        end
        if overlapCurrentRank[icore+1]["Slave"] > 0
            send_msg = vector[overlapCurrentRank[icore+1]["Slave"]]
            MPI.Send(send_msg, icore+1, 0, comm)
        end
        if overlapCurrentRank[icore+1]["Master"] > 0
            recv_msg = vector[overlapCurrentRank[icore+1]["Master"]]
            MPI.Recv!(recv_msg, 0, 0, comm)
            vector[overlapCurrentRank[icore+1]["Master"]]
        end
    end
    return vector
  
end
"""

function synch_slaves_to_master(comm, overlapnodes, vector, dof)
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
        if overlapnodes[rank+1][jcore]["Slave"] != []
            # check ut_create_overlap_map test
            # the function is clear there
            if dof == 1
                send_msg = vector[overlapnodes[rank+1][jcore]["Slave"]]
            else
                send_msg = vector[overlapnodes[rank+1][jcore]["Slave"], :]
            end
            MPI.Send(send_msg, jcore - 1, 0, comm)
        end
        # Receive
        if overlapnodes[rank+1][jcore]["Master"] != []
            if dof == 1
                recv_msg = similar(vector[overlapnodes[rank+1][jcore]["Master"]])
            else
                recv_msg = similar(vector[overlapnodes[rank+1][jcore]["Master"], :])
            end
            MPI.Recv!(recv_msg, jcore - 1, 0, comm)
            if typeof(recv_msg[1, 1]) == Bool
                continue
            end
            if dof == 1
                vector[overlapnodes[rank+1][jcore]["Master"]] += recv_msg
            else
                vector[overlapnodes[rank+1][jcore]["Master"], :] += recv_msg
            end
        end
    end
    return vector
end


function synch_master_to_slaves(comm, overlapnodes, vector, dof)

    ncores = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)

    if ncores == 1
        return vector
    end
    for jcore in 1:ncores
        if (rank + 1 == jcore)
            continue
        end
        if overlapnodes[rank+1][jcore]["Master"] != []
            # check ut_create_overlap_map test
            # the function is clear there
            if dof == 1
                send_msg = vector[overlapnodes[rank+1][jcore]["Master"]]
            else
                send_msg = vector[overlapnodes[rank+1][jcore]["Master"], :]
            end
            MPI.Send(send_msg, jcore - 1, 0, comm)
        end
        if overlapnodes[rank+1][jcore]["Slave"] != []
            if dof == 1
                recv_msg = similar(vector[overlapnodes[rank+1][jcore]["Slave"]])
            else
                recv_msg = similar(vector[overlapnodes[rank+1][jcore]["Slave"], :])
            end
            MPI.Recv!(recv_msg, jcore - 1, 0, comm)
            if dof == 1
                vector[overlapnodes[rank+1][jcore]["Slave"]] = recv_msg
            else
                vector[overlapnodes[rank+1][jcore]["Slave"], :] = recv_msg
            end
        end
    end
    return vector
end

function send_vectors_to_cores(comm, master, values, type)
    #tbd
    ncores = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    if type == String
        @error "Wrong type - String in function send_vectors_to_cores"
    end
    recv_msg = zeros(type, length(values), 1)
    if rank == master
        send_msg = zeros(type, length(values), 1)
        for i = 0:ncores-1
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

function synch_vector(comm, vector, topo)
    vector = get_vector(comm, vector, topo)
    vector = set_vector(comm, vector, topo)
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