import MPI

function send_single_value_from_vector(comm, master, values, type)
    size = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    if type == String
        @error "Wrong type String in function send_single_value_from_vector"
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
    return recv_msg
end
function send_string(comm, master, send_msg::String)
    recv_msg = ""
    recv_msg = MPI.bcast(send_msg, master, comm)
    if MPI.Comm_rank(comm) == master
        recv_msg = send_msg
    end
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

function set_overlap_information(comm, vector, topo)
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
                MPI.Recv!(vec[recvTopo[i+1][1]], i, 0, comm)
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