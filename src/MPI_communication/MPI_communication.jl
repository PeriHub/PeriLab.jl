import MPI
function init_data_fields(comm, globalLen, type)
    localLen = send_single_values(comm, 0, globalLen, type)
    local_data = zeros(type, localLen[1], 1)
    return local_data
end
function send_single_values(comm, master, values, type)
    size = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    recv_msg = zeros(type, 1, 1)
    if rank == master
        send_msg = zeros(type, 1, 1)
        for i = 0:size-1
            send_msg[1] = values[i+1]
            if i != master
                MPI.Send(send_msg, i, 0, comm)
            end
        end
        # +1 because the index of cores is zero based and julia matrices are one based
        recv_msg[master+1] = send_msg[master+1]
    else
        MPI.Recv!(recv_msg, master, 0, comm)
    end
    return recv_msg
end
