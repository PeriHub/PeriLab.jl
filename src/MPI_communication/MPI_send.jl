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
    receive_msg = ""
    receive_msg = MPI.bcast(send_msg, master, comm)
    if MPI.Comm_rank(comm) == master
        receive_msg = send_msg
    end
    return receive_msg
end