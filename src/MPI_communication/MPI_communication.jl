import MPI
function init_data_fields(comm, nfields, type)
    size = MPI.Comm_size(comm)

    for i in 1:size-1
        local_data = zeros(type, nfields, 1)
    end

    return local_data
end
function send_single_values(comm, master, values, type)
    size = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    recv_msg = zeros(type, 1, 1)
    send_msg = zeros(type, 1, 1)
    if rank == master
        for i = 0:size-1
            if i != master
                send_msg[1] = values[i+1]
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
function main()
    MPI.Init()
    comm = MPI.COMM_WORLD
    size = MPI.Comm_size(comm)
    type = Int32
    rank = 0

    if rank == 0
        data = [1, 2, 3, 4, 5]
        nn = [[1, 2, 3], [2, 3, 4], [3, 4, 5]]
        chunks = [[1, 2], [3, 4], [5]]
        lechunks = zeros(type, size, 1)
        for i in 1:size
            lechunks[i] = length(chunks[i])
        end
    end
    len_chunk = send_single_values(comm, 0, lechunks, type)
    println("Results on rank ", rank, ": ", len_chunk)

    #localData = init_data_fields(comm, len_chunk[1], "Int64")

    MPI.Barrier(comm)
    MPI.Finalize()
end
main()