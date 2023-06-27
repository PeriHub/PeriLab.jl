import MPI

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