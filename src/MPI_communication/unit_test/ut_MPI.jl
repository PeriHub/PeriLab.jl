import MPI
using Test
include("../MPI_communication.jl")
MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)

@testset "find_and_set_core_value_min and max" begin
    value = rank + 1
    value = find_and_set_core_value_min(comm, value)
    @test value == 1
    value = rank + 1
    value = find_and_set_core_value_max(comm, value)
    @test size == value
end
@testset "ut_send_value" begin
    if rank == 0
        send_msg = 100
    else
        send_msg = nothing
    end
    send_msg = send_value(comm, 0, send_msg)
    @test send_msg == 100
    if rank == 0
        send_msg = true
    else
        send_msg = nothing
    end
    send_msg = send_value(comm, 0, send_msg)
    @test send_msg
    if rank == 0
        send_msg = 100.5
    else
        send_msg = nothing
    end
    send_msg = send_value(comm, 0, send_msg)
    @test send_msg == 100.5
end
@testset "ut_send_vector_from_root_to_core_i" begin
    distribution = [[1, 2, 3], [1, 2, 3], [1, 2, 3]]
    if rank == 0
        send_msg = [2, 1, 5]
    else
        send_msg = nothing
    end
    recv_msg = [0, 0, 0]

    recv_msg = send_vector_from_root_to_core_i(comm, send_msg, recv_msg, distribution)
    @test recv_msg[1] == 2
    @test recv_msg[2] == 1
    @test recv_msg[3] == 5
    distribution = [[1, 2, 3], [3, 2, 1], [3, 2, 1]]
    recv_msg = send_vector_from_root_to_core_i(comm, send_msg, recv_msg, distribution)
    if rank != 0
        @test recv_msg[1] == 5
        @test recv_msg[2] == 1
        @test recv_msg[3] == 2
    end
end
MPI.Finalize()