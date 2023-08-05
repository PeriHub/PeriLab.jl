import MPI
using Test
include("../MPI_communication.jl")
MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)
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



MPI.Finalize()