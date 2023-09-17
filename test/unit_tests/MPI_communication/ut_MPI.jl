import MPI
using Test
include("../../../src/MPI_communication/MPI_communication.jl")
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

#if size == 3
include("../../../src/IO/mesh_data.jl")
include("../../../src/Support/data_manager.jl")
import .Read_Mesh
using .Data_manager
distribution = [[1, 2, 3], [2, 3, 4], [4, 1, 3]]
size = 3
ptc = [1, 2, 2, 3]
overlap_map = Read_Mesh.create_overlap_map(distribution, ptc, size)

overlap_map = Read_Mesh.get_local_overlap_map(overlap_map, distribution, size)
overlap_mapdof = 2
testDatamanager = Data_manager

if rank == 0
    testDatamanager.set_nmasters(1)
    testDatamanager.set_nslaves(2)
end
if rank == 1
    testDatamanager.set_nmasters(2)
    testDatamanager.set_nslaves(1)
end
if rank == 2
    testDatamanager.set_nmasters(1)
    testDatamanager.set_nslaves(2)
end
testDatamanager.set_dof(2)
A = testDatamanager.create_constant_node_field("A", Float32, 1)
B = testDatamanager.create_constant_node_field("B", Float32, 11)
C = testDatamanager.create_constant_node_field("C", Int64, 1)
D = testDatamanager.create_constant_node_field("D", Int64, 5)
E = testDatamanager.create_constant_node_field("E", Bool, 1)
if rank == 0
    A[1] = 1.4
end
if rank == 1
    A[1] = 3
    A[2] = 5
    B[1, 2] = 10
    B[1, 8] = 10.4
end
if rank == 2
    D[1, 4] = 3
    D[2, 1] = 3
end

A = synch_overlapnodes(comm, overlap_map, A)
if rank == 0
    @test A[1] == 1.4
    @test A[2] == 3
    @test A[3] == 5
end
if rank == 2
    @test A[1] == 3
    @test A[2] == 5
    @test A[3] == 0
end

if rank == 3
    @test A[1] == 0
    @test A[2] == 1.4
    @test A[3] == 5
end

#end
MPI.Finalize()