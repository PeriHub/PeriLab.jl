# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

import MPI
using Test
include("../../../src/MPI_communication/MPI_communication.jl")
MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
ncores = MPI.Comm_size(comm)

@testset "find_and_set_core_value_min and max" begin
    value = rank + 1
    value = find_and_set_core_value_min(comm, value)
    @test value == 1
    value = rank + 1
    value = find_and_set_core_value_max(comm, value)
    @test ncores == value
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

if ncores == 3
    include("../../../src/IO/mesh_data.jl")
    include("../../../src/Support/data_manager.jl")
    import .Read_Mesh
    using .Data_manager
    distribution = [[1, 2, 3], [2, 3, 4], [4, 1, 3]]
    ncores = 3
    ptc = [1, 2, 2, 3]
    overlap_map = Read_Mesh.create_overlap_map(distribution, ptc, ncores)

    overlap_map = Read_Mesh.get_local_overlap_map(overlap_map, distribution, ncores)

    test_Data_manager = Data_manager

    if rank == 0
        test_Data_manager.set_nmasters(1)
        test_Data_manager.set_nslaves(2)
    end
    if rank == 1
        test_Data_manager.set_nmasters(2)
        test_Data_manager.set_nslaves(1)
    end
    if rank == 2
        test_Data_manager.set_nmasters(1)
        test_Data_manager.set_nslaves(2)
    end
    test_Data_manager.set_dof(2)
    A = test_Data_manager.create_constant_node_field("A", Float64, 1)
    B = test_Data_manager.create_constant_node_field("B", Float64, 4)
    C = test_Data_manager.create_constant_node_field("C", Int64, 1)
    D = test_Data_manager.create_constant_node_field("D", Int64, 5)
    E = test_Data_manager.create_constant_node_field("E", Bool, 1)
    if rank == 0
        A[1] = 1.4
        A[2] = 3
        A[3] = 5
        B[3, 2] = -5
        B[1, 4] = 10.4
        C[1] = 1
        C[2] = 2
        C[3] = 3
        D[1, 1] = 2
        E[2] = true

    end
    if rank == 1

        A[1] = 1
        A[2] = 88
        A[3] = 1.6

        B[1, 2] = 10
        B[1, 4] = 10.4

        C[1] = 1
        C[2] = 2
        C[3] = 3
        D[1, 1] = 2
        E[1] = true
        E[2] = true
        E[3] = true
    end
    if rank == 2
        A[2] = -2.3
        A[3] = 1.1
        B[2, 4] = 3
        C[1] = 1
        C[2] = 2
        C[3] = 3
        D[1, 4] = 3
        D[2, 1] = 3
    end
    distribution = [[1, 2, 3], [2, 3, 4], [4, 1, 3]]

    # sammel ein und summiere -> zweite routine mit sende vom MAster an alle slave
    A[:] = synch_slaves_to_master(comm, overlap_map, A, 1)
    B[:] = synch_slaves_to_master(comm, overlap_map, B, 4)
    C[:] = synch_slaves_to_master(comm, overlap_map, C, 1)
    D[:] = synch_slaves_to_master(comm, overlap_map, D, 5)
    E[:] = synch_slaves_to_master(comm, overlap_map, E, 1)

    if rank == 0
        @testset "synch_slaves_to_master_rank_0" begin
            @test A[1] == Float64(-2.3 + 1.4)
            @test A[2] == 3
            @test A[3] == 5
            @test B[1, 4] == Float64(13.4)
            @test B[3, 2] == -5
            @test C[1] == 3
            @test C[2] == 2
            @test C[3] == 3
            @test D[1, 1] == 5
            @test E[1] == false
            @test E[2] == true
            @test E[3] == false

        end
    end
    if rank == 1
        @testset "synch_slaves_to_master_rank_1" begin
            @test A[1] == 3 + 1
            @test A[2] == Float64(88 + 5 + 1.1)
            @test A[3] == Float64(1.6)
            @test B[1, 2] == Float64(10)
            @test B[1, 4] == Float64(10.4)
            @test B[2, 2] == Float64(-5)
            @test C[1] == 3
            @test C[2] == 8
            @test C[3] == 3
            @test D[1, 1] == 2
            @test E[1] == true
            @test E[2] == true
            @test E[3] == true
        end
    end
    if rank == 2
        @testset "synch_slaves_to_master_rank_3" begin
            @test A[1] == Float64(1.6)
            @test A[2] == Float64(-2.3)
            @test A[3] == Float64(1.1)
            @test B[1, 1] == Float64(0.0)
            @test B[1, 2] == Float64(0.0)
            @test B[1, 3] == Float64(0.0)
            @test B[1, 4] == Float64(0.0)
            @test B[2, 3] == Float64(0.0)
            @test B[2, 4] == Float64(3.0)

            @test C[1] == 4
            @test C[2] == 2
            @test C[3] == 3

            @test D[1, 4] == 3
            @test D[2, 1] == 3

            @test E[1] == false
            @test E[2] == false
            @test E[3] == false
        end
    end
    MPI.Barrier(comm)

    A[:] = synch_master_to_slaves(comm, overlap_map, A, 1)
    B[:] = synch_master_to_slaves(comm, overlap_map, B, 4)
    C[:] = synch_master_to_slaves(comm, overlap_map, C, 1)
    D[:] = synch_master_to_slaves(comm, overlap_map, D, 5)
    E[:] = synch_master_to_slaves(comm, overlap_map, E, 1)
    if rank == 0
        @testset "synch_master_to_slaves_rank_0" begin
            @test A[1] == Float64(-0.9)
            @test A[2] == Float64(4)
            @test A[3] == Float64(94.1)
            @test B[1, 4] == Float64(13.4)
            @test B[2, 2] == Float64(10)
            @test B[2, 4] == Float64(10.4)
            @test B[3, 2] == Float64(-5.0)
            @test C[1] == 3
            @test C[2] == 3
            @test C[3] == 8
            @test D[1, 1] == 5
            @test D[2, 1] == 2
            @test E[1] == false
            @test E[2] == true
            @test E[3] == true

        end
    end
    if rank == 1
        @testset "synch_master_to_slaves_rank_1" begin
            @test A[1] == Float64(4.0)
            @test A[2] == Float64(94.1)
            @test A[3] == Float64(1.6)
            @test B[1, 2] == Float64(10)
            @test B[1, 4] == Float64(10.4)
            @test B[2, 2] == Float64(-5)
            @test C[1] == 3
            @test C[2] == 8
            @test C[3] == 4
            @test D[1, 1] == 2
            @test D[3, 3] == 0
            @test D[3, 4] == 3
            @test D[3, 5] == 0
            @test E[1] == true
            @test E[2] == true
            @test E[3] == false
        end
    end
    if rank == 2
        @testset "synch_master_to_slaves_rank_3" begin
            @test A[1] == Float64(1.6)
            @test A[2] == Float64(-0.9)
            @test A[3] == Float64(94.1)
            @test B[1, 1] == Float64(0.0)
            @test B[1, 2] == Float64(0.0)
            @test B[1, 3] == Float64(0.0)
            @test B[1, 4] == Float64(0.0)
            @test B[2, 3] == Float64(0.0)
            @test B[2, 4] == Float64(13.4)
            @test B[3, 2] == Float64(-5.0)
            @test C[1] == 4
            @test C[2] == 3
            @test C[3] == 8

            @test D[1, 4] == 3
            @test D[2, 1] == 5

            @test E[1] == false
            @test E[2] == false
            @test E[3] == true
        end
    end

end
MPI.Finalize()