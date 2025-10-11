# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

import MPI
using Test
using JSON3
using TimerOutputs
using PeriLab

function push_test!(dict::Dict, test::Bool, file::String, line::Int)
    push!(dict["tests"], test)
    push!(dict["line"], "$file:$line")
end

MPI.Init()

const to = TimerOutput()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
ncores = MPI.Comm_size(comm)

test_dict = Dict()

test = test_dict["find_and_set_core_value_min and max"] = Dict("tests" => [], "line" => [])
value = rank + 1
value = PeriLab.MPI_communication.find_and_set_core_value_min(comm, value)
push_test!(test, (value == 1), @__FILE__, @__LINE__)
value = rank + 1
value = PeriLab.MPI_communication.find_and_set_core_value_max(comm, value)
push_test!(test, (ncores == value), @__FILE__, @__LINE__)

test = test_dict["ut_broadcast_value"] = Dict("tests" => [], "line" => [])
if rank == 0
    send_msg = 100
else
    send_msg = nothing
end
send_msg = PeriLab.MPI_communication.broadcast_value(comm, send_msg)
push_test!(test, (send_msg == 100), @__FILE__, @__LINE__)
if rank == 0
    send_msg = true
else
    send_msg = nothing
end
send_msg = PeriLab.MPI_communication.broadcast_value(comm, send_msg)
push_test!(test, (send_msg), @__FILE__, @__LINE__)
if rank == 0
    send_msg = 100.5
else
    send_msg = nothing
end
send_msg = PeriLab.MPI_communication.broadcast_value(comm, send_msg)
push_test!(test, (send_msg == 100.5), @__FILE__, @__LINE__)

test = test_dict["ut_send_vector_from_root_to_core_i"] = Dict("tests" => [], "line" => [])
distribution = [[1, 2, 3], [1, 2, 3], [1, 2, 3]]
if rank == 0
    send_msg = [2, 1, 5]
else
    send_msg = nothing
end
recv_msg = [0, 0, 0]

recv_msg = PeriLab.MPI_communication.send_vector_from_root_to_core_i(comm, send_msg, recv_msg, distribution)
push_test!(test, (recv_msg[1] == 2), @__FILE__, @__LINE__)
push_test!(test, (recv_msg[2] == 1), @__FILE__, @__LINE__)
push_test!(test, (recv_msg[3] == 5), @__FILE__, @__LINE__)
distribution = [[1, 2, 3], [3, 2, 1], [3, 2, 1]]
recv_msg = PeriLab.MPI_communication.send_vector_from_root_to_core_i(comm, send_msg, recv_msg, distribution)
if rank != 0
    push_test!(test, (recv_msg[1] == 5), @__FILE__, @__LINE__)
    push_test!(test, (recv_msg[2] == 1), @__FILE__, @__LINE__)
    push_test!(test, (recv_msg[3] == 2), @__FILE__, @__LINE__)
end

push_test!(test,
           (isnothing(PeriLab.MPI_communication.send_single_value_from_vector(comm, 0, [1], String))),
           @__FILE__,
           @__LINE__)
if ncores == 3
    distribution = [[1, 2, 3], [2, 3, 4], [4, 1, 3]]
    ncores = 3
    dof = 2
    ptc = [1, 2, 2, 3]
    overlap_map = PeriLab.IO.create_overlap_map(distribution, ptc, ncores)

    overlap_map = PeriLab.IO.get_local_overlap_map(overlap_map, distribution, ncores)

    test_data_manager = PeriLab.Data_manager
    test_data_manager.initialize_data()
    test_data_manager.set_comm(comm)
    test_data_manager.create_constant_node_field("Block_Id", Int64, 1)

    if rank == 0
        test_data_manager.set_num_controller(1)
        test_data_manager.set_num_responder(2)
        block_Id = test_data_manager.get_field("Block_Id")
        block_Id .= 1
    end
    if rank == 1
        test_data_manager.set_num_controller(2)
        test_data_manager.set_num_responder(1)
        block_Id = test_data_manager.get_field("Block_Id")
        block_Id .= 2
    end
    if rank == 2
        test_data_manager.set_num_controller(1)
        test_data_manager.set_num_responder(2)
        block_Id = test_data_manager.get_field("Block_Id")
        block_Id .= 1
    end
    test_data_manager.set_block_name_list(["block_1", "block_2"])
    test_data_manager.set_block_id_list([1, 2])
    test_data_manager.set_dof(dof)
    A = test_data_manager.create_constant_node_field("A", Float64, 1)
    B = test_data_manager.create_constant_node_field("B", Float64, 4)
    C = test_data_manager.create_constant_node_field("C", Int64, 1)
    D = test_data_manager.create_constant_node_field("D", Int64, 5)
    E = test_data_manager.create_constant_node_field("E", Bool, 1)
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

    # sammel ein und summiere -> zweite routine mit sende vom Controller an alle responder
    A = PeriLab.MPI_communication.synch_responder_to_controller(comm, overlap_map, A, 1)
    B = PeriLab.MPI_communication.synch_responder_to_controller(comm, overlap_map, B, 4)
    C = PeriLab.MPI_communication.synch_responder_to_controller(comm, overlap_map, C, 1)
    D = PeriLab.MPI_communication.synch_responder_to_controller(comm, overlap_map, D, 5)
    E = PeriLab.MPI_communication.synch_responder_to_controller(comm, overlap_map, E, 1)

    if rank == 0
        test = test_dict["synch_responder_to_controller_rank_0"] = Dict("tests" => [],
                                                                        "line" => [])
        push_test!(test, (A[1] == Float64(-2.3 + 1.4)), @__FILE__, @__LINE__)
        push_test!(test, (A[2] == 3), @__FILE__, @__LINE__)
        push_test!(test, (A[3] == 5), @__FILE__, @__LINE__)
        push_test!(test, (B[1, 4] == Float64(13.4)), @__FILE__, @__LINE__)
        push_test!(test, (B[3, 2] == -5), @__FILE__, @__LINE__)
        push_test!(test, (C[1] == 3), @__FILE__, @__LINE__)
        push_test!(test, (C[2] == 2), @__FILE__, @__LINE__)
        push_test!(test, (C[3] == 3), @__FILE__, @__LINE__)
        push_test!(test, (D[1, 1] == 5), @__FILE__, @__LINE__)
        push_test!(test, (E[1] == false), @__FILE__, @__LINE__)
        push_test!(test, (E[2] == true), @__FILE__, @__LINE__)
        push_test!(test, (E[3] == false), @__FILE__, @__LINE__)
    end
    if rank == 1
        test = test_dict["synch_responder_to_controller_rank_1"] = Dict("tests" => [],
                                                                        "line" => [])
        push_test!(test, (A[1] == 3 + 1), @__FILE__, @__LINE__)
        push_test!(test, (A[2] == Float64(88 + 5 + 1.1)), @__FILE__, @__LINE__)
        push_test!(test, (A[3] == Float64(1.6)), @__FILE__, @__LINE__)
        push_test!(test, (B[1, 2] == Float64(10)), @__FILE__, @__LINE__)
        push_test!(test, (B[1, 4] == Float64(10.4)), @__FILE__, @__LINE__)
        push_test!(test, (B[2, 2] == Float64(-5)), @__FILE__, @__LINE__)
        push_test!(test, (C[1] == 3), @__FILE__, @__LINE__)
        push_test!(test, (C[2] == 8), @__FILE__, @__LINE__)
        push_test!(test, (C[3] == 3), @__FILE__, @__LINE__)
        push_test!(test, (D[1, 1] == 2), @__FILE__, @__LINE__)
        push_test!(test, (E[1] == true), @__FILE__, @__LINE__)
        push_test!(test, (E[2] == true), @__FILE__, @__LINE__)
        push_test!(test, (E[3] == true), @__FILE__, @__LINE__)
    end
    if rank == 2
        test = test_dict["synch_responder_to_controller_rank_2"] = Dict("tests" => [],
                                                                        "line" => [])
        push_test!(test, (A[1] == Float64(1.6)), @__FILE__, @__LINE__)
        push_test!(test, (A[2] == Float64(-2.3)), @__FILE__, @__LINE__)
        push_test!(test, (A[3] == Float64(1.1)), @__FILE__, @__LINE__)
        push_test!(test, (B[1, 1] == Float64(0.0)), @__FILE__, @__LINE__)
        push_test!(test, (B[1, 2] == Float64(0.0)), @__FILE__, @__LINE__)
        push_test!(test, (B[1, 3] == Float64(0.0)), @__FILE__, @__LINE__)
        push_test!(test, (B[1, 4] == Float64(0.0)), @__FILE__, @__LINE__)
        push_test!(test, (B[2, 3] == Float64(0.0)), @__FILE__, @__LINE__)
        push_test!(test, (B[2, 4] == Float64(3.0)), @__FILE__, @__LINE__)

        push_test!(test, (C[1] == 4), @__FILE__, @__LINE__)
        push_test!(test, (C[2] == 2), @__FILE__, @__LINE__)
        push_test!(test, (C[3] == 3), @__FILE__, @__LINE__)

        push_test!(test, (D[1, 4] == 3), @__FILE__, @__LINE__)
        push_test!(test, (D[2, 1] == 3), @__FILE__, @__LINE__)

        push_test!(test, (E[1] == false), @__FILE__, @__LINE__)
        push_test!(test, (E[2] == false), @__FILE__, @__LINE__)
        push_test!(test, (E[3] == false), @__FILE__, @__LINE__)
    end
    PeriLab.MPI_communication.barrier(comm)

    A = PeriLab.MPI_communication.synch_controller_to_responder(comm, overlap_map, A, 1)
    B = PeriLab.MPI_communication.synch_controller_to_responder(comm, overlap_map, B, 4)
    C = PeriLab.MPI_communication.synch_controller_to_responder(comm, overlap_map, C, 1)
    D = PeriLab.MPI_communication.synch_controller_to_responder(comm, overlap_map, D, 5)
    E = PeriLab.MPI_communication.synch_controller_to_responder(comm, overlap_map, E, 1)
    if rank == 0
        test = test_dict["synch_controller_to_responder_rank_0"] = Dict("tests" => [],
                                                                        "line" => [])
        push_test!(test, isapprox(A[1], Float64(-0.9)), @__FILE__, @__LINE__)
        push_test!(test, (A[2] == Float64(4)), @__FILE__, @__LINE__)
        push_test!(test, (A[3] == Float64(94.1)), @__FILE__, @__LINE__)
        push_test!(test, (B[1, 4] == Float64(13.4)), @__FILE__, @__LINE__)
        push_test!(test, (B[2, 2] == Float64(10)), @__FILE__, @__LINE__)
        push_test!(test, (B[2, 4] == Float64(10.4)), @__FILE__, @__LINE__)
        push_test!(test, (B[3, 2] == Float64(-5.0)), @__FILE__, @__LINE__)
        push_test!(test, (C[1] == 3), @__FILE__, @__LINE__)
        push_test!(test, (C[2] == 3), @__FILE__, @__LINE__)
        push_test!(test, (C[3] == 8), @__FILE__, @__LINE__)
        push_test!(test, (D[1, 1] == 5), @__FILE__, @__LINE__)
        push_test!(test, (D[2, 1] == 2), @__FILE__, @__LINE__)
        push_test!(test, (E[1] == false), @__FILE__, @__LINE__)
        push_test!(test, (E[2] == true), @__FILE__, @__LINE__)
        push_test!(test, (E[3] == true), @__FILE__, @__LINE__)
        test = test_dict["find_global_core_value!_0"] = Dict("tests" => [], "line" => [])
        push_test!(test,
                   (PeriLab.IO.find_global_core_value!(0, "Sum", 1, test_data_manager) == 3),
                   @__FILE__,
                   @__LINE__)
        push_test!(test,
                   (PeriLab.IO.find_global_core_value!(0, "Maximum", 1, test_data_manager) == 2),
                   @__FILE__,
                   @__LINE__)
        push_test!(test,
                   (PeriLab.IO.find_global_core_value!(0, "Minimum", 1, test_data_manager) == 0),
                   @__FILE__,
                   @__LINE__)
        push_test!(test,
                   (PeriLab.IO.find_global_core_value!(0, "Average", 1, test_data_manager) == 1),
                   @__FILE__,
                   @__LINE__)
    end
    if rank == 1
        test = test_dict["synch_controller_to_responder_rank_1"] = Dict("tests" => [],
                                                                        "line" => [])
        push_test!(test, (A[1] == Float64(4.0)), @__FILE__, @__LINE__)
        push_test!(test, (A[2] == Float64(94.1)), @__FILE__, @__LINE__)
        push_test!(test, (A[3] == Float64(1.6)), @__FILE__, @__LINE__)
        push_test!(test, (B[1, 2] == Float64(10)), @__FILE__, @__LINE__)
        push_test!(test, (B[1, 4] == Float64(10.4)), @__FILE__, @__LINE__)
        push_test!(test, (B[2, 2] == Float64(-5)), @__FILE__, @__LINE__)
        push_test!(test, (C[1] == 3), @__FILE__, @__LINE__)
        push_test!(test, (C[2] == 8), @__FILE__, @__LINE__)
        push_test!(test, (C[3] == 4), @__FILE__, @__LINE__)
        push_test!(test, (D[1, 1] == 2), @__FILE__, @__LINE__)
        push_test!(test, (D[3, 3] == 0), @__FILE__, @__LINE__)
        push_test!(test, (D[3, 4] == 3), @__FILE__, @__LINE__)
        push_test!(test, (D[3, 5] == 0), @__FILE__, @__LINE__)
        push_test!(test, (E[1] == true), @__FILE__, @__LINE__)
        push_test!(test, (E[2] == true), @__FILE__, @__LINE__)
        push_test!(test, (E[3] == false), @__FILE__, @__LINE__)
        test = test_dict["find_global_core_value!_1"] = Dict("tests" => [], "line" => [])
        push_test!(test,
                   (PeriLab.IO.find_global_core_value!(1, "Sum", 1, test_data_manager) == 3),
                   @__FILE__,
                   @__LINE__)
        push_test!(test,
                   (PeriLab.IO.find_global_core_value!(1, "Maximum", 1, test_data_manager) == 2),
                   @__FILE__,
                   @__LINE__)
        push_test!(test,
                   (PeriLab.IO.find_global_core_value!(1, "Minimum", 1, test_data_manager) == 0),
                   @__FILE__,
                   @__LINE__)
        push_test!(test,
                   (PeriLab.IO.find_global_core_value!(1, "Average", 1, test_data_manager) == 1),
                   @__FILE__,
                   @__LINE__)
    end
    if rank == 2
        test = test_dict["synch_controller_to_responder_rank_2"] = Dict("tests" => [],
                                                                        "line" => [])
        push_test!(test, (A[1] == Float64(1.6)), @__FILE__, @__LINE__)
        push_test!(test, isapprox(A[2], Float64(-0.9)), @__FILE__, @__LINE__)
        push_test!(test, (A[3] == Float64(94.1)), @__FILE__, @__LINE__)
        push_test!(test, (B[1, 1] == Float64(0.0)), @__FILE__, @__LINE__)
        push_test!(test, (B[1, 2] == Float64(0.0)), @__FILE__, @__LINE__)
        push_test!(test, (B[1, 3] == Float64(0.0)), @__FILE__, @__LINE__)
        push_test!(test, (B[1, 4] == Float64(0.0)), @__FILE__, @__LINE__)
        push_test!(test, (B[2, 3] == Float64(0.0)), @__FILE__, @__LINE__)
        push_test!(test, (B[2, 4] == Float64(13.4)), @__FILE__, @__LINE__)
        push_test!(test, (B[3, 2] == Float64(-5.0)), @__FILE__, @__LINE__)
        push_test!(test, (C[1] == 4), @__FILE__, @__LINE__)
        push_test!(test, (C[2] == 3), @__FILE__, @__LINE__)
        push_test!(test, (C[3] == 8), @__FILE__, @__LINE__)

        push_test!(test, (D[1, 4] == 3), @__FILE__, @__LINE__)
        push_test!(test, (D[2, 1] == 5), @__FILE__, @__LINE__)

        push_test!(test, (E[1] == false), @__FILE__, @__LINE__)
        push_test!(test, (E[2] == false), @__FILE__, @__LINE__)
        push_test!(test, (E[3] == true), @__FILE__, @__LINE__)
        test = test_dict["find_global_core_value!_2"] = Dict("tests" => [], "line" => [])
        push_test!(test,
                   (PeriLab.IO.find_global_core_value!(2, "Sum", 1, test_data_manager) == 3),
                   @__FILE__,
                   @__LINE__)
        push_test!(test,
                   (PeriLab.IO.find_global_core_value!(2, "Maximum", 1, test_data_manager) == 2),
                   @__FILE__,
                   @__LINE__)
        push_test!(test,
                   (PeriLab.IO.find_global_core_value!(2, "Minimum", 1, test_data_manager) == 0),
                   @__FILE__,
                   @__LINE__)
        push_test!(test,
                   (PeriLab.IO.find_global_core_value!(2, "Average", 1, test_data_manager) == 1),
                   @__FILE__,
                   @__LINE__)
    end
    nn = test_data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn .= 2
    h = test_data_manager.create_constant_node_field("Horizon", Float64, 1)
    nodes = test_data_manager.get_nnodes()
    h .= 5.0
    bf = test_data_manager.create_constant_bond_field("Bond Forces", Float64, dof)

    bdN, bdNP1 = test_data_manager.create_bond_field("Bond Damage", Float64, 1, 1)
    dbN,
    dbNP1 = test_data_manager.create_bond_field("Deformed Bond Geometry", Float64, dof,
                                                1)
    dbdN, dbdNP1 = test_data_manager.create_bond_field("Deformed Bond Length", Float64, 1)
    bg = test_data_manager.create_constant_bond_field("Bond Geometry", Float64, dof)
    bd = test_data_manager.create_constant_bond_field("Bond Length", Float64, 1, 1)
    for iID in 1:nodes
        for jID in 1:nn[iID]
            dbdNP1[iID][jID] = 1 + (-1)^iID * 0.1
        end
    end

    test_data_manager = PeriLab.Solver_control.Model_Factory.Material.Bondbased_Elastic.init_model(test_data_manager,
                                                     Vector{Int64}(1:nodes),
                                                     Dict("Bulk Modulus" => 1.0,
                                                          "Young's Modulus" => 1.0))
    test_data_manager = PeriLab.Solver_control.Model_Factory.Material.Bondbased_Elastic.compute_model(test_data_manager,
                                                        Vector{Int64}(1:nodes),
                                                        Dict("Bulk Modulus" => 1.0,
                                                             "Young's Modulus" => 1.0),
                                                        1,
                                                        0.0,
                                                        0.0,
                                                        to)

    bf = test_data_manager.get_field("Bond Forces")

    PeriLab.MPI_communication.synch_controller_bonds_to_responder(comm, overlap_map, bf, dof)

    if rank == 0
        test = test_dict["synch_controller_bonds_to_responder_rank_0"] = Dict("tests" => [],
                                                                              "line" => [])
        # push_test!(test, (bf[1] == Float64(-0.9)), @__FILE__, @__LINE__)
    end

    PeriLab.MPI_communication.synch_controller_bonds_to_responder_flattened(comm, overlap_map, bf, dof)
    if rank == 0
        test = test_dict["synch_controller_bonds_to_responder_flattened_rank_0"] = Dict("tests" => [],
                                                                                        "line" => [])
        # push_test!(test, (bf[1] == Float64(-0.9)), @__FILE__, @__LINE__)
    end

    solver_options = Dict("Models" => ["Material"])
    params = Dict("Blocks" => Dict("block_1" => Dict("Block ID" => 1,
                                                     "Material Model" => "Test 1"),
                                   "block_2" => Dict("Block ID" => 2,
                                                     "Material Model" => "Test 2")))
    PeriLab.IO.show_block_summary(solver_options, params, "", false, comm, test_data_manager)
    PeriLab.IO.show_block_summary(solver_options, params, "", true, comm, test_data_manager)
    PeriLab.IO.show_mpi_summary("", false, comm, test_data_manager)
    PeriLab.IO.show_mpi_summary("", true, comm, test_data_manager)
end

open("test_results_$rank.json", "w") do f
    JSON3.pretty(f, JSON3.write(test_dict))
end

MPI.Finalize()
