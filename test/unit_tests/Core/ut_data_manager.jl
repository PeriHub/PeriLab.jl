# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using MPI
using Test

#include("../../../src/PeriLab.jl")
#using .PeriLab

@testset "set_comm" begin
    # MPI.Init()
    comm = MPI.COMM_WORLD
    test_data_manager = PeriLab.Data_manager
    test_data_manager.initialize_data()
    test_data_manager.set_comm(comm)
    b = test_data_manager.get_comm()
    @test comm == b
    # MPI.Finalize()
end

@testset "add_and_get_models" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.initialize_data()
    @test length(keys(test_data_manager.get_active_models())) == 0
    test_data_manager.add_active_model("Test", Test)
    test_data_manager.add_active_model("PeriLab", PeriLab)
    test_list = test_data_manager.get_active_models()
    @test collect(keys(test_list)) == ["Test", "PeriLab"]
    @test test_list["Test"] == Test
    @test test_list["PeriLab"] == PeriLab
    test_data_manager.add_active_model("PeriLab", PeriLab)
    test_list = test_data_manager.get_active_models()
    @test collect(keys(test_list)) == ["Test", "PeriLab"]
    @test test_list["Test"] == Test
    @test test_list["PeriLab"] == PeriLab
    test_data_manager.initialize_data()
    test_data_manager.add_active_model("PeriLab", PeriLab)
    test_data_manager.add_active_model("Test", Test)
    test_list = test_data_manager.get_active_models()
    @test collect(keys(test_list)) == ["PeriLab", "Test"]
end

@testset "ut_set_get_accuracy_order" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.initialize_data()
    @test test_data_manager.get_accuracy_order() == 1
    test_data_manager.set_accuracy_order(2)
    @test test_data_manager.get_accuracy_order() == 2
    test_data_manager.set_accuracy_order(5)
    @test test_data_manager.get_accuracy_order() == 5
    @test isnothing(test_data_manager.set_accuracy_order(0))
end
@testset "ranks" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.set_rank(2)
    test_data_manager.set_max_rank(3)
    @test test_data_manager.get_rank() == 2
    @test test_data_manager.get_max_rank() == 3
    test_data_manager.set_rank(3)
    @test test_data_manager.get_rank() == 3
end

@testset "get_local_nodes" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.set_glob_to_loc(Dict{Int64,Int64}(1 => 1, 3 => 2, 2 => 3))
    @test test_data_manager.get_local_nodes([1, 2, 3]) == [1, 3, 2]
    @test test_data_manager.get_local_nodes([1]) == [1]
    @test test_data_manager.get_local_nodes([2]) == [3]
    @test test_data_manager.get_local_nodes([3]) == [2]
    @test test_data_manager.get_local_nodes([4]) == []
    @test test_data_manager.get_local_nodes([4, 2, 3]) == [3, 2]
    test_data_manager.set_glob_to_loc(Dict{Int64,Int64}(3 => 1, 2 => 2, 4 => 3))
    @test test_data_manager.get_local_nodes([1]) == []
    @test test_data_manager.get_local_nodes([4]) == [3]
    @test test_data_manager.get_local_nodes([1, 4]) == [3]
end

@testset "get_set_functions" begin
    test_data_manager = PeriLab.Data_manager
    for i in 1:20
        test_data_manager.set_dof(i)
        @test test_data_manager.get_dof() == i
        test_data_manager.set_num_controller(i)
        @test test_data_manager.get_nnodes() == i
    end
    test_data_manager.set_num_controller(97)
    test_data_manager.set_num_responder(5)

    @test test_data_manager.get_num_responder() == 5
    @test test_data_manager.get_nnodes() == 97
    @test test_data_manager.data["num_responder"] == 5
    @test test_data_manager.data["num_controller"] == 97
    @test test_data_manager.data["nnodes"] == 102
    nnodes = test_data_manager.get_nnodes()
    nnodes = 3
    @test test_data_manager.get_nnodes() == 97
    @test nnodes == 3
end
test_data_manager = PeriLab.Data_manager
num_controller = 3
num_responder = 2
test_data_manager.set_num_controller(num_controller)
test_data_manager.set_num_responder(num_responder)
nn = test_data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
nn[1] = 2
nn[2] = 3
nn[3] = 2
nn[4] = 2
nn[5] = 5

test_data_manager.create_constant_node_field("A", Float64, 1)
B = test_data_manager.create_node_field("B", Bool, 1)
C = test_data_manager.create_constant_node_field("C", Float64, 4)
C[1, 2] = 4
test_data_manager.create_node_field("D", Int64, 7)
test_data_manager.create_constant_bond_field("F", Float64, 1)
test_data_manager.create_bond_field("G", Bool, 1)
test_data_manager.create_constant_bond_field("H", Float64, 4)
test_data_manager.create_bond_field("I", Int64, 7)
testfield_keys = test_data_manager.get_all_field_keys()
@testset "create data fields -> get all fields" begin
    @test test_data_manager.get_nnodes() == num_controller
    @test B[1] == test_data_manager.get_field("B", "N")
    @test B[2] == test_data_manager.get_field("B", "NP1")
    @test C == test_data_manager.get_field("C")
    @test C == test_data_manager.get_field("C", "Constant")
    @test isnothing(test_data_manager.get_field("C", "N"))
    @test isnothing(test_data_manager.get_field("C", "non_existing"))
    @test "A" in testfield_keys
    @test ("AN" in testfield_keys) == false
    @test ("ANP1" in testfield_keys) == false
    @test ("B" in testfield_keys) == false
    @test "BN" in testfield_keys
    @test "BNP1" in testfield_keys
    @test "C" in testfield_keys
    @test ("CN" in testfield_keys) == false
    @test ("CNP1" in testfield_keys) == false
    @test ("D" in testfield_keys) == false
    @test "DN" in testfield_keys
    @test "DNP1" in testfield_keys
    @test ("E" in testfield_keys) == false
    @test ("EN" in testfield_keys) == false
    @test ("ENP1" in testfield_keys) == false
    @test "F" in testfield_keys
    @test ("FN" in testfield_keys) == false
    @test ("FNP1" in testfield_keys) == false
    @test ("G" in testfield_keys) == false
    @test "GN" in testfield_keys
    @test "GNP1" in testfield_keys
    @test "H" in testfield_keys
    @test ("HN" in testfield_keys) == false
    @test ("HNP1" in testfield_keys) == false
    @test ("I" in testfield_keys) == false
    @test "IN" in testfield_keys
    @test "INP1" in testfield_keys
end
@testset "ut_set_fem" begin
    @test test_data_manager.fem_active() == false
    test_data_manager.set_fem(true)
    @test test_data_manager.fem_active() == true
    test_data_manager.set_fem(false)
    @test test_data_manager.fem_active() == false
end

@testset "ut_number_of_elements" begin
    test_data_manager.set_num_elements(5)
    @test test_data_manager.get_num_elements() == 5
    test_data_manager.set_num_elements(10)
    @test test_data_manager.get_num_elements() == 10
    test_data_manager.set_num_elements(1)
    @test test_data_manager.get_num_elements() == 1
    test_data_manager.set_num_elements(0)
    @test test_data_manager.get_num_elements() == 0
    @test isnothing(test_data_manager.set_num_elements(-1))
end
@testset "ut_create_existing_field" begin
    field1, field2 = test_data_manager.create_node_field("D", Int64, 3)
    testfield_keys = test_data_manager.get_all_field_keys()
    @test "DN" in testfield_keys
    @test "DNP1" in testfield_keys
    # because field is not overwritten, the dof value stay
    @test size(field1) == (num_controller + num_responder, 7)
    @test size(field2) == (num_controller + num_responder, 7)
end

@testset "get_field" begin
    A = test_data_manager.get_field("A")
    @test typeof(A[1]) == Float64
    @test length(A) == test_data_manager.data["nnodes"] == num_controller + num_responder

    B = test_data_manager.get_field("B", "N")
    @test typeof(B[1]) == Bool
    @test length(B) == test_data_manager.data["nnodes"] == num_controller + num_responder

    C = test_data_manager.get_field("C")
    @test typeof(C[1, 1]) == Float64
    @test length(C[:, 1]) ==
          test_data_manager.data["nnodes"] ==
          num_controller + num_responder
    @test length(C[1, :]) == 4

    D = test_data_manager.get_field("D", "NP1")
    @test typeof(D[1, 1]) == Int64
    @test length(D[:, 1]) ==
          test_data_manager.data["nnodes"] ==
          num_controller + num_responder
    @test length(D[1, :]) == 7

    F = test_data_manager.get_field("F")
    @test typeof(F[1][1][1]) == Float64
    @test length(F[:, 1]) == num_controller + num_responder
    @test length(F[1]) == nn[1]
    @test length(F[2]) == nn[2]
    @test length(F[3]) == nn[3]
    G = test_data_manager.get_field("G", "N")
    @test typeof(G[1, 1][1]) == Bool
    @test length(G[:, 1]) == num_controller + num_responder

    H = test_data_manager.get_field("H")
    @test typeof(H[1][1][1][1]) == Float64
    @test length(H[1][:, 1]) == nn[1]
    @test length(H[1][1]) == 4
    @test length(H[:][:]) == num_controller + num_responder

    I = test_data_manager.get_field("I", "NP1")
    @test typeof(I[1][1][1]) == Int64
    for i in 1:(num_controller + num_responder)
        @test length(I[i]) == nn[i]
    end
    @test length(I[1][1]) == 7
    # @test length(I[:][:][:]) == num_controller + num_responder

end

@testset "ut_get_field_type" begin
    @test test_data_manager.get_field_type("A") == Float64
    @test test_data_manager.get_field_type("DN") == Int64
    @test test_data_manager.get_field_type("DNP1") == Int64
    @test test_data_manager.get_field_type("GN") == Bool
    @test isnothing(test_data_manager.get_field_type("not there"))
    @test isnothing(test_data_manager.get_field_type("D"))
end

@testset "ut_create_free_size_field" begin
    test = test_data_manager.create_constant_free_size_field("BMatrix", Float64, (50, 3))
    @test size(test) == (50, 3)
    @test test_data_manager.get_field_type("BMatrix") == Float64
    test2 = test_data_manager.get_field("BMatrix")
    @test test == test2
    test = test_data_manager.create_constant_free_size_field("BMatrix", Float64, (2, 3))
    @test size(test) == (50, 3)
    test = test_data_manager.create_constant_free_size_field("GN", Float64, (2, 3))
    @info test
    @test size(test) == (5,)
    test = test_data_manager.create_constant_node_field("BMatrix", Float64, 3)
    @test size(test) == (50, 3)
    test = test_data_manager.create_constant_free_size_field("Test_size", Float64,
                                                             (2, 3, 3))
    @test size(test) == (2, 3, 3)
    test = test_data_manager.create_constant_free_size_field("Test_size_2",
                                                             Float64,
                                                             (2, 3, 3, 4))
    @test size(test) == (2, 3, 3, 4)
    test = test_data_manager.create_constant_node_field("Test_size_3", Float64, 3,
                                                        VectorOrMatrix = "Matrix")
    @test size(test) == (5, 3, 3)
    test = test_data_manager.create_constant_node_field("Test_size_3", Float64, 3,
                                                        VectorOrMatrix = "Matrix")
    @test size(test) == (5, 3, 3)
    test,
    test2 = test_data_manager.create_free_size_field("Test_size_4", Float64,
                                                     (3, 3, 1, 3))
    @test size(test) == (3, 3, 1, 3)
    @test size(test2) == (3, 3, 1, 3)
    @test "Test_size_4N" in test_data_manager.get_all_field_keys()
    @test "Test_size_4NP1" in test_data_manager.get_all_field_keys()
    test,
    test2 = test_data_manager.create_node_field("Test_size_4", Float64, 3,
                                                VectorOrMatrix = "Matrix")
    @test size(test) == (3, 3, 1, 3)
    @test size(test2) == (3, 3, 1, 3)
    test = test_data_manager.create_constant_free_size_field("Int8Matrix", Int64, (50, 3))
    @test typeof(test[1]) == Int64
end

@testset "set_get_field" begin
    nn = test_data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1] = 2
    nn[2] = 3
    nn[3] = 2
    nn[4] = 2
    nn[5] = 5
    test = test_data_manager.create_constant_node_field("test", Float64, 1)
    @test test == test_data_manager.get_field("test")
    @test test_data_manager.create_constant_node_field("test", Float64, 1) ==
          test_data_manager.get_field("test")
    test = test_data_manager.create_constant_node_field("test2", Float64, 3)
    @test test == test_data_manager.get_field("test2")
    test1, test2 = test_data_manager.create_node_field("test3", Float64, 1)
    @test test1 == test_data_manager.get_field("test3", "N")
    @test test2 == test_data_manager.get_field("test3", "NP1")
    test1, test2 = test_data_manager.create_node_field("test4", Float64, 3)
    @test test1 == test_data_manager.get_field("test4", "N")
    @test test2 == test_data_manager.get_field("test4", "NP1")
    test = test_data_manager.create_constant_bond_field("test5", Float64, 1)
    @test test == test_data_manager.get_field("test5")
    test = test_data_manager.create_constant_bond_field("test6", Float64, 3)
    @test test == test_data_manager.get_field("test6")
    test1, test2 = test_data_manager.create_bond_field("test7", Float64, 1)
    @test test1 == test_data_manager.get_field("test7", "N")
    @test test2 == test_data_manager.get_field("test7", "NP1")
    test1, test2 = test_data_manager.create_bond_field("test8", Float64, 3)
    @test test1 == test_data_manager.get_field("test8", "N")
    @test test2 == test_data_manager.get_field("test8", "NP1")
    # testnewFloat = test_data_manager.create_constant_node_field("testnewFloat", Float16, 1)
    # @test typeof(testnewFloat[1]) == Float16
    # testnewInt = test_data_manager.create_constant_node_field("testnewInt", Int8, 1)
    # @test typeof(testnewInt[1]) == Int8

    testDoesnotExists = test_data_manager.get_field("does not exist", "NP1")
    @test isnothing(testDoesnotExists)
    testDoesnotExists = test_data_manager.get_field("does not exist")
    @test isnothing(testDoesnotExists)
end

@testset "Matrix" begin
    #Arrays
    test = test_data_manager.create_constant_node_field("test9", Float64, 2,
                                                        VectorOrMatrix = "Matrix")
    test[1, 1, 1] = 1.2
    test[1, 2, 1] = -1.2
    test[1, 1, 2] = 1.4
    test[1, 2, 2] = 1.2
    @test test == test_data_manager.get_field("test9")
    test = test_data_manager.create_constant_bond_field("test10", Float64, 3,
                                                        VectorOrMatrix = "Matrix")
    test[1][1, 1, 1] = 1.2
    test[2][1, 2, 1] = -1.2
    test[2][1, 1, 3] = 1.4
    test[2][1, 2, 2] = 1.2
    @test test == test_data_manager.get_field("test10")
    test1,
    test2 = test_data_manager.create_bond_field("test11", Float64, 6,
                                                VectorOrMatrix = "Matrix")
    @test test1 == test_data_manager.get_field("test11", "N")
    @test test2 == test_data_manager.get_field("test11", "NP1")
    test1,
    test2 = test_data_manager.create_node_field("test12", Float64, 3,
                                                VectorOrMatrix = "Matrix")
    @test test1 == test_data_manager.get_field("test12", "N")
    @test test2 == test_data_manager.get_field("test12", "NP1")
end

@testset "get_NP1_to_N_Dict" begin
    @test test_data_manager.data["NP1_to_N"]["B"][1] == "BN"
    @test test_data_manager.data["NP1_to_N"]["D"][1] == "DN"
    @test test_data_manager.data["NP1_to_N"]["G"][1] == "GN"
    @test test_data_manager.data["NP1_to_N"]["I"][1] == "IN"
    @test test_data_manager.data["NP1_to_N"]["B"][2] == "BNP1"
    @test test_data_manager.data["NP1_to_N"]["D"][2] == "DNP1"
    @test test_data_manager.data["NP1_to_N"]["G"][2] == "GNP1"
    @test test_data_manager.data["NP1_to_N"]["I"][2] == "INP1"
end
@testset "set_and_get_values" begin
    DN = test_data_manager.get_field("D", "N")
    DN[1, 3] = 10
    DNtest = test_data_manager.get_field("D", "N")
    @test DN[1, 3] == DNtest[1, 3]
end

bdn, bdnp1 = test_data_manager.create_bond_field("Bond Damage", Float64, 1, 1)
test_data_manager.create_constant_node_field("Active", Bool, 1, true)
@testset "switch_NP1_to_N" begin
    bmatrixN,
    bmatrixNP1 = test_data_manager.create_bond_field("Bmat", Float64, 2,
                                                     VectorOrMatrix = "Matrix")
    nmatrixN,
    nmatrixNP1 = test_data_manager.create_node_field("Nmat", Float64, 2,
                                                     VectorOrMatrix = "Matrix")
    DN = test_data_manager.get_field("D", "N")
    DNP1 = test_data_manager.get_field("D", "NP1")

    IN = test_data_manager.get_field("I", "N")
    INP1 = test_data_manager.get_field("I", "NP1")
    IN[2][1][3] = 0

    DNP1[2, 3] = 5
    @test DN[2, 3] == 0
    @test DNP1[2, 3] == 5
    # bonds
    INP1[2][1][3] = 5
    @test IN[2][1][3] == 0
    @test INP1[2][1][3] == 5
    nmatrixNP1[1, 1, 1] = 2
    nmatrixNP1[1, 1, 2] = 2
    nmatrixNP1[1, 2, 1] = 3
    nmatrixNP1[1, 2, 2] = 4
    bmatrixNP1[1][1, 1, 1] = 2
    bmatrixNP1[1][2, 1, 2] = 2
    bmatrixNP1[1][1, 2, 1] = 3
    bmatrixNP1[1][2, 2, 2] = 4
    # extra test, because Bond Damage is set to one, to avoid unneccessary operations
    bdn = test_data_manager.get_field("Bond Damage", "N")
    bdnp1 = test_data_manager.get_field("Bond Damage", "NP1")
    bdnp1[2][1] = 0.5
    @test bdn[2][1] == 1.0
    @test bdnp1[2][1] == 0.5

    test_data_manager.switch_NP1_to_N()

    bdn = test_data_manager.get_field("Bond Damage", "N")
    bdnp1 = test_data_manager.get_field("Bond Damage", "NP1")
    @test bdn[2][1] == 0.5
    @test bdn[2][2] == 1.0
    @test bdnp1[2][1] == 0.5
    @test bdnp1[2][2] == 1.0

    DN = test_data_manager.get_field("D", "N")
    DNP1 = test_data_manager.get_field("D", "NP1")
    nmatrixN = test_data_manager.get_field("Nmat", "N")
    nmatrixNP1 = test_data_manager.get_field("Nmat", "NP1")
    @test DN[2, 3] == 5
    @test nmatrixN[1, 1, 1] == 2
    @test nmatrixN[1, 1, 2] == 2
    @test nmatrixN[1, 2, 1] == 3
    @test nmatrixN[1, 2, 2] == 4
    @test nmatrixNP1[1, 1, 1] == 0.0
    @test nmatrixNP1[1, 1, 2] == 0.0
    @test nmatrixNP1[1, 2, 1] == 0.0
    @test nmatrixNP1[1, 2, 2] == 0.0
    bmatrixN = test_data_manager.get_field("Bmat", "N")
    bmatrixNP1 = test_data_manager.get_field("Bmat", "NP1")
    @test bmatrixN[1][1, 1, 1] == 2
    @test bmatrixN[1][2, 1, 2] == 2
    @test bmatrixN[1][1, 2, 1] == 3
    @test bmatrixN[1][2, 2, 2] == 4
    @test bmatrixNP1[1][1, 1, 1] == 0.0
    @test bmatrixNP1[1][2, 1, 2] == 0.0
    @test bmatrixNP1[1][1, 2, 1] == 0.0
    @test bmatrixNP1[1][2, 2, 2] == 0.0
    # dependency test

    @test DNP1[2, 3] == 0
    DNP1[2, 3] = 6
    @test DNP1[2, 3] == 6
    @test DN[2, 3] == 5

    IN = test_data_manager.get_field("I", "N")
    INP1 = test_data_manager.get_field("I", "NP1")
    @test INP1[2][1][3] == 0
    # bonds
    @test IN[2][1][3] == 5
end

@testset "ut_nodesets" begin
    @test test_data_manager.get_nnsets() == 0
    test_data_manager.set_nset("N1", [1, 2])
    @test test_data_manager.get_nnsets() == 1
    test_data_manager.set_nset("N2", [4, 5])
    @test test_data_manager.get_nnsets() == 2
    test_data_manager.set_nset("N3", [1, 12, 22])
    @test test_data_manager.get_nnsets() == 3
    nsets = test_data_manager.get_nsets()
    @test nsets["N1"] == [1, 2]
    @test nsets["N2"] == [4, 5]
    @test nsets["N3"] == [1, 12, 22]
    test_data_manager.set_nset("N3", [1, 12])
    @test nsets["N3"] == [1, 12]
end

@testset "ut_block_name_list" begin
    test_data_manager.set_block_name_list(Vector{String}())
    block_name_list = test_data_manager.get_block_name_list()
    @test length(block_name_list) == 0
    test_data_manager.set_block_name_list(["1", "2", "3", "4"])
    block_name_list = test_data_manager.get_block_name_list()
    @test length(block_name_list) == 4
    @test block_name_list == ["1", "2", "3", "4"]
end

@testset "ut_properties" begin
    test_data_manager.set_block_id_list([2, 3, 1])
    test_data_manager.init_properties()
    @test length(test_data_manager.data["properties"]) == 3
    @test isnothing(test_data_manager.get_property(1, "Material Model", "E"))
    test_data_manager.set_property(1, "Material Model", "E", 3)
    test_data_manager.get_property(1, "Material Model", "E")
    @test test_data_manager.get_property(1, "Material Model", "E") == 3
    test_data_manager.set_property(1, "Material Model", "C", "Hello Test")
    @test test_data_manager.get_property(1, "Material Model", "C") == "Hello Test"
    test_data_manager.set_property(2, "Material Model", "E", 1.1)
    @test test_data_manager.get_property(2, "Material Model", "E") == 1.1
    test_data_manager.set_property(2, "Thermal Model", "E", [3 1 2; 1 2 3; 1 3 4])
    @test test_data_manager.get_property(2, "Thermal Model", "E") == [3 1 2; 1 2 3; 1 3 4]
    test_data_manager.set_property(3, "Thermal Model", "Q", 23.1)
    @test test_data_manager.get_property(3, "Thermal Model", "Q") == 23.1
    test_data_manager.set_property(3, "Damage Model", "SS", 0.1)
    @test test_data_manager.get_property(3, "Damage Model", "SS") == 0.1
    test_data_manager.set_property(1, "Additive Model", "E", [1, 2, 3])
    @test test_data_manager.get_property(1, "Additive Model", "E") == [1, 2, 3]
    test_data_manager.set_property(2, "Additive Model", "Qd", true)
    @test test_data_manager.get_property(2, "Additive Model", "Qd") == true
    @test isnothing(test_data_manager.get_property(2, "Additive Model", "not there"))
    @test test_data_manager.get_properties(1, "Material Model") ==
          Dict("C" => "Hello Test", "E" => 3)
    @test test_data_manager.get_properties(1, "Thermal Model") == Dict()
    @test test_data_manager.get_properties(2, "Material Model") == Dict("E" => 1.1)
    @test test_data_manager.get_properties(2, "Thermal Model") ==
          Dict("E" => [3 1 2; 1 2 3; 1 3 4])
    @test test_data_manager.get_properties(1, "") == Dict()
    @test !test_data_manager.check_property(1, "This is not a property")
    @test isnothing(test_data_manager.get_property(1, "Thermal Model",
                                                   "This is not a property"))
    test_data_manager.set_properties("FEM", Dict("A" => 2, "C" => "Model"))
    @test test_data_manager.get_properties(1, "FEM") == Dict("A" => 2, "C" => "Model")
    @test test_data_manager.get_properties(2, "FEM") == Dict("A" => 2, "C" => "Model")
    @test test_data_manager.get_properties(3, "FEM") == Dict("A" => 2, "C" => "Model")
end

@testset "ut_get_and_set_inverse_nlist" begin
    # inv_nlist = test_data_manager.get_inverse_nlist()
    # @test typeof(inv_nlist) == Vector{Dict{Int64,Int64}}
    # @test length(inv_nlist) == 0
    test_data_manager.set_inverse_nlist([
                                            Dict{Int64,Int64}(1 => 2),
                                            Dict{Int64,Int64}(1 => 2, 2 => 1)
                                        ])
    inv_nlist = test_data_manager.get_inverse_nlist()
    @test typeof(inv_nlist) == Vector{Dict{Int64,Int64}}
    @test length(inv_nlist) == 2
    @test inv_nlist[1] == Dict{Int64,Int64}(1 => 2)
    @test inv_nlist[2] == Dict{Int64,Int64}(1 => 2, 2 => 1)
    test_data_manager.set_inverse_nlist([Dict{Int64,Int64}(1 => 2, 2 => 1, 9 => 2)])
    inv_nlist = test_data_manager.get_inverse_nlist()
    @test length(inv_nlist) == 1
    @test inv_nlist[1] == Dict{Int64,Int64}(1 => 2, 2 => 1, 9 => 2)
end

@testset "ut_rotation" begin
    rotation = test_data_manager.get_rotation()
    angles = test_data_manager.get_field("Angles")
    @test !rotation
    @test isnothing(angles)
    test_angles = test_data_manager.create_constant_node_field("Angles", Float64, 3)
    test_data_manager.set_rotation(true)
    rotation = test_data_manager.get_rotation()
    angles = test_data_manager.get_field("Angles")
    @test rotation
    @test angles == test_angles
    rotation = test_data_manager.get_element_rotation()
    angles = test_data_manager.get_field("Element Angles")
    @test !rotation
    @test isnothing(angles)
    test_angles = test_data_manager.create_constant_node_field("Element Angles", Float64, 3)# in code it has length number of elements * element integration points
    test_data_manager.set_element_rotation(true)
    rotation = test_data_manager.get_element_rotation()
    angles = test_data_manager.get_field("Element Angles")
    @test rotation
    @test angles == test_angles
end

@testset "ut_cancel" begin
    test_data_manager.set_cancel(true)
    @test test_data_manager.get_cancel()
    test_data_manager.set_cancel(false)
    @test !test_data_manager.get_cancel()
end

@testset "ut_crit_values_matrix" begin
    crit_values = zeros(Float64, (2, 2, 2))
    crit_values[1, 1, 1] = 1.0
    test_data_manager.set_crit_values_matrix(crit_values)
    @test test_data_manager.get_crit_values_matrix() == crit_values
end

@testset "ut_initialize_data" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.create_node_field("test4", Float64, 3)
    test_data_manager.initialize_data()

    @test test_data_manager.get_nnodes() == 0
    @test test_data_manager.get_dof() == 2
    @test length(test_data_manager.get_all_field_keys()) == 0
end

@testset "get_field_alloc" begin
    test_data_manager.initialize_data()
    test_data_manager.set_dof(2)
    num_controller = 3
    num_responder = 2
    test_data_manager.set_num_controller(num_controller)
    test_data_manager.set_num_responder(num_responder)
    nn = test_data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1] = 2
    nn[2] = 3
    nn[3] = 2
    nn[4] = 2
    nn[5] = 5

    alloc = 0
    alloc += @allocated test_data_manager.create_constant_node_field("A", Float64, 1)
    alloc += @allocated test_data_manager.create_node_field("B", Bool, 1)
    alloc += @allocated test_data_manager.create_constant_node_field("C", Float64, 4)
    alloc += @allocated test_data_manager.create_node_field("D", Int64, 7)
    alloc += @allocated test_data_manager.create_constant_bond_field("E", Float64, 1)
    alloc += @allocated test_data_manager.create_bond_field("G", Bool, 1)
    alloc += @allocated test_data_manager.create_constant_bond_field("H", Float64, 4)
    alloc += @allocated test_data_manager.create_bond_field("I", Int64, 7)
    alloc += @allocated test_data_manager.create_constant_free_size_field("J", Float64,
                                                                          (2, 3))
    alloc += @allocated test_data_manager.create_constant_node_field("K", Float64, 3)
    alloc += @allocated test_data_manager.create_constant_free_size_field("L",
                                                                          Float64,
                                                                          (2, 3, 3))
    alloc += @allocated test_data_manager.create_constant_free_size_field("M",
                                                                          Float64,
                                                                          (2, 3, 3, 4))
    alloc += @allocated test_data_manager.create_constant_node_field("N", Float64,
                                                                     3,
                                                                     VectorOrMatrix = "Matrix")
    alloc += @allocated test_data_manager.create_constant_node_field("O", Float64,
                                                                     3,
                                                                     VectorOrMatrix = "Matrix")
    alloc += @allocated test_data_manager.create_free_size_field("P", Float64, (3, 3, 1, 3))
    alloc += @allocated test_data_manager.create_node_field("Q", Float64, 3,
                                                            VectorOrMatrix = "Matrix")
    alloc += @allocated test_data_manager.create_constant_free_size_field("R", Int64,
                                                                          (50, 3))
    alloc += @allocated test_data_manager.create_constant_bond_field("S", Float64,
                                                                     3,
                                                                     VectorOrMatrix = "Matrix")
    alloc += @allocated test_data_manager.create_bond_field("T", Float64, 3,
                                                            VectorOrMatrix = "Matrix")

    @test alloc < 10947441 # 1.3684 MB

    alloc = 0
    alloc += @allocated dof = test_data_manager.get_dof()
    alloc += @allocated a = test_data_manager.get_field("A")
    alloc += @allocated b = test_data_manager.get_field("B", "NP1")
    alloc += @allocated c = test_data_manager.get_field("C")
    alloc += @allocated d = test_data_manager.get_field("D", "NP1")
    alloc += @allocated e = test_data_manager.get_field("E")
    alloc += @allocated g = test_data_manager.get_field("G", "NP1")
    alloc += @allocated h = test_data_manager.get_field("H")
    alloc += @allocated i = test_data_manager.get_field("I", "NP1")
    alloc += @allocated j = test_data_manager.get_field("J")
    alloc += @allocated k = test_data_manager.get_field("K")
    alloc += @allocated l = test_data_manager.get_field("L")
    alloc += @allocated m = test_data_manager.get_field("M")
    alloc += @allocated n = test_data_manager.get_field("N")
    alloc += @allocated o = test_data_manager.get_field("O")
    alloc += @allocated p = test_data_manager.get_field("P", "NP1")
    alloc += @allocated q = test_data_manager.get_field("Q", "NP1")
    alloc += @allocated r = test_data_manager.get_field("R")
    alloc += @allocated s = test_data_manager.get_field("S")
    alloc += @allocated t = test_data_manager.get_field("T", "NP1")

    @test alloc == 0
end
