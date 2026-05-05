# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using MPI
#using Test

#using PeriLab

@testset "set_comm" begin
    # MPI.Init()
    comm = MPI.COMM_WORLD
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_comm(comm)
    b = PeriLab.Data_Manager.get_comm()
    @test comm == b
    # MPI.Finalize()
end

@testset "add_and_get_models" begin
    PeriLab.Data_Manager.initialize_data()
    @test length(keys(PeriLab.Data_Manager.get_active_models())) == 0
    PeriLab.Data_Manager.add_active_model("Test", Test)
    PeriLab.Data_Manager.add_active_model("PeriLab", PeriLab)
    test_list = PeriLab.Data_Manager.get_active_models()
    @test collect(keys(test_list)) == ["Test", "PeriLab"]
    @test test_list["Test"] == Test
    @test test_list["PeriLab"] == PeriLab
    PeriLab.Data_Manager.add_active_model("PeriLab", PeriLab)
    test_list = PeriLab.Data_Manager.get_active_models()
    @test collect(keys(test_list)) == ["Test", "PeriLab"]
    @test test_list["Test"] == Test
    @test test_list["PeriLab"] == PeriLab
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.add_active_model("PeriLab", PeriLab)
    PeriLab.Data_Manager.add_active_model("Test", Test)
    test_list = PeriLab.Data_Manager.get_active_models()
    @test collect(keys(test_list)) == ["PeriLab", "Test"]
end

@testset "ut_set_get_accuracy_order" begin
    PeriLab.Data_Manager.initialize_data()
    @test PeriLab.Data_Manager.get_accuracy_order() == 1
    PeriLab.Data_Manager.set_accuracy_order(2)
    @test PeriLab.Data_Manager.get_accuracy_order() == 2
    PeriLab.Data_Manager.set_accuracy_order(5)
    @test PeriLab.Data_Manager.get_accuracy_order() == 5
    @test_logs (:error, "Accuracy order must be greater than zero.") PeriLab.Data_Manager.set_accuracy_order(0)
end
@testset "ranks" begin
    PeriLab.Data_Manager.set_rank(2)
    PeriLab.Data_Manager.set_max_rank(3)
    @test PeriLab.Data_Manager.get_rank() == 2
    @test PeriLab.Data_Manager.get_max_rank() == 3
    PeriLab.Data_Manager.set_rank(3)
    @test PeriLab.Data_Manager.get_rank() == 3
end

@testset "get_local_nodes" begin
    PeriLab.Data_Manager.set_glob_to_loc(Dict{Int64,Int64}(1 => 1, 3 => 2, 2 => 3))
    @test PeriLab.Data_Manager.get_local_nodes([1, 2, 3]) == [1, 3, 2]
    @test PeriLab.Data_Manager.get_local_nodes([1]) == [1]
    @test PeriLab.Data_Manager.get_local_nodes([2]) == [3]
    @test PeriLab.Data_Manager.get_local_nodes([3]) == [2]
    @test PeriLab.Data_Manager.get_local_nodes([4]) == []
    @test PeriLab.Data_Manager.get_local_nodes([4, 2, 3]) == [3, 2]
    PeriLab.Data_Manager.set_glob_to_loc(Dict{Int64,Int64}(3 => 1, 2 => 2, 4 => 3))
    @test PeriLab.Data_Manager.get_local_nodes([1]) == []
    @test PeriLab.Data_Manager.get_local_nodes([4]) == [3]
    @test PeriLab.Data_Manager.get_local_nodes([1, 4]) == [3]
end

@testset "get_set_functions" begin
    ref_dof = [0, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
    for i in 1:15
        if i == 2 || i == 3
            PeriLab.Data_Manager.set_dof(i)
        else
            @test_logs (:error, "Degree of freedom $i is not supported.") PeriLab.Data_Manager.set_dof(i)
        end

        @test PeriLab.Data_Manager.get_dof() == ref_dof[i]
        PeriLab.Data_Manager.set_num_controller(i)
        @test PeriLab.Data_Manager.get_nnodes() == i
    end
    PeriLab.Data_Manager.set_num_controller(97)
    PeriLab.Data_Manager.set_num_responder(5)

    @test PeriLab.Data_Manager.get_num_responder() == 5
    @test PeriLab.Data_Manager.get_nnodes() == 97
    @test PeriLab.Data_Manager.data["num_responder"] == 5
    @test PeriLab.Data_Manager.data["num_controller"] == 97
    @test PeriLab.Data_Manager.data["nnodes"] == 102
    nnodes = PeriLab.Data_Manager.get_nnodes()
    nnodes = 3
    @test PeriLab.Data_Manager.get_nnodes() == 97
    @test nnodes == 3
end
num_controller = 3
num_responder = 2
PeriLab.Data_Manager.set_num_controller(num_controller)
PeriLab.Data_Manager.set_num_responder(num_responder)
nn = PeriLab.Data_Manager.create_constant_node_scalar_field("Number of Neighbors", Int64)
nn[1] = 2
nn[2] = 3
nn[3] = 2
nn[4] = 2
nn[5] = 5

PeriLab.Data_Manager.create_constant_node_scalar_field("A", Float64)
B = PeriLab.Data_Manager.create_node_scalar_field("B", Bool)
C = PeriLab.Data_Manager.create_constant_node_vector_field("C", Float64, 4)
C[1, 2] = 4
PeriLab.Data_Manager.create_node_vector_field("D", Int64, 7)
PeriLab.Data_Manager.create_constant_bond_scalar_state("F", Float64)
PeriLab.Data_Manager.create_bond_scalar_state("G", Bool)
PeriLab.Data_Manager.create_constant_bond_vector_state("H", Float64, 4)
PeriLab.Data_Manager.create_bond_vector_state("I", Int64, 7)
testfield_keys = PeriLab.Data_Manager.get_all_field_keys()
@testset "create data fields -> get all fields" begin
    @test PeriLab.Data_Manager.get_nnodes() == num_controller
    @test B[1] == PeriLab.Data_Manager.get_field("B", "N")
    @test B[2] == PeriLab.Data_Manager.get_field("B", "NP1")
    @test C == PeriLab.Data_Manager.get_field("C")
    @test C == PeriLab.Data_Manager.get_field("C", "Constant")
    @test_logs (:error,
                "Field ''C'' does not exist. Check if it is initialized as constant.") PeriLab.Data_Manager.get_field("C",
                                                                                                                      "N")
    @test_logs (:error, "Time non_existing is not supported. Use 'constant', 'N', or 'NP1'") PeriLab.Data_Manager.get_field("C",
                                                                                                                            "non_existing")
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
    @test PeriLab.Data_Manager.data["NP1_to_N"]["B"].N == "BN"
    @test PeriLab.Data_Manager.data["NP1_to_N"]["B"].NP1 == "BNP1"
    PeriLab.Data_Manager.set_NP1_to_N("B", Bool)
    @test PeriLab.Data_Manager.data["NP1_to_N"]["B"].N == "BN"
    @test PeriLab.Data_Manager.data["NP1_to_N"]["B"].NP1 == "BNP1"

    @test !haskey(PeriLab.Data_Manager.data["NP1_to_N"], "test_set_NP1_to_N")

    PeriLab.Data_Manager.set_NP1_to_N("test_set_NP1_to_N", Int64)

    @test haskey(PeriLab.Data_Manager.data["NP1_to_N"], "test_set_NP1_to_N")
    @test PeriLab.Data_Manager.data["NP1_to_N"]["test_set_NP1_to_N"].N ==
          "test_set_NP1_to_NN"
    @test PeriLab.Data_Manager.data["NP1_to_N"]["test_set_NP1_to_N"].N ==
          "test_set_NP1_to_NNP1"
end
@testset "ut_set_fem" begin
    @test PeriLab.Data_Manager.fem_active() == false
    PeriLab.Data_Manager.set_fem(true)
    @test PeriLab.Data_Manager.fem_active() == true
    PeriLab.Data_Manager.set_fem(false)
    @test PeriLab.Data_Manager.fem_active() == false
end

@testset "ut_set_verbose" begin
    PeriLab.Data_Manager.set_verbose(false)
    @test PeriLab.Data_Manager.get_verbose() == false
    PeriLab.Data_Manager.set_verbose(true)
    @test PeriLab.Data_Manager.get_verbose() == true
end

@testset "ut_number_of_elements" begin
    PeriLab.Data_Manager.set_num_elements(5)
    @test PeriLab.Data_Manager.get_num_elements() == 5
    PeriLab.Data_Manager.set_num_elements(10)
    @test PeriLab.Data_Manager.get_num_elements() == 10
    PeriLab.Data_Manager.set_num_elements(1)
    @test PeriLab.Data_Manager.get_num_elements() == 1
    PeriLab.Data_Manager.set_num_elements(0)
    @test PeriLab.Data_Manager.get_num_elements() == 0
    @test_logs (:error, "Number of elements must be positive or zero.") PeriLab.Data_Manager.set_num_elements(-1)
end
@testset "ut_create_existing_field" begin
    field1, field2 = PeriLab.Data_Manager.create_node_vector_field("D", Int64, 3)
    testfield_keys = PeriLab.Data_Manager.get_all_field_keys()
    @test "DN" in testfield_keys
    @test "DNP1" in testfield_keys
    # because field is not overwritten, the dof value stay
    @test size(field1) == (num_controller + num_responder, 7)
    @test size(field2) == (num_controller + num_responder, 7)
end

@testset "get_field" begin
    A = PeriLab.Data_Manager.get_field("A")
    @test typeof(A[1]) == Float64
    @test length(A) == PeriLab.Data_Manager.data["nnodes"] == num_controller + num_responder

    B = PeriLab.Data_Manager.get_field("B", "N")
    @test typeof(B[1]) == Bool
    @test length(B) == PeriLab.Data_Manager.data["nnodes"] == num_controller + num_responder

    C = PeriLab.Data_Manager.get_field("C")
    @test typeof(C[1, 1]) == Float64
    @test length(C[:, 1]) ==
          PeriLab.Data_Manager.data["nnodes"] ==
          num_controller + num_responder
    @test length(C[1, :]) == 4

    D = PeriLab.Data_Manager.get_field("D", "NP1")
    @test typeof(D[1, 1]) == Int64
    @test length(D[:, 1]) ==
          PeriLab.Data_Manager.data["nnodes"] ==
          num_controller + num_responder
    @test length(D[1, :]) == 7

    F = PeriLab.Data_Manager.get_field("F")
    @test typeof(F[1][1][1]) == Float64
    @test length(F[:, 1]) == num_controller + num_responder
    @test length(F[1]) == nn[1]
    @test length(F[2]) == nn[2]
    @test length(F[3]) == nn[3]
    G = PeriLab.Data_Manager.get_field("G", "N")
    @test typeof(G[1, 1][1]) == Bool
    @test length(G[:, 1]) == num_controller + num_responder

    H = PeriLab.Data_Manager.get_field("H")
    @test typeof(H[1][1][1][1]) == Float64
    @test length(H[1][:, 1]) == nn[1]
    @test length(H[1][1]) == 4
    @test length(H[:][:]) == num_controller + num_responder

    I = PeriLab.Data_Manager.get_field("I", "NP1")
    @test typeof(I[1][1][1]) == Int64
    for i in 1:(num_controller + num_responder)
        @test length(I[i]) == nn[i]
    end
    @test length(I[1][1]) == 7
    # @test length(I[:][:][:]) == num_controller + num_responder

end

@testset "ut_get_field_type" begin
    @test PeriLab.Data_Manager.get_field_type("A") == Float64
    @test PeriLab.Data_Manager.get_field_type("DN") == Int64
    @test PeriLab.Data_Manager.get_field_type("DNP1") == Int64
    @test PeriLab.Data_Manager.get_field_type("GN") == Vector{Bool}
    @test_logs (:error, "Field ''not there'' does not exist.") PeriLab.Data_Manager.get_field_type("not there")
    @test_logs (:error, "Field ''D'' does not exist.") PeriLab.Data_Manager.get_field_type("D")
end

@testset "ut_create_free_size_field" begin
    test = PeriLab.Data_Manager.create_constant_free_size_field("BMatrix", Float64, (50, 3))
    @test size(test) == (50, 3)
    @test PeriLab.Data_Manager.get_field_type("BMatrix") == Float64
    test2 = PeriLab.Data_Manager.get_field("BMatrix")
    @test test == test2
    test = PeriLab.Data_Manager.create_constant_free_size_field("BMatrix", Float64, (2, 3))
    @test size(test) == (50, 3)
    test = PeriLab.Data_Manager.create_constant_free_size_field("GN", Float64, (2, 3))
    @info test
    @test size(test) == (5,)
    test = PeriLab.Data_Manager.create_constant_node_vector_field("BMatrix", Float64, 3)
    @test size(test) == (50, 3)
    test = PeriLab.Data_Manager.create_constant_free_size_field("Test_size", Float64,
                                                                (2, 3, 3))
    @test size(test) == (2, 3, 3)
    test = PeriLab.Data_Manager.create_constant_free_size_field("Test_size_2",
                                                                Float64,
                                                                (2, 3, 3, 4))
    @test size(test) == (2, 3, 3, 4)
    test = PeriLab.Data_Manager.create_constant_node_tensor_field("Test_size_3", Float64, 3)
    @test size(test) == (5, 3, 3)
    test = PeriLab.Data_Manager.create_constant_node_tensor_field("Test_size_3", Float64, 3)
    @test size(test) == (5, 3, 3)
    test,
    test2 = PeriLab.Data_Manager.create_free_size_field("Test_size_4", Float64,
                                                        (3, 3, 1, 3))
    @test size(test) == (3, 3, 1, 3)
    @test size(test2) == (3, 3, 1, 3)
    @test "Test_size_4N" in PeriLab.Data_Manager.get_all_field_keys()
    @test "Test_size_4NP1" in PeriLab.Data_Manager.get_all_field_keys()
    test,
    test2 = PeriLab.Data_Manager.create_node_tensor_field("Test_size_4", Float64, 3)
    @test size(test) == (3, 3, 1, 3)
    @test size(test2) == (3, 3, 1, 3)
    test = PeriLab.Data_Manager.create_constant_free_size_field("Int8Matrix", Int64,
                                                                (50, 3))
    @test typeof(test[1]) == Int64
end

@testset "set_get_field" begin
    nn = PeriLab.Data_Manager.create_constant_node_scalar_field("Number of Neighbors",
                                                                Int64)
    nn[1] = 2
    nn[2] = 3
    nn[3] = 2
    nn[4] = 2
    nn[5] = 5
    test = PeriLab.Data_Manager.create_constant_node_scalar_field("test", Float64)
    @test test == PeriLab.Data_Manager.get_field("test")
    @test PeriLab.Data_Manager.create_constant_node_scalar_field("test", Float64) ==
          PeriLab.Data_Manager.get_field("test")
    test = PeriLab.Data_Manager.create_constant_node_vector_field("test2", Float64, 3)
    @test test == PeriLab.Data_Manager.get_field("test2")
    test1, test2 = PeriLab.Data_Manager.create_node_scalar_field("test3", Float64)
    @test test1 == PeriLab.Data_Manager.get_field("test3", "N")
    @test test2 == PeriLab.Data_Manager.get_field("test3", "NP1")
    test1, test2 = PeriLab.Data_Manager.create_node_vector_field("test4", Float64, 3)
    @test test1 == PeriLab.Data_Manager.get_field("test4", "N")
    @test test2 == PeriLab.Data_Manager.get_field("test4", "NP1")
    test = PeriLab.Data_Manager.create_constant_bond_scalar_state("test5", Float64)
    @test test == PeriLab.Data_Manager.get_field("test5")
    test = PeriLab.Data_Manager.create_constant_node_vector_field("test6", Float64, 3)
    @test test == PeriLab.Data_Manager.get_field("test6")
    test1, test2 = PeriLab.Data_Manager.create_bond_scalar_state("test7", Float64)
    @test test1 == PeriLab.Data_Manager.get_field("test7", "N")
    @test test2 == PeriLab.Data_Manager.get_field("test7", "NP1")
    test1, test2 = PeriLab.Data_Manager.create_bond_vector_state("test8", Float64, 3)
    @test test1 == PeriLab.Data_Manager.get_field("test8", "N")
    @test test2 == PeriLab.Data_Manager.get_field("test8", "NP1")
    # testnewFloat = PeriLab.Data_Manager.create_constant_node_field("testnewFloat", Float16, 1)
    # @test typeof(testnewFloat[1]) == Float16
    # testnewInt = PeriLab.Data_Manager.create_constant_node_field("testnewInt", Int8, 1)
    # @test typeof(testnewInt[1]) == Int8

    @test_logs (:error,
                "Field ''does not exist'' does not exist. Check if it is initialized as constant.") PeriLab.Data_Manager.get_field("does not exist",
                                                                                                                                   "NP1")
    @test_logs (:error,
                "Field ''does not exist'' does not exist. \n - Check if it is initialized as non-constant. \n - Check if the model is not activated in the solver options, e.g. Pre Calculation Models: False") PeriLab.Data_Manager.get_field("does not exist")
end

@testset "Matrix" begin
    #Arrays
    test = PeriLab.Data_Manager.create_constant_node_tensor_field("test9", Float64, 2)
    test[1, 1, 1] = 1.2
    test[1, 2, 1] = -1.2
    test[1, 1, 2] = 1.4
    test[1, 2, 2] = 1.2
    @test test == PeriLab.Data_Manager.get_field("test9")
    test = PeriLab.Data_Manager.create_constant_bond_tensor_state("test10", Float64, 3)
    test[1][1, 1, 1] = 1.2
    test[2][1, 2, 1] = -1.2
    test[2][1, 1, 3] = 1.4
    test[2][1, 2, 2] = 1.2
    @test test == PeriLab.Data_Manager.get_field("test10")
    test1,
    test2 = PeriLab.Data_Manager.create_bond_tensor_state("test11", Float64, 6)
    @test test1 == PeriLab.Data_Manager.get_field("test11", "N")
    @test test2 == PeriLab.Data_Manager.get_field("test11", "NP1")
    test1,
    test2 = PeriLab.Data_Manager.create_node_tensor_field("test12", Float64, 3)
    @test test1 == PeriLab.Data_Manager.get_field("test12", "N")
    @test test2 == PeriLab.Data_Manager.get_field("test12", "NP1")
end

@testset "get_NP1_to_N_Dict" begin
    @test PeriLab.Data_Manager.data["NP1_to_N"]["B"].N == "BN"
    @test PeriLab.Data_Manager.data["NP1_to_N"]["D"].N == "DN"
    @test PeriLab.Data_Manager.data["NP1_to_N"]["G"].N == "GN"
    @test PeriLab.Data_Manager.data["NP1_to_N"]["I"].N == "IN"
    @test PeriLab.Data_Manager.data["NP1_to_N"]["B"].NP1 == "BNP1"
    @test PeriLab.Data_Manager.data["NP1_to_N"]["D"].NP1 == "DNP1"
    @test PeriLab.Data_Manager.data["NP1_to_N"]["G"].NP1 == "GNP1"
    @test PeriLab.Data_Manager.data["NP1_to_N"]["I"].NP1 == "INP1"
end
@testset "set_and_get_values" begin
    DN = PeriLab.Data_Manager.get_field("D", "N")
    DN[1, 3] = 10
    DNtest = PeriLab.Data_Manager.get_field("D", "N")
    @test DN[1, 3] == DNtest[1, 3]
end

bdn,
bdnp1 = PeriLab.Data_Manager.create_bond_scalar_state("Bond Damage", Float64;
                                                      default_value = 1)
PeriLab.Data_Manager.create_constant_node_scalar_field("Active", Bool; default_value = true)
@testset "switch_NP1_to_N" begin
    bmatrixN,
    bmatrixNP1 = PeriLab.Data_Manager.create_bond_tensor_state("Bmat", Float64, 2)
    nmatrixN,
    nmatrixNP1 = PeriLab.Data_Manager.create_node_tensor_field("Nmat", Float64, 2)
    DN = PeriLab.Data_Manager.get_field("D", "N")
    DNP1 = PeriLab.Data_Manager.get_field("D", "NP1")

    IN = PeriLab.Data_Manager.get_field("I", "N")
    INP1 = PeriLab.Data_Manager.get_field("I", "NP1")
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
    bdn = PeriLab.Data_Manager.get_field("Bond Damage", "N")
    bdnp1 = PeriLab.Data_Manager.get_field("Bond Damage", "NP1")
    bdnp1[2][1] = 0.5
    @test bdn[2][1] == 1.0
    @test bdnp1[2][1] == 0.5

    PeriLab.Data_Manager.switch_NP1_to_N()

    bdn = PeriLab.Data_Manager.get_field("Bond Damage", "N")
    bdnp1 = PeriLab.Data_Manager.get_field("Bond Damage", "NP1")
    @test bdn[2][1] == 0.5
    @test bdn[2][2] == 1.0
    @test bdnp1[2][1] == 0.5
    @test bdnp1[2][2] == 1.0

    DN = PeriLab.Data_Manager.get_field("D", "N")
    DNP1 = PeriLab.Data_Manager.get_field("D", "NP1")
    nmatrixN = PeriLab.Data_Manager.get_field("Nmat", "N")
    nmatrixNP1 = PeriLab.Data_Manager.get_field("Nmat", "NP1")
    @test DN[2, 3] == 5
    @test nmatrixN[1, 1, 1] == 2
    @test nmatrixN[1, 1, 2] == 2
    @test nmatrixN[1, 2, 1] == 3
    @test nmatrixN[1, 2, 2] == 4
    @test nmatrixNP1[1, 1, 1] == 0.0
    @test nmatrixNP1[1, 1, 2] == 0.0
    @test nmatrixNP1[1, 2, 1] == 0.0
    @test nmatrixNP1[1, 2, 2] == 0.0
    bmatrixN = PeriLab.Data_Manager.get_field("Bmat", "N")
    bmatrixNP1 = PeriLab.Data_Manager.get_field("Bmat", "NP1")
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

    IN = PeriLab.Data_Manager.get_field("I", "N")
    INP1 = PeriLab.Data_Manager.get_field("I", "NP1")
    @test INP1[2][1][3] == 0
    # bonds
    @test IN[2][1][3] == 5
end

@testset "ut_nodesets" begin
    @test PeriLab.Data_Manager.get_nnsets() == 0
    PeriLab.Data_Manager.set_nset("N1", [1, 2])
    @test PeriLab.Data_Manager.get_nnsets() == 1
    PeriLab.Data_Manager.set_nset("N2", [4, 5])
    @test PeriLab.Data_Manager.get_nnsets() == 2
    PeriLab.Data_Manager.set_nset("N3", [1, 12, 22])
    @test PeriLab.Data_Manager.get_nnsets() == 3
    nsets = PeriLab.Data_Manager.get_nsets()
    @test nsets["N1"] == [1, 2]
    @test nsets["N2"] == [4, 5]
    @test nsets["N3"] == [1, 12, 22]
    PeriLab.Data_Manager.set_nset("N3", [1, 12])
    @test nsets["N3"] == [1, 12]
end

@testset "ut_block_name_list" begin
    PeriLab.Data_Manager.set_block_name_list(Vector{String}())
    block_name_list = PeriLab.Data_Manager.get_block_name_list()
    @test length(block_name_list) == 0
    PeriLab.Data_Manager.set_block_name_list(["1", "2", "3", "4"])
    block_name_list = PeriLab.Data_Manager.get_block_name_list()
    @test length(block_name_list) == 4
    @test block_name_list == ["1", "2", "3", "4"]
end

@testset "ut_properties" begin
    PeriLab.Data_Manager.set_block_id_list([2, 3, 1])
    PeriLab.Data_Manager.init_properties()
    @test length(PeriLab.Data_Manager.data["properties"]) == 3
    @test isnothing(PeriLab.Data_Manager.get_property(1, "Material Model", "E"))
    PeriLab.Data_Manager.set_property(1, "Material Model", "E", 3)
    PeriLab.Data_Manager.get_property(1, "Material Model", "E")
    @test PeriLab.Data_Manager.get_property(1, "Material Model", "E") == 3
    PeriLab.Data_Manager.set_property(1, "Material Model", "C", "Hello Test")
    @test PeriLab.Data_Manager.get_property(1, "Material Model", "C") == "Hello Test"
    PeriLab.Data_Manager.set_property(2, "Material Model", "E", 1.1)
    @test PeriLab.Data_Manager.get_property(2, "Material Model", "E") == 1.1
    PeriLab.Data_Manager.set_property(2, "Thermal Model", "E", [3 1 2; 1 2 3; 1 3 4])
    @test PeriLab.Data_Manager.get_property(2, "Thermal Model", "E") ==
          [3 1 2; 1 2 3; 1 3 4]
    PeriLab.Data_Manager.set_property(3, "Thermal Model", "Q", 23.1)
    @test PeriLab.Data_Manager.get_property(3, "Thermal Model", "Q") == 23.1
    PeriLab.Data_Manager.set_property(3, "Damage Model", "SS", 0.1)
    @test PeriLab.Data_Manager.get_property(3, "Damage Model", "SS") == 0.1
    PeriLab.Data_Manager.set_property(1, "Additive Model", "E", [1, 2, 3])
    @test PeriLab.Data_Manager.get_property(1, "Additive Model", "E") == [1, 2, 3]
    PeriLab.Data_Manager.set_property(2, "Additive Model", "Qd", true)
    @test PeriLab.Data_Manager.get_property(2, "Additive Model", "Qd") == true
    @test isnothing(PeriLab.Data_Manager.get_property(2, "Additive Model", "not there"))
    @test PeriLab.Data_Manager.get_properties(1, "Material Model") ==
          Dict("C" => "Hello Test", "E" => 3)
    @test PeriLab.Data_Manager.get_properties(1, "Thermal Model") == Dict()
    @test PeriLab.Data_Manager.get_properties(2, "Material Model") == Dict("E" => 1.1)
    @test PeriLab.Data_Manager.get_properties(2, "Thermal Model") ==
          Dict("E" => [3 1 2; 1 2 3; 1 3 4])
    @test PeriLab.Data_Manager.get_properties(1, "") == Dict()
    @test !PeriLab.Data_Manager.check_property(1, "This is not a property")
    @test isnothing(PeriLab.Data_Manager.get_property(1, "Thermal Model",
                                                      "This is not a property"))
    PeriLab.Data_Manager.set_properties("FEM", Dict("A" => 2, "C" => "Model"))
    @test PeriLab.Data_Manager.get_properties(1, "FEM") == Dict("A" => 2, "C" => "Model")
    @test PeriLab.Data_Manager.get_properties(2, "FEM") == Dict("A" => 2, "C" => "Model")
    @test PeriLab.Data_Manager.get_properties(3, "FEM") == Dict("A" => 2, "C" => "Model")
end

@testset "ut_get_and_set_inverse_nlist" begin
    # inv_nlist = PeriLab.Data_Manager.get_inverse_nlist()
    # @test typeof(inv_nlist) == Vector{Dict{Int64,Int64}}
    # @test length(inv_nlist) == 0
    PeriLab.Data_Manager.set_inverse_nlist([
                                               Dict{Int64,Int64}(1 => 2),
                                               Dict{Int64,Int64}(1 => 2, 2 => 1)
                                           ])
    inv_nlist = PeriLab.Data_Manager.get_inverse_nlist()
    @test typeof(inv_nlist) == Vector{Dict{Int64,Int64}}
    @test length(inv_nlist) == 2
    @test inv_nlist[1] == Dict{Int64,Int64}(1 => 2)
    @test inv_nlist[2] == Dict{Int64,Int64}(1 => 2, 2 => 1)
    PeriLab.Data_Manager.set_inverse_nlist([Dict{Int64,Int64}(1 => 2, 2 => 1, 9 => 2)])
    inv_nlist = PeriLab.Data_Manager.get_inverse_nlist()
    @test length(inv_nlist) == 1
    @test inv_nlist[1] == Dict{Int64,Int64}(1 => 2, 2 => 1, 9 => 2)
end

@testset "ut_rotation" begin
    rotation = PeriLab.Data_Manager.get_rotation()
    @test !rotation
    @test_logs (:error,
                "Field ''Angles'' does not exist. \n - Check if it is initialized as non-constant. \n - Check if the model is not activated in the solver options, e.g. Pre Calculation Models: False") PeriLab.Data_Manager.get_field("Angles")
    test_angles = PeriLab.Data_Manager.create_constant_node_vector_field("Angles", Float64,
                                                                         3)
    PeriLab.Data_Manager.set_rotation(true)
    rotation = PeriLab.Data_Manager.get_rotation()
    angles = PeriLab.Data_Manager.get_field("Angles")
    @test rotation
    @test angles == test_angles
    rotation = PeriLab.Data_Manager.get_element_rotation()
    @test !rotation
    @test_logs (:error,
                "Field ''Element Angles'' does not exist. \n - Check if it is initialized as non-constant. \n - Check if the model is not activated in the solver options, e.g. Pre Calculation Models: False") PeriLab.Data_Manager.get_field("Element Angles")
    test_angles = PeriLab.Data_Manager.create_constant_node_vector_field("Element Angles",
                                                                         Float64, 3)# in code it has length number of elements * element integration points
    PeriLab.Data_Manager.set_element_rotation(true)
    rotation = PeriLab.Data_Manager.get_element_rotation()
    angles = PeriLab.Data_Manager.get_field("Element Angles")
    @test rotation
    @test angles == test_angles
end

@testset "ut_cancel" begin
    PeriLab.Data_Manager.set_cancel(true)
    @test PeriLab.Data_Manager.get_cancel()
    PeriLab.Data_Manager.set_cancel(false)
    @test !PeriLab.Data_Manager.get_cancel()
end

@testset "ut_crit_values_matrix" begin
    crit_values = zeros(Float64, (2, 2, 2))
    crit_values[1, 1, 1] = 1.0
    PeriLab.Data_Manager.set_crit_values_matrix(crit_values)
    @test PeriLab.Data_Manager.get_crit_values_matrix() == crit_values
end

@testset "ut_initialize_data" begin
    PeriLab.Data_Manager.create_node_vector_field("test4", Float64, 3)
    PeriLab.Data_Manager.initialize_data()

    @test PeriLab.Data_Manager.get_nnodes() == 0
    @test PeriLab.Data_Manager.get_dof() == 0
    @test length(PeriLab.Data_Manager.get_all_field_keys()) == 0
end

@testset "get_field_alloc" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_dof(2)
    num_controller = 3
    num_responder = 2
    PeriLab.Data_Manager.set_num_controller(num_controller)
    PeriLab.Data_Manager.set_num_responder(num_responder)
    nn = PeriLab.Data_Manager.create_constant_node_scalar_field("Number of Neighbors",
                                                                Int64)
    nn[1] = 2
    nn[2] = 3
    nn[3] = 2
    nn[4] = 2
    nn[5] = 5

    alloc = 0
    alloc += @allocated PeriLab.Data_Manager.create_constant_node_scalar_field("A", Float64)
    alloc += @allocated PeriLab.Data_Manager.create_node_scalar_field("B", Bool)
    alloc += @allocated PeriLab.Data_Manager.create_constant_node_vector_field("C", Float64,
                                                                               4)
    alloc += @allocated PeriLab.Data_Manager.create_node_vector_field("D", Int64, 7)
    alloc += @allocated PeriLab.Data_Manager.create_constant_bond_scalar_state("E", Float64)
    alloc += @allocated PeriLab.Data_Manager.create_bond_scalar_state("G", Bool)
    alloc += @allocated PeriLab.Data_Manager.create_constant_bond_vector_state("H", Float64,
                                                                               4)
    alloc += @allocated PeriLab.Data_Manager.create_bond_vector_state("I", Int64, 7)
    alloc += @allocated PeriLab.Data_Manager.create_constant_free_size_field("J", Float64,
                                                                             (2, 3))
    alloc += @allocated PeriLab.Data_Manager.create_constant_node_vector_field("K", Float64,
                                                                               3)
    alloc += @allocated PeriLab.Data_Manager.create_constant_free_size_field("L",
                                                                             Float64,
                                                                             (2, 3, 3))
    alloc += @allocated PeriLab.Data_Manager.create_constant_free_size_field("M",
                                                                             Float64,
                                                                             (2, 3, 3, 4))
    alloc += @allocated PeriLab.Data_Manager.create_constant_node_tensor_field("N", Float64,
                                                                               3)
    alloc += @allocated PeriLab.Data_Manager.create_constant_node_tensor_field("O", Float64,
                                                                               3)
    alloc += @allocated PeriLab.Data_Manager.create_free_size_field("P", Float64,
                                                                    (3, 3, 1, 3))
    alloc += @allocated PeriLab.Data_Manager.create_node_tensor_field("Q", Float64, 3)
    alloc += @allocated PeriLab.Data_Manager.create_constant_free_size_field("R", Int64,
                                                                             (50, 3))
    alloc += @allocated PeriLab.Data_Manager.create_constant_bond_tensor_state("S", Float64,
                                                                               3)
    alloc += @allocated PeriLab.Data_Manager.create_bond_tensor_state("T", Float64, 3)

    @test alloc < 10947441 # 1.3684 MB

    alloc = 0
    alloc += @allocated dof = PeriLab.Data_Manager.get_dof()
    alloc += @allocated a = PeriLab.Data_Manager.get_field("A")
    alloc += @allocated b = PeriLab.Data_Manager.get_field("B", "NP1")
    alloc += @allocated c = PeriLab.Data_Manager.get_field("C")
    alloc += @allocated d = PeriLab.Data_Manager.get_field("D", "NP1")
    alloc += @allocated e = PeriLab.Data_Manager.get_field("E")
    alloc += @allocated g = PeriLab.Data_Manager.get_field("G", "NP1")
    alloc += @allocated h = PeriLab.Data_Manager.get_field("H")
    alloc += @allocated i = PeriLab.Data_Manager.get_field("I", "NP1")
    alloc += @allocated j = PeriLab.Data_Manager.get_field("J")
    alloc += @allocated k = PeriLab.Data_Manager.get_field("K")
    alloc += @allocated l = PeriLab.Data_Manager.get_field("L")
    alloc += @allocated m = PeriLab.Data_Manager.get_field("M")
    alloc += @allocated n = PeriLab.Data_Manager.get_field("N")
    alloc += @allocated o = PeriLab.Data_Manager.get_field("O")
    alloc += @allocated p = PeriLab.Data_Manager.get_field("P", "NP1")
    alloc += @allocated q = PeriLab.Data_Manager.get_field("Q", "NP1")
    alloc += @allocated r = PeriLab.Data_Manager.get_field("R")
    alloc += @allocated s = PeriLab.Data_Manager.get_field("S")
    alloc += @allocated t = PeriLab.Data_Manager.get_field("T", "NP1")

    @test alloc == 0
end
