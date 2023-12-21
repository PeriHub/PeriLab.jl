# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../../src/Support/data_manager.jl")
using MPI
using Test
@testset "set_comm" begin
    # MPI.Init()
    comm = MPI.COMM_WORLD
    test_Data_manager = Data_manager
    test_Data_manager.set_comm(comm)
    b = test_Data_manager.get_comm()
    @test comm == b
    # MPI.Finalize()
end

@testset "ranks" begin
    test_Data_manager = Data_manager
    test_Data_manager.set_rank(2)
    test_Data_manager.set_max_rank(3)
    @test test_Data_manager.get_rank() == 2
    @test test_Data_manager.get_max_rank() == 3
    test_Data_manager.set_rank(3)
    @test test_Data_manager.get_rank() == 3
end

@testset "get_local_nodes" begin
    test_Data_manager = Data_manager
    test_Data_manager.set_glob_to_loc(Dict{Int64,Int64}(1 => 1, 3 => 2, 2 => 3))
    @test test_Data_manager.get_local_nodes([1, 2, 3]) == [1, 3, 2]
    @test test_Data_manager.get_local_nodes([1]) == [1]
    @test test_Data_manager.get_local_nodes([2]) == [3]
    @test test_Data_manager.get_local_nodes([3]) == [2]
    @test test_Data_manager.get_local_nodes([4]) == []
    @test test_Data_manager.get_local_nodes([4, 2, 3]) == [3, 2]
    test_Data_manager.set_glob_to_loc(Dict{Int64,Int64}(3 => 1, 2 => 2, 4 => 3))
    @test test_Data_manager.get_local_nodes([1]) == []
    @test test_Data_manager.get_local_nodes([4]) == [3]
    @test test_Data_manager.get_local_nodes([1, 4]) == [3]
end


@testset "get_set_functions" begin
    test_Data_manager = Data_manager
    for i in 1:20
        test_Data_manager.set_dof(i)
        @test test_Data_manager.get_dof() == i
        test_Data_manager.set_num_controller(i)
        @test test_Data_manager.get_nnodes() == i
    end
    test_Data_manager.set_num_controller(97)
    test_Data_manager.set_num_responder(5)

    @test test_Data_manager.get_num_responder() == 5
    @test test_Data_manager.get_nnodes() == 97
    @test test_Data_manager.num_responder == 5
    @test test_Data_manager.num_controller == 97
    @test test_Data_manager.nnodes == 102
    nnodes = test_Data_manager.get_nnodes()
    nnodes = 3
    @test test_Data_manager.get_nnodes() == 97
    @test nnodes == 3
end
test_Data_manager = Data_manager
num_controller = 3
num_responder = 2
test_Data_manager.set_num_controller(num_controller)
test_Data_manager.set_num_responder(num_responder)
nn = test_Data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
nn[1] = 2
nn[2] = 3
nn[3] = 2
nn[4] = 2
nn[5] = 5

test_Data_manager.create_constant_node_field("A", Float64, 1)
B = test_Data_manager.create_node_field("B", Bool, 1)
C = test_Data_manager.create_constant_node_field("C", Float64, 4)
C[1, 2] = 4
test_Data_manager.create_node_field("D", Int64, 7)
test_Data_manager.create_constant_bond_field("F", Float64, 1)
test_Data_manager.create_bond_field("G", Bool, 1)
test_Data_manager.create_constant_bond_field("H", Float64, 4)
test_Data_manager.create_bond_field("I", Int64, 7)
testfield_keys = test_Data_manager.get_all_field_keys()
@testset "create data fields -> get all fields" begin
    @test test_Data_manager.get_nnodes() == num_controller
    @test B[1] == test_Data_manager.get_field("BN")
    @test B[2] == test_Data_manager.get_field("B", "NP1")
    @test C == test_Data_manager.get_field("C")
    @test C == test_Data_manager.get_field("C", "CONSTANT")
    @test C == test_Data_manager.get_field("C", "Constant")
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
    @test test_Data_manager.fem_active() == false
    test_Data_manager.set_fem(true)
    @test test_Data_manager.fem_active() == true
    test_Data_manager.set_fem(false)
    @test test_Data_manager.fem_active() == false
end


@testset "ut_number_of_elements" begin
    test_Data_manager.set_num_elements(5)
    @test test_Data_manager.get_num_elements() == 5
    test_Data_manager.set_num_elements(10)
    @test test_Data_manager.get_num_elements() == 10
    test_Data_manager.set_num_elements(1)
    @test test_Data_manager.get_num_elements() == 1
    test_Data_manager.set_num_elements(0)
    @test test_Data_manager.get_num_elements() == 0
    @test isnothing(test_Data_manager.set_num_elements(-1))
end
@testset "ut_create_existing_field" begin
    field1, field2 = test_Data_manager.create_node_field("D", Int64, 3)
    testfield_keys = test_Data_manager.get_all_field_keys()
    @test "DN" in testfield_keys
    @test "DNP1" in testfield_keys
    # because field is not overwritten, the dof value stay
    @test size(field1) == (num_controller + num_responder, 7)
    @test size(field2) == (num_controller + num_responder, 7)
end

@testset "get_field" begin

    A = test_Data_manager.get_field("A")
    @test typeof(A[1]) == Float64
    @test length(A) == test_Data_manager.nnodes == num_controller + num_responder

    B = test_Data_manager.get_field("BN")
    @test typeof(B[1]) == Bool
    @test length(B) == test_Data_manager.nnodes == num_controller + num_responder

    C = test_Data_manager.get_field("C")
    @test typeof(C[1, 1]) == Float64
    @test length(C[:, 1]) == test_Data_manager.nnodes == num_controller + num_responder
    @test length(C[1, :]) == 4

    D = test_Data_manager.get_field("DNP1")
    @test typeof(D[1, 1]) == Int64
    @test length(D[:, 1]) == test_Data_manager.nnodes == num_controller + num_responder
    @test length(D[1, :]) == 7

    F = test_Data_manager.get_field("F")
    @test typeof(F[1, 1][1]) == Float64
    @test length(F[:, 1]) == num_controller + num_responder
    @test length(F[1]) == nn[1]
    @test length(F[2]) == nn[2]
    @test length(F[3]) == nn[3]
    G = test_Data_manager.get_field("GN")
    @test typeof(G[1, 1][1]) == Bool
    @test length(G[:, 1]) == num_controller + num_responder

    H = test_Data_manager.get_field("H")
    @test typeof(H[1][1, 1][1]) == Float64
    @test length(H[1][:, 1]) == nn[1]
    @test length(H[1][1, :]) == 4
    @test length(H[:][:, :]) == num_controller + num_responder

    I = test_Data_manager.get_field("INP1")
    @test typeof(I[1][1, 1]) == Int64
    for i in 1:num_controller+num_responder
        @test length(I[i][:, 1]) == nn[i]
    end
    @test length(I[1][1, :]) == 7
    @test length(I[:][:, :]) == num_controller + num_responder

end

@testset "ut_get_field_type" begin
    @test test_Data_manager.get_field_type("A") == Float64
    @test test_Data_manager.get_field_type("DN") == Int64
    @test test_Data_manager.get_field_type("DNP1") == Int64
    @test test_Data_manager.get_field_type("GN") == Bool
    @test isnothing(test_Data_manager.get_field_type("not there"))
    @test isnothing(test_Data_manager.get_field_type("D"))
end

@testset "ut_create_free_size_field" begin
    test = test_Data_manager.create_constant_free_size_field("BMatrix", Float64, (50, 3))
    @test size(test) == (50, 3)
    @test test_Data_manager.get_field_type("BMatrix") == Float64
    test2 = test_Data_manager.get_field("BMatrix")
    @test test == test2
    test = test_Data_manager.create_constant_free_size_field("BMatrix", Float64, (2, 3))
    @test size(test) == (50, 3)
    test = test_Data_manager.create_constant_free_size_field("GN", Float64, (2, 3))
    @test size(test) == (5,)
    test = test_Data_manager.create_constant_node_field("BMatrix", Float64, 3)
    @test size(test) == (50, 3)
    test = test_Data_manager.create_constant_free_size_field("Test_size", Float64, (2, 3, 3))
    @test size(test) == (2, 3, 3)
    test = test_Data_manager.create_constant_free_size_field("Test_size_2", Float64, (2, 3, 3, 4))
    @test size(test) == (2, 3, 3, 4)
    test = test_Data_manager.create_constant_node_field("Test_size_3", Float64, "Matrix", 3)
    @test size(test) == (5, 3, 3)
    test = test_Data_manager.create_constant_node_field("Test_size_3", Float64, "Matrix", 3)
    @test size(test) == (5, 3, 3)
    test, test2 = test_Data_manager.create_free_size_field("Test_size_4", Float64, (3, 3, 1, 3))
    @test size(test) == (3, 3, 1, 3)
    @test size(test2) == (3, 3, 1, 3)
    @test "Test_size_4N" in test_Data_manager.get_all_field_keys()
    @test "Test_size_4NP1" in test_Data_manager.get_all_field_keys()
    test, test2 = test_Data_manager.create_node_field("Test_size_4", Float64, "Matrix", 3)
    @test size(test) == (3, 3, 1, 3)
    @test size(test2) == (3, 3, 1, 3)
end

function create_constant_free_size_field(name::String, type::Type, dof::Tuple)
    if haskey(fields, vartype) == false
        fields[vartype] = Dict{String,Any}()
    end
    if name in get_all_field_keys()
        if size(get_field(name)) != dof
            @warn "Field $name exists already with different size. Predefined field is returned"
        end
        return get_field(name)
    end
    fields[vartype][name] = Array{type}(zeros(dof))
    field_types[name] = vartype
    field_array_type[name] = Dict("Type" => "Field", "Dof" => dof)
    return get_field(name)
end

@testset "set_get_field" begin
    nn = test_Data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1] = 2
    nn[2] = 3
    nn[3] = 2
    nn[4] = 2
    nn[5] = 5
    test = test_Data_manager.create_constant_node_field("test", Float64, 1)
    @test test == test_Data_manager.get_field("test")
    @test test_Data_manager.create_constant_node_field("test", Float64, 1) == test_Data_manager.get_field("test")
    test = test_Data_manager.create_constant_node_field("test2", Float64, 3)
    @test test == test_Data_manager.get_field("test2")
    test1, test2 = test_Data_manager.create_node_field("test3", Float64, 1)
    @test test1 == test_Data_manager.get_field("test3N")
    @test test2 == test_Data_manager.get_field("test3NP1")
    test1, test2 = test_Data_manager.create_node_field("test4", Float64, 3)
    @test test1 == test_Data_manager.get_field("test4N")
    @test test2 == test_Data_manager.get_field("test4", "NP1")
    test = test_Data_manager.create_constant_bond_field("test5", Float64, 1)
    @test test == test_Data_manager.get_field("test5")
    test = test_Data_manager.create_constant_bond_field("test6", Float64, 3)
    @test test == test_Data_manager.get_field("test6")
    test1, test2 = test_Data_manager.create_bond_field("test7", Float64, 1)
    @test test1 == test_Data_manager.get_field("test7", "N")
    @test test2 == test_Data_manager.get_field("test7", "NP1")
    test1, test2 = test_Data_manager.create_bond_field("test8", Float64, 3)
    @test test1 == test_Data_manager.get_field("test8", "N")
    @test test2 == test_Data_manager.get_field("test8", "NP1")
    testnewFloat = test_Data_manager.create_constant_node_field("testnewFloat", Float16, 1)
    @test typeof(testnewFloat[1]) == Float16
    testnewInt = test_Data_manager.create_constant_node_field("testnewInt", Int8, 1)
    @test typeof(testnewInt[1]) == Int8

    testDoesnotExists = test_Data_manager.get_field("does not exist", "NP1")
    @test isnothing(testDoesnotExists)
    testDoesnotExists = test_Data_manager.get_field("does not exist")
    @test isnothing(testDoesnotExists)
end

@testset "Matrix" begin
    #Arrays
    test = test_Data_manager.create_constant_node_field("test9", Float64, "Matrix", 2)
    test[1, 1, 1] = 1.2
    test[1, 2, 1] = -1.2
    test[1, 1, 2] = 1.4
    test[1, 2, 2] = 1.2
    @test test == test_Data_manager.get_field("test9")
    test = test_Data_manager.create_constant_bond_field("test10", Float64, "Matrix", 3)
    test[1][1, 1, 1] = 1.2
    test[2][1, 2, 1] = -1.2
    test[2][1, 1, 3] = 1.4
    test[2][1, 2, 2] = 1.2
    @test test == test_Data_manager.get_field("test10")
    test1, test2 = test_Data_manager.create_bond_field("test11", Float64, "Matrix", 6)
    @test test1 == test_Data_manager.get_field("test11", "N")
    @test test2 == test_Data_manager.get_field("test11", "NP1")
    test1, test2 = test_Data_manager.create_node_field("test12", Float64, "Matrix", 3)
    @test test1 == test_Data_manager.get_field("test12", "N")
    @test test2 == test_Data_manager.get_field("test12", "NP1")
end


testNP1NDict = test_Data_manager.get_NP1_to_N_Dict()

@testset "get_NP1_to_N_Dict" begin
    @test testNP1NDict["BNP1"] == "BN"
    @test testNP1NDict["DNP1"] == "DN"
    @test testNP1NDict["GNP1"] == "GN"
    @test testNP1NDict["INP1"] == "IN"
end
@testset "set_and_get_values" begin
    DN = test_Data_manager.get_field("DN")
    DN[1, 3] = 10
    DNtest = test_Data_manager.get_field("DN")
    @test DN[1, 3] == DNtest[1, 3]
end

DN = test_Data_manager.get_field("DN")
DNP1 = test_Data_manager.get_field("DNP1")

IN = test_Data_manager.get_field("IN")
INP1 = test_Data_manager.get_field("INP1")
bd = test_Data_manager.create_bond_field("Bond Damage", Float64, 1)
@testset "switch_NP1_to_N" begin
    bmatrixN, bmatrixNP1 = test_Data_manager.create_bond_field("Bmat", Float64, "Matrix", 2)
    nmatrixN, nmatrixNP1 = test_Data_manager.create_node_field("Nmat", Float64, "Matrix", 2)
    IN[2][1, 3] = 0

    DNP1[2, 3] = 5
    @test DN[2, 3] == 0
    @test DNP1[2, 3] == 5
    # bonds
    INP1[2][1, 3] = 5
    @test IN[2][1, 3] == 0
    @test INP1[2][1, 3] == 5
    nmatrixNP1[1, 1, 1] = 2
    nmatrixNP1[1, 1, 2] = 2
    nmatrixNP1[1, 2, 1] = 3
    nmatrixNP1[1, 2, 2] = 4
    bmatrixNP1[1][1, 1, 1] = 2
    bmatrixNP1[1][2, 1, 2] = 2
    bmatrixNP1[1][1, 2, 1] = 3
    bmatrixNP1[1][2, 2, 2] = 4
    # extra test, because Bond Damage is set to one, to avoid unneccessary operations
    bd = test_Data_manager.get_field("Bond Damage", "NP1")
    @test sum(maximum(bd)) == 0
    test_Data_manager.switch_NP1_to_N()

    @test DN[2, 3] == 5
    @test nmatrixN[1, 1, 1] == 2
    @test nmatrixN[1, 1, 2] == 2
    @test nmatrixN[1, 2, 1] == 3
    @test nmatrixN[1, 2, 2] == 4
    @test nmatrixNP1[1, 1, 1] == 0.0
    @test nmatrixNP1[1, 1, 2] == 0.0
    @test nmatrixNP1[1, 2, 1] == 0.0
    @test nmatrixNP1[1, 2, 2] == 0.0
    bmatrixN = test_Data_manager.get_field("Bmat", "N")
    bmatrixNP1 = test_Data_manager.get_field("Bmat", "NP1")
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

    @test INP1[2][1, 3] == 0
    # bonds
    @test IN[2][1, 3] == 5
    bd = test_Data_manager.get_field("Bond Damage", "NP1")
    for id in eachindex(bd)
        @test sum(bd[id]) == nn[id]
    end

end

@testset "ut_nodesets" begin
    @test test_Data_manager.get_nnsets() == 0
    test_Data_manager.set_nset("N1", [1, 2])
    @test test_Data_manager.get_nnsets() == 1
    test_Data_manager.set_nset("N2", [4, 5])
    @test test_Data_manager.get_nnsets() == 2
    test_Data_manager.set_nset("N3", [1, 12, 22])
    @test test_Data_manager.get_nnsets() == 3
    nsets = test_Data_manager.get_nsets()
    @test nsets["N1"] == [1, 2]
    @test nsets["N2"] == [4, 5]
    @test nsets["N3"] == [1, 12, 22]
    test_Data_manager.set_nset("N3", [1, 12])
    @test nsets["N3"] == [1, 12]
end

@testset "ut_blocklist" begin
    blocklist = test_Data_manager.get_block_list()
    @test length(blocklist) == 0
    test_Data_manager.set_block_list([1, 2, 3, 4, 4, 4, 1, 1, 1, 2, 2])
    blocklist = test_Data_manager.get_block_list()
    @test length(blocklist) == 4
    @test blocklist == [1, 2, 3, 4]
    test_Data_manager.set_block_list([4, 4, 2, 2, 1, 1])
    blocklist = test_Data_manager.get_block_list()
    @test length(blocklist) == 3
    @test blocklist == [1, 2, 4]
end

@testset "ut_properties" begin
    test_Data_manager.set_block_list([2, 3, 1, 1])
    test_Data_manager.init_property()
    @test length(test_Data_manager.properties) == 3
    @test isnothing(test_Data_manager.get_property(1, "Material Model", "E"))
    test_Data_manager.set_property(1, "Material Model", "E", 3)
    test_Data_manager.get_property(1, "Material Model", "E")
    @test test_Data_manager.get_property(1, "Material Model", "E") == 3
    test_Data_manager.set_property(1, "Material Model", "C", "Hello Test")
    @test test_Data_manager.get_property(1, "Material Model", "C") == "Hello Test"
    test_Data_manager.set_property(2, "Material Model", "E", 1.1)
    @test test_Data_manager.get_property(2, "Material Model", "E") == 1.1
    test_Data_manager.set_property(2, "Thermal Model", "E", [3 1 2; 1 2 3; 1 3 4])
    @test test_Data_manager.get_property(2, "Thermal Model", "E") == [3 1 2; 1 2 3; 1 3 4]
    test_Data_manager.set_property(3, "Thermal Model", "Q", 23.1)
    @test test_Data_manager.get_property(3, "Thermal Model", "Q") == 23.1
    test_Data_manager.set_property(3, "Damage Model", "SS", 0.1)
    @test test_Data_manager.get_property(3, "Damage Model", "SS") == 0.1
    test_Data_manager.set_property(1, "Additive Model", "E", [1, 2, 3])
    @test test_Data_manager.get_property(1, "Additive Model", "E") == [1, 2, 3]
    test_Data_manager.set_property(2, "Additive Model", "Qd", true)
    @test test_Data_manager.get_property(2, "Additive Model", "Qd") == true
    @test isnothing(test_Data_manager.get_property(2, "Additive Model", "not there"))
    @test test_Data_manager.get_properties(1, "Material Model") == Dict("C" => "Hello Test", "E" => 3)
    @test test_Data_manager.get_properties(1, "Thermal Model") == Dict()
    @test test_Data_manager.get_properties(2, "Material Model") == Dict("E" => 1.1)
    @test test_Data_manager.get_properties(2, "Thermal Model") == Dict("E" => [3 1 2; 1 2 3; 1 3 4])
    @test test_Data_manager.get_properties(1, "") == Dict()
    @test !test_Data_manager.check_property(1, "This is not a property")
    @test isnothing(test_Data_manager.get_property(1, "Thermal Model", "This is not a property"))
    test_Data_manager.set_properties("FEM", Dict("A" => 2, "C" => "Model"))
    @test test_Data_manager.get_properties(1, "FEM") == Dict("A" => 2, "C" => "Model")
    @test test_Data_manager.get_properties(2, "FEM") == Dict("A" => 2, "C" => "Model")
    @test test_Data_manager.get_properties(3, "FEM") == Dict("A" => 2, "C" => "Model")
end

@testset "get_physics_options" begin
    test_Data_manager = Data_manager
    physics_options = test_Data_manager.get_physics_options()
    @test physics_options["Deformed Bond Geometry"]
    @test !physics_options["Bond Associated Shape Tensor"]
    @test !physics_options["Bond Associated Deformation Gradient"]
    @test !physics_options["Deformation Gradient"]
    @test !physics_options["Shape Tensor"]
    test_Data_manager.physics_options["Deformed Bond Geometry"] = false
    physics_options = test_Data_manager.get_physics_options()
    @test !physics_options["Deformed Bond Geometry"]
    @test !physics_options["Bond Associated Shape Tensor"]
    @test !physics_options["Bond Associated Deformation Gradient"]
    @test !physics_options["Deformation Gradient"]
    @test !physics_options["Shape Tensor"]
    test_Data_manager.physics_options["Bond Associated Deformation Gradient"] = true
    physics_options = test_Data_manager.get_physics_options()
    @test physics_options["Deformed Bond Geometry"]
    @test physics_options["Bond Associated Shape Tensor"]
    @test physics_options["Bond Associated Deformation Gradient"]
    @test physics_options["Deformation Gradient"]
    @test physics_options["Shape Tensor"]
    test_Data_manager.physics_options["Deformed Bond Geometry"] = false
    test_Data_manager.physics_options["Shape Tensor"] = false
    test_Data_manager.physics_options["Bond Associated Shape Tensor"] = false
    test_Data_manager.physics_options["Deformation Gradient"] = false
    test_Data_manager.physics_options["Bond Associated Deformation Gradient"] = false
    test_Data_manager.physics_options["Deformation Gradient"] = true
    physics_options = test_Data_manager.get_physics_options()
    @test physics_options["Deformed Bond Geometry"]
    @test !physics_options["Bond Associated Shape Tensor"]
    @test !physics_options["Bond Associated Deformation Gradient"]
    @test physics_options["Deformation Gradient"]
    @test physics_options["Shape Tensor"]
    test_Data_manager.physics_options["Deformed Bond Geometry"] = false
    test_Data_manager.physics_options["Shape Tensor"] = false
    test_Data_manager.physics_options["Bond Associated Shape Tensor"] = false
    test_Data_manager.physics_options["Deformation Gradient"] = true
    test_Data_manager.physics_options["Bond Associated Deformation Gradient"] = true
    physics_options = test_Data_manager.get_physics_options()
    @test physics_options["Deformed Bond Geometry"]
    @test physics_options["Bond Associated Shape Tensor"]
    @test physics_options["Bond Associated Deformation Gradient"]
    @test physics_options["Deformation Gradient"]
    @test physics_options["Shape Tensor"]
end

@testset "ut_get_and_set_inverse_nlist" begin
    inv_nlist = test_Data_manager.get_inverse_nlist()
    @test typeof(inv_nlist) == Vector{Dict{Int64,Int64}}
    @test length(inv_nlist) == 0
    test_Data_manager.set_inverse_nlist([Dict{Int64,Int64}(1 => 2), Dict{Int64,Int64}(1 => 2, 2 => 1)])
    inv_nlist = test_Data_manager.get_inverse_nlist()
    @test typeof(inv_nlist) == Vector{Dict{Int64,Int64}}
    @test length(inv_nlist) == 2
    @test inv_nlist[1] == Dict{Int64,Int64}(1 => 2)
    @test inv_nlist[2] == Dict{Int64,Int64}(1 => 2, 2 => 1)
    test_Data_manager.set_inverse_nlist([Dict{Int64,Int64}(1 => 2, 2 => 1, 9 => 2)])
    inv_nlist = test_Data_manager.get_inverse_nlist()
    @test length(inv_nlist) == 1
    @test inv_nlist[1] == Dict{Int64,Int64}(1 => 2, 2 => 1, 9 => 2)
end

@testset "ut_rotation" begin
    rotation, angles = test_Data_manager.rotation_data()
    @test !rotation
    @test isnothing(angles)
    test_angles = test_Data_manager.create_constant_node_field("Angles", Float32, 3)
    rotation, angles = test_Data_manager.rotation_data()
    @test rotation
    @test angles == test_angles
    rotation, angles = test_Data_manager.rotation_data("Node")
    @test rotation
    @test angles == test_angles
    rotation, angles = test_Data_manager.rotation_data("Element")
    @test !rotation
    @test isnothing(angles)
    test_angles = test_Data_manager.create_constant_node_field("Element Angles", Float32, 3)# in code it has length number of elements * element integration points
    rotation, angles = test_Data_manager.rotation_data("Element")
    @test rotation
    @test angles == test_angles
    @test isnothing(test_Data_manager.rotation_data("Hello"))

end

