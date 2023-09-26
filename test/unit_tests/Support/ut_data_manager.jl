# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../../src/Support/data_manager.jl")
using MPI
using Test
@testset "set_comm" begin
    # MPI.Init()
    comm = MPI.COMM_WORLD
    testDatamanager = Data_manager
    testDatamanager.set_comm(comm)
    b = testDatamanager.get_comm()
    @test comm == b
    # MPI.Finalize()
end
@testset "get_local_nodes" begin
    testDatamanager = Data_manager
    testDatamanager.set_glob_to_loc(Dict{Int64,Int64}(1 => 1, 3 => 2, 2 => 3))
    @test testDatamanager.get_local_nodes([1, 2, 3]) == [1, 3, 2]
    @test testDatamanager.get_local_nodes([1]) == [1]
    @test testDatamanager.get_local_nodes([2]) == [3]
    @test testDatamanager.get_local_nodes([3]) == [2]
    @test testDatamanager.get_local_nodes([4]) == []
    @test testDatamanager.get_local_nodes([4, 2, 3]) == [3, 2]
    testDatamanager.set_glob_to_loc(Dict{Int64,Int64}(3 => 1, 2 => 2, 4 => 3))
    @test testDatamanager.get_local_nodes([1]) == []
    @test testDatamanager.get_local_nodes([4]) == [3]
end


@testset "get_set_functions" begin
    testDatamanager = Data_manager
    for i in 1:20
        testDatamanager.set_dof(i)
        @test testDatamanager.get_dof() == i
        testDatamanager.set_nmasters(i)
        @test testDatamanager.get_nnodes() == i
    end
    testDatamanager.set_nmasters(97)
    testDatamanager.set_nslaves(5)
    @test testDatamanager.get_nnodes() == 97
    @test testDatamanager.nslaves == 5
    @test testDatamanager.nmasters == 97
    @test testDatamanager.nnodes == 102
    nnodes = testDatamanager.get_nnodes()
    nnodes = 3
    @test testDatamanager.get_nnodes() == 97
    @test nnodes == 3
end
testDatamanager = Data_manager
nmaster = 3
nslave = 2
testDatamanager.set_nmasters(nmaster)
testDatamanager.set_nslaves(nslave)
nn = testDatamanager.create_constant_node_field("Number of Neighbors", Int64, 1)
nn[1] = 2
nn[2] = 3
nn[3] = 2
nn[4] = 2
nn[5] = 5

testDatamanager.create_constant_node_field("A", Float32, 1)
B = testDatamanager.create_node_field("B", Bool, 1)
testDatamanager.create_constant_node_field("C", Float32, 4)
testDatamanager.create_node_field("D", Int64, 7)
testDatamanager.create_constant_bond_field("F", Float32, 1)
testDatamanager.create_bond_field("G", Bool, 1)
testDatamanager.create_constant_bond_field("H", Float32, 4)
testDatamanager.create_bond_field("I", Int64, 7)
testfield_keys = testDatamanager.get_all_field_keys()
@testset "create data fields -> get all fields" begin
    @test testDatamanager.get_nnodes() == nmaster
    @test B[1] == testDatamanager.get_field("BN")
    @test B[2] == testDatamanager.get_field("B", "NP1")
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

@testset "create_existing_field" begin
    field1, field2 = testDatamanager.create_node_field("D", Int64, 3)
    testfield_keys = testDatamanager.get_all_field_keys()
    @test "DN" in testfield_keys
    @test "DNP1" in testfield_keys
    # because field is not overwritten, the dof value stay
    @test size(field1) == (nmaster + nslave, 7)
    @test size(field2) == (nmaster + nslave, 7)
end

@testset "get_field" begin

    A = testDatamanager.get_field("A")
    @test typeof(A[1]) == Float32
    @test length(A) == testDatamanager.nnodes == nmaster + nslave

    B = testDatamanager.get_field("BN")
    @test typeof(B[1]) == Bool
    @test length(B) == testDatamanager.nnodes == nmaster + nslave

    C = testDatamanager.get_field("C")
    @test typeof(C[1, 1]) == Float32
    @test length(C[:, 1]) == testDatamanager.nnodes == nmaster + nslave
    @test length(C[1, :]) == 4

    D = testDatamanager.get_field("DNP1")
    @test typeof(D[1, 1]) == Int64
    @test length(D[:, 1]) == testDatamanager.nnodes == nmaster + nslave
    @test length(D[1, :]) == 7

    F = testDatamanager.get_field("F")
    @test typeof(F[1, 1][1]) == Float32
    @test length(F[:, 1]) == nmaster # bonds are only of length nmaster
    @test length(F[1]) == nn[1]
    @test length(F[2]) == nn[2]
    @test length(F[3]) == nn[3]
    G = testDatamanager.get_field("GN")
    @test typeof(G[1, 1][1]) == Bool
    @test length(G[:, 1]) == nmaster

    H = testDatamanager.get_field("H")
    @test typeof(H[1][1, 1][1]) == Float32
    @test length(H[1][:, 1]) == nn[1]
    @test length(H[1][1, :]) == 4
    @test length(H[:][:, :]) == nmaster

    I = testDatamanager.get_field("INP1")
    @test typeof(I[1][1, 1]) == Int64
    for i in 1:nmaster
        @test length(I[i][:, 1]) == nn[i]
    end
    @test length(I[1][1, :]) == 7
    @test length(I[:][:, :]) == nmaster

end

@testset "set_get_field" begin
    nn = testDatamanager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1] = 2
    nn[2] = 3
    nn[3] = 2
    nn[4] = 2
    nn[5] = 5
    test = testDatamanager.create_constant_node_field("test", Float32, 1)
    @test test == testDatamanager.get_field("test")
    test = testDatamanager.create_constant_node_field("test2", Float32, 3)
    @test test == testDatamanager.get_field("test2")
    test1, test2 = testDatamanager.create_node_field("test3", Float32, 1)
    @test test1 == testDatamanager.get_field("test3N")
    @test test2 == testDatamanager.get_field("test3NP1")
    test1, test2 = testDatamanager.create_node_field("test4", Float32, 3)
    @test test1 == testDatamanager.get_field("test4N")
    @test test2 == testDatamanager.get_field("test4", "NP1")
    test = testDatamanager.create_constant_bond_field("test5", Float32, 1)
    @test test == testDatamanager.get_field("test5")
    test = testDatamanager.create_constant_bond_field("test6", Float32, 3)
    @test test == testDatamanager.get_field("test6")
    test1, test2 = testDatamanager.create_bond_field("test7", Float32, 1)
    @test test1 == testDatamanager.get_field("test7", "N")
    @test test2 == testDatamanager.get_field("test7", "NP1")
    test1, test2 = testDatamanager.create_bond_field("test8", Float32, 3)
    @test test1 == testDatamanager.get_field("test8", "N")
    @test test2 == testDatamanager.get_field("test8", "NP1")
end

testNP1NDict = testDatamanager.get_NP1_to_N_Dict()

@testset "get_NP1_to_N_Dict" begin
    @test testNP1NDict["BNP1"] == "BN"
    @test testNP1NDict["DNP1"] == "DN"
    @test testNP1NDict["GNP1"] == "GN"
    @test testNP1NDict["INP1"] == "IN"
end
@testset "set_and_get_values" begin
    DN = testDatamanager.get_field("DN")
    DN[1, 3] = 10
    DNtest = testDatamanager.get_field("DN")
    @test DN[1, 3] == DNtest[1, 3]
end

DN = testDatamanager.get_field("DN")
DNP1 = testDatamanager.get_field("DNP1")

IN = testDatamanager.get_field("IN")
INP1 = testDatamanager.get_field("INP1")
bd = testDatamanager.create_bond_field("Bond Damage", Float32, 1)
IN[2][1, 3] = 0
@testset "switch_NP1_to_N" begin
    DNP1[2, 3] = 5
    @test DN[2, 3] == 0
    @test DNP1[2, 3] == 5
    # bonds
    INP1[2][1, 3] = 5
    @test IN[2][1, 3] == 0
    @test INP1[2][1, 3] == 5
    # extra test, because Bond Damage is set to one, to avoid unneccessary operations
    bd = testDatamanager.get_field("Bond Damage", "NP1")
    @test sum(maximum(bd)) == 0
    testDatamanager.switch_NP1_to_N()
    @test DN[2, 3] == 5
    # dependency test

    @test DNP1[2, 3] == 0
    DNP1[2, 3] = 6
    @test DNP1[2, 3] == 6
    @test DN[2, 3] == 5

    @test INP1[2][1, 3] == 0
    # bonds
    @test IN[2][1, 3] == 5
    bd = testDatamanager.get_field("Bond Damage", "NP1")
    for id in eachindex(bd)
        @test sum(bd[id]) == nn[id]
    end

end

@testset "ut_nodesets" begin
    @test testDatamanager.get_nnsets() == 0
    testDatamanager.set_nset("N1", [1, 2])
    @test testDatamanager.get_nnsets() == 1
    testDatamanager.set_nset("N2", [4, 5])
    @test testDatamanager.get_nnsets() == 2
    testDatamanager.set_nset("N3", [1, 12, 22])
    @test testDatamanager.get_nnsets() == 3
    nsets = testDatamanager.get_nsets()
    @test nsets["N1"] == [1, 2]
    @test nsets["N2"] == [4, 5]
    @test nsets["N3"] == [1, 12, 22]
end

@testset "ut_blocklist" begin
    blocklist = testDatamanager.get_block_list()
    @test length(blocklist) == 0
    testDatamanager.set_block_list([1, 2, 3, 4, 4, 4, 1, 1, 1, 2, 2])
    blocklist = testDatamanager.get_block_list()
    @test length(blocklist) == 4
    @test blocklist == [1, 2, 3, 4]
    testDatamanager.set_block_list([4, 4, 2, 2, 1, 1])
    blocklist = testDatamanager.get_block_list()
    @test length(blocklist) == 3
    @test blocklist == [1, 2, 4]
end

@testset "ut_properties" begin
    testDatamanager.set_block_list([2, 3, 1, 1])
    testDatamanager.init_property()
    @test length(testDatamanager.properties) == 3
    testDatamanager.set_property(1, "Material Model", "E", 3)
    @test testDatamanager.get_property(1, "Material Model", "E") == 3
    testDatamanager.set_property(1, "Material Model", "C", "Hello Test")
    @test testDatamanager.get_property(1, "Material Model", "C") == "Hello Test"
    testDatamanager.set_property(2, "Material Model", "E", 1.1)
    @test testDatamanager.get_property(2, "Material Model", "E") == 1.1
    testDatamanager.set_property(2, "Thermal Model", "E", [3 1 2; 1 2 3; 1 3 4])
    @test testDatamanager.get_property(2, "Thermal Model", "E") == [3 1 2; 1 2 3; 1 3 4]
    testDatamanager.set_property(3, "Thermal Model", "Q", 23.1)
    @test testDatamanager.get_property(3, "Thermal Model", "Q") == 23.1
    testDatamanager.set_property(3, "Damage Model", "SS", 0.1)
    @test testDatamanager.get_property(3, "Damage Model", "SS") == 0.1
    testDatamanager.set_property(1, "Additive Model", "E", [1, 2, 3])
    @test testDatamanager.get_property(1, "Additive Model", "E") == [1, 2, 3]
    testDatamanager.set_property(2, "Additive Model", "Qd", true)
    @test testDatamanager.get_property(2, "Additive Model", "Qd") == true
    @test testDatamanager.get_property(2, "Additive Model", "not there") == Nothing
    @test testDatamanager.get_properties(1, "Material Model") == Dict("C" => "Hello Test", "E" => 3)
    @test testDatamanager.get_properties(1, "Thermal Model") == Dict()
    @test testDatamanager.get_properties(2, "Material Model") == Dict("E" => 1.1)
    @test testDatamanager.get_properties(2, "Thermal Model") == Dict("E" => [3 1 2; 1 2 3; 1 3 4])
end