include("../data_manager.jl")
using Test
import .Data_manager
@testset "get_set_functions" begin
    testDatamanager = Data_manager
    for i in 1:20
        testDatamanager.set_dof(i)
        @test testDatamanager.get_dof() == i
        testDatamanager.set_nnodes(i)
        @test testDatamanager.get_nnodes() == i
        testDatamanager.set_nbonds(i)
        @test testDatamanager.get_nbonds() == i
    end
    testDatamanager.set_nnodes(97)
    @test testDatamanager.get_nnodes() == 97
    nnodes = testDatamanager.get_nnodes()
    nnodes = 3
    @test testDatamanager.get_nnodes() == 97
    @test nnodes == 3
end
testDatamanager = Data_manager

testDatamanager.set_nnodes(5)
nn = testDatamanager.create_constant_node_field("Number of Neighbors", Int64, 1)
nn[1] = 2
nn[2] = 3
nn[3] = 2
nn[4] = 2
nn[5] = 5

testDatamanager.create_constant_node_field("A", Float32, 1)
testDatamanager.create_node_field("B", Bool, 1)
testDatamanager.create_constant_node_field("C", Float32, 4)
testDatamanager.create_node_field("D", Int64, 7)
testDatamanager.create_constant_bond_field("F", Float32, 1)
testDatamanager.create_bond_field("G", Bool, 1)
testDatamanager.create_constant_bond_field("H", Float32, 4)
testDatamanager.create_bond_field("I", Int64, 7)
testfield_keys = testDatamanager.get_all_field_keys()
@testset "create data fields -> get all fields" begin


    @test testDatamanager.get_nnodes() == 5

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

@testset "get_field" begin
    A = testDatamanager.get_field("A")
    @test typeof(A[1]) == Float32
    @test length(A) == testDatamanager.get_nnodes()

    B = testDatamanager.get_field("BN")
    @test typeof(B[1]) == Bool
    @test length(B) == testDatamanager.get_nnodes()

    C = testDatamanager.get_field("C")
    @test typeof(C[1, 1]) == Float32
    @test length(C[:, 1]) == testDatamanager.get_nnodes()
    @test length(C[1, :]) == 4

    D = testDatamanager.get_field("DNP1")
    @test typeof(D[1, 1]) == Int64
    @test length(D[:, 1]) == testDatamanager.get_nnodes()
    @test length(D[1, :]) == 7

    F = testDatamanager.get_field("F")
    @test typeof(F[1, 1][1]) == Float32
    @test length(F[:, 1]) == testDatamanager.get_nnodes()
    @test length(F[1]) == nn[1]
    @test length(F[2]) == nn[2]
    @test length(F[3]) == nn[3]
    @test length(F[4]) == nn[4]
    @test length(F[5]) == nn[5]
    G = testDatamanager.get_field("GN")
    @test typeof(G[1, 1][1]) == Bool
    @test length(G[:, 1]) == testDatamanager.get_nnodes()

    H = testDatamanager.get_field("H")
    @test typeof(H[1][1, 1][1]) == Float32
    @test length(H[1][:, 1]) == nn[1]
    @test length(H[1][1, :]) == 4
    @test length(H[:][:, :]) == testDatamanager.get_nnodes()

    I = testDatamanager.get_field("INP1")
    @test typeof(I[1][1, 1]) == Int64
    for i in 1:5
        @test length(I[i][:, 1]) == nn[i]
    end
    @test length(I[1][1, :]) == 7
    @test length(I[:][:, :]) == testDatamanager.get_nnodes()

end
