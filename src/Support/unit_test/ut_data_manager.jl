include("../data_manager.jl")
using Test
import Data_manager
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

end
@testset "create data fields" begin
    testDatamanager = Data_manager

    testDatamanager.set_nnodes(5)
    testDatamanager.create_constant_node_field(Float32, "A", 1)
    testDatamanager.create_node_field(Float32, "B", 1)
    testDatamanager.create_constant_node_field(Float32, "C", 4)
    testDatamanager.create_node_field(Int64, "D", 7)

end