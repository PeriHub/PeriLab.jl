
include("../../Support/data_manager.jl")
include("../IO.jl")
import .Data_manager
using Test
import .IO
using Exodus

testDatamanager = Data_manager
filename1 = "test1.e"
filename2 = "test2.e"
testDatamanager.set_nnodes(5)
testDatamanager.set_dof(2)
testDatamanager.create_constant_node_field("Coordinates", Float32, 2)
coordinates = testDatamanager.get_field("Coordinates")
testDatamanager.create_constant_node_field("Block_Id", Int64, 1)
block_Id = testDatamanager.get_field("Block_Id")
testDatamanager.create_node_field("Displacements", Float32, 2)
testDatamanager.create_node_field("Forces", Float32, 6)
params = Dict("Output" => Dict("Output1" => Dict("Output Filename" => filename1, "Output Variables" => Dict("Forces" => true)), "Output2" => Dict("Output Filename" => filename2, "Output Variables" => Dict("Displacements" => true, "Forces" => true))))
coordinates[1, 1] = 0
coordinates[1, 2] = 0
coordinates[2, 1] = 1
coordinates[2, 2] = 0
coordinates[3, 1] = 0
coordinates[3, 2] = 1
coordinates[4, 1] = 1
coordinates[4, 2] = 1
coordinates[5, 1] = 2
coordinates[5, 2] = 2
block_Id .+= 1
block_Id[end] = 2


@testset "ut_get_results_mapping" begin
    output, mapping = IO.get_results_mapping(params, testDatamanager)
    @test output[1] == ["Forcesxx", "Forcesxy", "Forcesxz", "Forcesyx", "Forcesyy", "Forcesyz"]
    @test output[2] == ["Displacementsx", "Displacementsy", "Forcesxx", "Forcesxy", "Forcesxz", "Forcesyx", "Forcesyy", "Forcesyz"]

    for i in 1:2
        for entry in keys(mapping[i])
            if occursin("Forces", entry)
                @test mapping[i][entry] == "ForcesNP1"
            else
                @test mapping[i][entry] == "DisplacementsNP1"
            end
        end
    end

end





