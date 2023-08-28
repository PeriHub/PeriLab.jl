
include("../../Support/data_manager.jl")
include("../IO.jl")
import .Data_manager
using Test
import .IO
using Exodus

testDatamanager = Data_manager
filename1 = "test1.e"
filename2 = "test2.e"
dof = 2
nnodes = 5
testDatamanager.set_nnodes(nnodes)
testDatamanager.set_dof(dof)
testDatamanager.create_constant_node_field("Coordinates", Float32, 2)
coordinates = testDatamanager.get_field("Coordinates")
testDatamanager.create_constant_node_field("Block_Id", Int64, 1)
block_Id = testDatamanager.get_field("Block_Id")
testDatamanager.create_node_field("Displacements", Float32, 2)
testDatamanager.create_node_field("Forces", Float32, 6)
params = Dict("Outputs" => Dict("Output1" => Dict("Output Filename" => filename1, "Output Variables" => Dict("Forces" => true)), "Output2" => Dict("Output Filename" => filename2, "Output Variables" => Dict("Displacements" => true, "Forces" => true))))
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
    output = IO.get_results_mapping(params, testDatamanager)
    @test sort(collect(keys(output[1]))) == ["Forcesxx", "Forcesxy", "Forcesxz", "Forcesyx", "Forcesyy", "Forcesyz"]
    @test sort(collect(keys(output[2]))) == ["Displacementsx", "Displacementsy", "Forcesxx", "Forcesxy", "Forcesxz", "Forcesyx", "Forcesyy", "Forcesyz"]
    for i in 1:2
        dofForce = 0
        dofDisp = 0
        for entry in keys(sort(output[i]))
            if occursin("Forces", entry)
                dofForce += 1
                @test output[i][entry][1] == "ForcesNP1"
                @test output[i][entry][2] == i
                @test output[i][entry][3] == dofForce
                @test output[i][entry][4] == Float32
            else
                dofDisp += 1
                @test output[i][entry][1] == "DisplacementsNP1"
                @test output[i][entry][2] == 1
                @test output[i][entry][3] == dofDisp
                @test output[i][entry][4] == Float32
            end
        end
    end
end

@testset "ut_init_write_results" begin
    exos, outputs = IO.init_write_results(params, testDatamanager)

    @test length(exos) == 2
    @test length(exos[1].nodal_var_name_dict) == 6
    entries = collect(keys(exos[1].nodal_var_name_dict))
    @test sort(entries) == ["Forcesxx", "Forcesxy", "Forcesxz", "Forcesyx", "Forcesyy", "Forcesyz"]
    @test length(exos[2].nodal_var_name_dict) == 8
    entries = collect(keys(exos[2].nodal_var_name_dict))
    @test sort(entries) == ["Displacementsx", "Displacementsy", "Forcesxx", "Forcesxy", "Forcesxz", "Forcesyx", "Forcesyy", "Forcesyz"]

    coords = vcat(transpose(coordinates))
    for exo in exos
        @test exo.init.num_dim == dof
        exo_coords = read_coordinates(exo)
        exo_nsets = read_sets(exo, NodeSet)
        @test coords == exo_coords
        @test exo_nsets == []
    end

    for i in 1:2
        dofForce = 0
        dofDisp = 0
        for entry in keys(sort(outputs[i]))
            if occursin("Forces", entry)
                dofForce += 1
                @test outputs[i][entry][1] == "ForcesNP1"
                @test outputs[i][entry][2] == i
                @test outputs[i][entry][3] == dofForce
                @test outputs[i][entry][4] == Float32
            else
                dofDisp += 1
                @test outputs[i][entry][1] == "DisplacementsNP1"
                @test outputs[i][entry][2] == 1
                @test outputs[i][entry][3] == dofDisp
                @test outputs[i][entry][4] == Float32
            end
        end
    end
    IO.close_files(exos)
    rm(filename1)
    rm(filename2)
end

