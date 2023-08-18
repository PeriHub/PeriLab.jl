include("../exodus_export.jl")
include("../../Support/data_manager.jl")
using Test
import .Data_manager
import .Write_Exodus_Results
using Exodus


@testset "ut_get_block_nodes" begin
    block_Id = [1, 1, 2, 2, 2, 3, 3, 3, 1, 3, 3]
    test = Write_Exodus_Results.get_block_nodes(block_Id, 1)
    @test test == [1 2 9]
    test = Write_Exodus_Results.get_block_nodes(block_Id, 2)
    @test test == [3 4 5]
    test = Write_Exodus_Results.get_block_nodes(block_Id, 3)
    @test test == [6 7 8 10 11]
end

if !isdir("tmp")
    mkdir("tmp")
end
filename = "./tmp/" * "test.e"
@testset "ut_create_result_file" begin
    nnodes = 4
    dof = 3
    exo = Write_Exodus_Results.create_result_file(filename, nnodes, dof, 1, 0)
    @test isfile(filename)
    @test exo.file_name == filename
    @test exo.init.num_dim == dof
    @test exo.init.num_nodes == nnodes
    @test exo.init.num_node_sets == 0
    @test exo.init.num_elem_blks == 1

    close(exo)
    rm(filename)
    nnodes = 300
    dof = 2
    exo = Write_Exodus_Results.create_result_file(filename, nnodes, dof, 3, 2)
    @test isfile(filename)
    @test exo.file_name == filename
    @test exo.init.num_dim == dof
    @test exo.init.num_nodes == nnodes
    @test exo.init.num_node_sets == 2
    @test exo.init.num_elem_blks == 3
    close(exo)
    rm(filename)
end

@testset "ut_init_results_in_exodus" begin
    nnodes = 5
    dof = 2
    testDatamanager = Data_manager
    testDatamanager.set_nnodes(nnodes)
    testDatamanager.set_dof(dof)
    coordinates = testDatamanager.create_constant_node_field("Coordinates", Float32, 2)
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
    testDatamanager.create_constant_node_field("Block_Id", Int64, 1)
    block_Id = testDatamanager.get_field("Block_Id")
    block_Id .+= 1
    block_Id[end] = 2
    outputs = ["Displacement", "ForcesNP1"]
    nsets = testDatamanager.get_nsets()
    coords = vcat(transpose(coordinates))

    exo = Write_Exodus_Results.create_result_file(filename, nnodes, dof, maximum(block_Id), 0)
    @test exo.init.num_dim == dof

    exo = Write_Exodus_Results.init_results_in_exodus(exo, outputs, coords, block_Id, nsets)
    @test length(exo.nodal_var_name_dict) == 2

    entries = collect(keys(exo.nodal_var_name_dict))

    @test entries[1] in outputs
    @test entries[2] in outputs

    exo_coords = read_coordinates(exo)
    exo_nsets = read_sets(exo, NodeSet)
    @test coords == exo_coords
    @test exo_nsets == []
    close(exo)

    rm(filename)
end
@testset "ut_paraview_specifics" begin
    @test Write_Exodus_Results.paraview_specifics(1) == "x"
    @test Write_Exodus_Results.paraview_specifics(2) == "y"
    @test Write_Exodus_Results.paraview_specifics(3) == "z"
end

@testset "ut_get_paraviewCoordinates" begin
    for i in 1:3
        @test Write_Exodus_Results.get_paraviewCoordinates(1, i) == "x"
        @test Write_Exodus_Results.get_paraviewCoordinates(2, i) == "y"
        @test Write_Exodus_Results.get_paraviewCoordinates(3, i) == "z"
    end
    for ref = 4:9
        for i in 1:3
            for j in 1:3
                @test Write_Exodus_Results.get_paraviewCoordinates((i - 1) * 3 + j, ref) == Write_Exodus_Results.paraview_specifics(i) * Write_Exodus_Results.paraview_specifics(j)
            end
        end
    end
end