include("../../../src/IO/exodus_export.jl")
include("../../../src/Support/data_manager.jl")
include("../../../src/Support/Parameters/parameter_handling.jl")
using Test
import .Write_Exodus_Results
using Exodus
using Pkg
@testset "ut_get_block_nodes" begin
    block_Id = [1, 1, 2, 2, 2, 3, 3, 3, 1, 3, 3, 4]
    test = Write_Exodus_Results.get_block_nodes(block_Id, 1)
    @test test == [1 2 9]
    test = Write_Exodus_Results.get_block_nodes(block_Id, 2)
    @test test == [3 4 5]
    test = Write_Exodus_Results.get_block_nodes(block_Id, 3)
    @test test == [6 7 8 10 11]
    test = Write_Exodus_Results.get_block_nodes(block_Id, 4)
    @test test[1] == 12
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

if !isdir("tmp")
    mkdir("tmp")
end
@testset "ut_create_result_file" begin
    filename = "./tmp/" * "test.e"
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
filename = "./tmp/" * "test_2.e"
nnodes = 5
dof = 2
testDatamanager = Data_manager
testDatamanager.set_nmasters(nnodes)
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
#outputs = ["Displacements", "Forces"]
testDatamanager.set_nset("Nset_1", [1, 2])
testDatamanager.set_nset("Nset_2", [5])

nsets = testDatamanager.get_nsets()
coords = vcat(transpose(coordinates))
outputs = Dict("Forcesxx" => ["ForcesNP1", 1, 1, Float32], "Forcesxy" => ["ForcesNP1", 1, 2, Float32], "Forcesxz" => ["ForcesNP1", 1, 3, Float32], "Forcesyx" => ["ForcesNP1", 1, 4, Float32], "Forcesyy" => ["ForcesNP1", 1, 5, Float32], "Forcesyz" => ["ForcesNP1", 1, 6, Float32], "Displacements" => ["DisplacementsNP1", 2, 1, Float32])
exo = Write_Exodus_Results.create_result_file(filename, nnodes, dof, maximum(block_Id), length(nsets))
exo = Write_Exodus_Results.init_results_in_exodus(exo, outputs, coords, block_Id, Vector{Int32}(1:maximum(block_Id)), nsets)
exos = []
push!(exos, exo)
exos[1] = Write_Exodus_Results.write_step_and_time(exos[1], 2, 2.2)
exos[1] = Write_Exodus_Results.write_step_and_time(exos[1], 3, 3.7)
exos[1] = Write_Exodus_Results.write_step_and_time(exos[1], 4, 4.7)
exos[1] = Write_Exodus_Results.write_step_and_time(exos[1], 5, 5.7)
exos[1] = Write_Exodus_Results.write_step_and_time(exos[1], 6, 6.7)

@testset "ut_init_results_in_exodus" begin
    @test exo.init.num_dim == dof
    @test length(exo.nodal_var_name_dict) == 7
    entries = collect(keys(exo.nodal_var_name_dict))
    ref = collect(keys(outputs))
    @test sort(entries) == sort(ref)
    exo_coords = read_coordinates(exo)
    exo_nsets = read_sets(exo, NodeSet)
    @test length(exo_nsets) == length(nsets)
    @test coords == exo_coords
    @warn "Info test deactivated"
    # @test ["PeriLab Version " * string(Pkg.project().version) * ", under BSD License", "Copyright (c) 2023, Christian Willberg, Jan-Timo Hesse", "compiled with Julia Version " * string(VERSION)] == read_info(exo)
    @test read_number_of_time_steps(exo) == 6
    @test read_time(exo, 2) == 2.2
    @test read_time(exo, 3) == 3.7
    @test read_time(exo, 4) == 4.7
    @test read_time(exo, 5) == 5.7
    @test read_time(exo, 6) == 6.7
    @test read_name(exo, Block, 1) == "Block_1"
    @test read_name(exo, Block, 2) == "Block_2"
end


testDatamanager.create_node_field("Forces", Float32, 6)
testDatamanager.create_node_field("Displacements", Float32, 1)
force = testDatamanager.get_field("Forces", "NP1")
disp = testDatamanager.get_field("Displacements", "NP1")
force[5, 1:6] .= 3.3
force[1:3, 6] .= 2.3
disp[1] = 3
disp[2] = 3.00001
disp[3] = 2.1
disp[4] = -1.8
disp[5] = 0

exo = Write_Exodus_Results.write_nodal_results_in_exodus(exo, 2, outputs, testDatamanager)

test_disp_step_zero = read_values(exo, NodalVariable, 1, 1, "Displacements")

test_disp_step_one = read_values(exo, NodalVariable, 2, 1, 1)
@testset "ut_write_results_in_exodus" begin
    @test test_disp_step_zero == zeros(5)
    for id in eachindex(test_disp_step_one)
        if disp[id] != 0
            @test test_disp_step_one[id] / disp[id] - 1 < 1e-8
        else
            @test test_disp_step_one[id] == disp[id]
        end
    end

    ftest = read_values(exo, NodalVariable, 2, 1, "Forcesxx")
    for id in 1:4
        @test ftest[id] == 0
    end
    @test ftest[5] / 3.3 - 1 < 1e-8
    ftest = read_values(exo, NodalVariable, 2, 1, "Forcesxy")
    for id in 1:4
        @test ftest[id] == 0
    end
    @test ftest[5] / 3.3 - 1 < 1e-8
    ftest = read_values(exo, NodalVariable, 2, 1, "Forcesxz")
    for id in 1:4
        @test ftest[id] == 0
    end
    @test ftest[5] / 3.3 - 1 < 1e-8
    ftest = read_values(exo, NodalVariable, 2, 1, "Forcesyx")
    for id in 1:4
        @test ftest[id] == 0
    end
    @test ftest[5] / 3.3 - 1 < 1e-8

    ftest = read_values(exo, NodalVariable, 2, 1, "Forcesyy")
    for id in 1:4
        @test ftest[id] == 0
    end
    @test ftest[5] / 3.3 - 1 < 1e-8

    ftest = read_values(exo, NodalVariable, 2, 1, "Forcesyz")
    @test ftest[1] / 2.3 - 1 < 1e-8
    @test ftest[2] / 2.3 - 1 < 1e-8
    @test ftest[3] / 2.3 - 1 < 1e-8
    @test ftest[4] == 0
    @test ftest[5] / 3.3 - 1 < 1e-8

end

@testset "ut_merge_exodus_file" begin
    Write_Exodus_Results.merge_exodus_file(exo.file_name)
end

close(exo)
rm(filename)
