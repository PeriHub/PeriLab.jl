include("../exodus_export.jl")
include("../../Support/data_manager.jl")
using Test
import .Data_manager
import .Write_Exodus_Results
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

#@testset "ut_init_results_in_exodus" begin
nnodes = 4
dof = 2
testDatamanager = Data_manager
testDatamanager.set_nnodes(nnodes)
testDatamanager.set_dof(dof)
coordinates = testDatamanager.create_constant_node_field("Coordinates", Float32, 2)
block_Id = testDatamanager.create_constant_node_field("Block_Id", Int64, 1)
block_Id .+= 1
testDatamanager.create_node_field("Displacements", Float32, dof)
testDatamanager.create_node_field("Forces", Float32, dof)
outputs = Dict(1 => "Displacements")
if dof == 2
    coords = vcat(coordinates[:, 1], coordinates[:, 2])
else
    coords = vcat(coordinates[:, 1], coordinates[:, 2], coordinates[:, 3])
end
exo = Write_Exodus_Results.create_result_file(filename, nnodes, dof, 1, 0)
@test exo.init.num_dim == dof
exo = init_results_in_exodus(exo, output, coords, block_Id, nsets)

close(exo)
rm(filename)
#end