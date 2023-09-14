include("../../../src/IO/mesh_data.jl")
include("../../../src/Support/data_manager.jl")
include("../../../src/Support/Parameters/parameter_handling.jl")
using Test
import .Read_Mesh
@testset "ut_neighbors" begin
    nlist = []

    for i in 1:4
        append!(nlist, [collect(1:3*i*i-2)])
    end

    lenNlist = Read_Mesh.get_number_of_neighbornodes(nlist)

    for i in 1:4
        @test lenNlist[i] == 3 * i * i - 2
    end

end

@testset "ut_glob_to_loc" begin

    distribution = [1, 2, 3, 4, 5]
    glob_to_loc = Read_Mesh.glob_to_loc(distribution)
    len = length(distribution)
    #check trivial case of global and local are identical
    for id in 1:len
        @test distribution[id] == glob_to_loc[id]
    end
    @test length(distribution) == length(glob_to_loc)
    # reverse -> glob_to_loc_to_glob
    distribution = [1, 4, 2, 5, 6]
    glob_to_loc = Read_Mesh.glob_to_loc(distribution)
    for id in 1:len
        @test distribution[id] == distribution[glob_to_loc[distribution[id]]]
    end

end

@testset "ut_define_nsets" begin

    numbers = [11, 12, 13, 44, 125]
    lenNumbers = length(numbers)
    filename = "test.txt"
    file = open(filename, "w")
    for number in numbers
        println(file, number)
    end
    close(file)
    params = Dict("Discretization" => Dict("Node Sets" => Dict("Nset_1" => "1 2 3 4 5 6 7", "Nset_2" => filename)))
    testDatamanager = Data_manager
    @test testDatamanager.get_nnsets() == 0
    Read_Mesh.define_nsets(params, testDatamanager)
    @test testDatamanager.get_nnsets() == 2
    nsets = testDatamanager.get_nsets()
    @test nsets["Nset_1"] == [1, 2, 3, 4, 5, 6, 7]
    @test nsets["Nset_2"] == [11, 12, 13, 44, 125]

    rm(filename)
end

@testset "get_bond_geometry" begin
    testDatamanager = Data_manager
    testDatamanager.set_nmasters(3)
    testDatamanager.set_dof(2)
    lenNlist = testDatamanager.create_constant_node_field("Number of Neighbors", Int64, 1)
    lenNlist[:] = [2, 2, 2]
    nlist = testDatamanager.create_constant_bond_field("Neighborhoodlist", Int64, 1)

    nlist[1] = [2, 3]
    nlist[2] = [1, 3]
    nlist[3] = [1, 2]
    coor = testDatamanager.create_constant_node_field("Coordinates", Float32, 2)
    coor[1, 1] = 0
    coor[1, 2] = 0
    coor[2, 1] = 1
    coor[2, 2] = 0
    coor[3, 1] = 0
    coor[3, 2] = 1
    Read_Mesh.get_bond_geometry(testDatamanager)
    bondgeom = testDatamanager.get_field("Bond Geometry")

    @test bondgeom[1][1, 1] == 1
    @test bondgeom[1][1, 2] == 0
    @test bondgeom[1][1, 3] == 1
    @test bondgeom[1][2, 1] == 0
    @test bondgeom[1][2, 2] == 1
    @test bondgeom[1][2, 3] == 1

    @test bondgeom[2][1, 1] == -1
    @test bondgeom[2][1, 2] == 0
    @test bondgeom[2][1, 3] == 1
    @test bondgeom[2][2, 1] == -1
    @test bondgeom[2][2, 2] == 1
    @test bondgeom[2][2, 3] / sqrt(2) - 1 < 1e-8

    @test bondgeom[3][1, 1] == 0
    @test bondgeom[3][1, 2] == -1
    @test bondgeom[3][1, 3] == 1
    @test bondgeom[3][2, 1] == 1
    @test bondgeom[3][2, 2] == -1
    @test bondgeom[3][2, 3] / sqrt(2) - 1 < 1e-8
end
