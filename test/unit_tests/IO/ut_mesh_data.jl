# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../../src/IO/mesh_data.jl")
include("../../../src/Support/data_manager.jl")
include("../../../src/Support/Parameters/parameter_handling.jl")
using Test
using .Read_Mesh
using .Data_manager
using DataFrames
@testset "ut_create_base_chunk" begin
    distribution, point_to_core = Read_Mesh.create_base_chunk(4, 1)
    @test length(distribution) == 1
    @test distribution[1] == Int64[1, 2, 3, 4]
    @test point_to_core == Int64[1, 1, 1, 1]
    distribution, point_to_core = Read_Mesh.create_base_chunk(4, 2)
    @test length(distribution) == 2
    @test distribution[1] == Int64[1, 2]
    @test distribution[2] == Int64[3, 4]
    @test point_to_core == Int64[1, 1, 2, 2]
    distribution, point_to_core = Read_Mesh.create_base_chunk(4, 3)
    @test length(distribution) == 3
    @test distribution[1] == Int64[1]
    @test distribution[2] == Int64[2]
    @test distribution[3] == Int64[3, 4]
    @test point_to_core == Int64[1, 2, 3, 3]
    distribution, point_to_core = Read_Mesh.create_base_chunk(4, 4)
    @test length(distribution) == 4
    @test distribution[1] == Int64[1]
    @test distribution[2] == Int64[2]
    @test distribution[3] == Int64[3]
    @test distribution[4] == Int64[4]
    @test point_to_core == Int64[1, 2, 3, 4]
    distribution, point_to_core = Read_Mesh.create_base_chunk(4, 5)
    @test isnothing(distribution)
    @test isnothing(point_to_core)

end


@testset "ut_check_mesh_elements" begin
    data = Dict(
        "x" => [1.0, 1.1, 3],
        "y" => [25, 30, 22],
        "volume" => [1, 1, 1],
        "block_id" => [1, 1, 1],
        "active" => [true, true, false]
    )

    df = DataFrame(data)
    meshInfoDict = Read_Mesh.check_mesh_elements(df, 2)
    @test haskey(meshInfoDict, "Coordinates")
    @test haskey(meshInfoDict, "Block_Id")
    @test haskey(meshInfoDict, "Volume")
    @test haskey(meshInfoDict, "active")
    @test meshInfoDict["Coordinates"]["Mesh ID"] == ["x", "y"]
    @test meshInfoDict["Coordinates"]["Type"] == Float64
    @test meshInfoDict["Block_Id"]["Mesh ID"] == ["block_id"]
    @test meshInfoDict["Block_Id"]["Type"] == Int64
    @test meshInfoDict["Volume"]["Mesh ID"] == ["volume"]
    @test meshInfoDict["Volume"]["Type"] == Int64
    @test meshInfoDict["active"]["Mesh ID"] == ["active"]
    @test meshInfoDict["active"]["Type"] == Bool
    data = Dict(
        "x" => [1, 1, 3],
        "y" => [25, 30, 22],
        "z" => [25, 30, 22],
        "volume" => [1.2, 0.8, 1],
        "block_id" => [1, 2, 1],
        "activex" => [true, true, false],
        "activey" => [true, true, false],
        "activez" => [true, true, false],
        "field" => [1.0, 3.3, 2.3]
    )
    df = DataFrame(data)
    meshInfoDict = Read_Mesh.check_mesh_elements(df, 3)

    @test meshInfoDict["Coordinates"]["Mesh ID"] == ["x", "y", "z"]
    @test meshInfoDict["Coordinates"]["Type"] == Int64
    @test meshInfoDict["Block_Id"]["Mesh ID"] == ["block_id"]
    @test meshInfoDict["Block_Id"]["Type"] == Int64
    @test meshInfoDict["Volume"]["Mesh ID"] == ["volume"]
    @test meshInfoDict["Volume"]["Type"] == Float64
    @test meshInfoDict["active"]["Mesh ID"] == ["activex", "activey", "activez"]
    @test meshInfoDict["active"]["Type"] == Bool
    @test meshInfoDict["field"]["Mesh ID"] == ["field"]
    @test meshInfoDict["field"]["Type"] == Float64

end
@testset "ut__init_overlap_map_" begin
    overlap_map = Read_Mesh._init_overlap_map_(1)
    @test sort(collect(keys(overlap_map))) == [1]
    @test sort(collect(keys(overlap_map[1]))) == []
    overlap_map = Read_Mesh._init_overlap_map_(2)
    @test sort(collect(keys(overlap_map))) == [1, 2]
    @test sort(collect(keys(overlap_map[1]))) == [2]
    @test sort(collect(keys(overlap_map[2]))) == [1]
    overlap_map = Read_Mesh._init_overlap_map_(3)
    @test sort(collect(keys(overlap_map))) == [1, 2, 3]
    @test sort(collect(keys(overlap_map[1]))) == [2, 3]
    @test sort(collect(keys(overlap_map[2]))) == [1, 3]
    @test sort(collect(keys(overlap_map[3]))) == [1, 2]
    overlap_map = Read_Mesh._init_overlap_map_(4)
    @test sort(collect(keys(overlap_map))) == [1, 2, 3, 4]
    @test sort(collect(keys(overlap_map[1]))) == [2, 3, 4]
    @test sort(collect(keys(overlap_map[2]))) == [1, 3, 4]
    @test sort(collect(keys(overlap_map[3]))) == [1, 2, 4]
    @test sort(collect(keys(overlap_map[4]))) == [1, 2, 3]
end

@testset "ut_create_overlap_map" begin
    distribution = [[1, 2, 3], [2, 3, 4], [4, 1, 3]]
    size = 3
    ptc = [1, 2, 2, 3]
    overlap_map = Read_Mesh.create_overlap_map(distribution, ptc, size)

    @test overlap_map[1][2]["Slave"] == overlap_map[2][1]["Master"]
    @test overlap_map[1][3]["Slave"] == overlap_map[3][1]["Master"]
    @test overlap_map[2][3]["Slave"] == overlap_map[3][2]["Master"]
    @test overlap_map[1][2]["Master"] == overlap_map[2][1]["Slave"]
    @test overlap_map[1][3]["Master"] == overlap_map[3][1]["Slave"]
    @test overlap_map[2][3]["Master"] == overlap_map[3][2]["Slave"]

    for i in 1:3
        for j in 1:3
            if i != j
                if overlap_map[i][j]["Slave"] != [] && overlap_map[i][j]["Master"] != []
                    @test overlap_map[i][j]["Slave"] != overlap_map[i][j]["Master"]
                end
            end
        end
    end
    distribution = [[1, 2, 3], [2, 3, 4], [4, 1, 3]]
    size = 3
    ptc = [1, 2, 2, 3]
    @test overlap_map[1][2]["Master"] == []
    @test overlap_map[1][2]["Slave"] == [2, 3]
    @test overlap_map[1][3]["Master"] == [1]
    @test overlap_map[1][3]["Slave"] == []
    @test overlap_map[2][3]["Master"] == [3]
    @test overlap_map[2][3]["Slave"] == [4]

end
@testset "ut_get_local_overlap_map" begin
    overlap_map = Read_Mesh._init_overlap_map_(3)
    distribution = [[1, 2, 3], [2, 3, 4], [4, 1, 3]]

    overlap_map[1][2]["Slave"] = []
    overlap_map[1][2]["Master"] = [2, 3]
    overlap_map[2][1]["Slave"] = [2, 3]
    overlap_map[2][1]["Master"] = []

    overlap_map[1][3]["Slave"] = [1]
    overlap_map[1][3]["Master"] = []
    overlap_map[3][1]["Slave"] = []
    overlap_map[3][1]["Master"] = [1]

    overlap_map[2][3]["Slave"] = [3]
    overlap_map[2][3]["Master"] = [4]
    overlap_map[3][2]["Slave"] = [4]
    overlap_map[3][2]["Master"] = [3]

    test_overlap_map = Read_Mesh.get_local_overlap_map(overlap_map, distribution, 1)
    @test test_overlap_map == overlap_map
    test_overlap_map = Read_Mesh.get_local_overlap_map(overlap_map, distribution, 3)

    @test sort(test_overlap_map[1][2]["Slave"]) == []
    @test sort(test_overlap_map[1][2]["Master"]) == [2, 3]
    @test sort(test_overlap_map[2][1]["Slave"]) == [1, 2]
    @test sort(test_overlap_map[2][1]["Master"]) == []
    @test sort(test_overlap_map[1][3]["Slave"]) == [1]
    @test sort(test_overlap_map[1][3]["Master"]) == []
    @test sort(test_overlap_map[3][1]["Slave"]) == []
    @test sort(test_overlap_map[3][1]["Master"]) == [2]
    @test sort(test_overlap_map[2][3]["Slave"]) == [2]
    @test sort(test_overlap_map[2][3]["Master"]) == [3]
    @test sort(test_overlap_map[3][2]["Slave"]) == [1]
    @test sort(test_overlap_map[3][2]["Master"]) == [3]
end

@testset "ut_neighbors" begin

    nlist = fill(Vector{Int64}([]), 4)
    for i in 1:4
        nlist[i] = Vector{Int64}(collect(1:3*i*i-2))
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
    println(file, "header: global_id")
    for number in numbers
        println(file, number)
    end
    close(file)
    params = Dict("Discretization" => Dict("Node Sets" => Dict("Nset_1" => "1 2 3 4 5 6 7", "Nset_2" => filename)))
    testDatamanager = Data_manager
    @test testDatamanager.get_nnsets() == 0
    Read_Mesh.define_nsets(params, "", testDatamanager)
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
    coor = testDatamanager.create_constant_node_field("Coordinates", Float64, 2)
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
