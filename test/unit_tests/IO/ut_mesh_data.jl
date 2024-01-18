# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../../src/IO/IO.jl")
include("../../../src/Support/data_manager.jl")
include("../../../src/Support/Parameters/parameter_handling.jl")
using Test
using TimerOutputs
using Reexport
@reexport using .Parameter_Handling
@reexport using .Data_manager
@reexport using .IO
using DataFrames
@testset "ut_read_mesh" begin
    params = Dict("Discretization" => Dict("Type" => "not supported"))
    @test isnothing(IO.read_mesh("./", params))

    params = Dict("Discretization" => Dict("Type" => "Text File"))
    @test isnothing(IO.read_mesh("./", params))
    path = "./unit_tests/IO/"

    data = IO.read_mesh(joinpath(path, "example_mesh.txt"), params)
    if isnothing(data)
        path = "./test/unit_tests/IO/"
        data = IO.read_mesh(joinpath(path, "example_mesh.txt"), params)
    end
    @test length(data[:, 1]) == 3
    @test data[!, "x"] == [-1.5, -0.5, 0.5]
    @test data[!, "y"] == [0.5, 0.2, 0.0]
    @test data[!, "block_id"] == [1, 1, 2]
    @test data[!, "volume"] == [0.1, 0.1, 0.1]

    @test isnothing(IO.read_external_topology("./"))
    data = IO.read_external_topology(joinpath(path, "example_FE_mesh.txt"))
    @test length(data[:, 1]) == 4
    @test collect(skipmissing(data[1, :])) == [1, 2, 3, 4]
    @test collect(skipmissing(data[2, :])) == [3, 4, 5, 6]
    @test collect(skipmissing(data[3, :])) == [2, 3, 4, 5, 6, 7]
    @test collect(skipmissing(data[4, :])) == [5, 6, 7, 8]

end

@testset "ut_check_for_duplicate_in_dataframe" begin
    path = "./unit_tests/IO/"
    params = Dict("Discretization" => Dict("Type" => "Text File"))
    data = IO.read_mesh(joinpath(path, "example_mesh.txt"), params)
    if isnothing(data)
        path = "./test/unit_tests/IO/"
        data = IO.read_mesh(joinpath(path, "example_mesh.txt"), params)
    end
    @test !(IO.check_for_duplicate_in_dataframe(data))
    data[1, :] = data[2, :]
    @test IO.check_for_duplicate_in_dataframe(data)
    data[1, :] = data[3, :]
    @test IO.check_for_duplicate_in_dataframe(data)
    data[2, :] = data[3, :]
    @test IO.check_for_duplicate_in_dataframe(data)
end


@testset "ut_create_consistent_neighborhoodlist" begin
    path = "./unit_tests/IO/"

    external_topology = IO.read_external_topology(joinpath(path, "example_FE_mesh.txt"))
    if isnothing(external_topology)
        path = "./test/unit_tests/IO/"
        external_topology = IO.read_external_topology(joinpath(path, "example_FE_mesh.txt"))
    end
    dof::Int64 = 2
    params = Dict()

    nlist::Vector{Vector{Int64}} = [[2, 3, 4, 11], [1, 3, 4], [1, 2, 22, 23], [4], [8], [9], [1, 6], [3, 2], [10]]

    nlist_test, topology, nodes_to_element = IO.create_consistent_neighborhoodlist(external_topology, params, nlist, dof)

    @test topology[1] == [1, 2, 3, 4]
    @test topology[2] == [3, 4, 5, 6]
    @test topology[3] == [2, 3, 4, 5, 6, 7]
    @test topology[4] == [5, 6, 7, 8]
    @test nodes_to_element == [[1], [1, 3], [1, 2, 3], [1, 2, 3], [2, 3, 4], [2, 3, 4], [3, 4], [4]]
    @test nlist == [[2, 3, 4], [1, 3, 4, 5, 6, 7], [1, 2, 4, 5, 6, 7], [1, 2, 3, 5, 6, 7], [3, 4, 6, 2, 7, 8], [3, 4, 5, 2, 7, 8], [2, 3, 4, 5, 6, 8], [5, 6, 7], [10]]
    params = Dict("Add Neighbor Search" => false)
    nlist = [[2, 3, 4, 11], [1, 3, 4], [1, 2, 22, 23], [4], [8], [9], [1, 6], [3, 2], [10]]
    nlist, topology, nodes_to_element = IO.create_consistent_neighborhoodlist(external_topology, params, nlist, dof)
    @test nlist == [[2, 3, 4], [1, 3, 4, 5, 6, 7], [1, 2, 4, 5, 6, 7], [1, 2, 3, 5, 6, 7], [3, 4, 6, 2, 7, 8], [3, 4, 5, 2, 7, 8], [2, 3, 4, 5, 6, 8], [5, 6, 7], [10]]

    params = Dict("Add Neighbor Search" => true)
    nlist = [[2, 3, 4, 11], [1, 3, 4], [1, 2, 22, 23], [4], [8], [9], [1, 6], [3, 2], [10]]
    nlist, topology, nodes_to_element = IO.create_consistent_neighborhoodlist(external_topology, params, nlist, dof)

    @test topology[1] == [1, 2, 3, 4]
    @test topology[2] == [3, 4, 5, 6]
    @test topology[3] == [2, 3, 4, 5, 6, 7]
    @test topology[4] == [5, 6, 7, 8]
    @test nodes_to_element == [[1], [1, 3], [1, 2, 3], [1, 2, 3], [2, 3, 4], [2, 3, 4], [3, 4], [4]]
    @test nlist == [[2, 3, 4, 11], [1, 3, 4, 5, 6, 7], [1, 2, 22, 23, 4, 5, 6, 7], [1, 2, 3, 5, 6, 7], [8, 3, 4, 6, 2, 7], [9, 3, 4, 5, 2, 7, 8], [1, 6, 2, 3, 4, 5, 8], [3, 2, 5, 6, 7], [10]]
end
@testset "ut_element_distribution" begin
    nnodes = 8
    topology = Vector([[1, 2, 3, 4], [2, 4, 5, 6], [7, 8, 5, 6]])

    ptc::Vector{Int64} = zeros(8)
    ptc[:] .= 1
    ranksize = 1
    element_distribution = IO.element_distribution(topology, ptc, ranksize)
    @test length(element_distribution) == 1
    @test element_distribution[1] == [1, 2, 3]
    ptc[5:8] .= 2
    ranksize = 2
    element_distribution = IO.element_distribution(topology, ptc, ranksize)
    @test length(element_distribution) == 2
    @test element_distribution[1] == [1]
    @test element_distribution[2] == [2, 3]

    topology = Vector([[7, 8, 5, 6], [1, 2, 3, 4], [2, 4, 5, 6]])
    element_distribution = IO.element_distribution(topology, ptc, ranksize)

    @test length(element_distribution) == 2
    @test element_distribution[1] == [2]
    @test element_distribution[2] == [3, 1]

end
@testset "ut_get_local_element_topology" begin

    test_Data_manager = Data_manager
    topology::Vector{Vector{Int64}} = [[1, 2, 3, 4], [3, 4, 2, 1]]
    distribution::Vector{Vector{Int64}} = [[2, 3, 4, 1], [1, 2, 3, 4], [5, 6, 3, 2, 8, 1, 4]]
    test_Data_manager = IO.get_local_element_topology(test_Data_manager, topology, distribution[1])
    topo = test_Data_manager.get_field("FE Topology")

    @test topo[1, 1] == 4
    @test topo[1, 2] == 1
    @test topo[1, 3] == 2
    @test topo[1, 4] == 3
    @test topo[2, 1] == 2
    @test topo[2, 2] == 3
    @test topo[2, 3] == 1
    @test topo[2, 4] == 4
    test_Data_manager = IO.get_local_element_topology(test_Data_manager, topology, distribution[2])
    topo = test_Data_manager.get_field("FE Topology")
    @test topo[1, 1] == 1
    @test topo[1, 2] == 2
    @test topo[1, 3] == 3
    @test topo[1, 4] == 4
    @test topo[2, 1] == 3
    @test topo[2, 2] == 4
    @test topo[2, 3] == 2
    @test topo[2, 4] == 1
    test_Data_manager = IO.get_local_element_topology(test_Data_manager, topology, distribution[3])
    topo = test_Data_manager.get_field("FE Topology")
    @test topo[1, 1] == 6
    @test topo[1, 2] == 4
    @test topo[1, 3] == 3
    @test topo[1, 4] == 7
    @test topo[2, 1] == 3
    @test topo[2, 2] == 7
    @test topo[2, 3] == 4
    @test topo[2, 4] == 6

    test_Data_manager = IO.get_local_element_topology(test_Data_manager, Vector([Vector{Int64}([])]), distribution[3])
    # nothing happens, because no field is initialized
    topo = test_Data_manager.get_field("FE Topology")
    @test topo[1, 1] == 6
    @test topo[1, 2] == 4
    @test topo[1, 3] == 3
    @test topo[1, 4] == 7
    @test topo[2, 1] == 3
    @test topo[2, 2] == 7
    @test topo[2, 3] == 4
    @test topo[2, 4] == 6

    topology = [[1, 2, 3, 4], [3, 4, 2, 1, 3]]

    @test isnothing(IO.get_local_element_topology(test_Data_manager, topology, distribution[3]))
end
@testset "ut_create_base_chunk" begin
    distribution, point_to_core = IO.create_base_chunk(4, 1)
    @test length(distribution) == 1
    @test distribution[1] == Int64[1, 2, 3, 4]
    @test point_to_core == Int64[1, 1, 1, 1]
    distribution, point_to_core = IO.create_base_chunk(4, 2)
    @test length(distribution) == 2
    @test distribution[1] == Int64[1, 2]
    @test distribution[2] == Int64[3, 4]
    @test point_to_core == Int64[1, 1, 2, 2]
    distribution, point_to_core = IO.create_base_chunk(4, 3)
    @test length(distribution) == 3
    @test distribution[1] == Int64[1]
    @test distribution[2] == Int64[2]
    @test distribution[3] == Int64[3, 4]
    @test point_to_core == Int64[1, 2, 3, 3]
    distribution, point_to_core = IO.create_base_chunk(4, 4)
    @test length(distribution) == 4
    @test distribution[1] == Int64[1]
    @test distribution[2] == Int64[2]
    @test distribution[3] == Int64[3]
    @test distribution[4] == Int64[4]
    @test point_to_core == Int64[1, 2, 3, 4]
    distribution, point_to_core = IO.create_base_chunk(4, 5)
    @test isnothing(distribution)
    @test isnothing(point_to_core)

end

@testset "ut_local_nodes_from_dict" begin
    glob_to_loc = Dict{Int64,Int64}(1 => 2, 2 => 4, 3 => 3, 4 => 1)
    global_nodes = Vector{Int64}(1:4)
    test = IO.local_nodes_from_dict(glob_to_loc, global_nodes)
    @test test == [2, 4, 3, 1]
end
@testset "ut_check_mesh_elements" begin
    data = Dict(
        "volume" => [1, 1, 1],
        "block_id" => [1, 1, 1]
    )
    df = DataFrame(data)
    @test isnothing(IO.check_mesh_elements(df, 2))
    data = Dict(
        "x" => [1.0, 1.1, 3],
        "y" => [25, 30, 22],
        "block_id" => [1, 1, 1]
    )
    df = DataFrame(data)
    @test isnothing(IO.check_mesh_elements(df, 2))
    data = Dict(
        "x" => [1.0, 1.1, 3],
        "y" => [25, 30, 22],
        "volume" => [1, 1, 1]
    )
    df = DataFrame(data)
    @test isnothing(IO.check_mesh_elements(df, 2))

    data = Dict(
        "x" => [1.0, 1.1, 3],
        "y" => [25, 30, 22],
        "volume" => [1, 1, 1],
        "block_id" => [1, 1, 1],
        "active" => [true, true, false]
    )

    df = DataFrame(data)
    meshInfoDict = IO.check_mesh_elements(df, 2)
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
        "active_x" => [true, true, false],
        "active_y" => [true, true, false],
        "active_z" => [true, true, false],
        "field" => [1.0, 3.3, 2.3]
    )
    df = DataFrame(data)
    meshInfoDict = IO.check_mesh_elements(df, 3)
    @test meshInfoDict["Coordinates"]["Mesh ID"] == ["x", "y", "z"]
    @test meshInfoDict["Coordinates"]["Type"] == Int64
    @test meshInfoDict["Block_Id"]["Mesh ID"] == ["block_id"]
    @test meshInfoDict["Block_Id"]["Type"] == Int64
    @test meshInfoDict["Volume"]["Mesh ID"] == ["volume"]
    @test meshInfoDict["Volume"]["Type"] == Float64
    @test meshInfoDict["active"]["Mesh ID"] == ["active_x", "active_y", "active_z"]
    @test meshInfoDict["active"]["Type"] == Bool
    @test meshInfoDict["field"]["Mesh ID"] == ["field"]
    @test meshInfoDict["field"]["Type"] == Float64

end
@testset "ut__init_overlap_map_" begin
    overlap_map = IO._init_overlap_map_(1)
    @test sort(collect(keys(overlap_map))) == [1]
    @test sort(collect(keys(overlap_map[1]))) == []
    overlap_map = IO._init_overlap_map_(2)
    @test sort(collect(keys(overlap_map))) == [1, 2]
    @test sort(collect(keys(overlap_map[1]))) == [2]
    @test sort(collect(keys(overlap_map[2]))) == [1]
    overlap_map = IO._init_overlap_map_(3)
    @test sort(collect(keys(overlap_map))) == [1, 2, 3]
    @test sort(collect(keys(overlap_map[1]))) == [2, 3]
    @test sort(collect(keys(overlap_map[2]))) == [1, 3]
    @test sort(collect(keys(overlap_map[3]))) == [1, 2]
    overlap_map = IO._init_overlap_map_(4)
    @test sort(collect(keys(overlap_map))) == [1, 2, 3, 4]
    @test sort(collect(keys(overlap_map[1]))) == [2, 3, 4]
    @test sort(collect(keys(overlap_map[2]))) == [1, 3, 4]
    @test sort(collect(keys(overlap_map[3]))) == [1, 2, 4]
    @test sort(collect(keys(overlap_map[4]))) == [1, 2, 3]
end

@testset "ut_create_overlap_map" begin
    distribution = Vector([Vector{Int64}([1, 2, 3]), Vector{Int64}([2, 3, 4]), Vector{Int64}([4, 1, 3])])
    size::Int64 = 3
    ptc = [1, 2, 2, 3]
    overlap_map = IO.create_overlap_map(distribution, ptc, size)

    @test overlap_map[1][2]["Responder"] == overlap_map[2][1]["Controller"]
    @test overlap_map[1][3]["Responder"] == overlap_map[3][1]["Controller"]
    @test overlap_map[2][3]["Responder"] == overlap_map[3][2]["Controller"]
    @test overlap_map[1][2]["Controller"] == overlap_map[2][1]["Responder"]
    @test overlap_map[1][3]["Controller"] == overlap_map[3][1]["Responder"]
    @test overlap_map[2][3]["Controller"] == overlap_map[3][2]["Responder"]

    for i in 1:3
        for j in 1:3
            if i != j
                if overlap_map[i][j]["Responder"] != [] && overlap_map[i][j]["Controller"] != []
                    @test overlap_map[i][j]["Responder"] != overlap_map[i][j]["Controller"]
                end
            end
        end
    end
    distribution = [[1, 2, 3], [2, 3, 4], [4, 1, 3]]
    size = 3
    ptc = [1, 2, 2, 3]
    @test overlap_map[1][2]["Controller"] == []
    @test overlap_map[1][2]["Responder"] == [2, 3]
    @test overlap_map[1][3]["Controller"] == [1]
    @test overlap_map[1][3]["Responder"] == []
    @test overlap_map[2][3]["Controller"] == [3]
    @test overlap_map[2][3]["Responder"] == [4]

end
@testset "ut_get_local_overlap_map" begin
    overlap_map = IO._init_overlap_map_(3)
    distribution = [[1, 2, 3], [2, 3, 4], [4, 1, 3]]

    overlap_map[1][2]["Responder"] = []
    overlap_map[1][2]["Controller"] = [2, 3]
    overlap_map[2][1]["Responder"] = [2, 3]
    overlap_map[2][1]["Controller"] = []

    overlap_map[1][3]["Responder"] = [1]
    overlap_map[1][3]["Controller"] = []
    overlap_map[3][1]["Responder"] = []
    overlap_map[3][1]["Controller"] = [1]

    overlap_map[2][3]["Responder"] = [3]
    overlap_map[2][3]["Controller"] = [4]
    overlap_map[3][2]["Responder"] = [4]
    overlap_map[3][2]["Controller"] = [3]

    test_overlap_map = IO.get_local_overlap_map(overlap_map, distribution, 1)
    @test test_overlap_map == overlap_map
    test_overlap_map = IO.get_local_overlap_map(overlap_map, distribution, 3)

    @test sort(test_overlap_map[1][2]["Responder"]) == []
    @test sort(test_overlap_map[1][2]["Controller"]) == [2, 3]
    @test sort(test_overlap_map[2][1]["Responder"]) == [1, 2]
    @test sort(test_overlap_map[2][1]["Controller"]) == []
    @test sort(test_overlap_map[1][3]["Responder"]) == [1]
    @test sort(test_overlap_map[1][3]["Controller"]) == []
    @test sort(test_overlap_map[3][1]["Responder"]) == []
    @test sort(test_overlap_map[3][1]["Controller"]) == [2]
    @test sort(test_overlap_map[2][3]["Responder"]) == [2]
    @test sort(test_overlap_map[2][3]["Controller"]) == [3]
    @test sort(test_overlap_map[3][2]["Responder"]) == [1]
    @test sort(test_overlap_map[3][2]["Controller"]) == [3]
end

@testset "ut_neighbors" begin

    nlist = fill(Vector{Int64}([]), 4)
    for i in 1:4
        nlist[i] = Vector{Int64}(collect(1:3*i*i-2))
    end

    lenNlist = IO.get_number_of_neighbornodes(nlist)

    for i in 1:4
        @test lenNlist[i] == 3 * i * i - 2
    end
    nlist[1] = []

    @test isnothing(IO.get_number_of_neighbornodes(nlist))
end

@testset "ut_glob_to_loc" begin

    distribution = [1, 2, 3, 4, 5]
    glob_to_loc = IO.glob_to_loc(distribution)
    len = length(distribution)
    #check trivial case of global and local are identical
    for id in 1:len
        @test distribution[id] == glob_to_loc[id]
    end
    @test length(distribution) == length(glob_to_loc)
    # reverse -> glob_to_loc_to_glob
    distribution = [1, 4, 2, 5, 6]
    glob_to_loc = IO.glob_to_loc(distribution)
    for id in 1:len
        @test distribution[id] == distribution[glob_to_loc[distribution[id]]]
    end

end

@testset "ut_define_nsets" begin

    nsets_predef = Dict{String,Any}("Nset_2" => [11, 12, 13, 44, 125], "Nset_1" => [1, 2, 3, 4, 5, 6, 7])

    test_Data_manager = Data_manager
    @test test_Data_manager.get_nnsets() == 0
    IO.define_nsets(nsets_predef, test_Data_manager)
    @test test_Data_manager.get_nnsets() == 2
    nsets = test_Data_manager.get_nsets()
    @test nsets["Nset_1"] == [1, 2, 3, 4, 5, 6, 7]
    @test nsets["Nset_2"] == [11, 12, 13, 44, 125]

end

@testset "get_bond_geometry" begin
    test_Data_manager = Data_manager
    test_Data_manager.set_num_controller(3)
    test_Data_manager.set_dof(2)
    lenNlist = test_Data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    lenNlist[:] = [2, 2, 2]
    nlist = test_Data_manager.create_constant_bond_field("Neighborhoodlist", Int64, 1)

    nlist[1] = [2, 3]
    nlist[2] = [1, 3]
    nlist[3] = [1, 2]
    coor = test_Data_manager.create_constant_node_field("Coordinates", Float64, 2)
    coor[1, 1] = 0
    coor[1, 2] = 0
    coor[2, 1] = 1
    coor[2, 2] = 0
    coor[3, 1] = 0
    coor[3, 2] = 1
    IO.get_bond_geometry(test_Data_manager)
    undeformed_bond = test_Data_manager.get_field("Bond Geometry")

    @test undeformed_bond[1][1, 1] == 1
    @test undeformed_bond[1][1, 2] == 0
    @test undeformed_bond[1][1, 3] == 1
    @test undeformed_bond[1][2, 1] == 0
    @test undeformed_bond[1][2, 2] == 1
    @test undeformed_bond[1][2, 3] == 1

    @test undeformed_bond[2][1, 1] == -1
    @test undeformed_bond[2][1, 2] == 0
    @test undeformed_bond[2][1, 3] == 1
    @test undeformed_bond[2][2, 1] == -1
    @test undeformed_bond[2][2, 2] == 1
    @test undeformed_bond[2][2, 3] / sqrt(2) - 1 < 1e-8

    @test undeformed_bond[3][1, 1] == 0
    @test undeformed_bond[3][1, 2] == -1
    @test undeformed_bond[3][1, 3] == 1
    @test undeformed_bond[3][2, 1] == 1
    @test undeformed_bond[3][2, 2] == -1
    @test undeformed_bond[3][2, 3] / sqrt(2) - 1 < 1e-8
end
