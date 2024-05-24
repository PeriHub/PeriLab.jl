# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../../../../src/Physics/Material/Material_Models/Bond_Associated_Correspondence.jl")
# include("../../../../../src/Core/data_manager.jl")
using Test
using TimerOutputs

const to = TimerOutput()

@testset "correspondence_name" begin
    @test Bond_Associated_Correspondence.correspondence_name() == "Correspondence Bond-Associated"
end
# @testset "compute_forces" begin
#     nodes = 3
#     test_Data_manager = Data_manager
#     test_Data_manager.clear_data_manager()
#     test_Data_manager.set_num_controller(nodes)
#     dof = 3
#     test_Data_manager.set_dof(dof)
#     nn = test_Data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
#     nn[1] = 2
#     nn[2] = 2
#     nn[3] = 2
#     nlist = test_Data_manager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
#     nlist[1] = [2, 3]
#     nlist[2] = [1, 3]
#     nlist[2] = [1, 2]
#     h = test_Data_manager.create_constant_node_field("Horizon", Float64, 1)

#     h[1:nodes] = 1:nodes
#     bf = test_Data_manager.create_constant_bond_field("Bond Forces", Float64, dof)

#     coor = test_Data_manager.create_constant_node_field("Coordinates", Float64, 2)
#     coor[1, 1] = 0
#     coor[1, 2] = 0
#     coor[2, 1] = 1
#     coor[2, 2] = 0
#     coor[3, 1] = 0
#     coor[3, 2] = 1
#     omega = test_Data_manager.create_constant_node_field("Influence Function", Float64, 2)
#     omega[1, 1] = 0.5
#     omega[1, 2] = 0.8
#     omega[2, 1] = 0.7
#     omega[2, 2] = 0.9
#     omega[3, 1] = 0.6
#     omega[3, 2] = 0.4
#     volume = test_Data_manager.create_constant_node_field("Volume", Float64, 1)
#     volume[1:nodes] = [1.0, 2.0, 3.0]
#     bw_volume = test_Data_manager.create_constant_node_field("Bond Weighted Volume", Float64, 2)
#     bw_volume[1, 1] = 0
#     bw_volume[1, 2] = 0
#     bw_volume[2, 1] = 0
#     bw_volume[2, 2] = 0
#     bw_volume[3, 1] = 0
#     bw_volume[3, 2] = 0

#     bdN, bdNP1 = test_Data_manager.create_bond_field("Bond Damage", Float64, 1)
#     dbN, dbNP1 = test_Data_manager.create_bond_field("Deformed Bond Geometry", Float64, dof)
#     dbdN, dbdNP1 = test_Data_manager.create_bond_field("Deformed Bond Length", Float64, 1)
#     bg = test_Data_manager.create_constant_bond_field("Bond Geometry", Float64, dof)
#     bd = test_Data_manager.create_constant_bond_field("Bond Length", Float64, 1)
#     for iID in 1:nodes
#         bdNP1[iID] .= 1
#         bd[iID] .= 1
#         dbdNP1[iID] .= 1 + (-1)^iID * 0.1
#         dbNP1[iID] .= 1
#     end

#     test_Data_manager = Bond_Associated_Correspondence.compute_forces(test_Data_manager, Vector{Int64}(1:nodes), Dict("Bulk Modulus" => 1.0, "Young's Modulus" => 1.0), 0.0, 0.0, to)

#     bf = test_Data_manager.get_field("Bond Forces")
#     @test isapprox(bf[1][1, 1], -0.31830988618379064)
#     @test isapprox(bf[1][1, 2], -0.31830988618379064)
#     @test isapprox(bf[1][1, 3], -0.31830988618379064)
#     @test isapprox(bf[1][2, 1], -0.31830988618379064)
#     @test isapprox(bf[1][2, 2], -0.31830988618379064)
#     @test isapprox(bf[1][2, 3], -0.31830988618379064)
#     @test isapprox(bf[2][1, 1], 0.016277210088943853)
#     @test isapprox(bf[2][1, 2], 0.016277210088943853)
#     @test isapprox(bf[2][1, 3], 0.016277210088943853)
#     @test isapprox(bf[2][2, 1], 0.016277210088943853)
#     @test isapprox(bf[2][2, 2], 0.016277210088943853)
#     @test isapprox(bf[2][2, 3], 0.016277210088943853)
#     @test isapprox(bf[2][3, 1], 0.016277210088943853)
#     @test isapprox(bf[2][3, 2], 0.016277210088943853)
#     @test isapprox(bf[2][3, 3], 0.016277210088943853)
# end
@testset "ut_find_local_neighbors" begin
    coordinates = zeros(Float64, 5, 2)
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

    nlist = [2, 3, 4, 5]

    bond_horizon::Float64 = 1

    @test Bond_Associated_Correspondence.find_local_neighbors(5, coordinates, nlist, bond_horizon) == []
    bond_horizon = 2.6
    @test Bond_Associated_Correspondence.find_local_neighbors(5, coordinates, nlist, bond_horizon) == [2, 3, 4]
    @test Bond_Associated_Correspondence.find_local_neighbors(2, coordinates, nlist, bond_horizon) == [3, 4, 5]
end

@testset "ut_compute_bond_associated_weighted_volume" begin
    # Initialize input arrays
    nodes = [1, 2, 3]
    nlist = [[2, 3], [1, 3], [1, 2]]
    coordinates = [0.0 0.0; 1.0 0.0; 0.0 1.0]
    bond_damage = [[1.0, 1.0], [0.0, 1.0], [1.0, 1.0]]
    omega = [[0.5, 0.8], [0.7, 0.9], [0.6, 0.4]]
    volume = [1.0, 2.0, 3.0]
    weighted_volume = [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
    bond_horizon = 2.2
    # Call the function under test
    weighted_volume = Bond_Associated_Correspondence.compute_bond_associated_weighted_volume(nodes, nlist, coordinates, bond_damage, omega, volume, bond_horizon, weighted_volume)
    # Check the expected output
    @test weighted_volume ≈ [[0.7058823529411765, 0.2941176470588235], [1.0, 0.0], [0.5714285714285715, 0.4285714285714286]]
end