# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../../../../src/Physics/Material/Material_Models/Bond_Associated_Correspondence.jl")
# include("../../../../../src/Support/data_manager.jl")
using Test

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

    @test Bond_Associated_Correspondence.find_local_neighbors(coordinates[5,:], coordinates[nlist[nlist.!=5],:], nlist[nlist.!=5], bond_horizon) == []
    bond_horizon = 2.6
    @test Bond_Associated_Correspondence.find_local_neighbors(coordinates[5,:], coordinates[nlist[nlist.!=5],:], nlist[nlist.!=5], bond_horizon) == [2, 3, 4]
    @test Bond_Associated_Correspondence.find_local_neighbors(coordinates[2,:], coordinates[nlist[nlist.!=2],:], nlist[nlist.!=2], bond_horizon) == [3, 4, 5]
end

@testset "ut_compute_bond_associated_weighted_volume"
    # Initialize input arrays
    nodes = [1, 2, 3]
    nlist = [ [2, 3], [1, 3], [1, 2] ]
    coordinates = [ 0.0 0.0; 1.0 0.0; 0.0 1.0 ]
    bond_damage = [ [1.0, 1.0], [0.0, 1.0], [1.0, 1.0] ]
    omega = [ [0.5, 0.8], [0.7, 0.9], [0.6, 0.4] ]
    volume = [1.0, 2.0, 3.0]
    weighted_volume =  [ [0.0, 0.0], [0.0, 0.0], [0.0, 0.0] ]
    bond_horizon = 2.2
    # Call the function under test
    weighted_volume = Bond_Associated_Correspondence.compute_bond_associated_weighted_volume(nodes, nlist, coordinates, bond_damage, omega, volume, bond_horizon, weighted_volume)
    # Check the expected output
    @test weighted_volume ≈ [[0.7058823529411765, 0.2941176470588235], [1.0, 0.0], [0.5714285714285715, 0.4285714285714286]]
end