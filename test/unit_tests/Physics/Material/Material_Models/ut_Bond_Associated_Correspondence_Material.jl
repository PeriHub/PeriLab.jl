# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../../../../src/Physics/Material/Material_Models/Bond_Associated_Correspondence_Material.jl")
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

@test find_local_neighbors(coordinates[5,:], coordinates[nlist[nlist.!=5],:], nlist[nlist.!=5], bond_horizon) == []
bond_horizon::Float64 = 2.6
@test find_local_neighbors(coordinates[5,:], coordinates[nlist[nlist.!=5],:], nlist[nlist.!=5], bond_horizon) == [2, 3, 4]
@test find_local_neighbors(coordinates[2,:], coordinates[nlist[nlist.!=2],:], nlist[nlist.!=2], bond_horizon) == [3, 4, 5]
end