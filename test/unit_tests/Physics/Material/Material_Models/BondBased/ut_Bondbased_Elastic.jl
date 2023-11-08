# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
include("../../../../../../src/Physics/Material/BondBased/Bondbased_Elastic.jl")
include("../../../../../../src/Support/data_manager.jl")

using .Bondbased_Elastic
using Test

@testset "material_name" begin
    @test Bondbased_Elastic.material_name() == "Bond-based Elastic"
end
@testset "compute_force" begin
    nodes = 2
    testDatamanager = Data_manager
    testDatamanager.set_nmasters(nodes)
    dof = 3
    testDatamanager.set_dof(dof)
    nn = testDatamanager.create_constant_node_field("Number of Neighbors", Int64, 1)
    h = testDatamanager.create_constant_node_field("Horizon", Float64, 1)

    h[1:nodes] = 1:nodes
    bf = testDatamanager.create_constant_bond_field("Bond Forces", Float64, dof)

    nn[1] = 2
    nn[2] = 3

    bdN, bdNP1 = testDatamanager.create_bond_field("Bond Damage", Float64, 1)
    dbN, dbNP1 = testDatamanager.create_bond_field("Deformed Bond Geometry", Float64, dof + 1)
    bg = testDatamanager.create_constant_bond_field("Bond Geometry", Float64, dof + 1)
    for iID in 1:nodes
        bdNP1[iID][:] .= 1
        bg[iID][:, end] .= 1
        dbNP1[iID][:, end] .= 1 + (-1)^iID * 0.1
        dbNP1[iID][:, 1:dof] .= 1
    end

    testDatamanager = Bondbased_Elastic.compute_force(testDatamanager, Vector{Int64}(1:nodes), Dict("Bulk Modulus" => 1.0), 0.0, 0.0)

    bf = testDatamanager.get_field("Bond Forces")
    @test isapprox(bf[1][1, 1], -0.31830988618379064)
    @test isapprox(bf[1][1, 2], -0.31830988618379064)
    @test isapprox(bf[1][1, 3], -0.31830988618379064)
    @test isapprox(bf[1][2, 1], -0.31830988618379064)
    @test isapprox(bf[1][2, 2], -0.31830988618379064)
    @test isapprox(bf[1][2, 3], -0.31830988618379064)
    @test isapprox(bf[2][1, 1], 0.016277210088943853)
    @test isapprox(bf[2][1, 2], 0.016277210088943853)
    @test isapprox(bf[2][1, 3], 0.016277210088943853)
    @test isapprox(bf[2][2, 1], 0.016277210088943853)
    @test isapprox(bf[2][2, 2], 0.016277210088943853)
    @test isapprox(bf[2][2, 3], 0.016277210088943853)
    @test isapprox(bf[2][3, 1], 0.016277210088943853)
    @test isapprox(bf[2][3, 2], 0.016277210088943853)
    @test isapprox(bf[2][3, 3], 0.016277210088943853)
end