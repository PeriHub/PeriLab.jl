# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
include("../../../../../../src/Physics/Material/BondBased/Bondbased_Elastic.jl")
include("../../../../../../src/Support/data_manager.jl")

using .Bondbased_Elastic
using Test
using TimerOutputs

const to = TimerOutput()

@testset "material_name" begin
    @test Bondbased_Elastic.material_name() == "Bond-based Elastic"
end
@testset "compute_forces" begin
    nodes = 2
    test_Data_manager = Data_manager
    test_Data_manager.set_num_controller(nodes)
    dof = 3
    test_Data_manager.set_dof(dof)
    nn = test_Data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    h = test_Data_manager.create_constant_node_field("Horizon", Float64, 1)

    h[1:nodes] = 1:nodes
    bf = test_Data_manager.create_constant_bond_field("Bond Forces", Float64, dof)

    nn[1] = 2
    nn[2] = 3

    bdN, bdNP1 = test_Data_manager.create_bond_field("Bond Damage", Float64, 1)
    dbN, dbNP1 = test_Data_manager.create_bond_field("Deformed Bond Geometry", Float64, dof + 1)
    bg = test_Data_manager.create_constant_bond_field("Bond Geometry", Float64, dof + 1)
    for iID in 1:nodes
        bdNP1[iID][:] .= 1
        bg[iID][:, end] .= 1
        dbNP1[iID][:, end] .= 1 + (-1)^iID * 0.1
        dbNP1[iID][:, 1:dof] .= 1
    end

    test_Data_manager = Bondbased_Elastic.compute_forces(test_Data_manager, Vector{Int64}(1:nodes), Dict("Bulk Modulus" => 1.0, "Young's Modulus" => 1.0), 0.0, 0.0, to)

    bf = test_Data_manager.get_field("Bond Forces")
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
    bf[1][:, :] .= 0
    bf[2][:, :] .= 0

    test_Data_manager = Bondbased_Elastic.compute_forces(test_Data_manager, Vector{Int64}(1:nodes), Dict("Bulk Modulus" => 1.0, "Young's Modulus" => 1.0, "Symmetry" => "here is something"), 0.0, 0.0, to)

    bf = test_Data_manager.get_field("Bond Forces")
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
    bf[1][:, :] .= 0
    bf[2][:, :] .= 0
    test_Data_manager = Bondbased_Elastic.compute_forces(test_Data_manager, Vector{Int64}(1:nodes), Dict("Bulk Modulus" => 1.0, "Young's Modulus" => 1.0, "Symmetry" => "plane strain"), 0.0, 0.0, to)

    bf = test_Data_manager.get_field("Bond Forces")
    @test isapprox(bf[1][1, 1], -0.16976527263135502)
    @test isapprox(bf[1][1, 2], -0.16976527263135502)
    @test isapprox(bf[1][2, 1], -0.16976527263135502)
    @test isapprox(bf[1][2, 2], -0.16976527263135502)
    @test isapprox(bf[2][1, 1], 0.01736235742820678)
    @test isapprox(bf[2][1, 2], 0.01736235742820678)
    @test isapprox(bf[2][2, 1], 0.01736235742820678)
    @test isapprox(bf[2][2, 2], 0.01736235742820678)

    bf[1][:, :] .= 0
    bf[2][:, :] .= 0
    test_Data_manager = Bondbased_Elastic.compute_forces(test_Data_manager, Vector{Int64}(1:nodes), Dict("Bulk Modulus" => 1.0, "Young's Modulus" => 1.0, "Symmetry" => "plane stress"), 0.0, 0.0, to)

    bf = test_Data_manager.get_field("Bond Forces")
    @test isapprox(bf[1][1, 1], -0.15915494309189532)
    @test isapprox(bf[1][1, 2], -0.15915494309189532)
    @test isapprox(bf[1][2, 1], -0.15915494309189532)
    @test isapprox(bf[1][2, 2], -0.15915494309189532)
    @test isapprox(bf[2][1, 1], 0.016277210088943856)
    @test isapprox(bf[2][1, 2], 0.016277210088943856)
    @test isapprox(bf[2][2, 1], 0.016277210088943856)
    @test isapprox(bf[2][2, 2], 0.016277210088943856)

    dbNP1[1][1, end] = 0
    @test isnothing(Bondbased_Elastic.compute_forces(test_Data_manager, Vector{Int64}(1:nodes), Dict("Bulk Modulus" => 1.0, "Young's Modulus" => 1.0, "Symmetry" => "plane stress"), 0.0, 0.0))
end