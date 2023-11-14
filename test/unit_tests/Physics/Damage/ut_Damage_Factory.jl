# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../../../src/Physics/Damage/Damage_Factory.jl")
include("../../../../src/Support/data_manager.jl")
using Test
using .Damage
@testset "damage_index" begin
    test_Data_manager = Data_manager
    test_Data_manager.set_num_controller(3)
    nn = test_Data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1] = 1
    nn[2] = 2
    nn[3] = 1
    damageN, damageNP1_test = test_Data_manager.create_node_field("Damage", Float64, 1)
    volume = test_Data_manager.create_constant_node_field("Volume", Float64, 1)
    nlist = test_Data_manager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
    bdN, bdNP1 = test_Data_manager.create_bond_field("Bond Damage", Float64, 1, 1)
    nlist[1] = [2]
    nlist[2] = [1, 3]
    nlist[3] = [1]
    volume[:] = [1.0, 2.0, 3.0]
    @test bdN[1][:] == [1.0]
    @test bdN[2][:] == [1.0, 1.0]
    @test bdN[3][:] == [1.0]
    @test bdNP1[1][:] == [1.0]
    @test bdNP1[2][:] == [1.0, 1.0]
    @test bdNP1[3][:] == [1.0]
    nodes = view(Vector(1:3), eachindex(Vector(1:3)))
    Damage.damage_index(test_Data_manager, nodes)
    @test damageNP1_test[1] == 0
    @test damageNP1_test[2] == 0
    @test damageNP1_test[3] == 0
    bdNP1[1][:] .= 0
    Damage.damage_index(test_Data_manager, nodes)
    @test damageNP1_test[1] == 1
    @test damageNP1_test[2] == 0
    @test damageNP1_test[3] == 0
    bdNP1[2][1] = 0
    Damage.damage_index(test_Data_manager, nodes)
    @test damageNP1_test[1] == 1
    @test damageNP1_test[2] == 0.25
    @test damageNP1_test[3] == 0
    bdNP1[2][1] = 1
    bdNP1[2][2] = 0
    Damage.damage_index(test_Data_manager, nodes)
    @test damageNP1_test[1] == 1
    @test damageNP1_test[2] == 0.75
    @test damageNP1_test[3] == 0
    bdNP1[2][1] = 0
    bdNP1[2][2] = 0
    Damage.damage_index(test_Data_manager, nodes)
    @test damageNP1_test[1] == 1
    @test damageNP1_test[2] == 1
    @test damageNP1_test[3] == 0
    bdNP1[3][:] .= 0
    Damage.damage_index(test_Data_manager, nodes)
    @test damageNP1_test[1] == 1
    @test damageNP1_test[2] == 1
    @test damageNP1_test[3] == 1
end
@testset "set_bond_damage" begin
    test_Data_manager = Data_manager
    test_Data_manager.set_num_controller(3)
    nn = test_Data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1] = 1
    nn[2] = 2
    nn[3] = 1
    bdN, bdNP1 = test_Data_manager.create_bond_field("Bond Damage", Float64, 1)

    for iID in 1:3
        for jID in 1:nn[iID]
            bdN[iID][jID] = jID * iID * iID + jID - iID
            @test bdNP1[iID][jID] == 0
        end
    end
    nodes = Vector{Int64}(1:3)
    Damage.set_bond_damage(test_Data_manager, nodes)
    bdN = test_Data_manager.get_field("Bond Damage", "N")
    bdNP1 = test_Data_manager.get_field("Bond Damage", "NP1")
    for iID in 1:3
        for jID in 1:nn[iID]
            @test bdNP1[iID][jID] == bdN[iID][jID]
            bdN[iID][jID] = 0
            @test bdNP1[iID][jID] == jID * iID * iID + jID - iID
            @test bdN[iID][jID] == 0
        end
    end
end

@testset "Damage_factory_exceptions" begin
    test_Data_manager = Data_manager
    @test isnothing(Damage.compute_damage_pre_calculation(test_Data_manager, Vector{Int64}(1:3), 1, Dict("Damage Model" => "not there"), 1, 0.0, 0.0))
    @test isnothing(Damage.compute_damage(test_Data_manager, Vector{Int64}(1:3), Dict("Damage Model" => "not there"), 1, 0.0, 0.0))
end
