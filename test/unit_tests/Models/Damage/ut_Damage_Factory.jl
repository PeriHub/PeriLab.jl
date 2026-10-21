# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

#using Test

@testset "init_fields" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_dof(3)
    PeriLab.Data_Manager.set_num_controller(4)
    nn = PeriLab.Data_Manager.create_constant_node_scalar_field("Number of Neighbors",
                                                                Int64)
    nn[1] = 1
    nn[2] = 2
    nn[3] = 1
    nn[4] = 2
    nlist = PeriLab.Data_Manager.create_constant_bond_scalar_state("Neighborhoodlist",
                                                                   Int64)
    nlist[1] = [2]
    nlist[2] = [1, 3]
    nlist[3] = [1]
    nlist[4] = [1, 3]
    PeriLab.Solver_Manager.Model_Factory.Damage.init_fields()
    field_keys = PeriLab.Data_Manager.get_all_field_keys()
    @test "DamageN" in field_keys
    @test "DamageNP1" in field_keys
end

@testset "damage_index" begin
    PeriLab.Data_Manager.set_num_controller(3)
    nn = PeriLab.Data_Manager.create_constant_node_scalar_field("Number of Neighbors",
                                                                Int64)
    nn[1] = 1
    nn[2] = 2
    nn[3] = 1
    damageN,
    damageNP1_test = PeriLab.Data_Manager.create_node_scalar_field("Damage",
                                                                   Float64)
    volume = PeriLab.Data_Manager.create_constant_node_scalar_field("Volume", Float64)
    nlist = PeriLab.Data_Manager.create_constant_bond_scalar_state("Neighborhoodlist",
                                                                   Int64)
    bdN,
    bdNP1 = PeriLab.Data_Manager.create_bond_scalar_state("Bond Damage", Float64;
                                                          default_value = 1)
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
    PeriLab.Solver_Manager.Model_Factory.Damage.damage_index(nodes)
    @test damageNP1_test[1] == 0
    @test damageNP1_test[2] == 0
    @test damageNP1_test[3] == 0
    bdNP1[1][:] .= 0
    PeriLab.Solver_Manager.Model_Factory.Damage.damage_index(nodes)
    @test damageNP1_test[1] == 1
    @test damageNP1_test[2] == 0
    @test damageNP1_test[3] == 0
    bdNP1[2][1] = 0
    PeriLab.Solver_Manager.Model_Factory.Damage.damage_index(nodes)
    @test damageNP1_test[1] == 1
    @test damageNP1_test[2] == 0.25
    @test damageNP1_test[3] == 0
    bdNP1[2][1] = 1
    bdNP1[2][2] = 0
    PeriLab.Solver_Manager.Model_Factory.Damage.damage_index(nodes)
    @test damageNP1_test[1] == 1
    @test damageNP1_test[2] == 0.75
    @test damageNP1_test[3] == 0
    bdNP1[2][1] = 0
    bdNP1[2][2] = 0
    PeriLab.Solver_Manager.Model_Factory.Damage.damage_index(nodes)
    @test damageNP1_test[1] == 1
    @test damageNP1_test[2] == 1
    @test damageNP1_test[3] == 0
    bdNP1[3][:] .= 0
    PeriLab.Solver_Manager.Model_Factory.Damage.damage_index(nodes)
    @test damageNP1_test[1] == 1
    @test damageNP1_test[2] == 1
    @test damageNP1_test[3] == 1
end
