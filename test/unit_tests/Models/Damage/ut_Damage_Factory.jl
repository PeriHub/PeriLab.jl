# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test

@testset "init_fields" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.initialize_data()
    test_data_manager.set_dof(3)
    test_data_manager.set_num_controller(4)
    nn = test_data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1] = 1
    nn[2] = 2
    nn[3] = 1
    nn[4] = 2
    nlist = test_data_manager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
    nlist[1] = [2]
    nlist[2] = [1, 3]
    nlist[3] = [1]
    nlist[4] = [1, 3]
    PeriLab.Solver_control.Model_Factory.Damage.init_fields(test_data_manager)
    field_keys = test_data_manager.get_all_field_keys()
    @test "DamageN" in field_keys
    @test "DamageNP1" in field_keys
end

@testset "damage_index" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.set_num_controller(3)
    nn = test_data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1] = 1
    nn[2] = 2
    nn[3] = 1
    damageN, damageNP1_test = test_data_manager.create_node_field("Damage", Float64, 1)
    volume = test_data_manager.create_constant_node_field("Volume", Float64, 1)
    nlist = test_data_manager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
    bdN, bdNP1 = test_data_manager.create_bond_field("Bond Damage", Float64, 1, 1)
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
    PeriLab.Solver_control.Model_Factory.Damage.damage_index(test_data_manager, nodes)
    @test damageNP1_test[1] == 0
    @test damageNP1_test[2] == 0
    @test damageNP1_test[3] == 0
    bdNP1[1][:] .= 0
    PeriLab.Solver_control.Model_Factory.Damage.damage_index(test_data_manager, nodes)
    @test damageNP1_test[1] == 1
    @test damageNP1_test[2] == 0
    @test damageNP1_test[3] == 0
    bdNP1[2][1] = 0
    PeriLab.Solver_control.Model_Factory.Damage.damage_index(test_data_manager, nodes)
    @test damageNP1_test[1] == 1
    @test damageNP1_test[2] == 0.25
    @test damageNP1_test[3] == 0
    bdNP1[2][1] = 1
    bdNP1[2][2] = 0
    PeriLab.Solver_control.Model_Factory.Damage.damage_index(test_data_manager, nodes)
    @test damageNP1_test[1] == 1
    @test damageNP1_test[2] == 0.75
    @test damageNP1_test[3] == 0
    bdNP1[2][1] = 0
    bdNP1[2][2] = 0
    PeriLab.Solver_control.Model_Factory.Damage.damage_index(test_data_manager, nodes)
    @test damageNP1_test[1] == 1
    @test damageNP1_test[2] == 1
    @test damageNP1_test[3] == 0
    bdNP1[3][:] .= 0
    PeriLab.Solver_control.Model_Factory.Damage.damage_index(test_data_manager, nodes)
    @test damageNP1_test[1] == 1
    @test damageNP1_test[2] == 1
    @test damageNP1_test[3] == 1
end
@testset "ut_Damage_factory_exceptions" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.data["properties"][1] = Dict("Damage Model" => Dict("Damage Model" => "not there"))
    @test isnothing(PeriLab.Solver_control.Model_Factory.Damage.init_model(test_data_manager, Vector{Int64}(1:3), 1))
end

@testset "ut_init_interface_crit_values" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.set_block_id_list([2, 3, 1])
    crit_values_matrix::Array{Float64,3} = fill(-1, (1, 1, 1))
    test_data_manager.set_crit_values_matrix(crit_values_matrix)
    damage_parameter = Dict("Critical Value" => 1.0,
                            "Interblock Damage" => Dict("Interblock Critical Value 1_2" => 0.2,
                                                        "Interblock Critical Value 2_3" => 0.3,
                                                        "Interblock Critical Value 2_1" => 0.4))
    PeriLab.Solver_control.Model_Factory.Damage.init_interface_crit_values(test_data_manager, damage_parameter, 1)
    @test test_data_manager.get_crit_values_matrix()[1, 2, 1] == 0.2
    @test test_data_manager.get_crit_values_matrix()[2, 3, 1] == 0.3
    @test test_data_manager.get_crit_values_matrix()[2, 1, 1] == 0.4
end
