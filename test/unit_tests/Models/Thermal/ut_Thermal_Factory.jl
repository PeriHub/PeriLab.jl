# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
using Test

@testset "init_fields" begin
    test_data_manager = PeriLab.Data_Manager
    test_data_manager.initialize_data()
    test_data_manager.set_dof(3)
    test_data_manager.set_num_controller(4)
    nn = test_data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1] = 2
    nn[2] = 3
    nn[3] = 1
    nn[4] = 2

    PeriLab.Solver_Manager.Model_Factory.Thermal.init_fields(test_data_manager)
    field_keys = test_data_manager.get_all_field_keys()
    @test "TemperatureN" in field_keys
    @test "TemperatureNP1" in field_keys
    @test "Heat FlowN" in field_keys
    @test "Heat FlowNP1" in field_keys
end
@testset "init_model" begin
    test_data_manager = PeriLab.Data_Manager
    test_data_manager.initialize_data()
    test_data_manager.set_dof(2)
    test_data_manager.set_num_controller(4)
    nn = test_data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1] = 2
    nn[2] = 3
    nn[3] = 1
    nn[4] = 2
    nlist = test_data_manager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
    nlist[1] = [2, 3]
    nlist[2] = [1, 3, 4]
    nlist[3] = [2]
    nlist[4] = [2, 3]
    bond_geometry = test_data_manager.create_constant_bond_field("Bond Geometry", Float64,
                                                                 2, 1)
    test_data_manager.set_block_id_list([1, 2])
    test_data_manager.init_properties()
    test_data_manager.set_properties(1,
                                     "Thermal Model",
                                     Dict("Thermal Model" => "Heat Transfer"))
    PeriLab.Solver_Manager.Model_Factory.Thermal.init_model(test_data_manager, [1], 1)
    test_data_manager.set_properties(2, "Thermal Model", Dict("Thermal Model" => "Missing"))
    @test isnothing(PeriLab.Solver_Manager.Model_Factory.Thermal.init_model(test_data_manager, [1], 2))
end
