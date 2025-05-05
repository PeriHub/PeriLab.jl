# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
using Test

include("../../../../src/Models/Thermal/Thermal_Factory.jl")
using .Thermal
# include("../../../../src/Core/Data_manager.jl")

@testset "init_fields" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.initialize_data()
    test_data_manager.set_dof(3)
    test_data_manager.set_num_controller(4)
    nn = test_data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1] = 2
    nn[2] = 3
    nn[3] = 1
    nn[4] = 2

    Thermal.init_fields(test_data_manager)
    field_keys = test_data_manager.get_all_field_keys()
    @test "TemperatureN" in field_keys
    @test "TemperatureNP1" in field_keys
    @test "Heat FlowN" in field_keys
    @test "Heat FlowNP1" in field_keys
end
@testset "init_model" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.initialize_data()
    test_data_manager.set_block_id_list([1, 2])
    test_data_manager.init_properties()
    test_data_manager.set_properties(1,
                                     "Thermal Model",
                                     Dict("Thermal Model" => "Heat Transfer"))
    Thermal.init_model(test_data_manager, [1], 1)
    test_data_manager.set_properties(2, "Thermal Model", Dict("Thermal Model" => "Missing"))
    @test isnothing(Thermal.init_model(test_data_manager, [1], 2))
end
