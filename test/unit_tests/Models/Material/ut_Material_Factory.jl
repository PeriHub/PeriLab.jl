# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../../../src/Models/Material/Material_Factory.jl")
# include("../../../../src/Core/data_manager.jl")
# include("../../../../src/Support/Parameters/parameter_handling.jl")
using Test
import .Material

@testset "ut_init_model_fields" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.initialize_data()
    test_data_manager.set_dof(3)
    test_data_manager.set_num_controller(4)
    test_data_manager.create_constant_node_field("Coordinates", Float64, 3)
    nn = test_data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1:4] = 1:4
    Material.init_fields(test_data_manager)
    field_keys = test_data_manager.get_all_field_keys()
    @test "ForcesN" in field_keys
    @test "ForcesNP1" in field_keys
    @test "Acceleration" in field_keys
    @test "VelocityN" in field_keys
    @test "VelocityNP1" in field_keys
    @test "Bond Forces" in field_keys
end
