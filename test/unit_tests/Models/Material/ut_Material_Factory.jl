# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

#using Test

@testset "ut_init_model_fields" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_dof(3)
    PeriLab.Data_Manager.set_num_controller(4)
    PeriLab.Data_Manager.create_constant_node_vector_field("Coordinates", Float64, 3)
    nn = PeriLab.Data_Manager.create_constant_node_scalar_field("Number of Neighbors",
                                                                Int64)
    nn[1:4] = 1:4
    PeriLab.Solver_Manager.Model_Factory.Material.init_fields()
    field_keys = PeriLab.Data_Manager.get_all_field_keys()
    @test "ForcesN" in field_keys
    @test "ForcesNP1" in field_keys
    @test "Acceleration" in field_keys
    @test "VelocityN" in field_keys
    @test "VelocityNP1" in field_keys
    @test "Bond Forces" in field_keys
end

@testset "ut_init_model" begin
    PeriLab.Data_Manager.set_block_id_list([2, 3, 1])
    PeriLab.Data_Manager.init_properties()
    PeriLab.Data_Manager.set_property(1, "Material Model", "E", 1)
    @test_logs (:error, "Block 1 has no material model defined.") PeriLab.Solver_Manager.Model_Factory.Material.init_model(Vector{Int64}(1:4),
                                                                                                                           1)
    @test_logs (:error, "Block 2 has no material model defined.") PeriLab.Solver_Manager.Model_Factory.Material.init_model(Vector{Int64}(1:4),
                                                                                                                           2)
end
