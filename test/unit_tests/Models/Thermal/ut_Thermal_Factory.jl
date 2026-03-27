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
    nn[1] = 2
    nn[2] = 3
    nn[3] = 1
    nn[4] = 2

    PeriLab.Solver_Manager.Model_Factory.Thermal.init_fields()
    field_keys = PeriLab.Data_Manager.get_all_field_keys()
    @test "TemperatureN" in field_keys
    @test "TemperatureNP1" in field_keys
    @test "Heat FlowN" in field_keys
    @test "Heat FlowNP1" in field_keys
end
@testset "init_model" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_dof(2)
    PeriLab.Data_Manager.set_num_controller(4)
    nn = PeriLab.Data_Manager.create_constant_node_scalar_field("Number of Neighbors",
                                                                Int64)
    nn[1] = 2
    nn[2] = 3
    nn[3] = 1
    nn[4] = 2
    nlist = PeriLab.Data_Manager.create_constant_bond_scalar_state("Neighborhoodlist",
                                                                   Int64)
    nlist[1] = [2, 3]
    nlist[2] = [1, 3, 4]
    nlist[3] = [2]
    nlist[4] = [2, 3]
    bond_geometry = PeriLab.Data_Manager.create_constant_bond_vector_state("Bond Geometry",
                                                                           Float64,
                                                                           2;
                                                                           default_value = 1)
    PeriLab.Data_Manager.set_block_id_list([1, 2])
    PeriLab.Data_Manager.init_properties()
    PeriLab.Data_Manager.set_properties(1,
                                        "Thermal Model",
                                        Dict("Thermal Model" => "Heat Transfer"))
    PeriLab.Solver_Manager.Model_Factory.Thermal.init_model([1], 1)
    PeriLab.Data_Manager.set_properties(2, "Thermal Model",
                                        Dict("Thermal Model" => "Missing"))
    @test_logs (:error, "No thermal model of name Missing exists.") PeriLab.Solver_Manager.Model_Factory.Thermal.init_model([1],
                                                                                                                            2)
end
