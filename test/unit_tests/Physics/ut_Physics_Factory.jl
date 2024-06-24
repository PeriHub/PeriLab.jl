# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../../src/Physics/Physics_Factory.jl")
# include("../../../src/PeriLab.jl")
# using .PeriLab
using Test
import .Physics
@testset "ut_get_block_model_definition" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.clear_data_manager()
    block_list = [1, 2, 3]
    test_data_manager.set_block_list(block_list)
    prop_keys = test_data_manager.init_property()
    params = Dict("Blocks" => Dict("block_1" => Dict("Material Model" => "a"), "block_2" => Dict("Material Model" => "c"), "block_3" => Dict("Material Model" => "a", "Damage Model" => "a", "Thermal Model" => "therm", "Additive Model" => "add")), "Physics" => Dict("Material Models" => Dict("a" => Dict("value" => 1), "c" => Dict("value" => [1 2], "value2" => 1)), "Damage Models" => Dict("a" => Dict("value" => 3), "c" => Dict("value" => [1 2], "value2" => 1)), "Thermal Models" => Dict("therm" => Dict("value" => "hot", "bool" => true)), "Additive Models" => Dict("add" => Dict("value" => "ad", "bool" => false))))

    for block in block_list
        Physics.get_block_model_definition(params, block, prop_keys, test_data_manager.set_properties)
    end
    @test test_data_manager.get_property(1, "Material Model", "value") == params["Physics"]["Material Models"]["a"]["value"]
    @test test_data_manager.get_property(2, "Material Model", "value") == params["Physics"]["Material Models"]["c"]["value"]
    @test test_data_manager.get_property(2, "Material Model", "value2") == params["Physics"]["Material Models"]["c"]["value2"]
    @test test_data_manager.get_property(3, "Material Model", "value") == params["Physics"]["Material Models"]["a"]["value"]
    @test test_data_manager.get_property(3, "Damage Model", "value") == params["Physics"]["Damage Models"]["a"]["value"]
    @test test_data_manager.get_property(3, "Thermal Model", "value") == params["Physics"]["Thermal Models"]["therm"]["value"]
    @test test_data_manager.get_property(3, "Thermal Model", "bool") == params["Physics"]["Thermal Models"]["therm"]["bool"]
    @test test_data_manager.get_property(3, "Additive Model", "value") == params["Physics"]["Additive Models"]["add"]["value"]
    @test test_data_manager.get_property(3, "Additive Model", "bool") == params["Physics"]["Additive Models"]["add"]["bool"]
end

@testset "ut_read_properties" begin
    test_data_manager_read_properties = PeriLab.Data_manager
    block_list = [1, 2, 3]
    test_data_manager_read_properties.set_block_list(block_list)

    params = Dict("Blocks" => Dict("block_1" => Dict("Material Model" => "a"), "block_2" => Dict("Material Model" => "c"), "block_3" => Dict("Material Model" => "a", "Damage Model" => "a", "Thermal Model" => "therm")), "Physics" => Dict("Material Models" => Dict("a" => Dict("value" => 1, "name" => "t4", "Material Model" => "Test"), "c" => Dict("value" => [1 2], "value2" => 1, "name" => "t4", "Material Model" => "Test")), "Damage Models" => Dict("a" => Dict("value" => 3, "name" => "t"), "c" => Dict("value" => [1 2], "value2" => 1, "name" => "t2")), "Thermal Models" => Dict("therm" => Dict("value" => "hot", "bool" => true, "name" => "t3"))))
    Physics.read_properties(params, test_data_manager_read_properties, false)


    @test test_data_manager_read_properties.get_property(1, "Material Model", "value") == params["Physics"]["Material Models"]["a"]["value"]
    @test test_data_manager_read_properties.get_property(2, "Material Model", "value") == params["Physics"]["Material Models"]["c"]["value"]
    @test test_data_manager_read_properties.get_property(2, "Material Model", "value2") == params["Physics"]["Material Models"]["c"]["value2"]
    @test test_data_manager_read_properties.get_property(3, "Material Model", "value") == params["Physics"]["Material Models"]["a"]["value"]
    @test test_data_manager_read_properties.get_property(3, "Damage Model", "value") == params["Physics"]["Damage Models"]["a"]["value"]
    @test test_data_manager_read_properties.get_property(3, "Thermal Model", "value") == params["Physics"]["Thermal Models"]["therm"]["value"]
    @test test_data_manager_read_properties.get_property(3, "Thermal Model", "bool") == params["Physics"]["Thermal Models"]["therm"]["bool"]

    Physics.read_properties(params, test_data_manager_read_properties, true)
    @test test_data_manager_read_properties.get_property(1, "Material Model", "value") == params["Physics"]["Material Models"]["a"]["value"]
    @test test_data_manager_read_properties.get_property(2, "Material Model", "value") == params["Physics"]["Material Models"]["c"]["value"]
    @test test_data_manager_read_properties.get_property(2, "Material Model", "value2") == params["Physics"]["Material Models"]["c"]["value2"]
    @test test_data_manager_read_properties.get_property(3, "Material Model", "value") == params["Physics"]["Material Models"]["a"]["value"]
    @test test_data_manager_read_properties.get_property(3, "Damage Model", "value") == params["Physics"]["Damage Models"]["a"]["value"]
    @test test_data_manager_read_properties.get_property(3, "Thermal Model", "value") == params["Physics"]["Thermal Models"]["therm"]["value"]
    @test test_data_manager_read_properties.get_property(3, "Thermal Model", "bool") == params["Physics"]["Thermal Models"]["therm"]["bool"]
end

@testset "init_pre_calculation" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.set_dof(3)
    test_data_manager.set_num_controller(4)
    nn = test_data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1] = 2
    nn[2] = 3
    nn[3] = 1
    nn[4] = 2
    horizon = test_data_manager.create_constant_node_field("Horizon", Float64, 1)
    options = Dict("Deformed Bond Geometry" => true, "Shape Tensor" => false, "Deformation Gradient" => false, "Bond Associated Shape Tensor" => false, "Bond Associated Deformation Gradient" => false)
    test_data_manager = Physics.init_pre_calculation(test_data_manager, options)

    @test "Deformed Bond GeometryN" in test_data_manager.get_all_field_keys()
    @test "Deformed Bond GeometryNP1" in test_data_manager.get_all_field_keys()
    @test !("Shape Tensor" in test_data_manager.get_all_field_keys())
    @test !("Deformation Gradient" in test_data_manager.get_all_field_keys())
    @test !("Bond Associated Shape Tensor" in test_data_manager.get_all_field_keys())
    @test !("Bond Associated Deformation Gradient" in test_data_manager.get_all_field_keys())

    options = Dict("Deformed Bond Geometry" => true, "Shape Tensor" => false, "Deformation Gradient" => true, "Bond Associated Shape Tensor" => false, "Bond Associated Deformation Gradient" => false)

    test_data_manager = Physics.init_pre_calculation(test_data_manager, options)
    @test "Deformed Bond GeometryN" in test_data_manager.get_all_field_keys()
    @test "Deformed Bond GeometryNP1" in test_data_manager.get_all_field_keys()
    @test "Shape Tensor" in test_data_manager.get_all_field_keys()
    @test "Inverse Shape Tensor" in test_data_manager.get_all_field_keys()
    @test "Deformation Gradient" in test_data_manager.get_all_field_keys()
    @test !("Bond Associated Shape Tensor" in test_data_manager.get_all_field_keys())
    @test !("Bond Associated Deformation Gradient" in test_data_manager.get_all_field_keys())

    options = Dict("Deformed Bond Geometry" => true, "Shape Tensor" => true, "Deformation Gradient" => true, "Bond Associated Shape Tensor" => true, "Bond Associated Deformation Gradient" => false)

    test_data_manager = Physics.init_pre_calculation(test_data_manager, options)
    @test "Deformed Bond GeometryN" in test_data_manager.get_all_field_keys()
    @test "Deformed Bond GeometryNP1" in test_data_manager.get_all_field_keys()
    @test "Shape Tensor" in test_data_manager.get_all_field_keys()
    @test "Inverse Shape Tensor" in test_data_manager.get_all_field_keys()
    @test "Deformation Gradient" in test_data_manager.get_all_field_keys()
    @test "Bond Associated Shape Tensor" in test_data_manager.get_all_field_keys()
    @test !("Bond Associated Deformation Gradient" in test_data_manager.get_all_field_keys())

    options = Dict("Deformed Bond Geometry" => true, "Shape Tensor" => false, "Deformation Gradient" => true, "Bond Associated Shape Tensor" => false, "Bond Associated Deformation Gradient" => true)

    test_data_manager = Physics.init_pre_calculation(test_data_manager, options)
    @test "Deformed Bond GeometryN" in test_data_manager.get_all_field_keys()
    @test "Deformed Bond GeometryNP1" in test_data_manager.get_all_field_keys()
    @test "Shape Tensor" in test_data_manager.get_all_field_keys()
    @test "Deformation Gradient" in test_data_manager.get_all_field_keys()
    @test "Bond Associated Shape Tensor" in test_data_manager.get_all_field_keys()
    @test "Bond Associated Deformation Gradient" in test_data_manager.get_all_field_keys()
end