# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../../src/Physics/Physics_Factory.jl")
include("../../../src/Support/data_manager.jl")
include("../../../src/Support/Parameters/parameter_handling.jl")
using Test
import .Physics
@testset "ut_get_block_model_definition" begin
    test_Data_manager = Data_manager
    block_list = [1, 2, 3]
    test_Data_manager.set_block_list(block_list)
    prop_keys = test_Data_manager.init_property()
    params = Dict("Blocks" => Dict("block_1" => Dict("Material Model" => "a"), "block_2" => Dict("Material Model" => "c"), "block_3" => Dict("Material Model" => "a", "Damage Model" => "a", "Thermal Model" => "therm", "Additive Model" => "add")), "Physics" => Dict("Material Models" => Dict("a" => Dict("value" => 1), "c" => Dict("value" => [1 2], "value2" => 1)), "Damage Models" => Dict("a" => Dict("value" => 3), "c" => Dict("value" => [1 2], "value2" => 1)), "Thermal Models" => Dict("therm" => Dict("value" => "hot", "bool" => true)), "Additive Models" => Dict("add" => Dict("value" => "ad", "bool" => false))))

    for block in block_list
        Physics.get_block_model_definition(params, block, prop_keys, test_Data_manager.set_properties)
    end
    @test test_Data_manager.get_property(1, "Material Model", "value") == params["Physics"]["Material Models"]["a"]["value"]
    @test test_Data_manager.get_property(2, "Material Model", "value") == params["Physics"]["Material Models"]["c"]["value"]
    @test test_Data_manager.get_property(2, "Material Model", "value2") == params["Physics"]["Material Models"]["c"]["value2"]
    @test test_Data_manager.get_property(3, "Material Model", "value") == params["Physics"]["Material Models"]["a"]["value"]
    @test test_Data_manager.get_property(3, "Damage Model", "value") == params["Physics"]["Damage Models"]["a"]["value"]
    @test test_Data_manager.get_property(3, "Thermal Model", "value") == params["Physics"]["Thermal Models"]["therm"]["value"]
    @test test_Data_manager.get_property(3, "Thermal Model", "bool") == params["Physics"]["Thermal Models"]["therm"]["bool"]
    @test test_Data_manager.get_property(3, "Additive Model", "value") == params["Physics"]["Additive Models"]["add"]["value"]
    @test test_Data_manager.get_property(3, "Additive Model", "bool") == params["Physics"]["Additive Models"]["add"]["bool"]
end

@testset "ut_read_properties" begin
    test_Data_manager_read_properties = Data_manager
    block_list = [1, 2, 3]
    test_Data_manager_read_properties.set_block_list(block_list)

    params = Dict("Blocks" => Dict("block_1" => Dict("Material Model" => "a"), "block_2" => Dict("Material Model" => "c"), "block_3" => Dict("Material Model" => "a", "Damage Model" => "a", "Thermal Model" => "therm")), "Physics" => Dict("Material Models" => Dict("a" => Dict("value" => 1, "name" => "t4"), "c" => Dict("value" => [1 2], "value2" => 1, "name" => "t4")), "Damage Models" => Dict("a" => Dict("value" => 3, "name" => "t"), "c" => Dict("value" => [1 2], "value2" => 1, "name" => "t2")), "Thermal Models" => Dict("therm" => Dict("value" => "hot", "bool" => true, "name" => "t3"))))
    Physics.read_properties(params, test_Data_manager_read_properties, false)


    @test test_Data_manager_read_properties.get_property(1, "Material Model", "value") == params["Physics"]["Material Models"]["a"]["value"]
    @test test_Data_manager_read_properties.get_property(2, "Material Model", "value") == params["Physics"]["Material Models"]["c"]["value"]
    @test test_Data_manager_read_properties.get_property(2, "Material Model", "value2") == params["Physics"]["Material Models"]["c"]["value2"]
    @test test_Data_manager_read_properties.get_property(3, "Material Model", "value") == params["Physics"]["Material Models"]["a"]["value"]
    @test test_Data_manager_read_properties.get_property(3, "Damage Model", "value") == params["Physics"]["Damage Models"]["a"]["value"]
    @test test_Data_manager_read_properties.get_property(3, "Thermal Model", "value") == params["Physics"]["Thermal Models"]["therm"]["value"]
    @test test_Data_manager_read_properties.get_property(3, "Thermal Model", "bool") == params["Physics"]["Thermal Models"]["therm"]["bool"]

    Physics.read_properties(params, test_Data_manager_read_properties, true)
    @test test_Data_manager_read_properties.get_property(1, "Material Model", "value") == params["Physics"]["Material Models"]["a"]["value"]
    @test test_Data_manager_read_properties.get_property(2, "Material Model", "value") == params["Physics"]["Material Models"]["c"]["value"]
    @test test_Data_manager_read_properties.get_property(2, "Material Model", "value2") == params["Physics"]["Material Models"]["c"]["value2"]
    @test test_Data_manager_read_properties.get_property(3, "Material Model", "value") == params["Physics"]["Material Models"]["a"]["value"]
    @test test_Data_manager_read_properties.get_property(3, "Damage Model", "value") == params["Physics"]["Damage Models"]["a"]["value"]
    @test test_Data_manager_read_properties.get_property(3, "Thermal Model", "value") == params["Physics"]["Thermal Models"]["therm"]["value"]
    @test test_Data_manager_read_properties.get_property(3, "Thermal Model", "bool") == params["Physics"]["Thermal Models"]["therm"]["bool"]
end

@testset "ut_init_material_model_fields" begin
    test_Data_manager = Data_manager
    test_Data_manager.set_dof(3)
    test_Data_manager.set_nmasters(4)
    test_Data_manager.create_constant_node_field("Coordinates", Float64, 3)
    nn = test_Data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1:4] = 1:4
    Physics.init_material_model_fields(test_Data_manager)
    fieldkeys = test_Data_manager.get_all_field_keys()
    @test "ForcesN" in fieldkeys
    @test "ForcesNP1" in fieldkeys
    @test "Acceleration" in fieldkeys
    @test "VelocityN" in fieldkeys
    @test "VelocityNP1" in fieldkeys
    @test "Bond Forces" in fieldkeys
end

@testset "init_damage_model_fields" begin
    test_Data_manager = Data_manager
    test_Data_manager.set_dof(3)
    test_Data_manager.set_nmasters(4)
    Physics.init_damage_model_fields(test_Data_manager)
    fieldkeys = test_Data_manager.get_all_field_keys()
    @test "DamageN" in fieldkeys
    @test "DamageNP1" in fieldkeys
end

@testset "init_thermal_model_fields" begin
    test_Data_manager = Data_manager
    test_Data_manager.set_dof(3)
    test_Data_manager.set_nmasters(4)
    nn = test_Data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1] = 2
    nn[2] = 3
    nn[3] = 1
    nn[4] = 2

    Physics.init_thermal_model_fields(test_Data_manager)
    fieldkeys = test_Data_manager.get_all_field_keys()
    @test "TemperatureN" in fieldkeys
    @test "TemperatureNP1" in fieldkeys
    @test "Heat FlowN" in fieldkeys
    @test "Heat FlowNP1" in fieldkeys
    @test "Bond Heat Flow" in fieldkeys

end

@testset "init_additive_model_fields" begin
    test_Data_manager = Data_manager
    test_Data_manager.set_dof(3)
    test_Data_manager.set_nmasters(4)
    nn = test_Data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1] = 1
    nn[2] = 2
    nn[3] = 1
    nn[4] = 2

    test_Data_manager.create_bond_field("Bond Damage", Float64, 1)
    test_Data_manager = Physics.init_additive_model_fields(test_Data_manager)
    fieldkeys = test_Data_manager.get_all_field_keys()
    @test "Active" in fieldkeys
    active = test_Data_manager.get_field("Active")
    for iID in 1:4
        @test active[iID] == false
    end
end
@testset "init_pre_calculation" begin
    test_Data_manager = Data_manager
    test_Data_manager.set_dof(3)
    test_Data_manager.set_nmasters(4)
    n = test_Data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    n .= 2
    options = Dict("Deformed Bond Geometry" => true, "Shape Tensor" => false, "Deformation Gradient" => false, "Bond Associated Shape Tensor" => false, "Bond Associated Deformation Gradient" => false)
    test_Data_manager = Physics.init_pre_calculation(test_Data_manager, options)

    @test "Deformed Bond GeometryN" in test_Data_manager.get_all_field_keys()
    @test "Deformed Bond GeometryNP1" in test_Data_manager.get_all_field_keys()
    @test !("Shape Tensor" in test_Data_manager.get_all_field_keys())
    @test !("Deformation Gradient" in test_Data_manager.get_all_field_keys())
    @test !("Bond Associated Shape Tensor" in test_Data_manager.get_all_field_keys())
    @test !("Bond Associated Deformation Gradient" in test_Data_manager.get_all_field_keys())

    options = Dict("Deformed Bond Geometry" => true, "Shape Tensor" => false, "Deformation Gradient" => true, "Bond Associated Shape Tensor" => false, "Bond Associated Deformation Gradient" => false)

    test_Data_manager = Physics.init_pre_calculation(test_Data_manager, options)
    @test "Deformed Bond GeometryN" in test_Data_manager.get_all_field_keys()
    @test "Deformed Bond GeometryNP1" in test_Data_manager.get_all_field_keys()
    @test !("Shape Tensor" in test_Data_manager.get_all_field_keys())
    @test !("Inverse Shape Tensor" in test_Data_manager.get_all_field_keys())
    @test "Deformation Gradient" in test_Data_manager.get_all_field_keys()
    @test !("Bond Associated Shape Tensor" in test_Data_manager.get_all_field_keys())
    @test !("Bond Associated Deformation Gradient" in test_Data_manager.get_all_field_keys())

    options = Dict("Deformed Bond Geometry" => true, "Shape Tensor" => true, "Deformation Gradient" => true, "Bond Associated Shape Tensor" => true, "Bond Associated Deformation Gradient" => false)

    test_Data_manager = Physics.init_pre_calculation(test_Data_manager, options)
    @test "Deformed Bond GeometryN" in test_Data_manager.get_all_field_keys()
    @test "Deformed Bond GeometryNP1" in test_Data_manager.get_all_field_keys()
    @test "Shape Tensor" in test_Data_manager.get_all_field_keys()
    @test "Inverse Shape Tensor" in test_Data_manager.get_all_field_keys()
    @test "Deformation Gradient" in test_Data_manager.get_all_field_keys()
    @test "Bond Associated Shape Tensor" in test_Data_manager.get_all_field_keys()
    @test !("Bond Associated Deformation Gradient" in test_Data_manager.get_all_field_keys())

    options = Dict("Deformed Bond Geometry" => true, "Shape Tensor" => true, "Deformation Gradient" => true, "Bond Associated Shape Tensor" => true, "Bond Associated Deformation Gradient" => true)

    test_Data_manager = Physics.init_pre_calculation(test_Data_manager, options)
    @test "Deformed Bond GeometryN" in test_Data_manager.get_all_field_keys()
    @test "Deformed Bond GeometryNP1" in test_Data_manager.get_all_field_keys()
    @test "Shape Tensor" in test_Data_manager.get_all_field_keys()
    @test "Deformation Gradient" in test_Data_manager.get_all_field_keys()
    @test "Bond Associated Shape Tensor" in test_Data_manager.get_all_field_keys()
    @test "Bond Associated Deformation Gradient" in test_Data_manager.get_all_field_keys()
end