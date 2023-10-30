# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../../src/Physics/Physics_Factory.jl")
include("../../../src/Support/data_manager.jl")
include("../../../src/Support/Parameters/parameter_handling.jl")
using Test
import .Physics
@testset "ut_get_block_model_definition" begin
    testDatamanager = Data_manager
    block_list = [1, 2, 3]
    testDatamanager.set_block_list(block_list)
    prop_keys = testDatamanager.init_property()
    params = Dict("Blocks" => Dict("block_1" => Dict("Material Model" => "a"), "block_2" => Dict("Material Model" => "c"), "block_3" => Dict("Material Model" => "a", "Damage Model" => "a", "Thermal Model" => "therm")), "Physics" => Dict("Material Models" => Dict("a" => Dict("value" => 1), "c" => Dict("value" => [1 2], "value2" => 1)), "Damage Models" => Dict("a" => Dict("value" => 3), "c" => Dict("value" => [1 2], "value2" => 1)), "Thermal Models" => Dict("therm" => Dict("value" => "hot", "bool" => true))))

    for block in block_list
        Physics.get_block_model_definition(params, block, prop_keys, testDatamanager.set_properties)
    end
    @test testDatamanager.get_property(1, "Material Model", "value") == params["Physics"]["Material Models"]["a"]["value"]
    @test testDatamanager.get_property(2, "Material Model", "value") == params["Physics"]["Material Models"]["c"]["value"]
    @test testDatamanager.get_property(2, "Material Model", "value2") == params["Physics"]["Material Models"]["c"]["value2"]
    @test testDatamanager.get_property(3, "Material Model", "value") == params["Physics"]["Material Models"]["a"]["value"]
    @test testDatamanager.get_property(3, "Damage Model", "value") == params["Physics"]["Damage Models"]["a"]["value"]
    @test testDatamanager.get_property(3, "Thermal Model", "value") == params["Physics"]["Thermal Models"]["therm"]["value"]
    @test testDatamanager.get_property(3, "Thermal Model", "bool") == params["Physics"]["Thermal Models"]["therm"]["bool"]
end

@testset "ut_read_properties" begin
    testDatamanager_read_properties = Data_manager
    block_list = [1, 2, 3]
    testDatamanager_read_properties.set_block_list(block_list)

    params = Dict("Blocks" => Dict("block_1" => Dict("Material Model" => "a"), "block_2" => Dict("Material Model" => "c"), "block_3" => Dict("Material Model" => "a", "Damage Model" => "a", "Thermal Model" => "therm")), "Physics" => Dict("Material Models" => Dict("a" => Dict("value" => 1, "name" => "t4"), "c" => Dict("value" => [1 2], "value2" => 1, "name" => "t4")), "Damage Models" => Dict("a" => Dict("value" => 3, "name" => "t"), "c" => Dict("value" => [1 2], "value2" => 1, "name" => "t2")), "Thermal Models" => Dict("therm" => Dict("value" => "hot", "bool" => true, "name" => "t3"))))

    Physics.read_properties(params, testDatamanager_read_properties)

    @test testDatamanager_read_properties.get_property(1, "Material Model", "value") == params["Physics"]["Material Models"]["a"]["value"]
    @test testDatamanager_read_properties.get_property(2, "Material Model", "value") == params["Physics"]["Material Models"]["c"]["value"]
    @test testDatamanager_read_properties.get_property(2, "Material Model", "value2") == params["Physics"]["Material Models"]["c"]["value2"]
    @test testDatamanager_read_properties.get_property(3, "Material Model", "value") == params["Physics"]["Material Models"]["a"]["value"]
    @test testDatamanager_read_properties.get_property(3, "Damage Model", "value") == params["Physics"]["Damage Models"]["a"]["value"]
    @test testDatamanager_read_properties.get_property(3, "Thermal Model", "value") == params["Physics"]["Thermal Models"]["therm"]["value"]
    @test testDatamanager_read_properties.get_property(3, "Thermal Model", "bool") == params["Physics"]["Thermal Models"]["therm"]["bool"]
end

@testset "ut_init_material_model_fields" begin
    testDatamanager = Data_manager
    testDatamanager.set_dof(3)
    testDatamanager.set_nmasters(4)
    testDatamanager.create_constant_node_field("Coordinates", Float64, 3)
    Physics.init_material_model_fields(testDatamanager)
    fieldkeys = testDatamanager.get_all_field_keys()
    @test "ForcesN" in fieldkeys
    @test "ForcesNP1" in fieldkeys
    @test "Acceleration" in fieldkeys
    @test "VelocityN" in fieldkeys
    @test "VelocityNP1" in fieldkeys
end

@testset "init_damage_model_fields" begin
    testDatamanager = Data_manager
    testDatamanager.set_dof(3)
    testDatamanager.set_nmasters(4)
    Physics.init_damage_model_fields(testDatamanager)
    fieldkeys = testDatamanager.get_all_field_keys()
    @test "DamageN" in fieldkeys
    @test "DamageNP1" in fieldkeys
end

@testset "init_thermal_model_fields" begin
    testDatamanager = Data_manager
    testDatamanager.set_dof(3)
    testDatamanager.set_nmasters(4)
    nn = testDatamanager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1] = 2
    nn[2] = 3
    nn[3] = 1
    nn[4] = 2

    Physics.init_thermal_model_fields(testDatamanager)
    fieldkeys = testDatamanager.get_all_field_keys()
    @test "TemperatureN" in fieldkeys
    @test "TemperatureNP1" in fieldkeys
    @test "Heat FlowN" in fieldkeys
    @test "Heat FlowNP1" in fieldkeys
    @test "Bond Heat Flow" in fieldkeys

end

@testset "init_additive_model_fields" begin
    testDatamanager = Data_manager
    testDatamanager.set_dof(3)
    testDatamanager.set_nmasters(4)
    active = testDatamanager.create_constant_node_field("Active", Bool, 1)
    nn = testDatamanager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1] = 1
    nn[2] = 2
    nn[3] = 1
    nn[4] = 2

    active[:] .= true
    testDatamanager.create_bond_field("Bond Damage", Float64, 1)
    testDatamanager = Physics.init_additive_model_fields(testDatamanager)
    fieldkeys = testDatamanager.get_all_field_keys()
    @test "Active" in fieldkeys
    active = testDatamanager.get_field("Active")
    for iID in 1:4
        @test active[iID] == false
    end
end
@testset "init_pre_calculation" begin
    testDatamanager = Data_manager
    testDatamanager.set_dof(3)
    testDatamanager.set_nmasters(4)
    n = testDatamanager.create_constant_node_field("Number of Neighbors", Int64, 1)
    n .= 2
    options = Dict("Deformed Bond Geometry" => true, "Shape Tensor" => false, "Deformation Gradient" => false, "Bond Associated Shape Tensor" => false, "Bond Associated Deformation Gradient" => false)
    testDatamanager = Physics.init_pre_calculation(testDatamanager, options)

    @test "Deformed Bond GeometryN" in testDatamanager.get_all_field_keys()
    @test "Deformed Bond GeometryNP1" in testDatamanager.get_all_field_keys()
    @test !("Shape Tensor" in testDatamanager.get_all_field_keys())
    @test !("Deformation Gradient" in testDatamanager.get_all_field_keys())
    @test !("Bond Associated Shape Tensor" in testDatamanager.get_all_field_keys())
    @test !("Bond Associated Deformation Gradient" in testDatamanager.get_all_field_keys())

    options = Dict("Deformed Bond Geometry" => true, "Shape Tensor" => false, "Deformation Gradient" => true, "Bond Associated Shape Tensor" => false, "Bond Associated Deformation Gradient" => false)

    testDatamanager = Physics.init_pre_calculation(testDatamanager, options)
    @test "Deformed Bond GeometryN" in testDatamanager.get_all_field_keys()
    @test "Deformed Bond GeometryNP1" in testDatamanager.get_all_field_keys()
    @test !("Shape Tensor" in testDatamanager.get_all_field_keys())
    @test !("Inverse Shape Tensor" in testDatamanager.get_all_field_keys())
    @test "Deformation Gradient" in testDatamanager.get_all_field_keys()
    @test !("Bond Associated Shape Tensor" in testDatamanager.get_all_field_keys())
    @test !("Bond Associated Deformation Gradient" in testDatamanager.get_all_field_keys())

    options = Dict("Deformed Bond Geometry" => true, "Shape Tensor" => true, "Deformation Gradient" => true, "Bond Associated Shape Tensor" => true, "Bond Associated Deformation Gradient" => false)

    testDatamanager = Physics.init_pre_calculation(testDatamanager, options)
    @test "Deformed Bond GeometryN" in testDatamanager.get_all_field_keys()
    @test "Deformed Bond GeometryNP1" in testDatamanager.get_all_field_keys()
    @test "Shape Tensor" in testDatamanager.get_all_field_keys()
    @test "Inverse Shape Tensor" in testDatamanager.get_all_field_keys()
    @test "Deformation Gradient" in testDatamanager.get_all_field_keys()
    @test "Bond Associated Shape Tensor" in testDatamanager.get_all_field_keys()
    @test !("Bond Associated Deformation Gradient" in testDatamanager.get_all_field_keys())

    options = Dict("Deformed Bond Geometry" => true, "Shape Tensor" => true, "Deformation Gradient" => true, "Bond Associated Shape Tensor" => true, "Bond Associated Deformation Gradient" => true)

    testDatamanager = Physics.init_pre_calculation(testDatamanager, options)
    @test "Deformed Bond GeometryN" in testDatamanager.get_all_field_keys()
    @test "Deformed Bond GeometryNP1" in testDatamanager.get_all_field_keys()
    @test "Shape Tensor" in testDatamanager.get_all_field_keys()
    @test "Deformation Gradient" in testDatamanager.get_all_field_keys()
    @test "Bond Associated Shape Tensor" in testDatamanager.get_all_field_keys()
    @test "Bond Associated Deformation Gradient" in testDatamanager.get_all_field_keys()
end