# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../../src/Models/Model_Factory.jl")
#include("../../../src/PeriLab.jl")
#using .PeriLab
using Test

@testset "ut_get_block_model_definition" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.initialize_data()
    block_list = [1, 2, 3]
    test_data_manager.set_block_list(block_list)
    prop_keys = test_data_manager.init_properties()
    params = Dict(
        "Blocks" => Dict(
            "block_1" => Dict("Material Model" => "a"),
            "block_2" => Dict("Material Model" => "c"),
            "block_3" => Dict(
                "Material Model" => "a",
                "Damage Model" => "a",
                "Thermal Model" => "therm",
                "Additive Model" => "add",
            ),
        ),
        "Models" => Dict(
            "Material Models" => Dict(
                "a" => Dict("value" => 1),
                "c" => Dict("value" => [1 2], "value2" => 1),
            ),
            "Damage Models" => Dict(
                "a" => Dict("value" => 3),
                "c" => Dict("value" => [1 2], "value2" => 1),
            ),
            "Thermal Models" => Dict("therm" => Dict("value" => "hot", "bool" => true)),
            "Additive Models" => Dict("add" => Dict("value" => "ad", "bool" => false)),
        ),
    )

    Model_Factory.get_block_model_definition(
        params,
        block_list,
        prop_keys,
        test_data_manager.set_properties,
    )

    @test test_data_manager.get_property(1, "Material Model", "value") ==
          params["Models"]["Material Models"]["a"]["value"]
    @test test_data_manager.get_property(2, "Material Model", "value") ==
          params["Models"]["Material Models"]["c"]["value"]
    @test test_data_manager.get_property(2, "Material Model", "value2") ==
          params["Models"]["Material Models"]["c"]["value2"]
    @test test_data_manager.get_property(3, "Material Model", "value") ==
          params["Models"]["Material Models"]["a"]["value"]
    @test test_data_manager.get_property(3, "Damage Model", "value") ==
          params["Models"]["Damage Models"]["a"]["value"]
    @test test_data_manager.get_property(3, "Thermal Model", "value") ==
          params["Models"]["Thermal Models"]["therm"]["value"]
    @test test_data_manager.get_property(3, "Thermal Model", "bool") ==
          params["Models"]["Thermal Models"]["therm"]["bool"]
    @test test_data_manager.get_property(3, "Additive Model", "value") ==
          params["Models"]["Additive Models"]["add"]["value"]
    @test test_data_manager.get_property(3, "Additive Model", "bool") ==
          params["Models"]["Additive Models"]["add"]["bool"]
end

@testset "ut_read_properties" begin
    test_data_manager_read_properties = PeriLab.Data_manager
    block_list = [1, 2, 3]
    test_data_manager_read_properties.set_block_list(block_list)

    params = Dict(
        "Blocks" => Dict(
            "block_1" => Dict("Material Model" => "a"),
            "block_2" => Dict("Material Model" => "c"),
            "block_3" => Dict(
                "Material Model" => "a",
                "Damage Model" => "a",
                "Thermal Model" => "therm",
            ),
        ),
        "Models" => Dict(
            "Material Models" => Dict(
                "a" => Dict("value" => 1, "name" => "t4", "Material Model" => "Test"),
                "c" => Dict(
                    "value" => [1 2],
                    "value2" => 1,
                    "name" => "t4",
                    "Material Model" => "Test",
                ),
            ),
            "Damage Models" => Dict(
                "a" => Dict("value" => 3, "name" => "t"),
                "c" => Dict("value" => [1 2], "value2" => 1, "name" => "t2"),
            ),
            "Thermal Models" => Dict(
                "therm" => Dict("value" => "hot", "bool" => true, "name" => "t3"),
            ),
        ),
    )
    Model_Factory.read_properties(params, test_data_manager_read_properties, false)


    @test test_data_manager_read_properties.get_property(1, "Material Model", "value") ==
          params["Models"]["Material Models"]["a"]["value"]
    @test test_data_manager_read_properties.get_property(2, "Material Model", "value") ==
          params["Models"]["Material Models"]["c"]["value"]
    @test test_data_manager_read_properties.get_property(2, "Material Model", "value2") ==
          params["Models"]["Material Models"]["c"]["value2"]
    @test test_data_manager_read_properties.get_property(3, "Material Model", "value") ==
          params["Models"]["Material Models"]["a"]["value"]
    @test test_data_manager_read_properties.get_property(3, "Damage Model", "value") ==
          params["Models"]["Damage Models"]["a"]["value"]
    @test test_data_manager_read_properties.get_property(3, "Thermal Model", "value") ==
          params["Models"]["Thermal Models"]["therm"]["value"]
    @test test_data_manager_read_properties.get_property(3, "Thermal Model", "bool") ==
          params["Models"]["Thermal Models"]["therm"]["bool"]

    Model_Factory.read_properties(params, test_data_manager_read_properties, true)
    @test test_data_manager_read_properties.get_property(1, "Material Model", "value") ==
          params["Models"]["Material Models"]["a"]["value"]
    @test test_data_manager_read_properties.get_property(2, "Material Model", "value") ==
          params["Models"]["Material Models"]["c"]["value"]
    @test test_data_manager_read_properties.get_property(2, "Material Model", "value2") ==
          params["Models"]["Material Models"]["c"]["value2"]
    @test test_data_manager_read_properties.get_property(3, "Material Model", "value") ==
          params["Models"]["Material Models"]["a"]["value"]
    @test test_data_manager_read_properties.get_property(3, "Damage Model", "value") ==
          params["Models"]["Damage Models"]["a"]["value"]
    @test test_data_manager_read_properties.get_property(3, "Thermal Model", "value") ==
          params["Models"]["Thermal Models"]["therm"]["value"]
    @test test_data_manager_read_properties.get_property(3, "Thermal Model", "bool") ==
          params["Models"]["Thermal Models"]["therm"]["bool"]
end