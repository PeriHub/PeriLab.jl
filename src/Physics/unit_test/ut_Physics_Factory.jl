include("../../Support/data_manager.jl")
include("../Physics_Factory.jl")
include("../../Support/Parameters/parameter_handling.jl")
using Test
import .Physics
import .Data_manager
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

    params = Dict("Blocks" => Dict("block_1" => Dict("Material Model" => "a"), "block_2" => Dict("Material Model" => "c"), "block_3" => Dict("Material Model" => "a", "Damage Model" => "a", "Thermal Model" => "therm")), "Physics" => Dict("Material Models" => Dict("a" => Dict("value" => 1), "c" => Dict("value" => [1 2], "value2" => 1)), "Damage Models" => Dict("a" => Dict("value" => 3), "c" => Dict("value" => [1 2], "value2" => 1)), "Thermal Models" => Dict("therm" => Dict("value" => "hot", "bool" => true))))

    Physics.read_properties(params, testDatamanager_read_properties)

    @test testDatamanager_read_properties.get_property(1, "Material Model", "value") == params["Physics"]["Material Models"]["a"]["value"]
    @test testDatamanager_read_properties.get_property(2, "Material Model", "value") == params["Physics"]["Material Models"]["c"]["value"]
    @test testDatamanager_read_properties.get_property(2, "Material Model", "value2") == params["Physics"]["Material Models"]["c"]["value2"]
    @test testDatamanager_read_properties.get_property(3, "Material Model", "value") == params["Physics"]["Material Models"]["a"]["value"]
    @test testDatamanager_read_properties.get_property(3, "Damage Model", "value") == params["Physics"]["Damage Models"]["a"]["value"]
    @test testDatamanager_read_properties.get_property(3, "Thermal Model", "value") == params["Physics"]["Thermal Models"]["therm"]["value"]
    @test testDatamanager_read_properties.get_property(3, "Thermal Model", "bool") == params["Physics"]["Thermal Models"]["therm"]["bool"]
end