include("../../Support/data_manager.jl")
include("../Physics_Factory.jl")

using Test
import .Physics
import .Data_manager
@testset "ut_get_block_model_definition" begin
    testDatamanager = Data_manager
    block_list = [1, 2, 3]
    testDatamanager.set_block_list(block_list)
    prop_keys = testDatamanager.init_property()
    params = Dict("Blocks" => Dict("block_1" => Dict("Material Models" => "a"), "block_2" => Dict("Material Models" => "c"), "block_3" => Dict("Material Models" => "a", "Damage Models" => "a", "Thermal Models" => "therm")), "Material Models" => Dict("a" => Dict("value" => 1), "c" => Dict("value" => [1 2], "value2" => 1)), "Damage Models" => Dict("a" => Dict("value" => 3), "c" => Dict("value" => [1 2], "value2" => 1)), "Thermal Models" => Dict("therm" => Dict("value" => "hot", "bool" => true)))

    for block in block_list
        Physics.get_block_model_definition(params, block, prop_keys, testDatamanager.set_properties)
    end
    @test testDatamanager.get_property(1, "Material Models", "value") == params["Material Models"]["a"]["value"]
    @test testDatamanager.get_property(2, "Material Models", "value") == params["Material Models"]["c"]["value"]
    @test testDatamanager.get_property(2, "Material Models", "value2") == params["Material Models"]["c"]["value2"]
    @test testDatamanager.get_property(3, "Material Models", "value") == params["Material Models"]["a"]["value"]
    @test testDatamanager.get_property(3, "Damage Models", "value") == params["Damage Models"]["a"]["value"]
    @test testDatamanager.get_property(3, "Thermal Models", "value") == params["Thermal Models"]["therm"]["value"]
    @test testDatamanager.get_property(3, "Thermal Models", "bool") == params["Thermal Models"]["therm"]["bool"]
end

@testset "ut_read_properties" begin
    testDatamanager_read_properties = Data_manager
    block_list = [1, 2, 3]
    testDatamanager_read_properties.set_block_list(block_list)

    params = Dict("Blocks" => Dict("block_1" => Dict("Material Models" => "a"), "block_2" => Dict("Material Models" => "c"), "block_3" => Dict("Material Models" => "a", "Damage Models" => "a", "Thermal Models" => "therm")), "Material Models" => Dict("a" => Dict("value" => 1), "c" => Dict("value" => [1 2], "value2" => 1)), "Damage Models" => Dict("a" => Dict("value" => 3), "c" => Dict("value" => [1 2], "value2" => 1)), "Thermal Models" => Dict("therm" => Dict("value" => "hot", "bool" => true)))

    Physics.read_properties(params, testDatamanager_read_properties)

    @test testDatamanager_read_properties.get_property(1, "Material Models", "value") == params["Material Models"]["a"]["value"]
    @test testDatamanager_read_properties.get_property(2, "Material Models", "value") == params["Material Models"]["c"]["value"]
    @test testDatamanager_read_properties.get_property(2, "Material Models", "value2") == params["Material Models"]["c"]["value2"]
    @test testDatamanager_read_properties.get_property(3, "Material Models", "value") == params["Material Models"]["a"]["value"]
    @test testDatamanager_read_properties.get_property(3, "Damage Models", "value") == params["Damage Models"]["a"]["value"]
    @test testDatamanager_read_properties.get_property(3, "Thermal Models", "value") == params["Thermal Models"]["therm"]["value"]
    @test testDatamanager_read_properties.get_property(3, "Thermal Models", "bool") == params["Thermal Models"]["therm"]["bool"]
end