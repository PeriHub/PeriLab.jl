# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

const UT_MF = PeriLab.Solver_Manager.Model_Factory

@testset "generic block model slot" begin
    PeriLab.Data_Manager.initialize_data()
    @test PeriLab.Data_Manager.get_block_model("Thermal Model", 1) === nothing
    PeriLab.Data_Manager.set_block_model("Thermal Model", 1, :x)
    @test PeriLab.Data_Manager.get_block_model("Thermal Model", 1) === :x
    @test PeriLab.Data_Manager.get_block_model("Thermal Model", 2) === nothing
    @test PeriLab.Data_Manager.get_block_model("Additive Model", 1) === nothing
    @test "Block Models" in PeriLab.Data_Manager.CHECKPOINT_EXCLUDED_KEYS
end

@testset "block_typed_model" begin
    blocks = Dict("block_1" => Dict("Block ID" => 1, "Density" => 1.0, "Horizon" => 1.0,
                                    "Damage Model" => "Dam"),
                  "block_2" => Dict("Block ID" => 2, "Density" => 1.0, "Horizon" => 1.0))
    models = Dict("Damage Models" => Dict("Dam" => Dict("Damage Model" => "Critical Stretch",
                                                        "Critical Value" => 0.1)))
    input = typed_input(Dict("Blocks" => blocks, "Models" => models))
    d = UT_MF.block_typed_model(input, "block_1", :damage_model, "Damage Model",
                                input.damages)
    @test d === input.damages["Dam"]
    @test UT_MF.block_typed_model(input, "block_2", :damage_model, "Damage Model",
                                  input.damages) === nothing
    @test_logs (:error,
                "Damage Model model with name Dam is defined in blocks, but missing in the Damage Models definition.") @test_throws PeriLab.PeriLabError UT_MF.block_typed_model(input,
                                                                                                                                                                              "block_1",
                                                                                                                                                                              :damage_model,
                                                                                                                                                                              "Damage Model",
                                                                                                                                                                              Dict{String,Any}())
end

@testset "typed_model helper" begin
    m = typed_model(:damage, Dict("Damage Model" => "Critical Stretch", "Critical Value" => 0.1);
                    name_key = "Damage Model")
    @test m isa PeriLab.ParameterSpec.WithBase
end
