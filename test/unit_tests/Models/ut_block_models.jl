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

@testset "additive dispatch reads the block model" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(3)
    AT = UT_MF.Additive.Additive_template
    p = AT.AdditiveTemplateParams()
    PeriLab.Data_Manager.set_block_model("Additive Model", 1, p)
    @test UT_MF.has_block_model(1, "Additive Model")
    @test UT_MF.block_model_parameters(1, "Additive Model") === p
    UT_MF.Additive.init_model([1, 2, 3], 1)
    UT_MF.Additive.compute_model([1, 2, 3], p, 1, 0.0, 1.0)
    UT_MF.Additive.fields_for_local_synchronization("Additive Model", 1)
    @test length(methods(AT.compute_model)) == 1
end

@testset "licensed additive models load without a license" begin
    withenv("LICENSE_SERVER_URL" => nothing, "PERIHUB_LICENSE_KEY" => nothing,
            "LICENSED_MODULES_CONFIG" => nothing, "LICENSED_MODULES_DIR" => nothing) do
        @test isempty(UT_MF.Additive.load_licensed_models())
    end
end

@testset "thermal decomposition on the typed interface" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(2)
    PeriLab.Data_Manager.set_dof(2)
    nn = PeriLab.Data_Manager.create_constant_node_scalar_field("Number of Neighbors", Int64)
    nn .= 1
    nlist = PeriLab.Data_Manager.create_constant_bond_scalar_state("Neighborhoodlist", Int64)
    nlist[1] = [2]
    nlist[2] = [1]
    PeriLab.Data_Manager.create_bond_scalar_state("Bond Damage", Float64; default_value = 1)
    PeriLab.Data_Manager.create_constant_node_scalar_field("Active", Bool;
                                                           default_value = true)
    _, temperature = PeriLab.Data_Manager.create_node_scalar_field("Temperature", Float64)
    temperature .= [10.0, 60.0]
    UT_MF.Degradation.init_fields()
    p = typed_model(:degradation, Dict("Degradation Model" => "Thermal Decomposition",
                                       "Decomposition Temperature" => 50);
                    name_key = "Degradation Model")
    PeriLab.Data_Manager.set_block_model("Degradation Model", 1, p)
    UT_MF.Degradation.init_model([1, 2], 1)
    UT_MF.Degradation.compute_model([1, 2], p, 1, 0.0, 1.0)
    @test PeriLab.Data_Manager.get_field("Active") == [true, false]
    bd = PeriLab.Data_Manager.get_bond_damage("NP1")
    @test bd[2][1] == 0.0 && bd[1][1] == 0.0
    @test length(methods(UT_MF.Degradation.Thermal_Decomposition.compute_model)) == 1
end
