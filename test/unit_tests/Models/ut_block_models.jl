# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

const UT_MF = PeriLab.Solver_Manager.Model_Factory

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
    PeriLab.Data_Manager.set_block_models(1, PeriLab.Data_Manager.BlockModels(additive = p))
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
    PeriLab.Data_Manager.set_block_models(1, PeriLab.Data_Manager.BlockModels(degradation = p))
    UT_MF.Degradation.init_model([1, 2], 1)
    UT_MF.Degradation.compute_model([1, 2], p, 1, 0.0, 1.0)
    @test PeriLab.Data_Manager.get_field("Active") == [true, false]
    bd = PeriLab.Data_Manager.get_bond_damage("NP1")
    @test bd[2][1] == 0.0 && bd[1][1] == 0.0
    @test length(methods(UT_MF.Degradation.Thermal_Decomposition.compute_model)) == 1
end

function ut_thermal_model(raw)
    return typed_model(:thermal, raw; name_key = "Thermal Model")
end

@testset "thermal composite runs every part" begin
    PeriLab.Data_Manager.initialize_data()
    th = ut_thermal_model(Dict("Thermal Model" => "Thermal Flow + Heat Transfer",
                               "Thermal Conductivity" => 1.0,
                               "Heat Transfer Coefficient" => 1.0,
                               "Environmental Temperature" => 30))
    parts = UT_MF.Thermal.model_parts(th.model)
    @test nameof.(typeof.(parts)) == (:ThermalFlowParams, :HeatTransferParams)
    th2 = ut_thermal_model(Dict("Thermal Model" => "Heat Transfer+Thermal Flow",
                                "Thermal Conductivity" => 1.0,
                                "Heat Transfer Coefficient" => 1.0,
                                "Environmental Temperature" => 30))
    @test nameof.(typeof.(UT_MF.Thermal.model_parts(th2.model))) ==
          (:HeatTransferParams, :ThermalFlowParams)
    single = ut_thermal_model(Dict("Thermal Model" => "Thermal Expansion",
                                   "Thermal Expansion Coefficient" => 1.0))
    @test UT_MF.Thermal.model_parts(single.model) == (single.model,)
end

@testset "environmental temperature expression" begin
    HT = UT_MF.Thermal.Heat_Transfer
    @test HT.environmental_temperature(30.0, 5.0) == 30.0
    @test HT.environmental_temperature("20+t", 5.0) == 25.0
end

@testset "thermal critical time step" begin
    PeriLab.Data_Manager.initialize_data()
    th = ut_thermal_model(Dict("Thermal Model" => "Heat Transfer",
                               "Thermal Conductivity" => 0.12,
                               "Heat Transfer Coefficient" => 1.0,
                               "Environmental Temperature" => 30))
    PeriLab.Data_Manager.set_block_models(1, PeriLab.Data_Manager.BlockModels(thermal = th))
    @test UT_MF.block_thermal_conductivity(1) == 0.12
    @test UT_MF.block_thermal_conductivity(2) === nothing
end

@testset "heat capacity from the typed blocks" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(3)
    input = typed_input(Dict("Blocks" => Dict("block_1" => Dict("Block ID" => 1,
                                                                "Density" => 1.0,
                                                                "Horizon" => 1.0,
                                                                "Specific Heat Capacity" => 5.0),
                                              "block_2" => Dict("Block ID" => 2,
                                                                "Density" => 1.0,
                                                                "Horizon" => 1.0))))
    hc = zeros(3)
    UT_MF.set_heat_capacity(input, Dict(1 => [1, 2]), hc)
    @test hc == [5.0, 5.0, 0.0]
    @test_logs (:error,
                "Specific Heat Capacity of block_2 is not defined") @test_throws PeriLab.PeriLabError UT_MF.set_heat_capacity(input,
                                                                                                                                Dict(2 => [3]),
                                                                                                                                hc)
end

@testset "read_properties stores pre-calculations" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(2)
    PeriLab.Data_Manager.set_dof(2)
    PeriLab.Data_Manager.set_block_name_list(["block_1", "block_2", "block_3"])
    PeriLab.Data_Manager.set_block_id_list([1, 2, 3])
    blocks = Dict("block_1" => Dict("Block ID" => 1, "Density" => 1.0, "Horizon" => 1.0),
                  "block_2" => Dict("Block ID" => 2, "Density" => 1.0, "Horizon" => 1.0,
                                    "Pre Calculation Model" => "Pre"),
                  "block_3" => Dict("Block ID" => 3, "Density" => 1.0, "Horizon" => 1.0,
                                    "Pre Calculation Model" => "Off"))
    models = Dict("Pre Calculation Global" => Dict("Shape Tensor" => true,
                                                   "Deformed Bond Geometry" => true),
                  "Pre Calculation Models" => Dict("Pre" => Dict("Deformation Gradient" => true),
                                                   "Off" => Dict("Shape Tensor" => false)))
    input = typed_input(Dict("Blocks" => blocks, "Models" => models))
    UT_MF.read_properties(input, false)
    get(b) = PeriLab.Data_Manager.get_block_models(b).pre_calculation
    @test get(1) == ["Deformed Bond Geometry", "Shape Tensor"]
    @test get(2) == ["Deformation Gradient"]          # replaces the global switches
    @test isempty(get(3))                             # nothing active
    @test UT_MF.has_block_model(1, "Pre Calculation Model")
    @test !UT_MF.has_block_model(3, "Pre Calculation Model")
    @test UT_MF.Pre_Calculation.order_pre_calculations(["Axis Symmetric", "Shape Tensor",
                                                        "Deformed Bond Geometry"]) ==
          ["Deformed Bond Geometry", "Shape Tensor", "Axis Symmetric"]
end

@testset "licensed additive models load only for additive decks" begin
    # a half-configured license aborts when licensed modules are loaded
    withenv("LICENSE_SERVER_URL" => "http://localhost:1", "PERIHUB_LICENSE_KEY" => nothing,
            "LICENSED_MODULES_CONFIG" => nothing, "LICENSED_MODULES_DIR" => nothing) do
        mechanical = Dict("PeriLab" => Dict("Models" => Dict("Material Models" => Dict())))
        @test isempty(UT_MF.Additive.load_licensed_models(mechanical))
        additive = Dict("PeriLab" => Dict("Models" => Dict("Additive Models" => Dict())))
        @test_throws PeriLab.PeriLabError UT_MF.Additive.load_licensed_models(additive)
    end
end
