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
    PeriLab.Data_Manager.create_constant_node_scalar_field("Volume", Float64;
                                                           default_value = 1)
    th = typed_model(:thermal,
                     Dict("Thermal Model" => "Heat Transfer",
                          "Heat Transfer Coefficient" => 1.0,
                          "Environmental Temperature" => 30);
                     name_key = "Thermal Model")
    PeriLab.Data_Manager.set_block_models(1, PeriLab.Data_Manager.BlockModels(thermal = th))
    PeriLab.Solver_Manager.Model_Factory.Thermal.init_model([1], 1)
    @test PeriLab.Data_Manager.has_key("Bond Face")
    @test PeriLab.Data_Manager.has_key("Extent")
end

@testset "licensed thermal models" begin
    T = PeriLab.Solver_Manager.Model_Factory.Thermal
    # no thermal models in the deck: nothing loaded, no license contact
    @test T.load_licensed_models(Dict("PeriLab" => Dict("Models" => Dict()))) == Any[]
    @test T.load_licensed_models(Dict("x" => 1)) == Any[]
    # local dev mode with a dummy licensed thermal module
    dir = mktempdir()
    write(joinpath(dir, "Dummy_thermal.jl"), """
        module Dummy_Thermal
        using ......ParameterSpec: @params, register_thermal
        @params struct DummyThermalParams
            value::Float64 = opt("Value"; default = 1.0)
        end
        __init__() = register_thermal("Dummy Licensed Thermal", DummyThermalParams)
        post_process(input, files) = (global called = files; nothing)
        end""")
    write(joinpath(dir, "config.yaml"),
          "features:\n  t:\n    modules:\n      - file: Dummy_thermal.jl\n        type: Thermal\n")
    # licensed_modules reloads a `.env` above the working directory (DotEnv), which would
    # undo the withenv below once an earlier test registered it; run from an empty directory
    saved = copy(T.licensed)
    cd(mktempdir()) do
        withenv("LICENSED_MODULES_CONFIG" => joinpath(dir, "config.yaml"),
                "LICENSED_MODULES_DIR" => dir) do
            empty!(T.licensed)
            deck = Dict("PeriLab" => Dict("Models" => Dict("Thermal Models" => Dict())))
            mods = T.load_licensed_models(deck)
            @test any(m -> m["Module Name"] == "Dummy_Thermal", mods)
            @test any(first.(PeriLab.ParameterSpec.registered_models(:thermal)) .==
                      "Dummy Licensed Thermal")
            T.post_process_licensed(nothing, ["a.e"])
            @test Base.invokelatest(getfield, T.Dummy_Thermal, :called) == ["a.e"]
        end
    end
    append!(empty!(T.licensed), saved)
end
