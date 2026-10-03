# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test

function ut_valid_params()
    return Dict{Any,Any}("PeriLab" => Dict{Any,Any}("Discretization" => Dict{Any,Any}("Type" => "Text File",
                                                                                     "Input Mesh File" => "m.txt"),
                                                    "Blocks" => Dict{Any,Any}("block_1" => Dict{Any,Any}("Block ID" => 1,
                                                                                                         "Density" => 1.0,
                                                                                                         "Horizon" => 1.0)),
                                                    "Models" => Dict{Any,Any}("Material Models" => Dict{Any,Any}("mat_1" => Dict{Any,Any}("Material Model" => "a"))),
                                                    "Solver" => Dict{Any,Any}("Initial Time" => 0.0,
                                                                              "Final Time" => 1.0,
                                                                              "Verlet" => Dict{Any,Any}())))
end

@testset "validate_yaml returns the deck dict unchanged" begin
    params = ut_valid_params()
    @test PeriLab.Parameter_Handling.validate_yaml(params) === params["PeriLab"]
end

@testset "validate_yaml aborts with all input errors" begin
    params = ut_valid_params()
    params["PeriLab"]["Blocks"]["block_1"]["Horizon"] = "1.0"
    params["PeriLab"]["Bocks"] = Dict{Any,Any}()
    @test_throws PeriLab.PeriLabError PeriLab.Parameter_Handling.validate_yaml(params)
end

@testset "no_strict flag" begin
    params = ut_valid_params()
    params["PeriLab"]["Solver"]["Unused Option"] = 1
    @test_throws PeriLab.PeriLabError PeriLab.Parameter_Handling.validate_yaml(params)
    @test PeriLab.Parameter_Handling.validate_yaml(params; no_strict = true) ===
          params["PeriLab"]
    params["PeriLab"]["Solver"]["final time"] = 2.0       # case-only typo: always an error
    @test_throws PeriLab.PeriLabError PeriLab.Parameter_Handling.validate_yaml(params;
                                                                              no_strict = true)
    @test PeriLab.parse_commandline(["--no_strict", "a.yaml"])["no_strict"] == true
    @test PeriLab.parse_commandline(["a.yaml"])["no_strict"] == false
end

@testset "legacy model validation still applies" begin
    params = ut_valid_params()
    params["PeriLab"]["Models"]["Material Models"]["mat_1"]["Material Model"] = 5
    @test_throws PeriLab.PeriLabError PeriLab.Parameter_Handling.validate_yaml(params)
end

@testset "validate_input returns the typed input" begin
    params = ut_valid_params()
    deck, input = PeriLab.Parameter_Handling.validate_input(params)
    @test deck === params["PeriLab"]
    @test input isa PeriLab.InputDeck.PeriLabInput
    @test input.sections.solver.verlet !== nothing
end

@testset "read_input_deck" begin
    dir = mktempdir()
    file = joinpath(dir, "deck.yaml")
    write(file,
          """
          PeriLab:
            Discretization:
              Type: Text File
              Input Mesh File: m.txt
            Blocks:
              block_1:
                Block ID: 1
                Density: 1.0
                Horizon: 1.0
            Models:
              Material Models:
                mat_1:
                  Material Model: a
            Solver:
              Initial Time: 0.0
              Final Time: 1.0
              Number of Steps: 4
              Verlet:
                Safety Factor: 0.5
          """)
    deck, input = PeriLab.IO.read_input_deck(file)
    @test deck["Solver"]["Number of Steps"] == 4
    @test input.sections.solver.number_of_steps === 4
    @test PeriLab.InputDeck.solver_steps(input) == [-1]
end
