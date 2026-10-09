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
                                                    "Models" => Dict{Any,Any}("Material Models" => Dict{Any,Any}("mat_1" => Dict{Any,Any}("Material Model" => "Bond-based Elastic"))),
                                                    "Solver" => Dict{Any,Any}("Initial Time" => 0.0,
                                                                              "Final Time" => 1.0,
                                                                              "Verlet" => Dict{Any,Any}())))
end

@testset "validate_input returns the typed input" begin
    params = ut_valid_params()
    deck, input = PeriLab.IO.validate_input(params)
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
                  Material Model: Bond-based Elastic
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


@testset "no_strict flag" begin
    params = ut_valid_params()
    params["PeriLab"]["Solver"]["Unused Option"] = 1
    @test_throws PeriLab.PeriLabError PeriLab.IO.validate_input(params)
    @test first(PeriLab.IO.validate_input(params; no_strict = true)) ===
          params["PeriLab"]
    params["PeriLab"]["Solver"]["final time"] = 2.0       # case-only typo: always an error
    @test_throws PeriLab.PeriLabError PeriLab.IO.validate_input(params;
                                                                              no_strict = true)
    @test PeriLab.parse_commandline(["--no_strict", "a.yaml"])["no_strict"] == true
    @test PeriLab.parse_commandline(["a.yaml"])["no_strict"] == false
end

@testset "material model names are validated" begin
    params = ut_valid_params()
    params["PeriLab"]["Models"]["Material Models"]["mat_1"]["Material Model"] = 5
    @test_throws PeriLab.PeriLabError PeriLab.IO.validate_input(params)
end

@testset "an unknown material key aborts with the typed error" begin
    params = ut_valid_params()
    params["PeriLab"]["Models"]["Material Models"]["mat_1"]["Youngs Modulus"] = 1.0
    @test_logs (:error, r"did you mean \"Young's Modulus\"") match_mode=:any @test_throws PeriLab.PeriLabError begin
        PeriLab.IO.validate_input(params)
    end
end

