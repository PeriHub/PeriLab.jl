# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const PS = PeriLab.ParameterSpec
const ID = PeriLab.InputDeck
function ut_st_solver(dict)
    ctx = PS.ParseContext()
    s = PS.parse_section(ID.SolverParams, dict, "Solver", ctx)
    @test isempty(ctx.errors)
    return s
end

@testset "active options and solver name" begin
    s = ut_st_solver(Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                      "Static" => Dict{String,Any}("m" => 3)))
    @test ID.active_options(s) === s.static
    @test ID.active_options(s) isa ID.SolverOptions
    @test ID.solver_name(s) == "Static"
    @test ID.active_options(s).safety_factor === 1.0
    names = [ID.solver_name(ut_st_solver(Dict{String,Any}("Initial Time" => 0.0,
                                                          "Final Time" => 1.0,
                                                          key => Dict{String,Any}())))
             for key in ("Verlet", "Static", "Linear Static Matrix Based",
                         "Verlet Matrix Based", "Newmark")]
    @test names ==
          ["Verlet", "Static", "Linear Static Matrix Based", "Verlet Matrix Based", "Newmark"]
end

@testset "number of steps given or not" begin
    s = ut_st_solver(Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                      "Static" => Dict{String,Any}()))
    @test s.number_of_steps === nothing
    s = ut_st_solver(Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                      "Number of Steps" => 7,
                                      "Static" => Dict{String,Any}("Fixed dt" => 0.1)))
    @test s.number_of_steps === 7 && s.static.fixed_dt === 0.1
end

@testset "start and end time depend on the current time" begin
    s = ut_st_solver(Dict{String,Any}("Initial Time" => 0.5, "Final Time" => 2.0,
                                      "Verlet" => Dict{String,Any}()))
    @test ID.start_time(s, 0.0) === 0.5
    @test ID.start_time(s, 1.0) === 1.0          # already past the initial time
    @test ID.end_time(s, 1.0) === 2.0
    s = ut_st_solver(Dict{String,Any}("Additional Time" => 3.0, "Step ID" => 2,
                                      "Verlet" => Dict{String,Any}()))
    @test ID.start_time(s, 1.5) === 1.5          # later multistep step
    @test ID.end_time(s, 1.5) === 4.5
    @test_throws PeriLab.PeriLabError ID.start_time(s, 0.0)
end

@testset "model options keep the old order" begin
    s = ut_st_solver(Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                      "Damage Models" => true, "Thermal Models" => true,
                                      "Additive Models" => true,
                                      "Verlet" => Dict{String,Any}()))
    @test ID.model_options(s) ==
          ["Additive", "Damage", "Pre_Calculation", "Thermal", "Material"]
end

function ut_st_deck(; multistep = nothing)
    deck = Dict{String,Any}("Discretization" => Dict{String,Any}("Type" => "Text File",
                                                                 "Input Mesh File" => "m.txt"),
                            "Blocks" => Dict{String,Any}("b" => Dict{String,Any}("Block ID" => 1,
                                                                                 "Density" => 1.0,
                                                                                 "Horizon" => 1.0)),
                            "Models" => Dict{String,Any}())
    if multistep === nothing
        deck["Solver"] = Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                          "Verlet" => Dict{String,Any}())
    else
        deck["Multistep Solver"] = multistep
    end
    input, ctx = ID.read_input(deck)
    @test isempty(ctx.errors)
    return input
end

@testset "solver steps are sorted" begin
    step(id) = Dict{String,Any}("Step ID" => id, "Initial Time" => 0.0, "Final Time" => 1.0,
                                "Verlet" => Dict{String,Any}())
    input = ut_st_deck(multistep = Dict{String,Any}("c" => step(3), "a" => step(1),
                                                    "b" => step(2)))
    @test ID.solver_steps(input) == [1, 2, 3]
    @test ID.solver_step(input, 2).step_id === 2
    @test_throws PeriLab.PeriLabError ID.solver_step(input, 9)
    input = ut_st_deck()
    @test ID.solver_steps(input) == [-1]
    @test ID.solver_step(input, -1) === input.sections.solver
end

@testset "Solver_Manager.init keeps its docstring" begin
    @test occursin("Initialize the solver", string(@doc PeriLab.Solver_Manager.init))
end
