# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const PS = PeriLab.ParameterSpec
const ID = PeriLab.InputDeck
const SM = PeriLab.Solver_Manager

function ut_minimal_deck()
    return Dict{String,Any}("Discretization" => Dict{String,Any}("Type" => "Text File",
                                                                 "Input Mesh File" => "mesh.txt"),
                            "Blocks" => Dict{String,Any}("block_1" => Dict{String,Any}("Block ID" => 1,
                                                                                       "Density" => 1.0,
                                                                                       "Horizon" => 1.0)),
                            "Models" => Dict{String,Any}(),
                            "Solver" => Dict{String,Any}("Initial Time" => 0.0,
                                                         "Final Time" => 1.0,
                                                         "Verlet" => Dict{String,Any}()))
end

function ut_solver(dict)
    ctx = PS.ParseContext()
    return PS.parse_section(ID.SolverParams, dict, "Solver", ctx), ctx
end

@testset "Verlet solver with defaults" begin
    s, ctx = ut_solver(Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1,
                                        "Verlet" => Dict{String,Any}("Safety Factor" => 0.9)))
    @test isempty(ctx.errors)
    @test s.final_time === 1.0 && s.number_of_steps === nothing && s.maximum_damage === Inf
    @test s.material_models && s.pre_calculation_models && !s.damage_models
    @test s.options isa SM.Verlet_Solver.VerletParams && s.options.safety_factor === 0.9 &&
          s.options.fixed_dt === -1.0
    @test s.options.numerical_damping === 0.0
    @test ID.solver_name(s) == "Verlet"
end

@testset "Static, matrix based, Newmark and model reduction options" begin
    s, ctx = ut_solver(Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                        "Static" => Dict{String,Any}("Fixed dt" => 0.1,
                                                                     "Residual scaling" => 70000)))
    @test isempty(ctx.errors)
    @test s.options.residual_scaling === 70000.0 && s.options.m === 15
    @test s.options.maximum_number_of_iterations === 100
    @test s.options.fixed_dt === 0.1
    s, ctx = ut_solver(Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                        "Verlet Matrix Based" => Dict{String,Any}("Model Reduction" => Dict{String,Any}("Type" => "Craig Bampton",
                                                                                                                     "Number of Modes" => 5))))
    @test isempty(ctx.errors)
    @test s.options.model_reduction.number_of_modes === 5
    @test s.options.model_reduction.material_point_region === true
    s, ctx = ut_solver(Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                        "Newmark" => Dict{String,Any}("Matrix Update" => true)))
    @test isempty(ctx.errors)
    @test s.options.matrix_update && s.options.newmark_delta === 0.5
    @test s.options.newmark_alpha === nothing
end

@testset "solver rules" begin
    s, ctx = ut_solver(Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0))
    @test ctx.errors[1].path == "Solver"
    @test ctx.errors[1].message ==
          "one solver is required: Linear Static Matrix Based, Newmark, Static, Verlet, Verlet Matrix Based"
    s, ctx = ut_solver(Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                        "Verlet" => Dict{String,Any}(),
                                        "Static" => Dict{String,Any}()))
    @test ctx.errors[1].message == "only one solver may be given, found: Static, Verlet"
    s, ctx = ut_solver(Dict{String,Any}("Initial Time" => 0.0, "Verlet" => Dict{String,Any}()))
    @test ctx.errors[1].message == "\"Final Time\" or \"Additional Time\" is required"
    s, ctx = ut_solver(Dict{String,Any}("Additional Time" => 1.0, "Step ID" => 2,
                                        "Verlet" => Dict{String,Any}()))
    @test isempty(ctx.errors) && s.step_id === 2
    s, ctx = ut_solver(Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                        "External" => Dict{String,Any}(),
                                        "Verlet" => Dict{String,Any}()))
    @test ctx.errors[1].path == "Solver.External"
    @test startswith(ctx.errors[1].message, "unknown key")
end

@testset "Model Reduction: false" begin
    s, ctx = ut_solver(Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                        "Verlet Matrix Based" => Dict{String,Any}("Model Reduction" => false)))
    @test isempty(ctx.errors)
    @test s.options.model_reduction === nothing
    # solvers that never read a key do not declare it
    s, ctx = ut_solver(Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                        "Newmark" => Dict{String,Any}("Model Reduction" => false)))
    @test ctx.errors[1].path == "Solver.Newmark.\"Model Reduction\""
    @test startswith(ctx.errors[1].message, "unknown key")
end

@testset "misspelled solver" begin
    s, ctx = ut_solver(Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                        "Verlte" => Dict{String,Any}()))
    @test s === nothing
    @test any(e -> e.path == "Solver.Verlte" &&
                   e.message == "unknown key — did you mean \"Verlet\"?", ctx.errors)
end

@testset "multistep steps with different solvers" begin
    deck = ut_minimal_deck()
    delete!(deck, "Solver")
    deck["Multistep Solver"] = Dict{String,Any}("Step 1" => Dict{String,Any}("Step ID" => 1,
                                                                              "Initial Time" => 0.0,
                                                                              "Final Time" => 1.0,
                                                                              "Verlet" => Dict{String,Any}()),
                                                "Step 2" => Dict{String,Any}("Step ID" => 2,
                                                                              "Additional Time" => 1.0,
                                                                              "Static" => Dict{String,Any}("m" => 4)))
    input, ctx = ID.read_input(deck)
    @test isempty(ctx.errors)
    @test ID.solver_step(input, 1).options isa SM.Verlet_Solver.VerletParams
    @test ID.solver_step(input, 2).options isa SM.Static_Solver.StaticParams
    @test ID.solver_step(input, 2).options.m === 4
end

@testset "Fallback names one solver" begin
    s, ctx = ut_solver(Dict{String,Any}("Final Time" => 1.0, "Static" => Dict{String,Any}(),
                                        "Fallback" => Dict{String,Any}("Verlet" => Dict{String,Any}("Safety Factor" => 0.9))))
    @test isempty(ctx.errors)
    @test s.fallback.options isa PeriLab.Solver_Manager.Verlet_Solver.VerletParams
    @test s.fallback.options.safety_factor == 0.9
    @test PeriLab.InputDeck.solver_name(s.fallback.options) == "Verlet"
    s, ctx = ut_solver(Dict{String,Any}("Final Time" => 1.0, "Static" => Dict{String,Any}()))
    @test isempty(ctx.errors) && s.fallback === nothing
end

@testset "Fallback errors" begin
    function fallback_errors(fallback)
        s, ctx = ut_solver(Dict{String,Any}("Final Time" => 1.0, "Static" => Dict{String,Any}(),
                                            "Fallback" => fallback))
        return join((e.path * ": " * e.message for e in ctx.errors), "\n")
    end
    @test occursin("one solver is required", fallback_errors(Dict{String,Any}()))
    @test occursin("Verlett", fallback_errors(Dict{String,Any}("Verlett" => Dict{String,Any}())))
    @test occursin("Fallback",
                   fallback_errors(Dict{String,Any}("Verlet" => Dict{String,Any}(),
                                                    "Fallback" => Dict{String,Any}("Static" => Dict{String,Any}()))))
end
