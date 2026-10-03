# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const PS = PeriLab.ParameterSpec
const ID = PeriLab.InputDeck

function ut_solver(dict)
    ctx = PS.ParseContext()
    return PS.parse_section(ID.SolverParams, dict, "Solver", ctx), ctx
end

@testset "Verlet solver with defaults" begin
    s, ctx = ut_solver(Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1,
                                        "Verlet" => Dict{String,Any}("Safety Factor" => 0.9)))
    @test isempty(ctx.errors)
    @test s.final_time === 1.0 && s.number_of_steps === 1 && s.maximum_damage === Inf
    @test s.material_models && s.pre_calculation_models && !s.damage_models
    @test s.verlet.safety_factor === 0.9 && s.verlet.fixed_dt === -1.0
    @test s.verlet.numerical_damping === 0.0
    @test s.static === nothing && s.newmark === nothing
end

@testset "Static, matrix based, Newmark and model reduction options" begin
    s, ctx = ut_solver(Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                        "Static" => Dict{String,Any}("NLSolve" => true,
                                                                     "Fixed dt" => 0.1,
                                                                     "Residual scaling" => 70000)))
    @test isempty(ctx.errors)
    @test s.static.residual_scaling === 70000.0 && s.static.m === 15
    @test s.static.maximum_number_of_iterations === 100 && s.static.nlsolve === true
    @test s.static.fixed_dt === 0.1
    s, ctx = ut_solver(Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                        "Verlet Matrix Based" => Dict{String,Any}("Model Reduction" => Dict{String,Any}("Type" => "Craig Bampton",
                                                                                                                     "Number of Modes" => 5))))
    @test isempty(ctx.errors)
    @test s.verlet_matrix_based.model_reduction.number_of_modes === 5
    @test s.verlet_matrix_based.model_reduction.material_point_region === true
    s, ctx = ut_solver(Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                        "Newmark" => Dict{String,Any}("Matrix Update" => true)))
    @test isempty(ctx.errors)
    @test s.newmark.matrix_update && s.newmark.newmark_delta === 0.5
    @test s.newmark.newmark_alpha === nothing
end

@testset "solver rules" begin
    s, ctx = ut_solver(Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0))
    @test ctx.errors[1].path == "Solver"
    @test ctx.errors[1].message ==
          "one solver is required: Verlet, Static, Linear Static Matrix Based, Verlet Matrix Based or Newmark"
    s, ctx = ut_solver(Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                        "Verlet" => Dict{String,Any}(),
                                        "Static" => Dict{String,Any}()))
    @test ctx.errors[1].message == "only one solver may be given, found: Verlet, Static"
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
    @test s.verlet_matrix_based.model_reduction === nothing
    s, ctx = ut_solver(Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                        "Newmark" => Dict{String,Any}("Model Reduction" => false)))
    @test isempty(ctx.errors)
    @test s.newmark.model_reduction === false
end
