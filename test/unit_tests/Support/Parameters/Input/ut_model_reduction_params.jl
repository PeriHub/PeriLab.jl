# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const PS = PeriLab.ParameterSpec
const ID = PeriLab.InputDeck

# Model_reduction.jl is included by Matrix_Verlet.jl
const UT_MR = PeriLab.Solver_Manager.Matrix_Verlet.Model_reduction

function ut_mr(dict)
    ctx = PS.ParseContext()
    p = PS.parse_section(ID.ModelReductionParams, dict, "MR", ctx)
    @test isempty(ctx.errors)
    return p
end

@testset "reduction blocks forms" begin
    @test UT_MR.parse_reduction_blocks(ut_mr(Dict{String,Any}("Type" => "Guyan",
                                                              "Reduction Blocks" => 2))) == [2]
    @test UT_MR.parse_reduction_blocks(ut_mr(Dict{String,Any}("Type" => "Guyan",
                                                              "Reduction Blocks" => "2, 3"))) ==
          [2, 3]
    @test UT_MR.parse_reduction_blocks(ut_mr(Dict{String,Any}("Type" => "Guyan"))) === nothing
end

@testset "Newmark model reduction is rejected clearly" begin
    ctx = PS.ParseContext()
    PS.parse_section(ID.SolverParams,
                     Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                      "Newmark" => Dict{String,Any}("Model Reduction" => true)),
                     "Solver", ctx)
    @test ctx.errors[1].path == "Solver.Newmark"
    @test ctx.errors[1].message ==
          "\"Model Reduction\" is not supported by the Newmark solver; use \"Verlet Matrix Based\" with a \"Model Reduction\" section"
    ctx = PS.ParseContext()
    PS.parse_section(ID.SolverParams,
                     Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                      "Newmark" => Dict{String,Any}("Model Reduction" => false)),
                     "Solver", ctx)
    @test isempty(ctx.errors)
end
