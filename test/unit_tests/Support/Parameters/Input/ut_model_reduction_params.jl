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
    p = PS.parse_section(UT_MR.ModelReductionParams, dict, "MR", ctx)
    @test isempty(ctx.errors)
    return p
end

@testset "reduction blocks forms" begin
    @test UT_MR.parse_reduction_blocks(ut_mr(Dict{String,Any}("Type" => "Static Condensation",
                                                              "Reduction Blocks" => 2))) ==
          [2]
    @test UT_MR.parse_reduction_blocks(ut_mr(Dict{String,Any}("Type" => "Static Condensation",
                                                              "Reduction Blocks" => "2, 3"))) ==
          [2, 3]
    @test UT_MR.parse_reduction_blocks(ut_mr(Dict{String,Any}("Type" => "Static Condensation"))) ===
          nothing
end

@testset "Newmark has no model reduction" begin
    ctx = PS.ParseContext()
    PS.parse_section(ID.SolverParams,
                     Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                      "Newmark" =>
                                          Dict{String,Any}("Model Reduction" => true)),
                     "Solver", ctx)
    @test ctx.errors[1].path == "Solver.Newmark.\"Model Reduction\""
    @test startswith(ctx.errors[1].message, "unknown key")
end

@testset "Craig Bampton Cascade options" begin
    p = ut_mr(Dict{String,Any}("Type" => "Craig Bampton Cascade", "Reduction Blocks" => 2,
                               "Number of Modes" => 0, "Maximum Frequency" => 1000,
                               "Number of Subregions" => 4,
                               "Reduced Matrix Cache" => "cache.jls"))
    @test p.number_of_modes == 0
    @test p.max_frequency == 1000.0
    @test p.number_of_subregions == 4
    @test p.reduced_matrix_cache == "cache.jls"

    ctx = PS.ParseContext()
    PS.parse_section(UT_MR.ModelReductionParams,
                     Dict{String,Any}("Type" => "Static Condensation",
                                      "Maximum Frequency" => 1000), "MR", ctx)
    @test ctx.errors[1].path == "MR.\"Maximum Frequency\""
end
