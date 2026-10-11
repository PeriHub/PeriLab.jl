# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test

module UTPluginSolver
import PeriLab
using PeriLab.ParameterSpec: @params
using PeriLab.InputDeck: SolverOptions, SolverParams
using PeriLab.Solver_Manager: TimeWindow, SolverRun, SolverContext

@params struct UTPluginParams <: SolverOptions
    steps::Int64 = opt("Steps"; default = 3, min = 1)
end

const CALLS = String[]

function init_solver(options::UTPluginParams, params::SolverParams, models::Vector{String},
                     bcs, block_nodes)
    push!(CALLS, "init")
    return TimeWindow(0.0, 1.0, 1.0 / options.steps, options.steps), (; models = models)
end

function run_solver(run::SolverRun{UTPluginParams}, ctx::SolverContext)
    push!(CALLS, "run $(run.time.nsteps)")
    return ctx.result_files
end
end

const UT_SM = PeriLab.Solver_Manager
const UT_BC = UT_SM.Boundary_Conditions.BoundaryCondition

function ut_plugin_params(dict)
    ctx = PeriLab.ParameterSpec.ParseContext()
    p = PeriLab.ParameterSpec.parse_section(PeriLab.InputDeck.SolverParams, dict, "Solver", ctx)
    return p, ctx
end

@testset "a registered solver runs through Solver_Manager" begin
    PeriLab.ParameterSpec.register_solver("UT Plugin", UTPluginSolver.UTPluginParams)
    try
        params, ctx = ut_plugin_params(Dict{String,Any}("Initial Time" => 0.0,
                                                        "Final Time" => 1.0,
                                                        "Calculate Strain" => true,
                                                        "UT Plugin" => Dict{String,Any}("Steps" => 4)))
        @test isempty(ctx.errors)
        @test PeriLab.InputDeck.solver_name(params) == "UT Plugin"
        calculation = PeriLab.InputDeck.calculation_options(params)
        @test calculation == PeriLab.InputDeck.CalculationOptions(false, false, true)
        run = UT_SM.create_run(params, ["Material"], ["Material", "Thermal"], calculation,
                               Dict{String,UT_BC}(), Dict{Int64,Vector{Int64}}())
        @test run isa UT_SM.SolverRun{UTPluginSolver.UTPluginParams}
        @test PeriLab.InputDeck.solver_name(run) == "UT Plugin"
        @test run.time == UT_SM.TimeWindow(0.0, 1.0, 0.25, 4)
        @test run.state.models == ["Material"] && run.all_models == ["Material", "Thermal"]
        sctx = UT_SM.SolverContext(Dict{Int64,Vector{Int64}}(), Dict{String,UT_BC}(),
                                   Dict{Int64,Dict{}}(), Dict[], (args...) -> nothing, true)
        @test UT_SM.solver(run, sctx) == Dict[]
        @test UTPluginSolver.CALLS == ["init", "run 4"]
    finally
        delete!(PeriLab.ParameterSpec.REGISTRY[:solver], "UT Plugin")
    end
end

@testset "with_steps" begin
    PeriLab.ParameterSpec.register_solver("UT Plugin", UTPluginSolver.UTPluginParams)
    try
        params, _ = ut_plugin_params(Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                                      "UT Plugin" => Dict{String,Any}()))
        run = UT_SM.create_run(params, ["Material"], ["Material"],
                               PeriLab.InputDeck.calculation_options(params),
                               Dict{String,UT_BC}(), Dict{Int64,Vector{Int64}}())
        short = UT_SM.with_steps(run, 10)
        @test short.time == UT_SM.TimeWindow(0.0, 1.0, run.time.dt, 10)
        @test short.params === run.params && short.state === run.state
        @test short.models === run.models && short.calculation == run.calculation
        @test short isa UT_SM.SolverRun{UTPluginSolver.UTPluginParams}
    finally
        delete!(PeriLab.ParameterSpec.REGISTRY[:solver], "UT Plugin")
    end
end

@testset "calculation options merge over steps" begin
    CO = PeriLab.InputDeck.CalculationOptions
    @test (CO(true, false, false) | CO(false, false, true)) == CO(true, false, true)
end
