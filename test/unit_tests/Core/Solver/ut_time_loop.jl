# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test

module UTLoopSolver
import PeriLab
using PeriLab.ParameterSpec: @params
using PeriLab.InputDeck: SolverOptions, SolverParams
using PeriLab.Solver_Manager: TimeWindow, SolverRun, SolverContext, StepClock, StepFailed
import PeriLab.Solver_Manager: prepare_solver, time_step!
const DM = PeriLab.Data_Manager

@params struct UTLoopParams <: SolverOptions
    steps::Int64 = opt("Steps"; default = 4, min = 1)
    fail_at::Int64 = opt("Fail At"; default = 0)
    cancel_at::Int64 = opt("Cancel At"; default = 0)
    damage_at::Int64 = opt("Damage At"; default = 0)
end

struct UTLoopWorkspace
    clocks::Vector{StepClock}
end

const WORKSPACES = UTLoopWorkspace[]

function init_solver(options::UTLoopParams, params::SolverParams, models, bcs, block_nodes)
    t0 = something(params.initial_time, 0.0)
    return TimeWindow(t0, 1.0, (1.0 - t0) / options.steps, options.steps), nothing
end

function prepare_solver(run::SolverRun{UTLoopParams}, ctx::SolverContext)
    ws = UTLoopWorkspace(StepClock[])
    push!(WORKSPACES, ws)
    return ws
end

function time_step!(ws::UTLoopWorkspace, run::SolverRun{UTLoopParams}, ctx::SolverContext,
                    clock::StepClock)
    push!(ws.clocks, clock)
    o = run.params.options
    DM.get_field("Displacements", "NP1") .= clock.idt
    clock.idt == o.damage_at && (DM.get_damage("NP1")[1] = 0.5)
    clock.idt == o.cancel_at && DM.set_cancel(true)
    clock.idt == o.fail_at && return StepFailed("UT failure")
    return nothing
end
end

const UTL_SM = PeriLab.Solver_Manager
const UTL_BC = UTL_SM.Boundary_Conditions.BoundaryCondition
const UTL_DM = PeriLab.Data_Manager

function utl_data()
    UTL_DM.initialize_data()
    UTL_DM.set_num_controller(3)
    UTL_DM.set_dof(2)
    UTL_DM.create_constant_node_scalar_field("Active", Bool; default_value = true)
    UTL_DM.create_node_vector_field("Displacements", Float64, 2)
    UTL_DM.create_node_scalar_field("Damage", Float64)
    UTL_DM.set_cancel(false)
    empty!(UTLoopSolver.WORKSPACES)
end

function utl_run(solver::Dict{String,Any}; models = ["Material"], extra = Dict{String,Any}())
    ctx = PeriLab.ParameterSpec.ParseContext()
    deck = merge(Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                  "UT Loop" => solver), extra)
    params = PeriLab.ParameterSpec.parse_section(PeriLab.InputDeck.SolverParams, deck,
                                                 "Solver", ctx)
    @test isempty(ctx.errors)
    return UTL_SM.create_run(params, models, models,
                             PeriLab.InputDeck.calculation_options(params),
                             Dict{String,UTL_BC}(), Dict{Int64,Vector{Int64}}())
end

function utl_ctx(writes; set_output_frequency = nsteps -> nothing)
    write = (rf, t, max_damage, outputs) -> (push!(writes, (t, max_damage)); rf)
    return UTL_SM.SolverContext(Dict{Int64,Vector{Int64}}(), Dict{String,UTL_BC}(),
                                Dict{Int64,Dict{}}(), Dict[], write, true,
                                set_output_frequency)
end

function utl_with_solver(f)
    PeriLab.ParameterSpec.register_solver("UT Loop", UTLoopSolver.UTLoopParams)
    try
        utl_data()
        f()
    finally
        delete!(PeriLab.ParameterSpec.REGISTRY[:solver], "UT Loop")
    end
end

@testset "time_loop runs nsteps + 1 steps" begin
    utl_with_solver() do
        writes = Tuple{Float64,Float64}[]
        result = UTL_SM.solver(utl_run(Dict{String,Any}("Steps" => 4)), utl_ctx(writes))
        @test result == Dict[]
        clocks = only(UTLoopSolver.WORKSPACES).clocks
        @test [c.idt for c in clocks] == 1:5
        @test [c.time for c in clocks] ≈ [0.0, 0.25, 0.5, 0.75, 1.0]
        @test [c.step_time for c in clocks] ≈ [0.0, 0.25, 0.5, 0.75, 1.0]
        @test all(c -> c.dt == 0.25, clocks)
        @test first.(writes) ≈ [0.0, 0.25, 0.5, 0.75, 1.0]
        @test all(w -> w[2] == 0.0, writes)
        @test UTL_DM.get_current_time() == 1.0
        @test !UTL_DM.get_cancel()
    end
end

@testset "time_loop stops on cancel after writing the step" begin
    utl_with_solver() do
        writes = Tuple{Float64,Float64}[]
        UTL_SM.solver(utl_run(Dict{String,Any}("Cancel At" => 2)), utl_ctx(writes))
        @test length(only(UTLoopSolver.WORKSPACES).clocks) == 2
        @test first.(writes) ≈ [0.0, 0.25]
        @test UTL_DM.get_current_time() == 1.0
    end
end

@testset "time_loop passes max damage and cancels at the damage limit" begin
    utl_with_solver() do
        writes = Tuple{Float64,Float64}[]
        run = utl_run(Dict{String,Any}("Damage At" => 3); models = ["Material", "Damage"],
                      extra = Dict{String,Any}("Maximum Damage" => 0.4))
        UTL_SM.solver(run, utl_ctx(writes))
        @test last.(writes) == [0.0, 0.0, 0.5]
        @test UTL_DM.get_cancel()
    end
end

@testset "time_loop ignores damage without damage models" begin
    utl_with_solver() do
        writes = Tuple{Float64,Float64}[]
        UTL_SM.solver(utl_run(Dict{String,Any}("Damage At" => 2)), utl_ctx(writes))
        @test all(w -> w[2] == 0.0, writes) && length(writes) == 5
    end
end

@testset "a failed step without fallback is discarded and cancels" begin
    utl_with_solver() do
        writes = Tuple{Float64,Float64}[]
        UTL_SM.solver(utl_run(Dict{String,Any}("Fail At" => 3)), utl_ctx(writes))
        @test length(only(UTLoopSolver.WORKSPACES).clocks) == 3
        @test first.(writes) ≈ [0.0, 0.25]
        @test all(UTL_DM.get_field("Displacements", "NP1") .== 0)
        @test all(UTL_DM.get_field("Displacements", "N") .== 2)
        @test UTL_DM.get_cancel()
    end
end

@testset "reset_NP1 restores the state after switch_NP1_to_N" begin
    UTL_DM.initialize_data()
    UTL_DM.set_num_controller(2)
    UTL_DM.set_dof(2)
    UTL_DM.create_constant_node_scalar_field("Active", Bool; default_value = true)
    UTL_DM.create_node_scalar_field("A", Float64)
    UTL_DM.get_field("A", "NP1") .= 3.0
    UTL_DM.switch_NP1_to_N()
    UTL_DM.get_field("A", "NP1") .= 7.0
    UTL_DM.reset_NP1()
    @test UTL_DM.get_field("A", "NP1") == [0.0, 0.0]
    @test UTL_DM.get_field("A", "N") == [3.0, 3.0]
end

function utl_fallback(steps; fail_at = 0)
    return Dict{String,Any}("Fallback" => Dict{String,Any}("UT Loop" => Dict{String,Any}("Steps" => steps,
                                                                                        "Fail At" => fail_at)))
end

@testset "a failed step continues with the fallback solver" begin
    utl_with_solver() do
        writes = Tuple{Float64,Float64}[]
        frequencies = Int64[]
        run = utl_run(Dict{String,Any}("Fail At" => 3); extra = utl_fallback(2))
        UTL_SM.solver(run, utl_ctx(writes; set_output_frequency = n -> push!(frequencies, n)))
        first_run, fallback = UTLoopSolver.WORKSPACES
        @test length(first_run.clocks) == 3
        @test [c.time for c in fallback.clocks] ≈ [0.5, 0.75, 1.0]
        @test fallback.clocks[1].step_time == 0.0
        @test first.(writes) ≈ [0.0, 0.25, 0.5, 0.75, 1.0]
        @test frequencies == [2]
        @test !UTL_DM.get_cancel()
        @test UTL_DM.get_current_time() == 1.0
    end
end

@testset "fallback after a failure in the first step" begin
    utl_with_solver() do
        writes = Tuple{Float64,Float64}[]
        UTL_SM.solver(utl_run(Dict{String,Any}("Fail At" => 1); extra = utl_fallback(2)),
                      utl_ctx(writes))
        @test [c.time for c in UTLoopSolver.WORKSPACES[2].clocks] ≈ [0.0, 0.5, 1.0]
        @test first.(writes) ≈ [0.0, 0.5, 1.0]
    end
end

@testset "a failing fallback cancels" begin
    utl_with_solver() do
        writes = Tuple{Float64,Float64}[]
        UTL_SM.solver(utl_run(Dict{String,Any}("Fail At" => 2); extra = utl_fallback(4; fail_at = 1)),
                      utl_ctx(writes))
        @test length(UTLoopSolver.WORKSPACES) == 2
        @test first.(writes) ≈ [0.0]
        @test UTL_DM.get_cancel()
    end
end

@testset "a dry run with a fallback does not crash" begin
    utl_with_solver() do
        writes = Tuple{Float64,Float64}[]
        run = UTL_SM.with_steps(utl_run(Dict{String,Any}("Fail At" => 2, "Steps" => 4);
                                        extra = utl_fallback(2)), 3)
        UTL_SM.solver(run, utl_ctx(writes))
        @test length(UTLoopSolver.WORKSPACES) == 2
        @test !UTL_DM.get_cancel()
    end
end
