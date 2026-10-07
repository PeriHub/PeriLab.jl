# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const PS = PeriLab.ParameterSpec
const ID = PeriLab.InputDeck
const PH = PeriLab.Parameter_Handling

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

const UT_ST_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", "..", ".."))
const UT_ST_SKIP = ("test/fullscale_tests/test_PD_Solid_Elastic/strain_xx_external.yaml",)

function ut_st_decks()
    decks = String[]
    for root in ("examples", "test"), (dir, _, files) in walkdir(joinpath(UT_ST_ROOT, root))
        append!(decks, [joinpath(dir, f) for f in files if endswith(f, ".yaml")])
    end
    return sort!(decks)
end

function ut_st_compare_all_decks(t)
    for file in ut_st_decks()
        relpath(file, UT_ST_ROOT) in UT_ST_SKIP && continue
        raw = PeriLab.IO.read_input(file)
        (raw isa AbstractDict && haskey(raw, "PeriLab")) || continue
        deck = raw["PeriLab"]
        input, ctx = ID.read_input(deck, dirname(file))
        deck_needs_license(ctx) && continue
        steps = PH.get_solver_steps(deck)
        @test ID.solver_steps(input) == steps
        for step in steps
            d = step == -1 ? deck["Solver"] : PH.get_solver_params(deck, step)
            s = ID.solver_step(input, step)
            options = ID.active_options(s)
            calc = PH.get_calculation_options(d)
            @test ID.solver_name(s) == PH.get_solver_name(d)
            @test something(s.number_of_steps, 1) == PH.get_nsteps(d)
            @test options.safety_factor == PH.get_safety_factor(d)
            @test options.fixed_dt == PH.get_fixed_dt(d)
            @test options.numerical_damping == PH.get_numerical_damping(d)
            @test s.maximum_damage == PH.get_max_damage(d)
            @test ID.model_options(s) == PH.get_model_options(d)
            @test s.calculate_cauchy == calc["Calculate Cauchy"]
            @test s.calculate_von_mises_stress == calc["Calculate von Mises stress"]
            @test s.calculate_strain == calc["Calculate Strain"]
            @test ID.end_time(s, t) == PH.get_final_time(d)
            if haskey(d, "Initial Time") || t != 0.0     # otherwise both abort
                @test ID.start_time(s, t) == PH.get_initial_time(d)
            end
        end
    end
end

@testset "typed solver options equal the Dict getters on every shipped deck" begin
    # the Dict time getters read the current time from Data_Manager
    dm = PeriLab.Data_Manager
    previous = haskey(dm.data, "current_time") ? dm.get_current_time() : nothing
    t = 0.0
    dm.set_current_time(t)
    try
        ut_st_compare_all_decks(t)
    finally
        previous === nothing ? delete!(dm.data, "current_time") :
        dm.set_current_time(previous)
    end
end

@testset "Solver_Manager.init keeps its docstring" begin
    @test occursin("Initialize the solver", string(@doc PeriLab.Solver_Manager.init))
end
