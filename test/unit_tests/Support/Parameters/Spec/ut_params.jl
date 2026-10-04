# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const PS = PeriLab.ParameterSpec

@enum UTMode FastMode SafeMode

PS.@params struct UTInner
    radius::Float64 = req("Radius"; min = 0)
    allow_contact::Bool = opt("Allow Contact"; default = false)
end

"""Docstring for UTOuter."""
PS.@params struct UTOuter
    horizon::Float64 = req("Horizon"; min = 0, quantity = :length,
                           description = "Neighborhood radius")
    steps::Int64 = opt("Number of Steps"; default = 10, min = 1)
    mode::UTMode = opt("Mode"; default = SafeMode)
    note::Union{Nothing,String} = opt("Note"; default = nothing)
    weights::Vector{Float64} = opt("Weights"; default = [1.0, 2.0])
    filter::UTInner = req("Filter")
    sets::Dict{String,UTInner} = opt("Sets"; default = Dict{String,UTInner}())
end

PS.@params struct UTDependentMat
    youngs_modulus::Dependent = req("Young's Modulus"; min = 0)
    poissons_ratio::Float64 = req("Poisson's Ratio"; min = -1, max = 0.5)
end

PS.@params struct UTDerived
    a::Float64 = req("A")
    twice_a::Float64 = opt("Twice A"; default = 0.0)
end
PS.derive(p::UTDerived) = UTDerived(p.a, 2 * p.a)

ut_errors_by_path(ctx) = Dict(e.path => e.message for e in ctx.errors)

@testset "@params generates a plain struct and its spec" begin
    @test PS.is_params(UTOuter)
    @test !PS.is_params(Float64)
    @test fieldnames(UTOuter) ==
          (:horizon, :steps, :mode, :note, :weights, :filter, :sets)
    @test isconcretetype(UTOuter)
    spec = PS.parameter_spec(UTOuter)
    @test [fs.alias for fs in spec] ==
          ["Horizon", "Number of Steps", "Mode", "Note", "Weights", "Filter", "Sets"]
    @test spec[1].required && spec[1].min === 0.0 && spec[1].quantity === :length
    @test spec[1].description == "Neighborhood radius"
    @test spec[2].default === 10
    @test spec[3].default === SafeMode
    @test occursin("Docstring for UTOuter", string(@doc UTOuter))
end

@testset "Dependent fields make the struct parametric" begin
    @test !isconcretetype(UTDependentMat)
    @test fieldtype(UTDependentMat{PS.Constant}, :youngs_modulus) === PS.Constant
    @test PS.parameter_spec(UTDependentMat)[1].type === PS.Dependent
end

@testset "build: valid input" begin
    dict = Dict{String,Any}("Horizon" => 1,
                            "Filter" => Dict{String,Any}("Radius" => 0.5),
                            "Mode" => "Fast Mode",
                            "Sets" => Dict{String,Any}("left" => Dict{String,Any}("Radius" => 1.0,
                                                                                  "Allow Contact" => true)))
    ctx = PS.ParseContext()
    p = PS.parse_section(UTOuter, dict, "Disc", ctx)
    @test isempty(ctx.errors)
    @test p isa UTOuter
    @test p.horizon === 1.0
    @test p.steps === 10
    @test p.mode === FastMode
    @test p.note === nothing
    @test p.weights == [1.0, 2.0]
    @test p.filter == UTInner(0.5, false)
    @test p.sets["left"] == UTInner(1.0, true)
    p2 = PS.parse_section(UTOuter, dict, "Disc", ctx)
    @test p2.weights !== p.weights        # defaults are not shared between instances
end

@testset "build: collects all errors" begin
    dict = Dict{String,Any}("Horizon" => -0.1,
                            "Number of Steps" => 0,
                            "Filter" => Dict{String,Any}("Radius" => "big", "Radus" => 1),
                            "Poissons Ratio" => 0.3,
                            "Globals" => Dict{String,Any}("anything" => 1))
    ctx = PS.ParseContext()
    @test PS.parse_section(UTOuter, dict, "Disc", ctx) === nothing
    msgs = ut_errors_by_path(ctx)
    @test msgs["Disc.Horizon"] == "-0.1 is below minimum 0"
    @test msgs["Disc.\"Number of Steps\""] == "0 is below minimum 1"
    @test msgs["Disc.Filter.Radius"] == "expected a number, got \"big\""
    @test msgs["Disc.Filter.Radus"] == "unknown key — did you mean \"Radius\"?"
    @test msgs["Disc.\"Poissons Ratio\""] == "unknown key"
    @test length(ctx.errors) == 5         # "Globals" is never reported
end

@testset "missing and empty values" begin
    ctx = PS.ParseContext()
    @test PS.parse_section(UTOuter, Dict{String,Any}("Horizon" => nothing), "Disc", ctx) ===
          nothing
    msgs = ut_errors_by_path(ctx)
    @test msgs["Disc.Horizon"] == "expected a number, got an empty value"
    @test msgs["Disc.Filter"] == "missing (required by UTOuter)"
end

@testset "Int field accepts integral float" begin
    ctx = PS.ParseContext()
    p = PS.parse_section(UTOuter,
                         Dict{String,Any}("Horizon" => 1.0, "Number of Steps" => 100.0,
                                          "Filter" => Dict{String,Any}("Radius" => 1)),
                         "Disc", ctx)
    @test isempty(ctx.errors)
    @test p.steps === 100
    @test p.filter.radius === 1.0
end

@testset "case variants are not accepted silently" begin
    ctx = PS.ParseContext()
    @test PS.parse_section(UTInner, Dict{String,Any}("radius" => 1.0), "F", ctx) === nothing
    msgs = ut_errors_by_path(ctx)
    @test msgs["F.radius"] == "unknown key — did you mean \"Radius\"?"
    @test msgs["F.Radius"] == "missing (required by UTInner)"
end

@testset "non-strict mode downgrades unknown keys to warnings" begin
    ctx = PS.ParseContext(strict = false)
    p = PS.parse_section(UTInner, Dict{String,Any}("Radius" => 1.0, "Colour" => "red"), "F",
                         ctx)
    @test p == UTInner(1.0, false)
    @test !PS.has_errors(ctx)
    @test ctx.errors[1].severity == :warning
    @test ctx.errors[1].path == "F.Colour"
end

@testset "derive runs after building" begin
    ctx = PS.ParseContext()
    p = PS.parse_section(UTDerived, Dict{String,Any}("A" => 2.0), "D", ctx)
    @test p.twice_a == 4.0
end

@testset "Dependent fields via parse_section" begin
    dir = mktempdir()
    write(joinpath(dir, "E.txt"), "header: Temperature Young's_Modulus\n0 200\n100 180\n")
    ctx = PS.ParseContext(directory = dir)
    pc = PS.parse_section(UTDependentMat,
                          Dict{String,Any}("Young's Modulus" => 210.0,
                                           "Poisson's Ratio" => 0.3), "M", ctx)
    @test pc isa UTDependentMat{PS.Constant}
    @test isconcretetype(typeof(pc))
    pt = PS.parse_section(UTDependentMat,
                          Dict{String,Any}("Young's Modulus" => "E.txt",
                                           "Poisson's Ratio" => 0.3), "M", ctx)
    @test pt isa UTDependentMat{PS.Table1D}
    @test isempty(ctx.errors)
    @test PS.parse_section(UTInner, Dict{String,Any}("Radius" => "E.txt"), "F", ctx) ===
          nothing
    @test occursin("does not support dependent values", ctx.errors[end].message)
end

function ut_definition_error(ex)
    try
        Core.eval(@__MODULE__, ex)
    catch e
        return e
    end
    return nothing
end

@testset "definition-time errors" begin
    cases = [(:(PS.@params struct UTBad1
                    x::Float64
                end),
              "UTBad1: field `x` must be written as `x::Type = req(\"YAML key\"; ...)` or `x::Type = opt(\"YAML key\"; default = ...)`"),
             (:(PS.@params struct UTBad2
                    x::Float64 = 3.0
                end),
              "UTBad2: field `x` must be written as"),
             (:(PS.@params struct UTBad3
                    x::Matrix{Float64} = req("X")
                end),
              "UTBad3.x: unsupported field type Matrix{Float64}"),
             (:(PS.@params struct UTBad4
                    x::Float64 = opt("X")
                end),
              "UTBad4.x: opt(\"X\") needs a default"),
             (:(PS.@params struct UTBad5
                    x::Float64 = opt("X"; default = -1.0, min = 0)
                end),
              "UTBad5.x: default -1.0 is invalid: -1.0 is below minimum 0"),
             (:(PS.@params struct UTBad6
                    m::UTMode = opt("Mode"; default = "Turbo")
                end),
              "UTBad6.m: default \"Turbo\" is invalid: \"Turbo\" is not one of: FastMode, SafeMode"),
             (:(PS.@params struct UTBad7
                    a::Float64 = req("X")
                    b::Float64 = req("X")
                end),
              "UTBad7.b: YAML key \"X\" is already used by field `a`"),
             (:(PS.@params struct UTBad8
                    s::String = req("S"; min = 0)
                end),
              "UTBad8.s: min/max are only allowed on numeric fields"),
             (:(PS.@params mutable struct UTBad9
                    x::Float64 = req("X")
                end),
              "@params structs must be immutable"),
             (:(PS.@params struct UTBad10
                    inner::UTDependentMat = req("Inner")
                end),
              "UTBad10.inner: nested section type UTDependentMat contains Dependent fields; this is not supported"),
             (:(PS.@params x = 1),
              "@params must be applied to a struct definition")]
    for (ex, expected) in cases
        e = ut_definition_error(ex)
        @test e isa PS.ParamsDefinitionError
        @test e !== nothing && occursin(expected, e.msg)
    end
end

PS.@params struct UTOptionalRadius
    radius::Float64 = opt("Radius"; default = 1.0)
end

@testset "case-only key mismatch is an error even in non-strict mode" begin
    ctx = PS.ParseContext(strict = false)
    PS.parse_section(UTOptionalRadius, Dict{String,Any}("radius" => 5.0), "F", ctx)
    @test PS.has_errors(ctx)
    @test ctx.errors[1].path == "F.radius"
    @test ctx.errors[1].message == "unknown key — did you mean \"Radius\"?"
    ctx2 = PS.ParseContext(strict = false)
    PS.parse_section(UTOptionalRadius, Dict{String,Any}("Colour" => 5.0), "F", ctx2)
    @test !PS.has_errors(ctx2)            # genuinely unknown keys stay warnings
end

@testset "optional Dependent fields" begin
    ex = :(PS.@params struct UTOptionalDependent
               youngs_modulus_x::Union{Nothing,Dependent} = opt("Young's Modulus X";
                                                                default = nothing, min = 0)
               poissons_ratio::Float64 = req("Poisson's Ratio")
           end)
    @test ut_definition_error(ex) === nothing
    T = getfield(@__MODULE__, :UTOptionalDependent)
    @test PS.parameter_spec(T)[1].type === Union{Nothing,PS.Dependent}
    dir = mktempdir()
    write(joinpath(dir, "Ex.txt"), "header: Temperature Young's_Modulus_X\n0 200\n100 180\n")
    ctx = PS.ParseContext(directory = dir)
    absent = PS.parse_section(T, Dict{String,Any}("Poisson's Ratio" => 0.3), "M", ctx)
    constant = PS.parse_section(T, Dict{String,Any}("Young's Modulus X" => 210.0,
                                                    "Poisson's Ratio" => 0.3), "M", ctx)
    table = PS.parse_section(T, Dict{String,Any}("Young's Modulus X" => "Ex.txt",
                                                 "Poisson's Ratio" => 0.3), "M", ctx)
    @test isempty(ctx.errors)
    @test absent.youngs_modulus_x === nothing && isconcretetype(typeof(absent))
    @test constant.youngs_modulus_x isa PS.Constant && isconcretetype(typeof(constant))
    @test table.youngs_modulus_x isa PS.Table1D && isconcretetype(typeof(table))
end

@testset "optional Dependent from a file: missing column and min" begin
    T = getfield(@__MODULE__, :UTOptionalDependent)
    dir = mktempdir()
    write(joinpath(dir, "bad.txt"), "header: Temperature Other\n0 1\n1 2\n")
    write(joinpath(dir, "neg.txt"), "header: Temperature Young's_Modulus_X\n0 -1\n1 2\n")
    ctx = PS.ParseContext(directory = dir)
    PS.parse_section(T, Dict{String,Any}("Young's Modulus X" => "bad.txt",
                                         "Poisson's Ratio" => 0.3), "M", ctx)
    @test occursin("has no column \"Young's_Modulus_X\"", ctx.errors[end].message)
    PS.parse_section(T, Dict{String,Any}("Young's Modulus X" => "neg.txt",
                                         "Poisson's Ratio" => 0.3), "M", ctx)
    @test ctx.errors[end].message == "-1.0 is below minimum 0"
end
