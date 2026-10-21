# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const PS = PeriLab.ParameterSpec

PS.@params struct UTUnions
    degree::Union{Int64,String} = req("Degree")
    value::Union{Float64,String} = req("Value")
    step::Union{Nothing,Int64,String} = opt("Step ID"; default = nothing)
    sets::Dict{String,Union{Int64,String}} = opt("Node Sets";
                                                 default = Dict{String,Union{Int64,String}}())
    flags::Dict{String,Bool} = opt("Flags"; default = Dict{String,Bool}())
end

PS.@params struct UTChecked
    low::Float64 = req("Low")
    high::Float64 = req("High")
end
function PS.check!(p::UTChecked, path::String, ctx::PS.ParseContext)
    p.low <= p.high || PS.add_error!(ctx, path, "\"Low\" must not exceed \"High\"")
    return nothing
end

PS.@params struct UTCheckedHolder
    items::Dict{String,UTChecked} = req("Items")
end

function ut_ext_conv(T, raw)
    ctx = PS.ParseContext()
    return PS.convert_value(T, raw, "k", ctx), ctx
end

@testset "scalar unions are supported" begin
    @test PS.supported_type(Union{Int64,String})
    @test PS.supported_type(Union{Nothing,Int64,String})
    @test PS.supported_type(Union{Float64,String})
    @test !PS.supported_type(Union{Int64,Float64})
    @test PS.supported_type(Dict{String,Union{Int64,String}})
    @test PS.supported_type(Dict{String,Bool})
end

@testset "scalar union conversion" begin
    @test ut_ext_conv(Union{Int64,String}, 5)[1] === 5
    @test ut_ext_conv(Union{Int64,String}, "1 1")[1] == "1 1"
    @test ut_ext_conv(Union{Int64,String}, 2.0)[1] === 2
    @test ut_ext_conv(Union{Float64,String}, -25)[1] === -25.0
    @test ut_ext_conv(Union{Float64,String}, "0*t")[1] == "0*t"
    @test ut_ext_conv(Union{Nothing,Int64,String}, nothing)[1] === nothing
    v, ctx = ut_ext_conv(Union{Int64,String}, 2.5)
    @test v === PS.FAILED
    @test ctx.errors[1].message == "expected an integer or text, got 2.5"
    v, ctx = ut_ext_conv(Union{Float64,String}, true)
    @test ctx.errors[1].message == "expected a number or text, got true"
end

@testset "dicts of scalars, Globals skipped in named sections" begin
    ctx = PS.ParseContext()
    p = PS.parse_section(UTUnions,
                         Dict{String,Any}("Degree" => "1 1", "Value" => 0,
                                          "Node Sets" => Dict{String,Any}("a" => 5,
                                                                          "b" => "2 3 4",
                                                                          "Globals" => Dict{String,Any}()),
                                          "Flags" => Dict{String,Any}("x" => true)),
                         "S", ctx)
    @test isempty(ctx.errors)
    @test p.degree == "1 1" && p.value === 0.0 && p.step === nothing
    @test p.sets == Dict{String,Union{Int64,String}}("a" => 5, "b" => "2 3 4")
    @test p.flags == Dict("x" => true)
    ctx = PS.ParseContext()
    PS.parse_section(UTUnions,
                     Dict{String,Any}("Degree" => 1, "Value" => 1,
                                      "Flags" => Dict{String,Any}("x" => "yes")), "S", ctx)
    @test ctx.errors[1].path == "S.Flags.x"
    @test ctx.errors[1].message == "expected true or false, got \"yes\""
end

@testset "check! runs after building, also for nested sections" begin
    ctx = PS.ParseContext()
    @test PS.parse_section(UTChecked, Dict{String,Any}("Low" => 2.0, "High" => 1.0), "C",
                           ctx) isa UTChecked
    @test ctx.errors[1].path == "C"
    @test ctx.errors[1].message == "\"Low\" must not exceed \"High\""
    ctx = PS.ParseContext()
    PS.parse_section(UTCheckedHolder,
                     Dict{String,Any}("Items" => Dict{String,Any}("a" => Dict{String,Any}("Low" => 3.0,
                                                                                         "High" => 1.0))),
                     "H", ctx)
    @test ctx.errors[1].path == "H.Items.a"
    ctx = PS.ParseContext()
    PS.parse_section(UTChecked, Dict{String,Any}("Low" => 1.0), "C", ctx)
    @test length(ctx.errors) == 1           # check! is not run on a struct that failed to build
end

abstract type UTAbstractOptions end

PS.@params struct UTConcreteOptions <: UTAbstractOptions
    factor::Float64 = opt("Factor"; default = 1.0)
end

PS.@params struct UTDependentOptions <: UTAbstractOptions
    value::Dependent = req("Value")
end

@testset "@params accepts a supertype" begin
    @test UTConcreteOptions <: UTAbstractOptions
    @test UTDependentOptions{PS.Constant} <: UTAbstractOptions
    ctx = PS.ParseContext()
    p = PS.parse_section(UTConcreteOptions, Dict{String,Any}("Factor" => 2), "O", ctx)
    @test isempty(ctx.errors) && p.factor === 2.0
end
