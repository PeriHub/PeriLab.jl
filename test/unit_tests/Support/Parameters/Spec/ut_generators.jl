# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const GPS = PeriLab.ParameterSpec

@enum UTGenMode UTGenFast UTGenSlow

GPS.@params struct UTGenInner
    depth::Int64 = req("Depth"; min = 1)
end

GPS.@params struct UTGenParams
    stiffness::Float64 = req("Stiffness"; min = 0, quantity = :stress,
                             description = "elastic stiffness")
    table::GPS.Dependent = opt("Table Value"; default = 1.0)
    alpha::Union{Float64,Vector{Float64}} = opt("Alpha"; default = 0.5, min = 0)
    kind::String = opt("Kind"; default = "a", allowed = ["a", "b"])
    mode::UTGenMode = opt("Mode"; default = UTGenFast)
    label::Union{Nothing,String} = opt("Label"; default = nothing)
    id_or_name::Union{Int64,String} = opt("Id"; default = 1)
    inner::Union{Nothing,UTGenInner} = opt("Inner"; default = nothing)
    named::Dict{String,UTGenInner} = opt("Named"; default = Dict{String,UTGenInner}())
end

@testset "type schemas" begin
    @test GPS.type_schema(Float64) == Dict("type" => "number")
    @test GPS.type_schema(Union{Nothing,Int64}) == Dict("type" => "integer")
    @test GPS.type_schema(Union{Int64,String})["type"] == ["integer", "string"]
    dep = GPS.type_schema(GPS.Dependent)["oneOf"]
    @test Dict("type" => "number") in dep
    @test any(s -> get(s, "type", "") == "string", dep)
    nl = GPS.type_schema(Union{Float64,Vector{Float64}})["oneOf"]
    @test Dict("type" => "number") in nl
    @test Dict("type" => "array", "items" => Dict("type" => "number")) in nl
    @test GPS.type_schema(Vector{Int64}) ==
          Dict("type" => "array", "items" => Dict("type" => "integer"))
end

@testset "struct schema" begin
    s = GPS.json_schema(UTGenParams)
    p = s["properties"]
    @test s["type"] == "object" && s["additionalProperties"] == false
    @test s["required"] == ["Stiffness"]
    @test p["Stiffness"]["minimum"] == 0.0
    @test occursin("elastic stiffness", p["Stiffness"]["description"])
    @test occursin("stress", p["Stiffness"]["description"])
    @test p["Kind"]["enum"] == ["a", "b"] && p["Kind"]["default"] == "a"
    @test all(b -> get(b, "minimum", 0.0) == 0.0 &&
                   get(get(b, "items", Dict()), "minimum", 0.0) == 0.0, p["Alpha"]["oneOf"])
    @test p["Inner"]["properties"]["Depth"]["minimum"] == 1.0
    @test p["Named"]["additionalProperties"]["required"] == ["Depth"]
    @test p["Globals"] == Dict("type" => "object")
    @test !haskey(p["Label"], "default")                 # `nothing` has no JSON default
end

@testset "enum fields are not strict in the schema" begin
    m = GPS.json_schema(UTGenParams)["properties"]["Mode"]
    @test m["type"] == "string"
    @test !haskey(m, "enum")
    @test occursin("UTGenFast", m["description"]) && occursin("UTGenSlow", m["description"])
    @test m["default"] == "UTGenFast"
end

@testset "parameter rows and markdown" begin
    rows = GPS.parameter_rows(UTGenParams)
    stiff = only(filter(r -> r.key == "Stiffness", rows))
    @test stiff.type == "number" && stiff.required == "required" && stiff.range == "≥ 0"
    @test stiff.quantity == "stress" && stiff.description == "elastic stiffness"
    @test only(filter(r -> r.key == "Table Value", rows)).type == "number or data file"
    @test only(filter(r -> r.key == "Kind", rows)).range == "one of: a, b"
    @test only(filter(r -> r.key == "Label", rows)).default == "—"
    @test first.(GPS.nested_params(UTGenParams)) == ["Inner", "Named"]
    md = GPS.markdown_table(rows)
    @test startswith(md, "| YAML key | Type | Required | Default | Range | Quantity | Description |")
    @test occursin("| Stiffness | number | required |", md)
end

@testset "registered models" begin
    names = first.(GPS.registered_models(:material))
    @test "Correspondence Plastic" in names
    @test issorted(names)
    @test all(T -> GPS.is_params(T), last.(GPS.registered_models(:material)))
end
