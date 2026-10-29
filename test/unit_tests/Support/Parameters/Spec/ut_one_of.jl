# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const PS = PeriLab.ParameterSpec

abstract type UTDrive end

PS.@params struct UTCar <: UTDrive
    speed::Float64 = opt("Speed"; default = 1.0, min = 0)
end

PS.@params struct UTBoat <: UTDrive
    depth::Float64 = req("Depth")
end

PS.@params struct UTTrip
    distance::Float64 = req("Distance")
    drive::UTDrive = one_of(:ut_drive; description = "how to travel")
end

PS.@params struct UTNoDrive
    x::Float64 = opt("X"; default = 0.0)
end

function ut_trip(dict)
    ctx = PS.ParseContext()
    return PS.parse_section(UTTrip, dict, "Trip", ctx), ctx
end

function with_ut_drives(f)
    PS.register_one_of_type!(:ut_drive, UTDrive)
    PS.register_model!(:ut_drive, "Car", UTCar)
    PS.register_model!(:ut_drive, "Boat", UTBoat)
    try
        f()
    finally
        delete!(PS.REGISTRY, :ut_drive)
        delete!(PS.ONE_OF_TYPES, :ut_drive)
    end
end

ut_messages(ctx) = [(e.path, e.message) for e in ctx.errors]

@testset "one_of selects the registered sub-section" begin
    with_ut_drives() do
        t, ctx = ut_trip(Dict{String,Any}("Distance" => 3.0,
                                          "Car" => Dict{String,Any}("Speed" => 2.0)))
        @test isempty(ctx.errors)
        @test t.drive isa UTCar && t.drive.speed === 2.0 && t.distance === 3.0
        t, ctx = ut_trip(Dict{String,Any}("Distance" => 3.0,
                                          "Boat" => Dict{String,Any}("Depth" => 1.5)))
        @test isempty(ctx.errors) && t.drive isa UTBoat && t.drive.depth === 1.5
        @test PS.one_of_keys(UTTrip) == Set(["Boat", "Car"])
        @test PS.aliases(UTTrip) == Set(["Distance"])
    end
end

@testset "one_of: none, several, misspelled" begin
    with_ut_drives() do
        t, ctx = ut_trip(Dict{String,Any}("Distance" => 3.0))
        @test t === nothing
        @test ut_messages(ctx) == [("Trip", "one ut drive is required: Boat, Car")]
        t, ctx = ut_trip(Dict{String,Any}("Distance" => 3.0, "Car" => Dict{String,Any}(),
                                          "Boat" => Dict{String,Any}("Depth" => 1.0)))
        @test ut_messages(ctx) == [("Trip", "only one ut drive may be given, found: Boat, Car")]
        t, ctx = ut_trip(Dict{String,Any}("Distance" => 3.0, "Cart" => Dict{String,Any}()))
        @test ("Trip", "one ut drive is required: Boat, Car") in ut_messages(ctx)
        @test ("Trip.Cart", "unknown key — did you mean \"Car\"?") in ut_messages(ctx)
    end
end

@testset "case variant is an error" begin
    with_ut_drives() do
        t, ctx = ut_trip(Dict{String,Any}("Distance" => 3.0, "car" => Dict{String,Any}()))
        @test t === nothing
        @test ("Trip", "one ut drive is required: Boat, Car") in ut_messages(ctx)
        @test ("Trip.car", "unknown key — did you mean \"Car\"?") in ut_messages(ctx)
    end
end

@testset "empty sub-section means defaults" begin
    with_ut_drives() do
        t, ctx = ut_trip(Dict{String,Any}("Distance" => 3.0, "Car" => nothing))
        @test isempty(ctx.errors) && t.drive isa UTCar && t.drive.speed === 1.0
    end
end

@testset "non-section value" begin
    with_ut_drives() do
        t, ctx = ut_trip(Dict{String,Any}("Distance" => 3.0, "Car" => true))
        @test t === nothing
        @test ut_messages(ctx) ==
              [("Trip.Car", "expected a section of `key: value` entries, got true")]
    end
end

@testset "one_of: nested errors" begin
    with_ut_drives() do
        t, ctx = ut_trip(Dict{String,Any}("Distance" => 3.0, "Boat" => Dict{String,Any}()))
        @test ut_messages(ctx) == [("Trip.Boat.Depth", "missing (required by UTBoat)")]
        t, ctx = ut_trip(Dict{String,Any}("Distance" => 3.0,
                                          "Car" => Dict{String,Any}("Sped" => 1.0)))
        @test ut_messages(ctx) == [("Trip.Car.Sped", "unknown key — did you mean \"Speed\"?")]
    end
end

@testset "one_of: unavailable entry" begin
    with_ut_drives() do
        PS.register_unavailable!(:ut_drive, "Plane")
        t, ctx = ut_trip(Dict{String,Any}("Distance" => 3.0, "Plane" => Dict{String,Any}()))
        @test ut_messages(ctx) ==
              [("Trip.Plane", "ut drive \"Plane\" requires a license that is not available")]
    end
end

@testset "one_of: definition errors" begin
    e = try
        Core.eval(@__MODULE__,
                  :(PS.@params struct UTBadOneOf1
                        x::Float64 = one_of(:ut_drive)
                    end))
        nothing
    catch err
        err
    end
    @test e isa PS.ParamsDefinitionError
    @test occursin("UTBadOneOf1.x: one_of needs an abstract type, got Float64", e.msg)
    e = try
        Core.eval(@__MODULE__,
                  :(PS.@params struct UTBadOneOf2
                        a::UTDrive = one_of(:ut_drive)
                        b::UTDrive = one_of(:ut_drive)
                    end))
        nothing
    catch err
        err
    end
    @test e isa PS.ParamsDefinitionError
    @test occursin("UTBadOneOf2.b: only one one_of field per struct", e.msg)
    with_ut_drives() do
        e = try
            PS.register_model!(:ut_drive, "Nothing", UTNoDrive)
            nothing
        catch err
            err
        end
        @test e isa PS.ParamsDefinitionError
        @test e.msg == "ut_drive model \"Nothing\": UTNoDrive must be a subtype of UTDrive"
    end
    @test_throws PS.ParamsDefinitionError PS.register_one_of_type!(:ut_drive, UTCar)
end

@testset "one_of: schema, rows, nested sections" begin
    with_ut_drives() do
        s = PS.json_schema(UTTrip)
        @test s["required"] == ["Distance"]
        # an empty sub-section (null) means the model's defaults
        @test s["properties"]["Car"] == merge(PS.json_schema(UTCar),
                                              Dict("type" => ["object", "null"]))
        @test s["properties"]["Boat"] == merge(PS.json_schema(UTBoat),
                                               Dict("type" => ["object", "null"]))
        @test s["oneOf"] == Any[Dict{String,Any}("required" => ["Boat"]),
                                Dict{String,Any}("required" => ["Car"])]
        @test !haskey(s["properties"], "")
        rows = PS.parameter_rows(UTTrip)
        @test rows[2].key == "<ut drive>" && rows[2].type == "one of: Boat, Car"
        @test rows[2].required == "required" && rows[2].description == "how to travel"
        @test PS.nested_params(UTTrip) == ["Boat" => UTBoat, "Car" => UTCar]
        @test occursin("| <ut drive> | one of: Boat, Car | required |",
                       PS.markdown_table(rows))
    end
end

@testset "one_of: template" begin
    with_ut_drives() do
        io = IOBuffer()
        PeriLab.InputDeck._template(io, UTTrip, "  ")
        text = String(take!(io))
        @test occursin("  # exactly one of: Boat, Car\n  Boat:\n    Depth: <number>", text)
    end
end
