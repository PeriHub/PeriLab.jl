# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const PS = PeriLab.ParameterSpec

@enum UTHardening LinearHardening ExponentialHardening
@enum UTSymmetry PlaneStress PlaneStrain Full3D
PS.enum_aliases(::Type{UTSymmetry}) = Dict{String,UTSymmetry}("3D" => Full3D)

function ut_conv(T, raw; directory = "", alias = "")
    ctx = PS.ParseContext(directory = directory)
    v = PS.convert_value(T, raw, "k", ctx; alias = alias)
    return v, ctx
end

function ut_conv_error(T, raw)
    v, ctx = ut_conv(T, raw)
    @test v === PS.FAILED
    @test length(ctx.errors) == 1
    @test ctx.errors[1].path == "k"
    return ctx.errors[1].message
end

@testset "Float64" begin
    @test ut_conv(Float64, 2)[1] === 2.0
    @test ut_conv(Float64, 2.5)[1] === 2.5
    @test ut_conv_error(Float64, "abc") == "expected a number, got \"abc\""
    @test ut_conv_error(Float64, "E.txt") ==
          "expected a number, got file path \"E.txt\" — this parameter does not support dependent values"
    @test ut_conv_error(Float64, true) == "expected a number, got true"
    @test ut_conv_error(Float64, nothing) == "expected a number, got an empty value"
end

@testset "Int64" begin
    @test ut_conv(Int64, 3)[1] === 3
    @test ut_conv(Int64, 100.0)[1] === 100
    @test ut_conv_error(Int64, 2.5) == "expected an integer, got 2.5"
    @test ut_conv_error(Int64, "3") == "expected an integer, got \"3\""
end

@testset "Bool and String" begin
    @test ut_conv(Bool, true)[1] === true
    @test ut_conv_error(Bool, "yes") == "expected true or false, got \"yes\""
    @test ut_conv(String, "abc")[1] == "abc"
    @test ut_conv_error(String, 5) == "expected text, got 5"
end

@testset "Enum" begin
    @test ut_conv(UTHardening, "LinearHardening")[1] === LinearHardening
    @test ut_conv(UTHardening, "Linear Hardening")[1] === LinearHardening
    @test ut_conv(UTHardening, "exponential_hardening")[1] === ExponentialHardening
    @test ut_conv(UTHardening, ExponentialHardening)[1] === ExponentialHardening
    @test ut_conv(UTSymmetry, "3D")[1] === Full3D
    @test ut_conv(UTSymmetry, "plane stress")[1] === PlaneStress
    @test ut_conv_error(UTHardening, "Cubic") ==
          "\"Cubic\" is not one of: LinearHardening, ExponentialHardening"
    @test ut_conv_error(UTSymmetry, 3) ==
          "3 is not one of: PlaneStress, PlaneStrain, Full3D, 3D"
end

@testset "Vectors" begin
    v = ut_conv(Vector{Float64}, [1, 2.5])[1]
    @test v == [1.0, 2.5] && v isa Vector{Float64}
    @test ut_conv(Vector{Int64}, [1, 2])[1] == [1, 2]
    @test ut_conv(Vector{String}, ["a", "b"])[1] == ["a", "b"]
    v, ctx = ut_conv(Vector{Float64}, [1, "a"])
    @test v === PS.FAILED
    @test ctx.errors[1].path == "k[2]"
    @test ut_conv_error(Vector{Float64}, 3) == "expected a list, got 3"
end

@testset "Union{Nothing,T}" begin
    @test ut_conv(Union{Nothing,Float64}, nothing)[1] === nothing
    @test ut_conv(Union{Nothing,Float64}, 2)[1] === 2.0
    @test ut_conv_error(Union{Nothing,Float64}, "x") == "expected a number, got \"x\""
end

@testset "Dependent" begin
    @test ut_conv(PS.Dependent, 3)[1] === PS.Constant(3.0)
    dir = mktempdir()
    write(joinpath(dir, "E.txt"), "header: Temperature Young's_Modulus\n0 200\n100 180\n")
    t, ctx = ut_conv(PS.Dependent, "E.txt"; directory = dir, alias = "Young's Modulus")
    @test !PS.has_errors(ctx)
    @test t isa PS.Table1D
    @test t.y == [200.0, 180.0]
    v, ctx = ut_conv(PS.Dependent, "missing.txt"; directory = dir, alias = "Young's Modulus")
    @test v === PS.FAILED
    @test occursin("not found", ctx.errors[1].message)
    @test ut_conv_error(PS.Dependent, true) ==
          "expected a number or a data file path, got true"
end

@testset "req / opt / FieldSpec" begin
    d = PS.req("Horizon"; min = 0, quantity = :length, description = "radius")
    @test d.alias == "Horizon" && d.required && d.min === 0.0 && d.max === nothing
    @test d.quantity === :length && d.description == "radius"
    o = PS.opt("Steps"; default = 10, allowed = [10, 20])
    @test !o.required && o.default == 10 && o.allowed == Any[10, 20]
    @test PS.opt("X").default === PS.NO_DEFAULT
    fs = PS.FieldSpec(:steps, Int64, o)
    @test fs.name === :steps && fs.type === Int64 && fs.alias == "Steps" && fs.default == 10
end

@testset "check_constraints!" begin
    ctx = PS.ParseContext()
    fs = PS.FieldSpec(:h, Float64, PS.req("Horizon"; min = 0, max = 10))
    @test PS.check_constraints!(fs, 5.0, "h", ctx)
    @test !PS.check_constraints!(fs, -0.1, "h", ctx)
    @test ctx.errors[end].message == "-0.1 is below minimum 0"
    @test !PS.check_constraints!(fs, 11.0, "h", ctx)
    @test ctx.errors[end].message == "11.0 is above maximum 10"
    vec_fs = PS.FieldSpec(:v, Vector{Float64}, PS.req("V"; min = 0))
    @test !PS.check_constraints!(vec_fs, [1.0, -2.0], "v", ctx)
    @test ctx.errors[end].message == "-2.0 is below minimum 0"
    dir = mktempdir()
    write(joinpath(dir, "E.txt"), "header: Temperature Young's_Modulus\n0 200\n100 180\n")
    table = PS.read_table(joinpath(dir, "E.txt"), "Young's Modulus", "e", ctx)
    dep_fs = PS.FieldSpec(:e, PS.Dependent, PS.req("Young's Modulus"; max = 190))
    @test !PS.check_constraints!(dep_fs, table, "e", ctx)
    @test ctx.errors[end].message == "200.0 is above maximum 190"
    @test PS.check_constraints!(dep_fs, PS.Constant(100.0), "e", ctx)
    type_fs = PS.FieldSpec(:t, String, PS.req("Type"; allowed = ["Exodus", "CSV"]))
    @test PS.check_constraints!(type_fs, "CSV", "t", ctx)
    @test !PS.check_constraints!(type_fs, "VTK", "t", ctx)
    @test ctx.errors[end].message == "\"VTK\" is not one of: \"Exodus\", \"CSV\""
    union_fs = PS.FieldSpec(:u, Union{Nothing,Float64}, PS.opt("U"; default = nothing, min = 0))
    @test PS.check_constraints!(union_fs, nothing, "u", ctx)
end
