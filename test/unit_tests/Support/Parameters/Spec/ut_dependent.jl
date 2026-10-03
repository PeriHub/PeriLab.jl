# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const PS = PeriLab.ParameterSpec

const UT_TABLE_DIR = mktempdir()
function ut_write(name, content)
    path = joinpath(UT_TABLE_DIR, name)
    write(path, content)
    return path
end

const UT_E_T = """
# Young's modulus and Poisson's ratio over temperature
header: Temperature Young's_Modulus Poisson's_Ratio
0.0 200.0 0.30
100.0 180.0 0.30
200.0 150.0 0.31
"""

ut_table(alias) = PS.read_table(ut_write("E_T_$(hash(alias)).txt", UT_E_T), alias, "p",
                                PS.ParseContext())

@testset "Constant" begin
    c = PS.Constant(5.0)
    @test PS.value(c, 1) == 5.0
    @test PS.value(c, 99) == 5.0
    @test (@inferred PS.value(c, 3)) == 5.0
end

@testset "read_table" begin
    file = ut_write("E_T.txt", UT_E_T)
    ctx = PS.ParseContext()
    t = PS.read_table(file, "Young's Modulus", "M.\"Young's Modulus\"", ctx)
    @test !PS.has_errors(ctx)
    @test t isa PS.Table1D
    @test t.field_name == "Temperature"
    @test t.x == [0.0, 100.0, 200.0]
    @test t.y == [200.0, 180.0, 150.0]
    @test t.source == file
    @test !t.bound[]
    nu = PS.read_table(file, "Poisson's Ratio", "p", ctx)
    @test nu.y == [0.30, 0.30, 0.31]
end

@testset "read_table tolerates CRLF, tabs and blank lines" begin
    file = ut_write("crlf.txt",
                    "# comment\r\nheader:\tTemperature Young's_Modulus\r\n\r\n0.0\t200.0\r\n100.0  180.0\r\n")
    ctx = PS.ParseContext()
    t = PS.read_table(file, "Young's Modulus", "p", ctx)
    @test !PS.has_errors(ctx)
    @test t.x == [0.0, 100.0]
    @test t.y == [200.0, 180.0]
end

@testset "read_table errors" begin
    cases = [("missing.txt", nothing, "data file"),
             ("no_header.txt", "0.0 1.0\n1.0 2.0\n", "line 1: expected 'header:"),
             ("no_column.txt", "header: Temperature Density\n0 1\n1 2\n",
              "has no column \"Young's_Modulus\""),
             ("bad_value.txt", "header: Temperature Young's_Modulus\n0 1\n1 abc\n",
              "line 3: non-numeric value"),
             ("bad_count.txt", "header: Temperature Young's_Modulus\n0 1\n1\n",
              "line 3: expected 2 values, got 1"),
             ("unsorted.txt", "header: Temperature Young's_Modulus\n10 1\n0 2\n",
              "must be strictly increasing"),
             ("one_row.txt", "header: Temperature Young's_Modulus\n0 1\n",
              "at least 2 data rows")]
    for (name, content, expected) in cases
        file = content === nothing ? joinpath(UT_TABLE_DIR, name) : ut_write(name, content)
        ctx = PS.ParseContext()
        @test PS.read_table(file, "Young's Modulus", "p", ctx) === nothing
        @test length(ctx.errors) == 1
        @test occursin(expected, ctx.errors[1].message)
        @test ctx.errors[1].path == "p"
    end
end

@testset "Table1D value and binding" begin
    t = ut_table("Young's Modulus")
    @test_throws ArgumentError PS.value(t, 1)
    temperature = [0.0, 100.0, 200.0, -50.0, 500.0]
    PS.bind_table!(t, temperature)
    @test t.bound[]
    @test PS.value(t, 1) ≈ 200.0
    @test PS.value(t, 2) ≈ 180.0
    @test PS.value(t, 3) ≈ 150.0
    @test t.warn[]
    @test PS.value(t, 4) ≈ 200.0          # below range: nearest boundary value
    @test !t.warn[]                        # warned once, never again
    @test PS.value(t, 5) ≈ 150.0          # above range: nearest boundary value
    temperature[1] = 100.0                 # bound by reference: sees field updates
    @test PS.value(t, 1) ≈ 180.0
    @test (@inferred PS.value(t, 2)) ≈ 180.0
end

@testset "combine" begin
    @test PS.combine(+, PS.Constant(1.0), PS.Constant(2.0)) == PS.Constant(3.0)
    E = ut_table("Young's Modulus")
    doubled = PS.combine(*, E, PS.Constant(2.0))
    @test doubled isa PS.Table1D
    @test doubled.x == E.x
    @test doubled.y == [400.0, 360.0, 300.0]
    @test doubled.field_name == "Temperature"
    @test !doubled.bound[]
    @test PS.combine(-, PS.Constant(1000.0), E).y == [800.0, 820.0, 850.0]
    nu = ut_table("Poisson's Ratio")
    G = PS.combine((e, n) -> e / (2 * (1 + n)), E, nu)
    @test G.x == [0.0, 100.0, 200.0]
    @test G.y ≈ [200.0 / 2.6, 180.0 / 2.6, 150.0 / 2.62]
    other = PS.read_table(ut_write("E_D.txt", "header: Damage Young's_Modulus\n0 1\n1 2\n"),
                          "Young's Modulus", "p", PS.ParseContext())
    @test_throws ArgumentError PS.combine(+, E, other)
end

@testset "read_table tolerates a UTF-8 byte order mark" begin
    file = ut_write("bom.txt", "﻿header: Temperature Young's_Modulus\n0 200\n100 180\n")
    ctx = PS.ParseContext()
    t = PS.read_table(file, "Young's Modulus", "p", ctx)
    @test isempty(ctx.errors)
    @test t.y == [200.0, 180.0]
end
