# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const PS = PeriLab.ParameterSpec

@testset "join_path" begin
    @test PS.join_path("", "Blocks") == "Blocks"
    @test PS.join_path("Blocks", "block_1") == "Blocks.block_1"
    @test PS.join_path("Models", "Material Models") == "Models.\"Material Models\""
    @test PS.join_path("M", "Young's Modulus") == "M.\"Young's Modulus\""
end

@testset "levenshtein" begin
    @test PS.levenshtein("kitten", "sitting") == 3
    @test PS.levenshtein("", "abc") == 3
    @test PS.levenshtein("same", "same") == 0
end

@testset "suggest" begin
    candidates = ["Young's Modulus", "Poisson's Ratio", "Density"]
    @test PS.suggest("Poissons Ratio", candidates) == "Poisson's Ratio"
    @test PS.suggest("young's modulus", candidates) == "Young's Modulus"
    @test PS.suggest("Youngs Modulos", candidates) == "Young's Modulus"
    @test PS.suggest("Horizon", candidates) === nothing
    @test PS.suggest("x", String[]) === nothing
    @test PS.unknown_key_message("Densty", candidates) ==
          "unknown key — did you mean \"Density\"?"
    @test PS.unknown_key_message("Horizon", candidates) == "unknown key"
end

@testset "ParseContext and report!" begin
    ctx = PS.ParseContext(directory = "/tmp", strict = true)
    @test ctx.directory == "/tmp"
    @test ctx.strict
    @test !PS.has_errors(ctx)
    PS.add_warning!(ctx, "a", "only a warning")
    @test !PS.has_errors(ctx)
    @test PS.report!(ctx) === nothing
    PS.add_error!(ctx, "Blocks.block_1.Horizon", "-0.1 is below minimum 0")
    PS.add_error!(ctx, "x", "missing")
    @test PS.has_errors(ctx)
    errors = filter(e -> e.severity == :error, ctx.errors)
    @test PS.format_errors(errors) ==
          "Input errors (2):\n  Blocks.block_1.Horizon: -0.1 is below minimum 0\n  x: missing"
    @test_throws PeriLab.PeriLabExceptions.PeriLabError PS.report!(ctx)
end

@testset "strict_mode" begin
    @test PS.strict_mode(Dict{String,Any}()) == true
    @test PS.strict_mode(Dict{String,Any}("Strict Validation" => false)) == false
    @test PS.strict_mode(Dict{String,Any}("Strict Validation" => true);
                         no_strict_flag = true) == false
    @test_throws PeriLab.PeriLabExceptions.PeriLabError PS.strict_mode(Dict{String,Any}("Strict Validation" => "no"))
end

@testset "ParamsDefinitionError" begin
    e = PS.ParamsDefinitionError("X.y: bad")
    @test sprint(showerror, e) == "ParamsDefinitionError: X.y: bad"
end

@testset "unknown-key errors mention how to downgrade them" begin
    errors = [PS.InputError("Solver.Foo", "unknown key", :error),
              PS.InputError("x", "missing", :error)]
    @test PS.format_errors(errors) ==
          "Input errors (2):\n  Solver.Foo: unknown key\n  x: missing\n" *
          "Unknown keys can be reported as warnings instead: run with --no_strict or set \"Strict Validation: false\"."
end
