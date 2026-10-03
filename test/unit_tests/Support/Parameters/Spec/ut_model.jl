# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const PS = PeriLab.ParameterSpec

PS.@params struct UTElastic
    youngs_modulus::Dependent = req("Young's Modulus"; min = 0)
    poissons_ratio::Float64 = req("Poisson's Ratio")
end

PS.@params struct UTPlastic
    youngs_modulus::Dependent = req("Young's Modulus"; min = 0)
    yield_stress::Float64 = req("Yield Stress"; min = 0)
end

PS.@params struct UTConflicting
    youngs_modulus::Float64 = req("Young's Modulus")
end

PS.register_model!(:ut_material, "UT Elastic", UTElastic)
PS.register_model!(:ut_material, "UT Plastic", UTPlastic)
PS.register_model!(:ut_material, "UT Conflicting", UTConflicting)
PS.register_unavailable!(:ut_material, "UT Licensed")

const UT_KEY = "Material Model"
const UT_PATH = "Models.\"Material Models\".Steel"

function ut_parse(dict; strict = true)
    ctx = PS.ParseContext(strict = strict)
    return PS.parse_model(:ut_material, dict, UT_PATH, ctx; name_key = UT_KEY), ctx
end

@testset "registry" begin
    @test PS.lookup_model(:ut_material, "UT Elastic") === UTElastic
    @test PS.lookup_model(:ut_material, "nope") === nothing
    @test PS.lookup_model(:no_such_category, "UT Elastic") === nothing
    @test PS.lookup_model(:ut_material, "UT Licensed") isa PS.UnavailableModel
    @test PS.registered_names(:ut_material) == ["UT Conflicting", "UT Elastic", "UT Plastic"]
    PS.register_model!(:ut_material, "UT Elastic", UTElastic)      # same type again: fine
    e = try
        PS.register_model!(:ut_material, "UT Elastic", UTPlastic)
    catch err
        err
    end
    @test e isa PS.ParamsDefinitionError
    @test e.msg == "ut_material model \"UT Elastic\" is already registered by UTElastic"
    @test_throws PS.ParamsDefinitionError PS.register_model!(:ut_material, "X", Float64)
    PS.register_unavailable!(:ut_material, "UT Elastic")           # never hides a real model
    @test PS.lookup_model(:ut_material, "UT Elastic") === UTElastic
    PS.register_unavailable!(:ut_material, "UT Later")
    PS.register_model!(:ut_material, "UT Later", UTPlastic)        # real replaces stub
    @test PS.lookup_model(:ut_material, "UT Later") === UTPlastic
end

@testset "no model" begin
    m, ctx = ut_parse(nothing)
    @test m === PS.NoModel()
    @test isempty(ctx.errors)
end

@testset "single model" begin
    m, ctx = ut_parse(Dict{String,Any}(UT_KEY => "UT Elastic", "Young's Modulus" => 210.0,
                                       "Poisson's Ratio" => 0.3))
    @test isempty(ctx.errors)
    @test m isa UTElastic{PS.Constant}
end

@testset "composite model" begin
    m, ctx = ut_parse(Dict{String,Any}(UT_KEY => "UT Elastic + UT Plastic",
                                       "Young's Modulus" => 210.0,
                                       "Poisson's Ratio" => 0.3, "Yield Stress" => 5.0))
    @test isempty(ctx.errors)
    @test m isa PS.Composite{Tuple{UTElastic{PS.Constant},UTPlastic{PS.Constant}}}
    @test isconcretetype(typeof(m))
    @test m.parts[2].yield_stress == 5.0
    @test m.parts[1].youngs_modulus.value == m.parts[2].youngs_modulus.value == 210.0
end

@testset "composite: a key is unknown only if no part declares it" begin
    m, ctx = ut_parse(Dict{String,Any}(UT_KEY => "UT Elastic + UT Plastic",
                                       "Young's Modulus" => 210.0,
                                       "Poisson's Ratio" => 0.3, "Yeild Stress" => 5.0))
    @test m === nothing
    msgs = Dict(e.path => e.message for e in ctx.errors)
    @test msgs["$UT_PATH.\"Yield Stress\""] == "missing (required by UT Plastic)"
    @test msgs["$UT_PATH.\"Yeild Stress\""] == "unknown key — did you mean \"Yield Stress\"?"
    @test length(ctx.errors) == 2
end

@testset "model name errors" begin
    name_path = "$UT_PATH.\"Material Model\""
    cases = [("UT Elastc", "model \"UT Elastc\" not found — did you mean \"UT Elastic\"?"),
             ("Completely Different",
              "model \"Completely Different\" not found; it may require a licensed module"),
             ("UT Licensed",
              "model \"UT Licensed\" requires a license that is not available"),
             ("UT Elastic + ", "empty model name in \"UT Elastic + \""),
             (42, "expected a model name, got 42")]
    for (name, expected) in cases
        m, ctx = ut_parse(Dict{String,Any}(UT_KEY => name, "Young's Modulus" => 1.0,
                                           "Poisson's Ratio" => 0.3))
        @test m === nothing
        @test ctx.errors[1].path == name_path
        @test ctx.errors[1].message == expected
    end
    m, ctx = ut_parse(Dict{String,Any}("Young's Modulus" => 1.0))
    @test m === nothing
    @test ctx.errors[1].path == name_path
    @test ctx.errors[1].message == "missing (names the model to use)"
end

@testset "alias type conflict between combined models" begin
    m, ctx = ut_parse(Dict{String,Any}(UT_KEY => "UT Elastic + UT Conflicting",
                                       "Young's Modulus" => 1.0, "Poisson's Ratio" => 0.3))
    @test m === nothing
    @test ctx.errors[1].path == "$UT_PATH.\"Young's Modulus\""
    @test ctx.errors[1].message ==
          "declared as Dependent by \"UT Elastic\" but as Float64 by \"UT Conflicting\"; models combined with + must agree"
end
