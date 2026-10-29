# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const PS = PeriLab.ParameterSpec

# A module written the way a module author would, loaded at runtime like a
# licensed module (registration in __init__).
const UT_AUTHOR_MODULE = raw"""
module UTAuthorMaterial
using PeriLab.ParameterSpec: @params, register_model!, value, combine, Constant
import PeriLab.ParameterSpec: derive

@params struct UTLinearElastic
    youngs_modulus::Dependent = req("Young's Modulus"; min = 0, quantity = :stress)
    poissons_ratio::Float64 = req("Poisson's Ratio"; min = -1, max = 0.5)
    shear_modulus::Dependent = opt("Shear Modulus"; default = 0.0, quantity = :stress)
end

function derive(p::UTLinearElastic)
    G = combine((E, nu) -> E / (2 * (1 + nu)), p.youngs_modulus, Constant(p.poissons_ratio))
    return UTLinearElastic(p.youngs_modulus, p.poissons_ratio, G)
end

function sum_shear(p::UTLinearElastic, nodes)
    G = p.shear_modulus
    s = 0.0
    for iID in nodes
        s += value(G, iID)
    end
    return s
end

function __init__()
    register_model!(:ut_author, "UT Linear Elastic", UTLinearElastic)
end
end
"""
Base.include_string(@__MODULE__, UT_AUTHOR_MODULE, "UTAuthorMaterial.jl")

const UT_E2E_DIR = mktempdir()
write(joinpath(UT_E2E_DIR, "E_T.txt"),
      "header: Temperature Young's_Modulus\n0.0 200.0\n100.0 180.0\n")

function ut_author_model(youngs_modulus)
    ctx = PS.ParseContext(directory = UT_E2E_DIR)
    m = PS.parse_model(:ut_author,
                       Dict{String,Any}("Material Model" => "UT Linear Elastic",
                                        "Young's Modulus" => youngs_modulus,
                                        "Poisson's Ratio" => 0.3),
                       "Models.\"Material Models\".A", ctx; name_key = "Material Model")
    return m, ctx
end

@testset "runtime-loaded module registers in __init__" begin
    @test PS.lookup_model(:ut_author, "UT Linear Elastic") === UTAuthorMaterial.UTLinearElastic
end

@testset "constant parameters: derive, type stability, no allocations" begin
    m, ctx = ut_author_model(260.0)
    @test isempty(ctx.errors)
    @test m.shear_modulus isa PS.Constant
    @test m.shear_modulus.value ≈ 100.0
    @test isconcretetype(typeof(m))
    @test (@inferred UTAuthorMaterial.sum_shear(m, 1:2)) ≈ 200.0
    UTAuthorMaterial.sum_shear(m, 1:2)
    @test (@allocated UTAuthorMaterial.sum_shear(m, 1:2)) == 0
end

@testset "table parameters: derive, bind, evaluate" begin
    m, ctx = ut_author_model("E_T.txt")
    @test isempty(ctx.errors)
    @test m.youngs_modulus isa PS.Table1D
    @test m.shear_modulus isa PS.Table1D
    temperature = [0.0, 100.0]
    PS.bind_dependents!(m, name -> name == "Temperature" ? temperature : nothing,
                        "Models.\"Material Models\".A", ctx)
    @test isempty(ctx.errors)
    @test m.youngs_modulus.bound[] && m.shear_modulus.bound[]
    @test (@inferred UTAuthorMaterial.sum_shear(m, 1:2)) ≈ 200.0 / 2.6 + 180.0 / 2.6
end

@testset "binding errors" begin
    m, _ = ut_author_model("E_T.txt")
    ctx = PS.ParseContext()
    PS.bind_dependents!(m, name -> nothing, "B", ctx)
    @test ctx.errors[1].path == "B.\"Young's Modulus\""
    @test ctx.errors[1].message ==
          "field \"Temperature\" required by $(joinpath(UT_E2E_DIR, "E_T.txt")) does not exist"
    m2, _ = ut_author_model("E_T.txt")
    ctx2 = PS.ParseContext()
    PS.bind_dependents!(m2, name -> [1, 2], "B", ctx2)
    @test ctx2.errors[1].message ==
          "field \"Temperature\" must be a per-node Vector{Float64}, got Vector{Int64}"
end

@testset "binding walks composites and named entries" begin
    a, _ = ut_author_model("E_T.txt")
    b, _ = ut_author_model("E_T.txt")
    c, _ = ut_author_model("E_T.txt")
    temperature = [0.0, 100.0]
    lookup = name -> temperature
    ctx = PS.ParseContext()
    PS.bind_dependents!(PS.Composite((a, b)), lookup, "C", ctx)
    PS.bind_dependents!(Dict("block_1" => c), lookup, "D", ctx)
    @test isempty(ctx.errors)
    @test a.youngs_modulus.bound[] && b.youngs_modulus.bound[] && c.youngs_modulus.bound[]
end

# Licensed modules are loaded by `Base.include` from inside a running function;
# their generated methods are newer than the caller's world.
function ut_load_and_parse_in_one_function()
    Base.include_string(@__MODULE__, """
                        module UTWorldAgeMaterial
                        using PeriLab.ParameterSpec: @params, register_model!
                        @params struct UTWorldAge
                            stiffness::Dependent = req("Stiffness"; min = 0)
                        end
                        __init__() = register_model!(:ut_world_age, "UT World Age", UTWorldAge)
                        end
                        """, "UTWorldAgeMaterial.jl")
    ctx = PS.ParseContext()
    m = PS.parse_model(:ut_world_age,
                       Dict{String,Any}("Material Model" => "UT World Age", "Stiffness" => 3.0),
                       "M", ctx; name_key = "Material Model")
    PS.bind_dependents!(m, name -> nothing, "M", ctx)
    section = PS.parse_section(typeof(m), Dict{String,Any}("Stiffness" => 4.0), "S", ctx)
    return m, section, ctx
end

@testset "models loaded at runtime can be parsed in the same call" begin
    m, section, ctx = ut_load_and_parse_in_one_function()
    @test isempty(ctx.errors)
    @test m.stiffness.value == 3.0
    @test section.stiffness.value == 4.0
end
