using Test
include("../../../../src/Physics/Material/material_basis.jl")

@testset "get_all_elastic_moduli" begin
    parameter = Dict("Bulk Modulus" => 0, "Young's Modulus" => 0, "Shear Modulus" => 0, "Poisson's Ratio" => 0)
    test = get_all_elastic_moduli(Dict{String,Float32}())

    @test sort(collect(keys(test))) == sort(collect(keys(parameter)))
    test = get_all_elastic_moduli(Dict{String,Float32}("Bulk Modulus" => 100, "Young's Modulus" => 10))
    @test sort(collect(keys(test))) == sort(collect(keys(parameter)))
    test = get_all_elastic_moduli(Dict{String,Float32}("Bulk Modulus" => 1, "Shear Modulus" => 10))
    @test sort(collect(keys(test))) == sort(collect(keys(parameter)))
    test = get_all_elastic_moduli(Dict{String,Float32}("Bulk Modulus" => 1, "Shear Modulus" => 10, "Poisson's Ratio" => 0.2))
    @test sort(collect(keys(test))) == sort(collect(keys(parameter)))

    println()
    test = get_all_elastic_moduli(Dict{String,Float32}("Bulk Modulus" => 10, "Shear Modulus" => 10))
    @test test["Young's Modulus"] == Float32(22.5)
    @test test["Poisson's Ratio"] == Float32(0.25)
    @test test["Bulk Modulus"] == 10
    @test test["Shear Modulus"] == 10

    test = get_all_elastic_moduli(Dict{String,Float32}("Bulk Modulus" => 5, "Shear Modulus" => 1.25))
    @test test["Young's Modulus"] / 3.4615384615384617 - 1 < 1e-7
    @test test["Poisson's Ratio"] / 0.45454545454545453 - 1 < 1e-7
    @test test["Bulk Modulus"] == 5
    @test test["Shear Modulus"] == Float32(1.25)
    test = get_all_elastic_moduli(Dict{String,Float32}("Bulk Modulus" => 5, "Young's Modulus" => 1.25))
    @test test["Shear Modulus"] / 4.2857142857142855e-1 - 1 < 1e-7
    @test test["Poisson's Ratio"] / 0.4583333333333333 - 1 < 1e-7

    test = get_all_elastic_moduli(Dict{String,Float32}("Poisson's Ratio" => 0.45, "Shear Modulus" => 1.25))
    @test test["Young's Modulus"] / 3.625e+0 - 1 < 1e-8
    @test test["Bulk Modulus"] / 1.2083333333333336e+1 - 1 < 1e-7
    @test test["Poisson's Ratio"] == Float32(0.45)
    @test test["Shear Modulus"] == Float32(1.25)
    test = get_all_elastic_moduli(Dict{String,Float32}("Young's Modulus" => 5, "Poisson's Ratio" => 0.125))
    @test test["Bulk Modulus"] / 2.2222222222222223e+0 - 1 < 1e-7
    @test test["Shear Modulus"] / 2.2222222222222223e+0 - 1 < 1e-7
    @test test["Poisson's Ratio"] == Float32(0.125)
    @test test["Young's Modulus"] == 5
end