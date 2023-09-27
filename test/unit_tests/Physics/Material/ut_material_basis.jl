# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
include("../../../../src/Physics/Material/material_basis.jl")

@testset "get_all_elastic_moduli" begin
    parameter = Dict("Bulk Modulus" => 0, "Young's Modulus" => 0, "Shear Modulus" => 0, "Poisson's Ratio" => 0)
    test = get_all_elastic_moduli(Dict{String,Any}())

    @test test == Dict{String,Any}()
    test = get_all_elastic_moduli(Dict{String,Any}("Bulk Modulus" => 100, "Young's Modulus" => 10))
    @test sort(collect(keys(test))) == sort(collect(keys(parameter)))
    test = get_all_elastic_moduli(Dict{String,Any}("Bulk Modulus" => 1, "Shear Modulus" => 10))
    @test sort(collect(keys(test))) == sort(collect(keys(parameter)))
    test = get_all_elastic_moduli(Dict{String,Any}("Bulk Modulus" => 1, "Shear Modulus" => 10, "Poisson's Ratio" => 0.2))
    @test sort(collect(keys(test))) == sort(collect(keys(parameter)))

    println()
    test = get_all_elastic_moduli(Dict{String,Any}("Bulk Modulus" => 10, "Shear Modulus" => 10))
    @test test["Young's Modulus"] == Float32(22.5)
    @test test["Poisson's Ratio"] == Float32(0.25)
    @test test["Bulk Modulus"] == 10
    @test test["Shear Modulus"] == 10

    test = get_all_elastic_moduli(Dict{String,Any}("Bulk Modulus" => 5, "Shear Modulus" => 1.25))
    @test test["Young's Modulus"] / 3.4615384615384617 - 1 < 1e-7
    @test test["Poisson's Ratio"] / 0.45454545454545453 - 1 < 1e-7
    @test test["Bulk Modulus"] == 5
    @test test["Shear Modulus"] == Float32(1.25)
    test = get_all_elastic_moduli(Dict{String,Any}("Bulk Modulus" => 5, "Young's Modulus" => 1.25))
    @test test["Shear Modulus"] / 4.2857142857142855e-1 - 1 < 1e-7
    @test test["Poisson's Ratio"] / 0.4583333333333333 - 1 < 1e-7

    test = get_all_elastic_moduli(Dict{String,Any}("Poisson's Ratio" => 0.45, "Shear Modulus" => 1.25))
    @test test["Young's Modulus"] / 3.625e+0 - 1 < 1e-8
    @test test["Bulk Modulus"] / 1.2083333333333336e+1 - 1 < 1e-7
    @test test["Poisson's Ratio"] == Float32(0.45)
    @test test["Shear Modulus"] == Float32(1.25)
    test = get_all_elastic_moduli(Dict{String,Any}("Young's Modulus" => 5, "Poisson's Ratio" => 0.125))
    @test test["Bulk Modulus"] / 2.2222222222222223e+0 - 1 < 1e-7
    @test test["Shear Modulus"] / 2.2222222222222223e+0 - 1 < 1e-7
    @test test["Poisson's Ratio"] == Float32(0.125)
    @test test["Young's Modulus"] == 5
end
"""
@testset "get_Hook_matrix" begin
    parameter = get_all_elastic_moduli(Dict{String,Float32}("Bulk Modulus" => 5, "Shear Modulus" => 1.25))

    symmetry = "isotropic"
    E = parameter["Young's Modulus"]
    nu = parameter["Poisson's Ratio"]
    temp = 1 / ((1 + nu) * (1 - 2 * nu))
    C = get_Hook_matrix(parameter, symmetry, 3)
    for iID in 1:3
        @test C[iID, iID] / (E * (1 - nu) * temp) - 1 < 1e-7
        @test C[iID+3, iID+3] == parameter["Shear Modulus"]
        for jID in 1:3
            if iID != jID
                @test C[iID, jID] / (E * nu * temp) - 1 < 1e-7
            end
        end
    end

    symmetry = "isotropic_plain_strain"
    C2D = get_Hook_matrix(parameter, symmetry, 2)
    for iID in 1:2
        @test C2D[iID, iID] / (E * (1 - nu) * temp) - 1 < 1e-7
        for jID in 1:2
            if iID != jID
                @test C2D[iID, jID] / (E * nu * temp) - 1 < 1e-7
            end
        end
    end
    @test C2D[3, 3] == parameter["Shear Modulus"]
    symmetry = "isotropic_plain_stress"
    C2D_test = zeros(3, 3)
    Cinv = inv(C)
    C2D_test[1:2, 1:2] = Cinv[1:2, 1:2]
    C2D_test[3, 3] = Cinv[6, 6]
    C2D_test = inv(C2D_test)
    C = get_Hook_matrix(parameter, symmetry, 2)
    for iID in 1:3
        for jID in 1:3
            if C2D_test[iID, jID] != 0
                @test C[iID, jID] / C2D_test[iID, jID] - 1 < 1e-7
            end
        end
    end

    parameter = Dict{String,Float32}()
    for iID in 1:6
        for jID in 1:6
            parameter["C"*string(iID)*string(jID)] = iID * jID + jID
        end
    end
    symmetry = "anisotropic"
    C = get_Hook_matrix(parameter, symmetry, 3)
    for iID in 1:6
        for jID in 1:6
            @test C[iID, jID] == C[jID, iID]
            if jID >= iID
                @test C[iID, jID] == parameter["C"*string(iID)*string(jID)]
            end
        end
    end
    symmetry = "anisotropic_plain_strain"
    C = get_Hook_matrix(parameter, symmetry, 2)
    for iID in 1:2
        for jID in 1:2
            @test C[iID, jID] == C[jID, iID]
            if jID >= iID
                @test C[iID, jID] == parameter["C"*string(iID)*string(jID)]
            end
        end
    end
    @test C[3, 3] == parameter["C55"]
    @test C[1, 3] == parameter["C15"]
    @test C[2, 3] == parameter["C25"]
    @test C[3, 1] == parameter["C15"]
    @test C[3, 2] == parameter["C25"]

    symmetry = "anisotropic_plain_stress"

end
"""
