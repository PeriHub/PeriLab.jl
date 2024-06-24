# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
# include("../../../../src/Physics/Material/material_basis.jl")

@testset "ut_flaw_function" begin

    stress::Float64 = 5.3
    @test flaw_function(Dict(), Vector{Float64}([1, 2]), stress) == stress

    @test isnothing(flaw_function(Dict("Flaw Function" => Dict()), Vector{Float64}([1, 2]), stress))
    @test isnothing(flaw_function(Dict("Flaw Function" => Dict("Active" => false)), Vector{Float64}([1, 2]), stress))
    @test flaw_function(Dict("Flaw Function" => Dict("Active" => false, "Function" => "Pre-defined")), Vector{Float64}([1, 2]), stress) == stress
    @test isnothing(flaw_function(Dict("Flaw Function" => Dict("Active" => true, "Function" => "Pre-defined", "Flaw Location X" => 1.1, "Flaw Location Y" => 1.1, "Flaw Magnitude" => 1.3, "Flaw Size" => 0.2)), Vector{Float64}([1, 2]), stress))
    @test isnothing(flaw_function(Dict("Flaw Function" => Dict("Active" => true, "Function" => "Pre-defined", "Flaw Location X" => 1.1, "Flaw Location Y" => 1.1, "Flaw Magnitude" => -1.3, "Flaw Size" => 0.2)), Vector{Float64}([1, 2]), stress))

    @test isapprox(flaw_function(Dict("Flaw Function" => Dict("Active" => true, "Function" => "Pre-defined", "Flaw Location X" => 1.1, "Flaw Location Y" => 1.1, "Flaw Magnitude" => 0.3, "Flaw Size" => 0.2)), Vector{Float64}([1, 2]), stress), 5.29999999)

    @test isapprox(flaw_function(Dict("Flaw Function" => Dict("Active" => true, "Function" => "Pre-defined", "Flaw Location X" => 1.1, "Flaw Location Y" => 1.1, "Flaw Location Z" => 2.1, "Flaw Magnitude" => 0.3, "Flaw Size" => 0.2)), Vector{Float64}([1, 2, 3]), stress), 5.29999999)


    #  @test flaw_function(Dict("Flaw Function" => Dict("Active" => true, "Function" => "x*x")), Vector{Float64}([1, 2]), stress) == 1
end
@testset "check_symmetry" begin
    @test isnothing(check_symmetry(Dict(), 2))
    @test isnothing(check_symmetry(Dict(), 3))
    @test isnothing(check_symmetry(Dict("Symmetry" => "a"), 2))
    @test isnothing(check_symmetry(Dict("Symmetry" => "plane"), 2))
    @test isnothing(check_symmetry(Dict("Symmetry" => "stress"), 2))
    @test isnothing(check_symmetry(Dict("Symmetry" => "strain"), 2))
    check_symmetry(Dict("Symmetry" => "plane stress"), 2)
    check_symmetry(Dict("Symmetry" => "plane strain"), 2)
    check_symmetry(Dict("Symmetry" => "3D"), 3)
end

@testset "get_symmetry" begin
    @test get_symmetry(Dict()) == "3D"
    @test get_symmetry(Dict("Symmetry" => "iso plane stress")) == "plane stress"
    @test get_symmetry(Dict("Symmetry" => "iso plane stress")) == "plane stress"
    @test get_symmetry(Dict("Symmetry" => "iso Plane Stress")) == "plane stress"

    @test get_symmetry(Dict("Symmetry" => "plane strain")) == "plane strain"
    @test get_symmetry(Dict("Symmetry" => "plane Strain")) == "plane strain"
    @test get_symmetry(Dict("Symmetry" => "PLANE strain")) == "plane strain"
    @test get_symmetry(Dict("Symmetry" => "plan strain")) == "3D"
end

@testset "get_all_elastic_moduli" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.clear_data_manager()
    test_data_manager.set_num_controller(3)
    ref_parameter = Dict("Bulk Modulus" => 0, "Computed" => true, "Young's Modulus" => 0, "Shear Modulus" => 0, "Poisson's Ratio" => 0)
    test = get_all_elastic_moduli(test_data_manager, Dict{String,Any}())
    @test isnothing(test)

    parameter = Dict{String,Any}("Bulk Modulus" => 1000, "Young's Modulus" => 10)
    get_all_elastic_moduli(test_data_manager, parameter)
    @test sort(collect(keys(parameter))) == sort(collect(keys(ref_parameter)))

    parameter = Dict{String,Any}("Bulk Modulus" => 1, "Shear Modulus" => 10)
    get_all_elastic_moduli(test_data_manager, parameter)
    @test sort(collect(keys(parameter))) == sort(collect(keys(ref_parameter)))

    parameter = Dict{String,Any}("Bulk Modulus" => 1, "Shear Modulus" => 10, "Poisson's Ratio" => 0.2)
    get_all_elastic_moduli(test_data_manager, parameter)
    @test sort(collect(keys(parameter))) == sort(collect(keys(ref_parameter)))

    parameter = Dict{String,Any}("Bulk Modulus" => 10, "Shear Modulus" => 10)
    get_all_elastic_moduli(test_data_manager, parameter)
    @test parameter["Young's Modulus"] == Float64(22.5)
    @test parameter["Poisson's Ratio"] == Float64(0.125)
    @test parameter["Bulk Modulus"] == 10
    @test parameter["Shear Modulus"] == 10

    parameter = Dict{String,Any}("Bulk Modulus" => 5, "Shear Modulus" => 1.25)
    get_all_elastic_moduli(test_data_manager, parameter)
    @test parameter["Young's Modulus"] / 3.4615384615384617 - 1 < 1e-7
    @test parameter["Poisson's Ratio"] / 0.45454545454545453 - 1 < 1e-7
    @test parameter["Bulk Modulus"] == 5
    @test parameter["Shear Modulus"] == Float64(1.25)

    parameter = Dict{String,Any}("Bulk Modulus" => 5, "Young's Modulus" => 1.25)
    get_all_elastic_moduli(test_data_manager, parameter)
    @test parameter["Shear Modulus"] / 4.2857142857142855e-1 - 1 < 1e-7
    @test parameter["Poisson's Ratio"] / 0.4583333333333333 - 1 < 1e-7

    parameter = Dict{String,Any}("Poisson's Ratio" => 0.45, "Shear Modulus" => 1.25)
    get_all_elastic_moduli(test_data_manager, parameter)
    @test parameter["Young's Modulus"] / 3.625e+0 - 1 < 1e-8
    @test parameter["Bulk Modulus"] / 1.2083333333333336e+1 - 1 < 1e-7
    @test parameter["Poisson's Ratio"] == Float64(0.45)
    @test parameter["Shear Modulus"] == Float64(1.25)

    parameter = Dict{String,Any}("Young's Modulus" => 5, "Poisson's Ratio" => 0.125)
    get_all_elastic_moduli(test_data_manager, parameter)
    @test parameter["Bulk Modulus"] / 2.2222222222222223e+0 - 1 < 1e-7
    @test parameter["Shear Modulus"] / 2.2222222222222223e+0 - 1 < 1e-7
    @test parameter["Poisson's Ratio"] == Float64(0.125)
    @test parameter["Young's Modulus"] == 5

    test_data_manager.create_constant_node_field("Bulk_Modulus", Float64, 1, 10)
    parameter = Dict{String,Any}("Shear Modulus" => 10)
    get_all_elastic_moduli(test_data_manager, parameter)
    @test parameter["Young's Modulus"] == [22.5, 22.5, 22.5]
    @test parameter["Poisson's Ratio"] == [0.125, 0.125, 0.125]
    @test parameter["Bulk Modulus"] == [10, 10, 10]
    @test parameter["Shear Modulus"] == [10, 10, 10]
end

@testset "get_Hooke_matrix" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.clear_data_manager()
    parameter = Dict{String,Any}("Bulk Modulus" => 5, "Shear Modulus" => 1.25, "Poisson's Ratio" => 0.2)
    get_all_elastic_moduli(test_data_manager, parameter)

    symmetry = "isotropic"
    E = parameter["Young's Modulus"]
    nu = parameter["Poisson's Ratio"]
    temp = 1 / ((1 + nu) * (1 - 2 * nu))
    C = get_Hooke_matrix(parameter, symmetry, 3)
    for iID in 1:3
        @test C[iID, iID] / (E * (1 - nu) * temp) - 1 < 1e-7
        @test C[iID+3, iID+3] == parameter["Shear Modulus"]
        for jID in 1:3
            if iID != jID
                @test C[iID, jID] / (E * nu * temp) - 1 < 1e-7
            end
        end
    end

    symmetry = "isotropic plane strain"
    C2D = get_Hooke_matrix(parameter, symmetry, 2)
    for iID in 1:2
        @test C2D[iID, iID] / (E * (1 - nu) * temp) - 1 < 1e-7
        for jID in 1:2
            if iID != jID
                @test C2D[iID, jID] / (E * nu * temp) - 1 < 1e-7
            end
        end
    end
    @test C2D[3, 3] == parameter["Shear Modulus"]

    symmetry = "missing"
    C2D = get_Hooke_matrix(parameter, symmetry, 2)
    for iID in 1:2
        @test C2D[iID, iID] / (E * (1 - nu) * temp) - 1 < 1e-7
        for jID in 1:2
            if iID != jID
                @test C2D[iID, jID] / (E * nu * temp) - 1 < 1e-7
            end
        end
    end
    @test C2D[3, 3] == parameter["Shear Modulus"]

    symmetry = "isotropic plane stress"
    C2D_test = zeros(3, 3)
    Cinv = inv(C)
    C2D_test[1:2, 1:2] = Cinv[1:2, 1:2]
    C2D_test[3, 3] = Cinv[6, 6]
    C2D_test = inv(C2D_test)
    C = get_Hooke_matrix(parameter, symmetry, 2)
    for iID in 1:3
        for jID in 1:3
            if C2D_test[iID, jID] != 0
                @test C[iID, jID] / C2D_test[iID, jID] - 1 < 1e-7
            end
        end
    end

    for iID in 1:6
        for jID in 1:6
            parameter["C"*string(iID)*string(jID)] = iID * jID + jID
        end
    end

    symmetry = "isotropic missing"
    @test isnothing(get_Hooke_matrix(parameter, symmetry, 2))

    symmetry = "anisotropic"
    C = get_Hooke_matrix(parameter, symmetry, 3)
    for iID in 1:6
        for jID in 1:6
            @test C[iID, jID] == C[jID, iID]
            if jID >= iID
                @test C[iID, jID] == parameter["C"*string(iID)*string(jID)]
            end
        end
    end
    symmetry = "anisotropic plane strain"
    C = get_Hooke_matrix(parameter, symmetry, 2)
    for iID in 1:2
        for jID in 1:2
            @test C[iID, jID] == C[jID, iID]
            if jID >= iID
                @test C[iID, jID] == parameter["C"*string(iID)*string(jID)]
            end
        end
    end
    @test C[3, 3] == parameter["C66"]
    @test C[1, 3] == parameter["C16"]
    @test C[2, 3] == parameter["C26"]
    @test C[3, 1] == parameter["C16"]
    @test C[3, 2] == parameter["C26"]

    symmetry = "anisotropic plane stress"
    C = get_Hooke_matrix(parameter, symmetry, 2)
    for iID in 1:3
        for jID in 1:3
            if C2D_test[iID, jID] != 0
                @test C[iID, jID] / C2D_test[iID, jID] - 1 < 1e-7
            end
        end
    end

    parameter = Dict{String,Any}("C11" => 2.0, "C12" => 3.0, "C13" => 4.0, "C22" => 5.0, "C33" => 7.0)
    @test isnothing(get_Hooke_matrix(parameter, symmetry, 2))

    symmetry = "anisotropic missing"
    @test isnothing(get_Hooke_matrix(parameter, symmetry, 2))

    symmetry = "missing"
    parameter = Dict{String,Any}("Missing" => 5)
    @test isnothing(get_Hooke_matrix(parameter, symmetry, 2))

end

@testset "ut_matrix_to_voigt" begin
    matrix = [1 2; 3 4]
    voigt = matrix_to_voigt(matrix)
    @test voigt[1] == 1
    @test voigt[2] == 4
    @test voigt[3] == 2.5
    matrix = [1 2 3; 4 5 6; 7 8 9]
    voigt = matrix_to_voigt(matrix)
    @test voigt[1] == 1
    @test voigt[2] == 5
    @test voigt[3] == 9
    @test voigt[4] == 7
    @test voigt[5] == 5
    @test voigt[6] == 3
    matrix = [1 2 3 3; 4 5 6 3; 7 8 9 3]
    @test isnothing(matrix_to_voigt(matrix))
end
@testset "ut_voigt_to_matrix" begin
    @test isnothing(voigt_to_matrix([1, 2.2]))
end