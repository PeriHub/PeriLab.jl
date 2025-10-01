# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
#include("../../../../src/PeriLab.jl")
#using .PeriLab

@testset "ut_init_local_damping_due_to_damage" begin
    test_data_manager = PeriLab.Data_Manager
    test_data_manager.initialize_data()
    test_data_manager.set_num_controller(3)
    nn = test_data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn .= 2
    @test isnothing(PeriLab.Solver_Manager.Material_Basis.init_local_damping_due_to_damage(test_data_manager,
                                                                                           collect(1:2),
                                                                                           Dict(),
                                                                                           Dict("Local Damping" => Dict())))
    @test isnothing(PeriLab.Solver_Manager.Material_Basis.init_local_damping_due_to_damage(test_data_manager,
                                                                                           collect(1:2),
                                                                                           Dict(),
                                                                                           Dict("Local Damping" => Dict("Representative Young's modulus" => 0))))
end

@testset "ut_apply_pointwise_E" begin
    nodes = 2:3

    bond_force=[[ones(2), ones(2)], [ones(2), ones(2)], [ones(2), ones(2)]]
    E = 3.3
    bond_force[3][2][2] = 3
    PeriLab.Solver_Manager.Material_Basis.apply_pointwise_E(nodes, E, bond_force)
    @test bond_force[1][1][1] == 1
    @test bond_force[1][1][2] == 1
    @test bond_force[1][2][1] == 1
    @test bond_force[1][2][2] == 1
    @test bond_force[2][2][1] == E
    @test bond_force[2][2][2] == E
    @test bond_force[2][2][1] == E
    @test bond_force[2][2][2] == E
    @test bond_force[3][1][1] == E
    @test bond_force[3][1][2] == E
    @test bond_force[3][2][1] == E
    @test bond_force[3][2][2] == 3 * E

    bond_force=[[ones(2), ones(2)], [ones(2), ones(2)], [ones(2), ones(2), ones(2)]]
    E = zeros(4)
    PeriLab.Solver_Manager.Material_Basis.apply_pointwise_E(nodes, E, bond_force)
    @test bond_force[1][1][1] == 1
    @test bond_force[1][1][2] == 1
    @test bond_force[1][2][1] == 1
    @test bond_force[1][2][2] == 1
    @test bond_force[2][1][1] == 0
    @test bond_force[2][1][2] == 0
    @test bond_force[2][2][1] == 0
    @test bond_force[2][2][2] == 0
    @test bond_force[3][1][1] == 0
    @test bond_force[3][1][2] == 0
    @test bond_force[3][2][1] == 0
    @test bond_force[3][2][2] == 0
    @test bond_force[3][3][1] == 0
    @test bond_force[3][3][2] == 0
end

@testset "ut_distribute_forces" begin
    nodes = [1, 2]
    dof = 2
    nlist = [[2, 3], [1, 3]]
    nBonds = fill(dof, 2)
    bond_force = [[fill(1.0, dof) for j in 1:n] for n in nBonds]
    volume = [1.0, 1.0, 1.0]
    bond_damage = [fill(0.5, n) for n in nBonds]
    force_densities = zeros(3, 2)

    expected_force_densities = copy(force_densities)
    for iID in nodes
        expected_force_densities[iID,
        :] .+= transpose(sum(bond_damage[iID] .*
                                                           mapreduce(permutedims, vcat,
                                                                     bond_force[iID]) .*
                                                           volume[nlist[iID]],
                                                           dims = 1))
        expected_force_densities[nlist[iID],
        :] .-= bond_damage[iID] .*
                                                    mapreduce(permutedims, vcat,
                                                              bond_force[iID]) .*
                                                    volume[iID]
    end

    PeriLab.Solver_Manager.Material_Basis.distribute_forces!(force_densities, nodes, nlist,
                                                             bond_force, volume,
                                                             bond_damage)
    @test force_densities ≈ expected_force_densities
end
@testset "ut_flaw_function" begin
    stress::Float64 = 5.3
    @test PeriLab.Solver_Manager.Material_Basis.flaw_function(Dict(),
                                                              Vector{Float64}([1, 2]),
                                                              stress) == stress

    @test isnothing(PeriLab.Solver_Manager.Material_Basis.flaw_function(Dict("Flaw Function" => Dict()),
                                                                        Vector{Float64}([
                                                                                            1,
                                                                                            2
                                                                                        ]),
                                                                        stress))
    @test isnothing(PeriLab.Solver_Manager.Material_Basis.flaw_function(Dict("Flaw Function" => Dict("Active" => false)),
                                                                        Vector{Float64}([
                                                                                            1,
                                                                                            2
                                                                                        ]),
                                                                        stress))
    @test PeriLab.Solver_Manager.Material_Basis.flaw_function(Dict("Flaw Function" => Dict("Active" => false,
                                                                                           "Function" => "Pre-defined")),
                                                              Vector{Float64}([1, 2]),
                                                              stress) == stress
    @test isnothing(PeriLab.Solver_Manager.Material_Basis.flaw_function(Dict("Flaw Function" => Dict("Active" => true,
                                                                                                     "Function" => "Pre-defined",
                                                                                                     "Flaw Location X" => 1.1,
                                                                                                     "Flaw Location Y" => 1.1,
                                                                                                     "Flaw Magnitude" => 1.3,
                                                                                                     "Flaw Size" => 0.2)),
                                                                        Vector{Float64}([
                                                                                            1,
                                                                                            2
                                                                                        ]),
                                                                        stress))
    @test isnothing(PeriLab.Solver_Manager.Material_Basis.flaw_function(Dict("Flaw Function" => Dict("Active" => true,
                                                                                                     "Function" => "Pre-defined",
                                                                                                     "Flaw Location X" => 1.1,
                                                                                                     "Flaw Location Y" => 1.1,
                                                                                                     "Flaw Magnitude" => -1.3,
                                                                                                     "Flaw Size" => 0.2)),
                                                                        Vector{Float64}([
                                                                                            1,
                                                                                            2
                                                                                        ]),
                                                                        stress))

    @test isapprox(PeriLab.Solver_Manager.Material_Basis.flaw_function(Dict("Flaw Function" => Dict("Active" => true,
                                                                                                    "Function" => "Pre-defined",
                                                                                                    "Flaw Location X" => 1.1,
                                                                                                    "Flaw Location Y" => 1.1,
                                                                                                    "Flaw Magnitude" => 0.3,
                                                                                                    "Flaw Size" => 0.2)),
                                                                       Vector{Float64}([
                                                                                           1,
                                                                                           2
                                                                                       ]),
                                                                       stress),
                   5.29999999)

    @test isapprox(PeriLab.Solver_Manager.Material_Basis.flaw_function(Dict("Flaw Function" => Dict("Active" => true,
                                                                                                    "Function" => "Pre-defined",
                                                                                                    "Flaw Location X" => 1.1,
                                                                                                    "Flaw Location Y" => 1.1,
                                                                                                    "Flaw Location Z" => 2.1,
                                                                                                    "Flaw Magnitude" => 0.3,
                                                                                                    "Flaw Size" => 0.2)),
                                                                       Vector{Float64}([
                                                                                           1,
                                                                                           2,
                                                                                           3
                                                                                       ]),
                                                                       stress),
                   5.29999999)

    #  @test PeriLab.Solver_Manager.Material_Basis.flaw_function(Dict("Flaw Function" => Dict("Active" => true, "Function" => "x*x")), Vector{Float64}([1, 2]), stress) == 1
end
@testset "check_symmetry" begin
    @test isnothing(PeriLab.Solver_Manager.Material_Basis.check_symmetry(Dict(), 2))
    @test isnothing(PeriLab.Solver_Manager.Material_Basis.check_symmetry(Dict(), 3))
    @test isnothing(PeriLab.Solver_Manager.Material_Basis.check_symmetry(Dict("Symmetry" => "a"),
                                                                         2))
    @test isnothing(PeriLab.Solver_Manager.Material_Basis.check_symmetry(Dict("Symmetry" => "plane"),
                                                                         2))
    @test isnothing(PeriLab.Solver_Manager.Material_Basis.check_symmetry(Dict("Symmetry" => "stress"),
                                                                         2))
    @test isnothing(PeriLab.Solver_Manager.Material_Basis.check_symmetry(Dict("Symmetry" => "strain"),
                                                                         2))
    PeriLab.Solver_Manager.Material_Basis.check_symmetry(Dict("Symmetry" => "plane stress"),
                                                         2)
    PeriLab.Solver_Manager.Material_Basis.check_symmetry(Dict("Symmetry" => "plane strain"),
                                                         2)
    PeriLab.Solver_Manager.Material_Basis.check_symmetry(Dict("Symmetry" => "3D"), 3)
end

@testset "get_symmetry" begin
    @test PeriLab.Solver_Manager.Material_Basis.get_symmetry(Dict()) == "3D"
    @test PeriLab.Solver_Manager.Material_Basis.get_symmetry(Dict("Symmetry" => "iso plane stress")) ==
          "plane stress"
    @test PeriLab.Solver_Manager.Material_Basis.get_symmetry(Dict("Symmetry" => "iso plane stress")) ==
          "plane stress"
    @test PeriLab.Solver_Manager.Material_Basis.get_symmetry(Dict("Symmetry" => "iso Plane Stress")) ==
          "plane stress"

    @test PeriLab.Solver_Manager.Material_Basis.get_symmetry(Dict("Symmetry" => "plane strain")) ==
          "plane strain"
    @test PeriLab.Solver_Manager.Material_Basis.get_symmetry(Dict("Symmetry" => "plane Strain")) ==
          "plane strain"
    @test PeriLab.Solver_Manager.Material_Basis.get_symmetry(Dict("Symmetry" => "PLANE strain")) ==
          "plane strain"
    @test PeriLab.Solver_Manager.Material_Basis.get_symmetry(Dict("Symmetry" => "plan strain")) ==
          "3D"
end

@testset "get_all_elastic_moduli" begin
    test_data_manager = PeriLab.Data_Manager
    test_data_manager.initialize_data()
    test_data_manager.set_num_controller(3)
    test_data_manager.set_dof(2)
    ref_parameter = Dict("Material Model" => "PD Solid Elastic",
                         "Bulk Modulus" => 0,
                         "Computed" => true,
                         "Young's Modulus" => 0,
                         "Shear Modulus" => 0,
                         "Poisson's Ratio" => 0,
                         "Symmetry" => "isotropic")
    test = PeriLab.Solver_Manager.Material_Basis.get_all_elastic_moduli(test_data_manager,
                                                                        Dict{String,Any}("Material Model" => "PD Solid Elastic"))
    @test isnothing(test)

    parameter = Dict{String,Any}("Material Model" => "PD Solid Elastic",
                                 "Bulk Modulus" => 1000,
                                 "Young's Modulus" => 10)
    PeriLab.Solver_Manager.Material_Basis.get_all_elastic_moduli(test_data_manager,
                                                                 parameter)
    @test sort(collect(keys(parameter))) == sort(collect(keys(ref_parameter)))

    parameter = Dict{String,Any}("Material Model" => "PD Solid Elastic",
                                 "Bulk Modulus" => 1,
                                 "Shear Modulus" => 10)
    PeriLab.Solver_Manager.Material_Basis.get_all_elastic_moduli(test_data_manager,
                                                                 parameter)
    @test sort(collect(keys(parameter))) == sort(collect(keys(ref_parameter)))

    parameter = Dict{String,Any}("Material Model" => "PD Solid Elastic",
                                 "Bulk Modulus" => 1,
                                 "Shear Modulus" => 10,
                                 "Poisson's Ratio" => 0.2)
    PeriLab.Solver_Manager.Material_Basis.get_all_elastic_moduli(test_data_manager,
                                                                 parameter)
    @test sort(collect(keys(parameter))) == sort(collect(keys(ref_parameter)))

    parameter = Dict{String,Any}("Material Model" => "PD Solid Elastic",
                                 "Bulk Modulus" => 10)
    @test isnothing(PeriLab.Solver_Manager.Material_Basis.get_all_elastic_moduli(test_data_manager,
                                                                                 parameter))

    parameter = Dict{String,Any}("Material Model" => "PD Solid Elastic",
                                 "Bulk Modulus" => 10,
                                 "Shear Modulus" => 10)
    PeriLab.Solver_Manager.Material_Basis.get_all_elastic_moduli(test_data_manager,
                                                                 parameter)
    @test parameter["Young's Modulus"] == Float64(22.5)
    @test parameter["Poisson's Ratio"] == Float64(0.125)
    @test parameter["Bulk Modulus"] == 10
    @test parameter["Shear Modulus"] == 10

    parameter = Dict{String,Any}("Material Model" => "PD Solid Elastic",
                                 "Bulk Modulus" => 5,
                                 "Shear Modulus" => 1.25)
    PeriLab.Solver_Manager.Material_Basis.get_all_elastic_moduli(test_data_manager,
                                                                 parameter)
    @test parameter["Young's Modulus"] / 3.4615384615384617 - 1 < 1e-7
    @test parameter["Poisson's Ratio"] / 0.45454545454545453 - 1 < 1e-7
    @test parameter["Bulk Modulus"] == 5
    @test parameter["Shear Modulus"] == Float64(1.25)

    parameter = Dict{String,Any}("Material Model" => "PD Solid Elastic",
                                 "Bulk Modulus" => 5,
                                 "Young's Modulus" => 1.25)
    PeriLab.Solver_Manager.Material_Basis.get_all_elastic_moduli(test_data_manager,
                                                                 parameter)
    @test parameter["Shear Modulus"] / 4.2857142857142855e-1 - 1 < 1e-7
    @test parameter["Poisson's Ratio"] / 0.4583333333333333 - 1 < 1e-7

    parameter = Dict{String,Any}("Material Model" => "PD Solid Elastic",
                                 "Poisson's Ratio" => 0.45,
                                 "Shear Modulus" => 1.25)
    PeriLab.Solver_Manager.Material_Basis.get_all_elastic_moduli(test_data_manager,
                                                                 parameter)
    @test parameter["Young's Modulus"] / 3.625e+0 - 1 < 1e-8
    @test parameter["Bulk Modulus"] / 1.2083333333333336e+1 - 1 < 1e-7
    @test parameter["Poisson's Ratio"] == Float64(0.45)
    @test parameter["Shear Modulus"] == Float64(1.25)

    parameter = Dict{String,Any}("Material Model" => "PD Solid Elastic",
                                 "Young's Modulus" => 5,
                                 "Poisson's Ratio" => 0.125)
    PeriLab.Solver_Manager.Material_Basis.get_all_elastic_moduli(test_data_manager,
                                                                 parameter)
    @test parameter["Bulk Modulus"] / 2.2222222222222223e+0 - 1 < 1e-7
    @test parameter["Shear Modulus"] / 2.2222222222222223e+0 - 1 < 1e-7
    @test parameter["Poisson's Ratio"] == Float64(0.125)
    @test parameter["Young's Modulus"] == 5

    parameter = Dict{String,Any}("Material Model" => "Bond-based Elastic",
                                 "Young's Modulus" => 5)
    PeriLab.Solver_Manager.Material_Basis.get_all_elastic_moduli(test_data_manager,
                                                                 parameter)
    @test parameter["Bulk Modulus"] == 5
    @test parameter["Shear Modulus"] == 1.875
    @test parameter["Poisson's Ratio"] == Float64(1 / 3)
    @test parameter["Young's Modulus"] == 5

    parameter = Dict{String,Any}("Material Model" => "Bond-based Elastic",
                                 "Young's Modulus" => 5,
                                 "Poisson's Ratio" => 0.125)
    PeriLab.Solver_Manager.Material_Basis.get_all_elastic_moduli(test_data_manager,
                                                                 parameter)
    @test parameter["Bulk Modulus"] == 5
    @test parameter["Shear Modulus"] == 1.875
    @test parameter["Poisson's Ratio"] == Float64(1 / 3)
    @test parameter["Young's Modulus"] == 5

    test_data_manager.create_constant_node_field("Bulk_Modulus", Float64, 1, 10)
    parameter = Dict{String,Any}("Material Model" => "PD Solid Elastic",
                                 "Shear Modulus" => 10)
    PeriLab.Solver_Manager.Material_Basis.get_all_elastic_moduli(test_data_manager,
                                                                 parameter)
    @test parameter["Young's Modulus"] == [22.5, 22.5, 22.5]
    @test parameter["Poisson's Ratio"] == [0.125, 0.125, 0.125]
    @test parameter["Bulk Modulus"] == [10, 10, 10]
    @test parameter["Shear Modulus"] == [10, 10, 10]

    parameter = Dict{String,Any}("Material Model" => "Unified Bond-based Elastic",
                                 "Young's Modulus" => 5,
                                 "Poisson's Ratio" => 0.125)
    PeriLab.Solver_Manager.Material_Basis.get_all_elastic_moduli(test_data_manager,
                                                                 parameter)

    @test parameter["Young's Modulus"] == [22.5, 22.5, 22.5]
    @test parameter["Poisson's Ratio"] == [0.125, 0.125, 0.125]
    @test parameter["Bulk Modulus"] == [10, 10, 10]
    @test parameter["Shear Modulus"] == [10, 10, 10]

    parameter = Dict{String,Any}("Material Model" => "PD Solid Elastic",
                                 "Symmetry" => "Anisotropic",
                                 "C11" => 5)
    @test isnothing(PeriLab.Solver_Manager.Material_Basis.get_all_elastic_moduli(test_data_manager,
                                                                                 parameter))

    parameter = Dict{String,Any}("Material Model" => "PD Solid Elastic",
                                 "Symmetry" => "Orthotropic",
                                 "Young's Modulus X" => 5)
    @test isnothing(PeriLab.Solver_Manager.Material_Basis.get_all_elastic_moduli(test_data_manager,
                                                                                 parameter))
end

@testset "get_Hooke_matrix" begin
    test_data_manager = PeriLab.Data_Manager
    test_data_manager.initialize_data()
    parameter = Dict{String,Any}("Material Model" => "PD Solid Elastic",
                                 "Bulk Modulus" => 5,
                                 "Shear Modulus" => 1.25,
                                 "Poisson's Ratio" => 0.2,
                                 "Compute_Hook" => false)
    PeriLab.Solver_Manager.Material_Basis.get_all_elastic_moduli(test_data_manager,
                                                                 parameter)

    symmetry = "isotropic"
    E = parameter["Young's Modulus"]
    nu = parameter["Poisson's Ratio"]
    temp = 1 / ((1 + nu) * (1 - 2 * nu))
    C = PeriLab.Solver_Manager.Material_Basis.get_Hooke_matrix(test_data_manager, parameter,
                                                               symmetry, 3)
    for iID in 1:3
        @test isapprox(C[iID, iID], E * (1 - nu) * temp)
        @test C[iID+3, iID+3] == (1 - 2 * nu) * temp * E
        for jID in 1:3
            if iID != jID
                @test isapprox(C[iID, jID], E * nu * temp)
            end
        end
    end

    symmetry = "isotropic plane strain"
    C2D = PeriLab.Solver_Manager.Material_Basis.get_Hooke_matrix(test_data_manager,
                                                                 parameter, symmetry, 2)
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
    C2D = PeriLab.Solver_Manager.Material_Basis.get_Hooke_matrix(test_data_manager,
                                                                 parameter, symmetry, 2)
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
    C = PeriLab.Solver_Manager.Material_Basis.get_Hooke_matrix(test_data_manager, parameter,
                                                               symmetry, 2)
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
    @test isnothing(PeriLab.Solver_Manager.Material_Basis.get_Hooke_matrix(test_data_manager,
                                                                           parameter,
                                                                           symmetry, 2))

    symmetry = "anisotropic"
    C = PeriLab.Solver_Manager.Material_Basis.get_Hooke_matrix(test_data_manager, parameter,
                                                               symmetry, 3)
    for iID in 1:6
        for jID in 1:6
            @test C[iID, jID] == C[jID, iID]
            if jID >= iID
                @test C[iID, jID] == parameter["C"*string(iID)*string(jID)]
            end
        end
    end

    symmetry = "anisotropic plane strain"
    C = PeriLab.Solver_Manager.Material_Basis.get_Hooke_matrix(test_data_manager, parameter,
                                                               symmetry, 2)
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

    # symmetry = "anisotropic plane stress"
    # C = PeriLab.Solver_Manager.Material_Basis.get_Hooke_matrix(test_data_manager, parameter, symmetry, 2)
    # for iID = 1:3
    #     for jID = 1:3
    #         if C2D_test[iID, jID] != 0
    #             @test C[iID, jID] / C2D_test[iID, jID] - 1 < 1e-7
    #         end
    #     end
    # end

    #TODO: Check above

    symmetry = "Orthotropic"
    parameter = Dict{String,Any}("Material Model" => "PD Solid Elastic",
                                 "Young's Modulus X" => 5,
                                 "Young's Modulus Y" => 6,
                                 "Young's Modulus Z" => 7,
                                 "Poisson's Ratio XY" => 0.1,
                                 "Poisson's Ratio YZ" => 0.2,
                                 "Poisson's Ratio XZ" => 0.3,
                                 "Shear Modulus XY" => 1,
                                 "Shear Modulus YZ" => 2,
                                 "Shear Modulus XZ" => 3,
                                 "Compute_Hook" => true)
    C = PeriLab.Solver_Manager.Material_Basis.get_Hooke_matrix(test_data_manager, parameter,
                                                               symmetry, 3)
    @test C[1, 1] == 5.9692770078477215
    @test C[1, 2] == 1.2773417932876945
    @test C[1, 3] == 2.8051427617298383
    @test C[2, 1] == 1.277341793287694
    @test C[2, 2] == 6.567039572549674
    @test C[2, 3] == 2.068792786775756
    @test C[3, 1] == 2.8051427617298383
    @test C[3, 2] == 2.068792786775756
    @test C[3, 3] == 8.660878276840876
    @test C[4, 4] == 4
    @test C[5, 5] == 6
    @test C[6, 6] == 2
    E = 7000
    nu = 0.3
    G = E / (2 * (1 + nu))
    symmetry = "Orthotropic"
    parameter = Dict{String,Any}("Material Model" => "PD Solid Elastic",
                                 "Young's Modulus X" => E,
                                 "Young's Modulus Y" => E,
                                 "Young's Modulus Z" => E,
                                 "Poisson's Ratio XY" => nu,
                                 "Poisson's Ratio YZ" => nu,
                                 "Poisson's Ratio XZ" => nu,
                                 "Shear Modulus XY" => G,
                                 "Shear Modulus YZ" => G,
                                 "Shear Modulus XZ" => G,
                                 "Compute_Hook" => true)
    C = PeriLab.Solver_Manager.Material_Basis.get_Hooke_matrix(test_data_manager, parameter,
                                                               symmetry, 3)

    symmetry = "isotropic"
    parameter = Dict{String,Any}("Material Model" => "PD Solid Elastic",
                                 "Young's Modulus" => E,
                                 "Poisson's Ratio" => nu,
                                 "Shear Modulus" => G,
                                 "Compute_Hook" => true)
    @test C == PeriLab.Solver_Manager.Material_Basis.get_Hooke_matrix(test_data_manager,
                                                                 parameter, symmetry, 3)
end

@testset "ut_compute_Piola_Kirchhoff_stress" begin
    stress = [1.0 0.0; 0.0 1.0]
    deformation_gradient = [2.0 0.0; 0.0 2.0]
    expected_result = [2.0 0.0; 0.0 2.0]
    result = zeros(2, 2)
    PeriLab.Solver_Manager.Material_Basis.compute_Piola_Kirchhoff_stress!(result, stress,
                                                                          deformation_gradient)
    @test isapprox(result, expected_result)
end

@testset "ut_matrix_to_voigt" begin
    matrix = Matrix{Float64}([1 2; 3 4])
    voigt = matrix_to_voigt(matrix)
    @test voigt[1] == 1
    @test voigt[2] == 4
    @test voigt[3] == 2.5
    matrix = Matrix{Float64}([1 2 3; 4 5 6; 7 8 9])
    voigt = matrix_to_voigt(matrix)
    @test voigt[1] == 1
    @test voigt[2] == 5
    @test voigt[3] == 9
    @test voigt[4] == 7
    @test voigt[5] == 5
    @test voigt[6] == 3
    matrix = Matrix{Float64}([1 2 3 3; 4 5 6 3; 7 8 9 3])
    @test isnothing(matrix_to_voigt(matrix))
end

@testset "ut_voigt_to_matrix" begin
    @test isnothing(voigt_to_matrix([1, 2.2]))
end
