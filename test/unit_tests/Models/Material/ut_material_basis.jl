# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

#using Test
#include("../../../../src/PeriLab.jl")
#using .PeriLab

@testset "ut_init_local_damping_due_to_damage" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(3)
    nn = PeriLab.Data_Manager.create_constant_node_scalar_field("Number of Neighbors",
                                                                Int64)
    nn .= 2
    @test_logs (:error,
                "Representative Young's modulus is missing.") @test_throws PeriLab.PeriLabError begin
        PeriLab.Solver_Manager.Material_Basis.init_local_damping_due_to_damage(collect(1:2),
                                                                               "3D",
                                                                               Dict("Local Damping" =>
                                                                                        Dict()))
    end
    @test_logs (:error, "Damping coefficient is missing.") @test_throws PeriLab.PeriLabError begin
        PeriLab.Solver_Manager.Material_Basis.init_local_damping_due_to_damage(collect(1:2),
                                                                               "3D",
                                                                               Dict("Local Damping" =>
                                                                                        Dict("Representative Young's modulus" =>
                                                                                                 0)))
    end
end

@testset "ut_apply_pointwise_E" begin
    nodes = 2:3

    bond_force = [[ones(2), ones(2)], [ones(2), ones(2)], [ones(2), ones(2)]]
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

    bond_force = [[ones(2), ones(2)], [ones(2), ones(2)], [ones(2), ones(2), ones(2)]]
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
@testset "ut_compute_Piola_Kirchhoff_stress" begin
    stress = [1.0 0.0; 0.0 1.0]
    deformation_gradient = [2.0 0.0; 0.0 2.0]
    expected_result = [2.0 0.0; 0.0 2.0]
    result = zeros(2, 2)
    PeriLab.Solver_Manager.Material_Basis.compute_Piola_Kirchhoff_stress!(result, stress,
                                                                          deformation_gradient)
    @test isapprox(result, expected_result)
end
