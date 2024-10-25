# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
include("../../../../src/Models/Pre_calculation/pre_bond_associated_correspondence.jl")
#include("../../../../src/PeriLab.jl")
#using .PeriLab


@testset "ut_calculate_Q" begin

    accuracy_order = 1
    dof = 2
    undeformed_bond = [1.0, 2.0]
    horizon = 3.0
    expected_Q = [1.0 / 3.0, 2.0 / 3.0]
    Q = zeros(2)
    Q = Pre_Bond_Associated_Correspondence.calculate_Q(
        accuracy_order,
        dof,
        undeformed_bond,
        horizon,
        Q,
    )

    @test isapprox(Q, expected_Q)

    accuracy_order = 2
    dof = 2
    undeformed_bond = [1.0, 2.0]
    horizon = 3.0
    expected_Q = [
        0.3333333333333333,
        0.6666666666666666,
        0.1111111111111111,
        0.2222222222222222,
        0.4444444444444444,
    ]
    Q = zeros(5)
    Q = Pre_Bond_Associated_Correspondence.calculate_Q(
        accuracy_order,
        dof,
        undeformed_bond,
        horizon,
        Q,
    )

    @test isapprox(Q, expected_Q)

    accuracy_order = 2
    dof = 3
    undeformed_bond = [1.0, 2.0, 5]
    horizon = 3.0
    expected_Q = [
        0.3333333333333333,
        0.6666666666666666,
        1.6666666666666667,
        0.1111111111111111,
        0.2222222222222222,
        0.5555555555555556,
        0.4444444444444444,
        1.1111111111111112,
        2.777777777777778,
    ]
    Q = zeros(9)
    Q = Pre_Bond_Associated_Correspondence.calculate_Q(
        accuracy_order,
        dof,
        undeformed_bond,
        horizon,
        Q,
    )

    @test isapprox(Q, expected_Q)
    Q = zeros(3)
    undeformed_bond = [1.0, 0.0, 0]
    @test Pre_Bond_Associated_Correspondence.calculate_Q(1, dof, undeformed_bond, 1.0, Q) ==
          [1, 0, 0]
    undeformed_bond = [0.0, 1.0, 0]
    @test Pre_Bond_Associated_Correspondence.calculate_Q(1, dof, undeformed_bond, 1.0, Q) ==
          [0, 1, 0]
    undeformed_bond = [0.0, 0.0, 1.0]
    @test Pre_Bond_Associated_Correspondence.calculate_Q(1, dof, undeformed_bond, 1.0, Q) ==
          [0, 0, 1]
end


@testset "ut_compute_weighted_volume" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.initialize_data()
    test_data_manager.set_num_controller(4)
    test_data_manager.set_dof(3)

    volume = test_data_manager.create_constant_node_field("Volume", Float64, 1)
    volume[:] = [1.0, 2.0, 3.0, 4.0]

    nn = test_data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1] = 3
    nn[2:4] .= 1
    nodes = [1]

    nlist = test_data_manager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
    nlist[1] = [2, 3, 4]
    nlist[2] = [1]
    nlist[3] = [1]
    nlist[4] = [1]

    bond_damage = test_data_manager.create_constant_bond_field("Bond Damage", Float64, 1)
    bond_damage[1][:] .= 1
    omega = test_data_manager.create_constant_bond_field("Influence Function", Float64, 1)
    omega[1][:] = [1.0, 0.8, 0.6]
    weighted_volume =
        test_data_manager.create_constant_node_field("Weighted Volume", Float64, 1)


    result = Pre_Bond_Associated_Correspondence.compute_weighted_volume(
        nodes,
        nlist,
        volume,
        bond_damage,
        omega,
        weighted_volume,
    )

    expected_weighted_volume = sum(bond_damage[1][:] .* omega[1][:] .* volume[nlist[1]])
    @test isapprox(result[1], expected_weighted_volume)

end

@testset "ut_compute_Lagrangian_gradient_weights" begin
    nodes = [1, 2]
    dof = 2
    accuracy_order = 2
    nlist = [[2, 3, 5, 6, 7], [1, 3, 4, 5, 6, 7]]# tbd
    horizon = [1.0, 1.0]
    bond_damage = [[0.9, 0.8, 1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]]
    omega = [[1.0, 1.0, 1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]]
    undeformed_bond = [
        [1.0 0.0; 0.0 1.0; 1.0 2.0; 0.2 1.0; 2.0 1.0],
        [-6.0 8.5; 6.0 1.5; 6.0 2.5; -0.9 0.8; 1.0 -0.9; 0.0 -0.9],
    ]
    volume = [1.0, 2.0, 1.0, 1.0, 1.0, 4.0, 1.0, 1.0]
    gradient_weights = [
        [0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0],
        [0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0],
    ]
    Q = zeros(5)
    Minv = zeros(5, 5)
    M = zeros(5, 5)

    Pre_Bond_Associated_Correspondence.compute_Lagrangian_gradient_weights(
        nodes,
        dof,
        accuracy_order,
        volume,
        nlist,
        horizon,
        bond_damage,
        omega,
        Q,
        M,
        Minv,
        undeformed_bond,
        gradient_weights,
    )

    @test isapprox(gradient_weights[1][1, :], [0.5000000000000014, -0.24999999999999992])
    @test isapprox(gradient_weights[1][2, :], [-2.5000000000000817, -1.000000000000037])
    @test isapprox(gradient_weights[1][3, :], [7.105427357601002e-15, -0.5])
    @test isapprox(gradient_weights[1][4, :], [0.6944444444444704, 0.6944444444444562])
    @test isapprox(gradient_weights[1][5, :], [-0.2777777777777821, 0.22222222222221788])

    @test isapprox(gradient_weights[2][1, :], [0.0008338678006687417, 0.010443292032181972])
    @test isapprox(gradient_weights[2][2, :], [-0.11142760668950946, -0.06982123424402786])
    @test isapprox(gradient_weights[2][3, :], [0.09771932233446412, 0.06511412698788455])
    @test isapprox(gradient_weights[2][4, :], [-0.36478006840748456, -0.17392108050722965])
    @test isapprox(gradient_weights[2][5, :], [0.189737712841881, -0.016406644181643948])
end
