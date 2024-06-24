# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../../../../src/Physics/Material/Material_Models/Bond_Associated_Correspondence.jl")
# include("../../../../../src/Core/data_manager.jl")
using Test
using TimerOutputs

const to = TimerOutput()

@testset "correspondence_name" begin
    @test Bond_Associated_Correspondence.correspondence_name() == "Correspondence Bond-Associated"
end
@testset "ut_calculate_Q" begin

    accuracy_order = 1
    dof = 2
    undeformed_bond = [1.0, 2.0]
    horizon = 3.0
    expected_Q = [1.0 / 3.0, 2.0 / 3.0]

    Q = Bond_Associated_Correspondence.calculate_Q(accuracy_order, dof, undeformed_bond, horizon)

    @test isapprox(Q, expected_Q)

    accuracy_order = 2
    dof = 2
    undeformed_bond = [1.0, 2.0]
    horizon = 3.0
    expected_Q = [0.3333333333333333, 0.6666666666666666, 0.1111111111111111, 0.2222222222222222, 0.4444444444444444]

    Q = Bond_Associated_Correspondence.calculate_Q(accuracy_order, dof, undeformed_bond, horizon)

    @test isapprox(Q, expected_Q)

    accuracy_order = 2
    dof = 3
    undeformed_bond = [1.0, 2.0, 5]
    horizon = 3.0
    expected_Q = [0.3333333333333333, 0.6666666666666666, 1.6666666666666667, 0.1111111111111111, 0.2222222222222222, 0.5555555555555556, 0.4444444444444444, 1.1111111111111112, 2.777777777777778]

    Q = Bond_Associated_Correspondence.calculate_Q(accuracy_order, dof, undeformed_bond, horizon)

    @test isapprox(Q, expected_Q)

    undeformed_bond = [1.0, 0.0, 0]
    @test Bond_Associated_Correspondence.calculate_Q(1, dof, undeformed_bond, 1.0) == [1, 0, 0]
    undeformed_bond = [0.0, 1.0, 0]
    @test Bond_Associated_Correspondence.calculate_Q(1, dof, undeformed_bond, 1.0) == [0, 1, 0]
    undeformed_bond = [0.0, 0.0, 1.0]
    @test Bond_Associated_Correspondence.calculate_Q(1, dof, undeformed_bond, 1.0) == [0, 0, 1]
end
@testset "ut_compute_Lagrangian_gradient_weights" begin
    nodes = [1, 2]
    dof = 3
    qdim = 9
    accuracy_order = 2
    nlist = [[2, 3, 5, 6, 7], [1, 3, 4, 5, 6, 7]]# tbd
    horizon = [1.0, 1.0]
    bond_damage = [[0.9, 0.8, 1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]]
    omega = [[1.0, 1.0, 1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]]
    undeformed_bond = [[1.0 0.0 0.0; 0.0 1.0 0.0; 1.0 2.0 5; 0.0 1.0 1.0; 2.0 1.0 1.0], [-6.0 8.5 2; 6.0 1.5 2; 6.0 2.5 5; -0.9 0.8 0.1; 1.0 -0.9 0.2; 0.0 -0.9 0.2]]
    volume = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    gradient_weights = [[0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0], [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]]

    test = Bond_Associated_Correspondence.compute_Lagrangian_gradient_weights(nodes, dof, qdim, accuracy_order, volume, nlist, horizon, bond_damage, omega, undeformed_bond, gradient_weights)

    @test isapprox(test[1][:], [2.28125, 0.5, 0.5, -0.25, -0.625, -0.125, -0.052109066205533794, -0.4375, 1.455078125, 0.5, 0.15625, -2.125, -0.125, 1.875, 0.0625])
    @test isapprox(test[2][:], [0.0012369601821692333, -0.03382323766941053, -0.006481215557232245, -0.20526883154875758, 0.2312955695281001, -0.2774482376609969, -0.006243107506445966, -0.018203513750354894, -0.0023387447599407807, 0.38834658871731265, 0.06936211266128167, -0.376723321917216, -0.025025064685230802, 0.007650389750137876, -0.014050789610079661, 0.6916462673154423, 0.7237356806277135, 0.3279228446136444])
end

@testset "ut_find_local_neighbors" begin
    coordinates = zeros(Float64, 5, 2)
    coordinates[1, 1] = 0
    coordinates[1, 2] = 0
    coordinates[2, 1] = 1
    coordinates[2, 2] = 0
    coordinates[3, 1] = 0
    coordinates[3, 2] = 1
    coordinates[4, 1] = 1
    coordinates[4, 2] = 1
    coordinates[5, 1] = 2
    coordinates[5, 2] = 2

    nlist = [2, 3, 4, 5]

    bond_horizon::Float64 = 1

    @test Bond_Associated_Correspondence.find_local_neighbors(5, coordinates, nlist, bond_horizon) == []
    bond_horizon = 2.6
    @test Bond_Associated_Correspondence.find_local_neighbors(5, coordinates, nlist, bond_horizon) == [2, 3, 4]
    @test Bond_Associated_Correspondence.find_local_neighbors(2, coordinates, nlist, bond_horizon) == [3, 4, 5]
end

@testset "ut_compute_bond_associated_weighted_volume" begin
    # Initialize input arrays
    nodes = [1, 2, 3]
    nlist = [[2, 3], [1, 3], [1, 2]]
    coordinates = [0.0 0.0; 1.0 0.0; 0.0 1.0]
    bond_damage = [[1.0, 1.0], [0.0, 1.0], [1.0, 1.0]]
    omega = [[0.5, 0.8], [0.7, 0.9], [0.6, 0.4]]
    volume = [1.0, 2.0, 3.0]
    weighted_volume = [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
    bond_horizon = 2.2
    # Call the function under test
    weighted_volume = Bond_Associated_Correspondence.compute_bond_associated_weighted_volume(nodes, nlist, coordinates, bond_damage, omega, volume, bond_horizon, weighted_volume)
    # Check the expected output
    @test weighted_volume ≈ [[0.7058823529411765, 0.2941176470588235], [1.0, 0.0], [0.5714285714285715, 0.4285714285714286]]
end