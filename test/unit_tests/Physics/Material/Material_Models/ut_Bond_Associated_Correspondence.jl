# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../../../../src/Physics/Material/Material_Models/Bond_Associated_Correspondence.jl")
# include("../../../../../src/Core/data_manager.jl")
using Test
#include("../../../../../src/PeriLab.jl")
#using .PeriLab


@testset "ut_correspondence_name" begin
    @test Bond_Associated_Correspondence.correspondence_name() == "Correspondence Bond-Associated"
end

@testset "ut_compute_bond_strain" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.clear_data_manager()
    test_data_manager.set_num_controller(2)
    test_data_manager.set_dof(2)

    nn = test_data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1] = 3
    nn[2] = 1
    nodes = [1, 2]

    nlist = test_data_manager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
    nlist[1] = [2, 3]
    nlist[2] = [1]


    deformation_gradient = test_data_manager.create_constant_bond_field("Deformation Gradient", Float64, "Matrix", 2)
    strain = test_data_manager.create_constant_bond_field("Strain", Float64, "Matrix", 2)

    deformation_gradient[1][1, :, :] = [1 0; 0 1]
    deformation_gradient[1][2, :, :] = [1 1; 1 1]

    deformation_gradient[2][1, :, :] = [0 1; 0 1]

    strain = Bond_Associated_Correspondence.compute_bond_strain(nodes, nlist, deformation_gradient, strain)

    @test strain[1][1, :, :] == [0 0; 0 0]
    @test strain[1][2, :, :] == [0.5 1.0; 1.0 0.5]
    @test strain[2][1, :, :] == [-0.5 0; 0 0.5]

end

@testset "ut_update_Green_Langrange_strain" begin
    dt = 0.1
    deformation_gradient = [1.0 0.2 0.0; 0.1 1.0 0.3; 0.0 0.1 1.0]
    deformation_gradient_dot = [0.05 0.01 0.0; 0.02 0.05 0.01; 0.0 0.02 0.05]
    expected_strain = 0.5 * dt * (deformation_gradient * deformation_gradient_dot + (deformation_gradient * deformation_gradient_dot)')
    computed_strain = Bond_Associated_Correspondence.update_Green_Langrange_strain(dt, deformation_gradient, deformation_gradient_dot)
    @test isapprox(computed_strain, expected_strain)
    deformation_gradient = zeros(3, 3)
    deformation_gradient_dot = zeros(3, 3)
    expected_strain = zeros(3, 3)
    computed_strain = Bond_Associated_Correspondence.update_Green_Langrange_strain(dt, deformation_gradient, deformation_gradient_dot)
    @test computed_strain == expected_strain
end
@testset "ut_compute_weighted_volume" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.clear_data_manager()
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
    weighted_volume = test_data_manager.create_constant_node_field("Weighted Volume", Float64, 1)


    result = Bond_Associated_Correspondence.compute_weighted_volume(nodes, nlist, volume, bond_damage, omega, weighted_volume)

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
    undeformed_bond = [[1.0 0.0; 0.0 1.0; 1.0 2.0; 0.2 1.0; 2.0 1.0], [-6.0 8.5; 6.0 1.5; 6.0 2.5; -0.9 0.8; 1.0 -0.9; 0.0 -0.9]]
    volume = [1.0, 2.0, 1.0, 1.0, 1.0, 4.0, 1.0, 1.0]
    gradient_weights = [[0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0], [0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0]]

    gradient_weights = Bond_Associated_Correspondence.compute_Lagrangian_gradient_weights(nodes, dof, accuracy_order, volume, nlist, horizon, bond_damage, omega, undeformed_bond, gradient_weights)

    @test isapprox(gradient_weights[1][1, :], [0.5555555555555571, -0.2777777777777777])
    @test isapprox(gradient_weights[1][2, :], [-3.125000000000102, -1.2500000000000462])
    @test isapprox(gradient_weights[1][3, :], [7.105427357601002e-15, -0.5])
    @test isapprox(gradient_weights[1][4, :], [0.6944444444444704, 0.6944444444444562])
    @test isapprox(gradient_weights[1][5, :], [-0.2777777777777821, 0.22222222222221788])

    @test isapprox(gradient_weights[2][1, :], [0.0008338678006687417, 0.010443292032181972])
    @test isapprox(gradient_weights[2][2, :], [-0.11142760668950946, -0.06982123424402786])
    @test isapprox(gradient_weights[2][3, :], [0.09771932233446412, 0.06511412698788455])
    @test isapprox(gradient_weights[2][4, :], [-0.36478006840748456, -0.17392108050722965])
    @test isapprox(gradient_weights[2][5, :], [0.189737712841881, -0.016406644181643948])
end
@testset "ut_init_Bond-Associated" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.clear_data_manager()
    test_data_manager.set_num_controller(2)
    test_data_manager.set_dof(3)
    nn = test_data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn .= 1
    nodes = Vector{Int64}(1:2)

    @test isnothing(Bond_Associated_Correspondence.init_material_model(test_data_manager, nodes, Dict()))

    material_parameter = Dict{String,Any}("Symmetry" => "isotropic")
    test_data_manager = Bond_Associated_Correspondence.init_material_model(test_data_manager, nodes, material_parameter)
    @test haskey(material_parameter, "Accuracy Order")
    @test material_parameter["Accuracy Order"] == 1
    @test "Bond StrainN" in test_data_manager.get_all_field_keys()
    @test "Bond StrainNP1" in test_data_manager.get_all_field_keys()
    @test "Bond Cauchy StressN" in test_data_manager.get_all_field_keys()
    @test "Bond Cauchy StressNP1" in test_data_manager.get_all_field_keys()
    @test "Bond Strain Increment" in test_data_manager.get_all_field_keys()
    @test "Weighted Volume" in test_data_manager.get_all_field_keys()
    @test "Lagrangian Gradient Weights" in test_data_manager.get_all_field_keys()
    @test "Integral Nodal Stress" in test_data_manager.get_all_field_keys()
    material_parameter = Dict("Symmetry" => "isotropic", "Accuracy Order" => 2)
    test_data_manager = Bond_Associated_Correspondence.init_material_model(test_data_manager, nodes, material_parameter)

    @test material_parameter["Accuracy Order"] == 2

    @test isnothing(Bond_Associated_Correspondence.init_material_model(test_data_manager, nodes, Dict("Symmetry" => "isotropic", "Accuracy Order" => 0)))

end
@testset "ut_compute_stress_integral" begin
    test_data_manager = PeriLab.Data_manager
    dof = 2
    test_data_manager.clear_data_manager()
    test_data_manager.set_num_controller(2)
    test_data_manager.set_dof(dof)
    nn = test_data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn .= 1
    nlist = test_data_manager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
    nlist[1][:] = [2]
    nlist[2][:] = [1]
    nodes = Vector{Int64}(1:2)
    omega = test_data_manager.create_constant_bond_field("Influence Function", Float64, 1)

    omega[1][:] = [1.0]
    omega[2][:] = [2.0]
    bond_damage = test_data_manager.create_constant_bond_field("Bond Damage", Float64, 1)
    bond_damage[1][:] = [1.0]
    bond_damage[2][:] = [1.0]
    volume = test_data_manager.create_constant_node_field("Volume", Float64, 1)
    volume .= 1
    weighted_volume = test_data_manager.create_constant_node_field("Weighted Volume", Float64, 1)
    weighted_volume .= 1

    bond_geometry = test_data_manager.create_constant_bond_field("Bond Geometry", Float64, dof)
    bond_geometry[1][:] = [1.0, 0.0]
    bond_geometry[2][:] = [-1.0, 0.0]

    bond_length = test_data_manager.create_constant_bond_field("Bond Length", Float64, 1)
    bond_length[1][:] = [1.0]
    bond_length[2][:] = [1.0]

    deformation_gradient = test_data_manager.create_constant_bond_field("Bond Deformation Gradient", Float64, "Matrix", dof)
    deformation_gradient[1][1, :, :] = [1.0 0.0; 0.0 1.0]
    deformation_gradient[2][1, :, :] = [1.0 0.0; 0.0 1.0]

    bond_stresses = test_data_manager.create_constant_bond_field("Bond Cauchy Stress", Float64, "Matrix", dof)
    bond_stresses[1][1, :, :] = [1.0 0.0; 0.0 1.0]
    bond_stresses[2][1, :, :] = [1.0 0.0; 0.0 1.0]

    stress_integral = test_data_manager.create_constant_node_field("Stress Integral", Float64, "Matrix", dof)

    stress_integral = Bond_Associated_Correspondence.compute_stress_integral(nodes, dof, nlist, omega, bond_damage, volume, weighted_volume, bond_geometry, bond_length, deformation_gradient, bond_stresses, stress_integral)

    expected_stress_integral = [
        [0.0 0.0; 0.0 1.0],
        [0.0 0.0; 0.0 2.0]
    ]
    @test isapprox(stress_integral[1, :, :], expected_stress_integral[1][:, :])
    @test isapprox(stress_integral[2, :, :], expected_stress_integral[2][:, :])
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
    @test isapprox(weighted_volume, [[0.7058823529411765, 0.2941176470588235], [1.0, 0.0], [0.5714285714285715, 0.4285714285714286]])
end