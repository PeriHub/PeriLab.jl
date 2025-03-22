# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../../../../../src/Models/Material/Material_Models/Correspondence/Bond_Associated_Correspondence.jl")
using Test
# include("../../../../../../src/PeriLab.jl")
# using .PeriLab

@testset "ut_compute_bond_strain" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.initialize_data()
    test_data_manager.set_num_controller(2)
    test_data_manager.set_dof(2)

    nn = test_data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1] = 1
    nn[2] = 1
    nodes = [1, 2]

    nlist = test_data_manager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
    nlist[1] = [2]
    nlist[2] = [1]

    deformation_gradient = test_data_manager.create_constant_bond_field("Deformation Gradient",
                                                                        Float64,
                                                                        "Matrix",
                                                                        2)
    strain, strainN = test_data_manager.create_bond_field("Strain", Float64, "Matrix", 2)
    strain_inc = test_data_manager.create_constant_bond_field("Strain Increment",
                                                              Float64,
                                                              "Matrix",
                                                              2)
    deformation_gradient[1][1, :, :] = [1 0; 0 1]
    deformation_gradient[2][1, :, :] = [0 1; 0 1]

    strainN[1][1, :, :] = [-1 -1; -1 -1]
    strainN[2][1, :, :] = [0.5 0; 0 0.5]
    Bond_Associated_Correspondence.compute_bond_strain(nodes,
                                                       nlist,
                                                       deformation_gradient,
                                                       strain,
                                                       strainN,
                                                       strain_inc)

    @test strain[1][1, :, :] == [0 0; 0 0]
    @test strain[2][1, :, :] == [-0.5 0; 0 0.5]
    @test strain_inc[1][1, :, :] == [1 1; 1 1]
    @test strain_inc[2][1, :, :] == [-1.0 0.0; 0.0 0.0]
end

@testset "ut_update_Green_Langrange_strain" begin
    dt = 0.1
    deformation_gradient = [1.0 0.2 0.0; 0.1 1.0 0.3; 0.0 0.1 1.0]
    deformation_gradient_dot = [0.05 0.01 0.0; 0.02 0.05 0.01; 0.0 0.02 0.05]
    expected_strain = 0.5 .* dt .* (deformation_gradient * deformation_gradient_dot +
                       (deformation_gradient * deformation_gradient_dot)')
    computed_strain = zeros(3, 3)
    Bond_Associated_Correspondence.update_Green_Langrange_strain(dt,
                                                                 deformation_gradient,
                                                                 deformation_gradient_dot,
                                                                 computed_strain)

    # for i in 1:3
    #     for j in 1:3
    #         @test isapprox(computed_strain[i, j], expected_strain[i, j])
    #     end
    # end
    @test computed_strain == expected_strain
    deformation_gradient = zeros(3, 3)
    deformation_gradient_dot = zeros(3, 3)
    expected_strain = zeros(3, 3)
    computed_strain = zeros(3, 3)
    Bond_Associated_Correspondence.update_Green_Langrange_strain(dt,
                                                                 deformation_gradient,
                                                                 deformation_gradient_dot,
                                                                 computed_strain)
    @test computed_strain == expected_strain
end
@testset "ut_init_Bond-Associated" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.initialize_data()
    test_data_manager.set_num_controller(2)
    test_data_manager.set_dof(3)
    nn = test_data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[:] .= 1
    nlist = test_data_manager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
    nlist[1] = [2]
    nlist[2] = [1]
    nodes = Vector{Int64}(1:2)
    test_data_manager.create_constant_node_field("Volume", Float64, 1)
    test_data_manager.create_constant_bond_field("Influence Function", Float64, 1)
    test_data_manager.create_bond_field("Bond Damage", Float64, 1)

    @test isnothing(Bond_Associated_Correspondence.init_model(test_data_manager, nodes,
                                                              Dict()))

    material_parameter = Dict{String,Any}("Symmetry" => "isotropic")
    test_data_manager = Bond_Associated_Correspondence.init_model(test_data_manager,
                                                                  nodes,
                                                                  material_parameter)

    @test test_data_manager.get_accuracy_order() == 1

    material_parameter = Dict("Symmetry" => "isotropic", "Accuracy Order" => 2)
    test_data_manager = Bond_Associated_Correspondence.init_model(test_data_manager,
                                                                  nodes,
                                                                  material_parameter)

    @test test_data_manager.get_accuracy_order() == 2
end
@testset "ut_compute_stress_integral" begin
    test_data_manager = PeriLab.Data_manager
    dof = 2
    test_data_manager.initialize_data()
    test_data_manager.set_num_controller(2)
    test_data_manager.set_dof(dof)
    nn = test_data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn = test_data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[:] .= 1
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
    weighted_volume = test_data_manager.create_constant_node_field("Weighted Volume",
                                                                   Float64, 1)
    weighted_volume .= 1

    bond_geometry = test_data_manager.create_constant_bond_field("Bond Geometry", Float64,
                                                                 dof)
    bond_geometry[1][1] = [1.0, 0.0]
    bond_geometry[2][1] = [-1.0, 0.0]

    bond_length = test_data_manager.create_constant_bond_field("Bond Length", Float64, 1, 1)

    deformation_gradient = test_data_manager.create_constant_bond_field("Bond Deformation Gradient",
                                                                        Float64,
                                                                        "Matrix",
                                                                        dof)
    deformation_gradient[1][1, :, :] = [1 0; 0 1]
    deformation_gradient[2][1, :, :] = [1 0; 0 1]

    bond_stresses = test_data_manager.create_constant_bond_field("Bond Cauchy Stress",
                                                                 Float64,
                                                                 "Matrix",
                                                                 dof)
    bond_stresses[1][1, :, :] = [1.0 0.0; 0.0 1.0]
    bond_stresses[2][1, :, :] = [1.0 0.0; 0.0 1.0]

    stress_integral = test_data_manager.create_constant_node_field("Stress Integral",
                                                                   Float64,
                                                                   "Matrix",
                                                                   dof)

    stress_integral = Bond_Associated_Correspondence.compute_stress_integral(nodes,
                                                                             dof,
                                                                             nlist,
                                                                             omega,
                                                                             bond_damage,
                                                                             volume,
                                                                             weighted_volume,
                                                                             bond_geometry,
                                                                             bond_length,
                                                                             deformation_gradient,
                                                                             bond_stresses,
                                                                             stress_integral)

    expected_stress_integral = [[0.0 0.0; 0.0 1.0], [0.0 0.0; 0.0 2.0]]
    @test isapprox(stress_integral[1, :, :], expected_stress_integral[1][:, :])
    @test isapprox(stress_integral[2, :, :], expected_stress_integral[2][:, :])
end
