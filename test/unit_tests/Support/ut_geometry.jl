# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
#include("../../../src/PeriLab.jl")
#using .PeriLab
@testset "ut_undeformed_bond" begin
    nnodes = 4
    dof = 2
    test_Data_manager = PeriLab.Data_manager
    test_Data_manager.clear_data_manager()
    test_Data_manager.set_num_controller(nnodes)
    test_Data_manager.set_num_responder(0)
    test_Data_manager.set_dof(dof)
    nn = test_Data_manager.create_constant_node_field("Number of Neighbors", Int32, 1)
    nn .= [2, 2, 2, 1]
    coor = test_Data_manager.create_constant_node_field("Coordinates", Float64, 2)
    undeformed_bond = test_Data_manager.create_constant_bond_field("Bond Geometry", Float64, 2)
    undeformed_bond_length = test_Data_manager.create_constant_bond_field("Bond Length", Float64, 1)
    nlist = test_Data_manager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
    nlist[1] = [2, 3]
    nlist[2] = [1, 3]
    nlist[3] = [1, 2]
    nlist[4] = [2]

    coor[1, 1] = 0
    coor[1, 2] = 0
    coor[2, 1] = 0.5
    coor[2, 2] = 0.5
    coor[3, 1] = 1
    coor[3, 2] = 0
    coor[4, 1] = 0
    coor[4, 2] = 1

    u_bond, u_bond_length = PeriLab.IO.Geometry.calculate_bond_length(1, coor, nlist[1])
    @test u_bond[1, 1] == 0.5
    @test u_bond[1, 2] == 0.5
    @test isapprox(u_bond_length[1], sqrt(0.5))
    @test u_bond[2, 1] == 1
    @test u_bond[2, 2] == 0
    @test u_bond_length[2] == 1

    undeformed_bond, undeformed_bond_length = PeriLab.IO.Geometry.bond_geometry(Vector(1:nnodes), nlist, coor, undeformed_bond, undeformed_bond_length)

    @test undeformed_bond[1][1, 1] == 0.5
    @test undeformed_bond[1][1, 2] == 0.5
    @test undeformed_bond_length[1][1] / sqrt(0.5) - 1 < 1e-8
    @test undeformed_bond[1][2, 1] == 1
    @test undeformed_bond[1][2, 2] == 0
    @test undeformed_bond_length[1][2] == 1

    @test undeformed_bond[2][1, 1] == -0.5
    @test undeformed_bond[2][1, 2] == -0.5
    @test undeformed_bond_length[2][1] / sqrt(0.5) - 1 < 1e-8
    @test undeformed_bond[2][2, 1] == 0.5
    @test undeformed_bond[2][2, 2] == -0.5
    @test undeformed_bond_length[2][2] / sqrt(1.25) - 1 < 1e-8

    @test undeformed_bond[3][1, 1] == -1
    @test undeformed_bond[3][1, 2] == 0
    @test undeformed_bond_length[3][1] == 1
    @test undeformed_bond[3][2, 1] == -0.5
    @test undeformed_bond[3][2, 2] == 0.5
    @test undeformed_bond_length[3][2] / sqrt(1.25) - 1 < 1e-8

    @test undeformed_bond[4][1, 1] == 0.5
    @test undeformed_bond[4][1, 2] == -0.5
    @test undeformed_bond_length[4][1] / sqrt(1.25) - 1 < 1e-8
    undeformed_bond, undeformed_bond_length = PeriLab.IO.Geometry.bond_geometry(Vector(1:nnodes), nlist, coor, undeformed_bond, undeformed_bond_length)
    # test if a sum exists or not
    @test undeformed_bond[1][1, 1] == 0.5
    @test undeformed_bond[1][1, 2] == 0.5
    @test undeformed_bond_length[1][1] / sqrt(0.5) - 1 < 1e-8
    @test undeformed_bond[1][2, 1] == 1
    @test undeformed_bond[1][2, 2] == 0
    @test undeformed_bond_length[1][2] == 1

    @test undeformed_bond[2][1, 1] == -0.5
    @test undeformed_bond[2][1, 2] == -0.5
    @test undeformed_bond_length[2][1] / sqrt(0.5) - 1 < 1e-8
    @test undeformed_bond[2][2, 1] == 0.5
    @test undeformed_bond[2][2, 2] == -0.5
    @test undeformed_bond_length[2][2] / sqrt(1.25) - 1 < 1e-8

    @test undeformed_bond[3][1, 1] == -1
    @test undeformed_bond[3][1, 2] == 0
    @test undeformed_bond_length[3][1] == 1
    @test undeformed_bond[3][2, 1] == -0.5
    @test undeformed_bond[3][2, 2] == 0.5
    @test undeformed_bond_length[3][2] / sqrt(1.25) - 1 < 1e-8

    @test undeformed_bond[4][1, 1] == 0.5
    @test undeformed_bond[4][1, 2] == -0.5
    @test undeformed_bond_length[4][1] / sqrt(1.25) - 1 < 1e-8



    coor[:, :] .= 0
    undeformed_bond[1][:, :] .= 0
    undeformed_bond[2][:, :] .= 0
    undeformed_bond[3][:, :] .= 0
    undeformed_bond[4][:, :] .= 0

    undeformed_bond = PeriLab.IO.Geometry.bond_geometry(Vector(1:nnodes), nlist, coor, undeformed_bond, undeformed_bond_length)
    @test isnothing(undeformed_bond)
end
@testset "ut_shape_tensor_and_deformation_gradient" begin
    nnodes = 4
    dof = 2
    nodes = Vector{Int64}(1:nnodes)
    test_Data_manager = PeriLab.Data_manager
    test_Data_manager.clear_data_manager()
    test_Data_manager.set_num_controller(nnodes)
    test_Data_manager.set_dof(dof)
    nn = test_Data_manager.create_constant_node_field("Number of Neighbors", Int32, 1)
    nn .= [3, 3, 3, 3]

    coor = test_Data_manager.create_constant_node_field("Coordinates", Float64, 2)
    undeformed_bond = test_Data_manager.create_constant_bond_field("Bond Geometry", Float64, dof)
    undeformed_bond_length = test_Data_manager.create_constant_bond_field("Bond Length", Float64, 1)
    # does not have to be NP1 for testing
    deformed_bond = test_Data_manager.create_constant_bond_field("Deformed Bond Geometry", Float64, dof)
    deformed_bond_length = test_Data_manager.create_constant_bond_field("Deformed Bond Length", Float64, 1)

    nlist = test_Data_manager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
    volume = test_Data_manager.create_constant_node_field("Volume", Float64, 1)
    omega = test_Data_manager.create_constant_bond_field("Influence Function", Float64, 1)

    bond_damage = test_Data_manager.create_constant_bond_field("Bond Damage", Float64, 1)
    shape_tensor = test_Data_manager.create_constant_node_field("Shape Tensor", Float64, "Matrix", dof)
    inverse_shape_tensor = test_Data_manager.create_constant_node_field("Inverse Shape Tensor", Float64, "Matrix", dof)
    deformation_gradient = test_Data_manager.create_constant_node_field("Deformation Gradient", Float64, "Matrix", dof)
    omega[1][:] .= 1
    omega[2][:] .= 1
    omega[3][:] .= 1
    omega[4][:] .= 1
    bond_damage[1][:] .= 1
    bond_damage[2][:] .= 1
    bond_damage[3][:] .= 1
    bond_damage[4][:] .= 1
    volume[:] .= 1
    nlist[1] = [2, 3, 4]
    nlist[2] = [1, 3, 4]
    nlist[3] = [1, 2, 4]
    nlist[4] = [1, 2, 3]
    coor[1, 1] = 0
    coor[1, 2] = 0
    coor[2, 1] = 1
    coor[2, 2] = 0
    coor[3, 1] = 0
    coor[3, 2] = 0.5
    coor[4, 1] = 1
    coor[4, 2] = 0.5

    undeformed_bond, undeformed_bond_length = PeriLab.IO.Geometry.bond_geometry(Vector(1:nnodes), nlist, coor, undeformed_bond, undeformed_bond_length)
    shape_tensor, inverse_shape_tensor = PeriLab.IO.Geometry.shape_tensor(view(nodes, eachindex(nodes)), dof, nlist, volume, omega, bond_damage, undeformed_bond, shape_tensor, inverse_shape_tensor)

    deformed_coor = copy(coor)

    deformation_gradient = PeriLab.IO.Geometry.deformation_gradient(view(nodes, eachindex(nodes)), dof, nlist, volume, omega, bond_damage, undeformed_bond, undeformed_bond, inverse_shape_tensor, deformation_gradient)
    for i in 1:nnodes
        @test deformation_gradient[i, 1, 1] - 1 < 1e-7
        @test deformation_gradient[i, 1, 2] < 1e-7
        @test deformation_gradient[i, 2, 1] < 1e-7
        @test deformation_gradient[i, 2, 2] - 1 < 1e-7
    end
    deformed_coor[1, 1] = -0.25
    deformed_coor[3, 1] = -0.25
    deformed_coor[2, 1] = 0.25
    deformed_coor[4, 1] = 0.25

    deformed_bond, deformed_bond_length = PeriLab.IO.Geometry.bond_geometry(Vector(1:nnodes), nlist, deformed_coor, deformed_bond, deformed_bond_length)
    deformation_gradient = PeriLab.IO.Geometry.deformation_gradient(view(nodes, eachindex(nodes)), dof, nlist, volume, omega, bond_damage, deformed_bond, undeformed_bond, inverse_shape_tensor, deformation_gradient)

    for i in 1:nnodes
        for j in 1:nn[i]
            test_vector = deformation_gradient[i, :, :] * undeformed_bond[i][j, 1:dof] - deformed_bond[i][j, 1:dof]
            for k in 1:dof
                @test abs(test_vector[k]) < 1e-7
            end
        end
    end

    deformed_coor = copy(coor)
    deformed_coor[3, 2] = 1.5
    deformed_coor[4, 2] = 1.5

    deformed_bond, deformed_bond_length = PeriLab.IO.Geometry.bond_geometry(Vector(1:nnodes), nlist, deformed_coor, deformed_bond, deformed_bond_length)
    deformation_gradient = PeriLab.IO.Geometry.deformation_gradient(view(nodes, eachindex(nodes)), dof, nlist, volume, omega, bond_damage, deformed_bond, undeformed_bond, inverse_shape_tensor, deformation_gradient)
    for i in 1:nnodes
        for j in nn[i]
            test_vector = deformation_gradient[i, :, :] * undeformed_bond[i][j, 1:dof] - deformed_bond[i][j, 1:dof]
            for k in 1:dof
                @test abs(test_vector[k]) < 1e-7
            end
        end
    end

    deformed_coor[1, 1] = 0.0
    deformed_coor[1, 2] = 0
    deformed_coor[2, 1] = 1
    deformed_coor[2, 2] = 0
    deformed_coor[3, 1] = 0.5
    deformed_coor[3, 2] = 0.5
    deformed_coor[4, 1] = 1.5
    deformed_coor[4, 2] = 0.5

    deformed_bond, deformed_bond_length = PeriLab.IO.Geometry.bond_geometry(Vector(1:nnodes), nlist, deformed_coor, deformed_bond, deformed_bond_length)
    deformation_gradient = PeriLab.IO.Geometry.deformation_gradient(view(nodes, eachindex(nodes)), dof, nlist, volume, omega, bond_damage, deformed_bond, undeformed_bond, inverse_shape_tensor, deformation_gradient)
    for i in 1:nnodes
        for j in nn[i]
            test_vector = deformation_gradient[i, :, :] * undeformed_bond[i][j, 1:dof] - deformed_bond[i][j, 1:dof]
            for k in 1:dof
                @test abs(test_vector[k]) < 1e-7
            end
        end
    end
    bond_damage[1][:] .= 0
    bond_damage[2][:] .= 0
    bond_damage[3][:] .= 0
    bond_damage[4][:] .= 0

    shape_tensor, inverse_shape_tensor = PeriLab.IO.Geometry.shape_tensor(view(nodes, eachindex(nodes)), dof, nlist, volume, omega, bond_damage, undeformed_bond, shape_tensor, inverse_shape_tensor)
    @test isnothing(shape_tensor)
    @test isnothing(inverse_shape_tensor)
end

@testset "ut_strain" begin
    test_Data_manager = PeriLab.Data_manager
    dof = test_Data_manager.get_dof()
    nnodes = 4
    nodes = Vector{Int64}(1:nnodes)
    deformation_gradient = test_Data_manager.create_constant_node_field("Deformation Gradient", Float64, "Matrix", dof)
    strain = test_Data_manager.create_constant_node_field("Strain", Float64, "Matrix", dof)
    nlist = test_Data_manager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
    volume = test_Data_manager.create_constant_node_field("Volume", Float64, 1)
    omega = test_Data_manager.create_constant_bond_field("Influence Function", Float64, 1)
    bond_damage = test_Data_manager.create_constant_bond_field("Bond Damage", Float64, 1)
    undeformed_bond = test_Data_manager.create_constant_bond_field("Bond Geometry", Float64, dof)
    undeformed_bond_length = test_Data_manager.create_constant_bond_field("Bond Length", Float64, 1)
    shape_tensor = test_Data_manager.create_constant_node_field("Shape Tensor", Float64, "Matrix", dof)
    inverse_shape_tensor = test_Data_manager.create_constant_node_field("Inverse Shape Tensor", Float64, "Matrix", dof)


    deformation_gradient = PeriLab.IO.Geometry.deformation_gradient(view(nodes, eachindex(nodes)), dof, nlist, volume, omega, bond_damage, undeformed_bond, undeformed_bond, inverse_shape_tensor, deformation_gradient)
    strain = PeriLab.IO.Geometry.strain(view(nodes, eachindex(nodes)), deformation_gradient, strain)
    deformation_gradient = PeriLab.IO.Geometry.deformation_gradient(view(nodes, eachindex(nodes)), dof, nlist, volume, omega, bond_damage, undeformed_bond, undeformed_bond, inverse_shape_tensor, deformation_gradient)
    strain = PeriLab.IO.Geometry.strain(view(nodes, eachindex(nodes)), deformation_gradient, strain) - strain

    for i in 1:nnodes
        @test strain[i, 1, 1] == 0
        @test strain[i, 2, 1] == 0
        @test strain[i, 1, 2] == 0
        @test strain[i, 2, 2] == 0
    end
    deformation_gradient_3D = test_Data_manager.create_constant_node_field("Deformation Gradient 3D", Float64, "Matrix", 3)
    deformation_gradient_3D[1, 1, 1] = 2.0
    deformation_gradient_3D[1, 1, 2] = 1.0
    deformation_gradient_3D[1, 1, 3] = 2.0
    deformation_gradient_3D[1, 2, 1] = 2.0
    deformation_gradient_3D[1, 2, 2] = 1.0
    deformation_gradient_3D[1, 2, 3] = 2.3
    deformation_gradient_3D[1, 3, 1] = 2.0
    deformation_gradient_3D[1, 3, 2] = -1.0
    deformation_gradient_3D[1, 3, 3] = 3.0
    strain_3D = test_Data_manager.create_constant_node_field("Strain_3D", Float64, "Matrix", 3)
    strain_3D = PeriLab.IO.Geometry.strain(view(nodes, eachindex(nodes)), deformation_gradient_3D, strain_3D)
    identity = zeros(3, 3)
    identity[1, 1] = 1
    identity[2, 2] = 1
    identity[3, 3] = 1
    test = 0.5 * (transpose(deformation_gradient_3D[1, :, :]) * deformation_gradient_3D[1, :, :] - identity)
    for i in 1:dof
        for j in 1:dof
            @test strain_3D[1, i, j] == test[i, j]
        end
    end
end
@testset "ut_rotation_tensor" begin
    rot = PeriLab.IO.Geometry.rotation_tensor(fill(Float64(0), (1)))
    @test rot[1, 1] == 1
    @test rot[1, 2] == 0
    @test rot[2, 1] == 0
    @test rot[2, 2] == 1
    rot = PeriLab.IO.Geometry.rotation_tensor(fill(Float64(0), (3)))
    @test rot[1, 1] == 1
    @test rot[1, 2] == 0
    @test rot[1, 3] == 0
    @test rot[2, 1] == 0
    @test rot[2, 2] == 1
    @test rot[2, 3] == 0
    @test rot[3, 1] == 0
    @test rot[3, 2] == 0
    @test rot[3, 3] == 1
    rot = PeriLab.IO.Geometry.rotation_tensor(fill(Float64(90), (1)))
    @test rot[1, 1] < 1e-10
    @test rot[1, 2] == -1
    @test rot[2, 1] == 1
    @test rot[2, 2] < 1e-10
end

@testset "ut_compute_weighted_deformation_gradient" begin
    nnodes = [1, 2]
    dof = 3
    nlist = [[2], [1]]
    volume = [0.1, 0.2]
    gradient_weight = [0.5 0.5 0.5; 0.5 0.5 0.5]
    displacement = [0.0 0.0 0.0; 1.0 1.0 1.0]
    velocity = [0.0 0.0 0.0; 1.0 1.0 1.0]
    deformation_gradient = zeros(Float64, length(nnodes), dof, dof)
    deformation_gradient_dot = zeros(Float64, length(nnodes), dof, dof)

    deformation_gradient, deformation_gradient_dot = PeriLab.IO.Geometry.compute_weighted_deformation_gradient(nnodes, dof, nlist, volume, gradient_weight, displacement, velocity, deformation_gradient, deformation_gradient_dot)
    println()
    @test deformation_gradient[1, :, :] == [1.1 0.1 0.1; 0.1 1.1 0.1; 0.1 0.1 1.1]
    @test deformation_gradient[2, :, :] == [0.95 -0.05 -0.05; -0.05 0.95 -0.05; -0.05 -0.05 0.95]
    @test deformation_gradient_dot[1, :, :] == [0.1 0.1 0.1; 0.1 0.1 0.1; 0.1 0.1 0.1]
    @test deformation_gradient_dot[2, :, :] == [-0.05 -0.05 -0.05; -0.05 -0.05 -0.05; -0.05 -0.05 -0.05]
end