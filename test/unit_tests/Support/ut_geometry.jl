# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
include("../../../src/Support/data_manager.jl")
include("../../../src/Support/geometry.jl")
import .Geometry

@testset "ut_bond_geometry" begin
    nnodes = 4
    dof = 2
    test_Data_manager = Data_manager
    test_Data_manager.set_nmasters(nnodes)
    test_Data_manager.set_dof(dof)
    nn = test_Data_manager.create_constant_node_field("Number of Neighbors", Int32, 1)
    nn[:] = [2, 2, 2, 1]
    coor = test_Data_manager.create_constant_node_field("Coordinates", Float64, 2)
    bond_geometry = test_Data_manager.create_constant_bond_field("Bond Geometry", Float64, 3)
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

    bond_geometry = Geometry.bond_geometry(Vector(1:nnodes), dof, nlist, coor, bond_geometry)

    @test bond_geometry[1][1, 1] == 0.5
    @test bond_geometry[1][1, 2] == 0.5
    @test bond_geometry[1][1, 3] / sqrt(0.5) - 1 < 1e-8
    @test bond_geometry[1][2, 1] == 1
    @test bond_geometry[1][2, 2] == 0
    @test bond_geometry[1][2, 3] == 1

    @test bond_geometry[2][1, 1] == -0.5
    @test bond_geometry[2][1, 2] == -0.5
    @test bond_geometry[2][1, 3] / sqrt(0.5) - 1 < 1e-8
    @test bond_geometry[2][2, 1] == 0.5
    @test bond_geometry[2][2, 2] == -0.5
    @test bond_geometry[2][2, 3] / sqrt(1.25) - 1 < 1e-8

    @test bond_geometry[3][1, 1] == -1
    @test bond_geometry[3][1, 2] == 0
    @test bond_geometry[3][1, 3] == 1
    @test bond_geometry[3][2, 1] == -0.5
    @test bond_geometry[3][2, 2] == 0.5
    @test bond_geometry[3][2, 3] / sqrt(1.25) - 1 < 1e-8

    @test bond_geometry[4][1, 1] == 0.5
    @test bond_geometry[4][1, 2] == -0.5
    @test bond_geometry[4][1, 3] / sqrt(1.25) - 1 < 1e-8
    coor[:, :] .= 0
    bond_geometry[1][:, :] .= 0
    bond_geometry[2][:, :] .= 0
    bond_geometry[3][:, :] .= 0
    bond_geometry[4][:, :] .= 0

    bond_geometry = Geometry.bond_geometry(Vector(1:nnodes), dof, nlist, coor, bond_geometry)
    println()
    @test isnothing(bond_geometry)
end
@testset "ut_shapeTensorAnddeformation_gradient" begin
    nnodes = 4
    dof = 2
    nodes = Vector{Int64}(1:nnodes)
    test_Data_manager = Data_manager
    test_Data_manager.set_nmasters(nnodes)
    test_Data_manager.set_dof(dof)
    nn = test_Data_manager.create_constant_node_field("Number of Neighbors", Int32, 1)
    nn[:] = [3, 3, 3, 3]
    delete!(test_Data_manager.fields[Int64], "Neighborhoodlist")
    delete!(test_Data_manager.field_types, "Neighborhoodlist")
    delete!(test_Data_manager.fields[Float64], "Bond Geometry")
    delete!(test_Data_manager.field_types, "Bond Geometry")
    coor = test_Data_manager.create_constant_node_field("Coordinates", Float64, 2)
    bond_geometry = test_Data_manager.create_constant_bond_field("Bond Geometry", Float64, dof + 1)
    # does not have to be NP1 for testing
    deformed_bond = test_Data_manager.create_constant_bond_field("Deformed Bond Geometry", Float64, dof + 1)

    nlist = test_Data_manager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
    volume = test_Data_manager.create_constant_node_field("Volume", Float64, 1)
    omega = test_Data_manager.create_constant_bond_field("Influence Function", Float64, 1)

    bond_damage = test_Data_manager.create_constant_bond_field("Bond Damage", Float64, 1)
    shapeTensor = test_Data_manager.create_constant_node_field("Shape Tensor", Float64, "Matrix", dof)
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

    bond_geometry = Geometry.bond_geometry(Vector(1:nnodes), dof, nlist, coor, bond_geometry)
    shapeTensor, inverse_shape_tensor = Geometry.shape_tensor(view(nodes, eachindex(nodes)), dof, nlist, volume, omega, bond_damage, bond_geometry, shapeTensor, inverse_shape_tensor)

    deformed_coor = copy(coor)

    deformation_gradient = Geometry.deformation_gradient(view(nodes, eachindex(nodes)), dof, nlist, volume, omega, bond_damage, bond_geometry, bond_geometry, inverse_shape_tensor, deformation_gradient)
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

    deformed_bond = Geometry.bond_geometry(Vector(1:nnodes), dof, nlist, deformed_coor, deformed_bond)
    deformation_gradient = Geometry.deformation_gradient(view(nodes, eachindex(nodes)), dof, nlist, volume, omega, bond_damage, deformed_bond, bond_geometry, inverse_shape_tensor, deformation_gradient)

    for i in 1:nnodes
        for j in 1:nn[i]
            test_vector = deformation_gradient[i, :, :] * bond_geometry[i][j, 1:dof] - deformed_bond[i][j, 1:dof]
            for k in 1:dof
                @test abs(test_vector[k]) < 1e-7
            end
        end
    end

    deformed_coor = copy(coor)
    deformed_coor[3, 2] = 1.5
    deformed_coor[4, 2] = 1.5

    deformed_bond = Geometry.bond_geometry(Vector(1:nnodes), dof, nlist, deformed_coor, deformed_bond)
    deformation_gradient = Geometry.deformation_gradient(view(nodes, eachindex(nodes)), dof, nlist, volume, omega, bond_damage, deformed_bond, bond_geometry, inverse_shape_tensor, deformation_gradient)
    for i in 1:nnodes
        for j in nn[i]
            test_vector = deformation_gradient[i, :, :] * bond_geometry[i][j, 1:dof] - deformed_bond[i][j, 1:dof]
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

    deformed_bond = Geometry.bond_geometry(Vector(1:nnodes), dof, nlist, deformed_coor, deformed_bond)
    deformation_gradient = Geometry.deformation_gradient(view(nodes, eachindex(nodes)), dof, nlist, volume, omega, bond_damage, deformed_bond, bond_geometry, inverse_shape_tensor, deformation_gradient)
    for i in 1:nnodes
        for j in nn[i]
            test_vector = deformation_gradient[i, :, :] * bond_geometry[i][j, 1:dof] - deformed_bond[i][j, 1:dof]
            for k in 1:dof
                @test abs(test_vector[k]) < 1e-7
            end
        end
    end
    bond_damage[1][:] .= 0
    bond_damage[2][:] .= 0
    bond_damage[3][:] .= 0
    bond_damage[4][:] .= 0

    shapeTensor, inverse_shape_tensor = Geometry.shape_tensor(view(nodes, eachindex(nodes)), dof, nlist, volume, omega, bond_damage, bond_geometry, shapeTensor, inverse_shape_tensor)
    @test isnothing(shapeTensor)
    @test isnothing(inverse_shape_tensor)
end

@testset "ut_strain" begin
    test_Data_manager = Data_manager
    dof = test_Data_manager.get_dof()
    nnodes = 4
    nodes = Vector{Int64}(1:nnodes)
    deformation_gradient = test_Data_manager.create_constant_node_field("Deformation Gradient", Float64, "Matrix", dof)
    strainInc = test_Data_manager.create_constant_node_field("Strain", Float64, "Matrix", dof)
    nlist = test_Data_manager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
    volume = test_Data_manager.create_constant_node_field("Volume", Float64, 1)
    omega = test_Data_manager.create_constant_bond_field("Influence Function", Float64, 1)
    bond_damage = test_Data_manager.create_constant_bond_field("Bond Damage", Float64, 1)
    bond_geometry = test_Data_manager.create_constant_bond_field("Bond Geometry", Float64, dof + 1)
    shapeTensor = test_Data_manager.create_constant_node_field("Shape Tensor", Float64, "Matrix", dof)
    inverse_shape_tensor = test_Data_manager.create_constant_node_field("Inverse Shape Tensor", Float64, "Matrix", dof)


    deformation_gradient = Geometry.deformation_gradient(view(nodes, eachindex(nodes)), dof, nlist, volume, omega, bond_damage, bond_geometry, bond_geometry, inverse_shape_tensor, deformation_gradient)
    strainInc = Geometry.strain(view(nodes, eachindex(nodes)), deformation_gradient, strainInc)
    deformation_gradient = Geometry.deformation_gradient(view(nodes, eachindex(nodes)), dof, nlist, volume, omega, bond_damage, bond_geometry, bond_geometry, inverse_shape_tensor, deformation_gradient)
    strainInc = Geometry.strain(view(nodes, eachindex(nodes)), deformation_gradient, strainInc) - strainInc
    strainInc = Geometry.strain(view(nodes, eachindex(nodes)), deformation_gradient, strainInc) - strainInc


    for i in 1:nnodes
        @test strainInc[i, 1, 1] == 0
        @test strainInc[i, 2, 1] == 0
        @test strainInc[i, 1, 2] == 0
        @test strainInc[i, 2, 2] == 0
    end

    deformation_gradient = zeros(4, 3, 3)
    strainInc = zeros(4, 3, 3)


    deformation_gradient[1, 1, 1] = 2.0
    deformation_gradient[1, 1, 2] = 1.0
    deformation_gradient[1, 1, 3] = 2.0
    deformation_gradient[1, 2, 1] = 2.0
    deformation_gradient[1, 2, 2] = 1.0
    deformation_gradient[1, 2, 3] = 2.3
    deformation_gradient[1, 3, 1] = 2.0
    deformation_gradient[1, 3, 2] = -1.0
    deformation_gradient[1, 3, 3] = 3.0

    strainInc = Geometry.strain(view(nodes, eachindex(nodes)), deformation_gradient, strainInc)
    identity = zeros(3, 3)
    identity[1, 1] = 1
    identity[2, 2] = 1
    identity[3, 3] = 1
    test = 0.5 * (transpose(deformation_gradient[1, :, :]) * deformation_gradient[1, :, :] - identity)
    for i in 1:dof
        for j in 1:dof
            @test strainInc[1, i, j] == test[i, j]
        end
    end
end
@testset "ut_rotation_tensor" begin
    rot = Geometry.rotation_tensor(fill(Float64(0), (1)))
    @test rot[1, 1] == 1
    @test rot[1, 2] == 0
    @test rot[1, 3] == 0
    @test rot[2, 1] == 0
    @test rot[2, 2] == 1
    @test rot[2, 3] == 0
    @test rot[3, 1] == 0
    @test rot[3, 2] == 0
    @test rot[3, 3] == 1
    rot = Geometry.rotation_tensor(fill(Float64(0), (3)))
    @test rot[1, 1] == 1
    @test rot[1, 2] == 0
    @test rot[1, 3] == 0
    @test rot[2, 1] == 0
    @test rot[2, 2] == 1
    @test rot[2, 3] == 0
    @test rot[3, 1] == 0
    @test rot[3, 2] == 0
    @test rot[3, 3] == 1
    rot = Geometry.rotation_tensor(fill(Float64(90), (1)))
    @test rot[1, 1] < 1e-10
    @test rot[1, 2] == -1
    @test rot[1, 3] < 1e-10
    @test rot[2, 1] == 1
    @test rot[2, 2] < 1e-10
    @test rot[2, 3] < 1e-10
    @test rot[3, 1] < 1e-10
    @test rot[3, 2] < 1e-10
    @test rot[3, 3] == 1
end

