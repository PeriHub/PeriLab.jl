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
    testDatamanager = Data_manager
    testDatamanager.set_nmasters(nnodes)
    testDatamanager.set_dof(dof)
    nn = testDatamanager.create_constant_node_field("Number of Neighbors", Int32, 1)
    nn[:] = [2, 2, 2, 1]
    coor = testDatamanager.create_constant_node_field("Coordinates", Float32, 2)
    bondGeom = testDatamanager.create_constant_bond_field("Bond Geometry", Float32, 3)
    nlist = testDatamanager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
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

    bondGeom = Geometry.bond_geometry(Vector(1:nnodes), dof, nlist, coor, bondGeom)

    @test bondGeom[1][1, 1] == 0.5
    @test bondGeom[1][1, 2] == 0.5
    @test bondGeom[1][1, 3] / sqrt(0.5) - 1 < 1e-8
    @test bondGeom[1][2, 1] == 1
    @test bondGeom[1][2, 2] == 0
    @test bondGeom[1][2, 3] == 1

    @test bondGeom[2][1, 1] == -0.5
    @test bondGeom[2][1, 2] == -0.5
    @test bondGeom[2][1, 3] / sqrt(0.5) - 1 < 1e-8
    @test bondGeom[2][2, 1] == 0.5
    @test bondGeom[2][2, 2] == -0.5
    @test bondGeom[2][2, 3] / sqrt(1.25) - 1 < 1e-8

    @test bondGeom[3][1, 1] == -1
    @test bondGeom[3][1, 2] == 0
    @test bondGeom[3][1, 3] == 1
    @test bondGeom[3][2, 1] == -0.5
    @test bondGeom[3][2, 2] == 0.5
    @test bondGeom[3][2, 3] / sqrt(1.25) - 1 < 1e-8

    @test bondGeom[4][1, 1] == 0.5
    @test bondGeom[4][1, 2] == -0.5
    @test bondGeom[4][1, 3] / sqrt(1.25) - 1 < 1e-8

end
@testset "ut_shapeTenorAndDefGrad" begin
    nnodes = 4
    dof = 2
    testDatamanager = Data_manager
    testDatamanager.set_nmasters(nnodes)
    testDatamanager.set_dof(dof)
    nn = testDatamanager.create_constant_node_field("Number of Neighbors", Int32, 1)
    nn[:] = [3, 3, 3, 3]
    delete!(testDatamanager.fields[Int64], "Neighborhoodlist")
    delete!(testDatamanager.field_names, "Neighborhoodlist")
    delete!(testDatamanager.fields[Float32], "Bond Geometry")
    delete!(testDatamanager.field_names, "Bond Geometry")
    coor = testDatamanager.create_constant_node_field("Coordinates", Float32, 2)
    bondGeom = testDatamanager.create_constant_bond_field("Bond Geometry", Float32, dof + 1)
    # does not have to be NP1 for testing
    deformedBond = testDatamanager.create_constant_bond_field("Deformed Bond Geometry", Float32, dof + 1)

    nlist = testDatamanager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
    volume = testDatamanager.create_constant_node_field("Volume", Float32, 1)
    omega = testDatamanager.create_constant_bond_field("Influence Function", Float32, 1)

    bondDamage = testDatamanager.create_constant_bond_field("Bond Damage", Float32, 1)
    shapeTensor = testDatamanager.create_constant_node_field("Shape Tensor", Float32, "Matrix", dof)
    invShapeTensor = testDatamanager.create_constant_node_field("Inverse Shape Tensor", Float32, "Matrix", dof)
    defGrad, defGradNP1 = testDatamanager.create_node_field("Deformation Gradient", Float32, "Matrix", dof)
    omega[1][:] .= 1
    omega[2][:] .= 1
    omega[3][:] .= 1
    omega[4][:] .= 1
    bondDamage[1][:] .= 1
    bondDamage[2][:] .= 1
    bondDamage[3][:] .= 1
    bondDamage[4][:] .= 1
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

    bondGeom = Geometry.bond_geometry(Vector(1:nnodes), dof, nlist, coor, bondGeom)
    shapeTensor, invShapeTensor = Geometry.shape_tensor(Vector(1:nnodes), dof, nlist, volume, omega, bondDamage, bondGeom, shapeTensor, invShapeTensor)

    defcoor = copy(coor)

    defGrad = Geometry.deformation_gradient(Vector(1:nnodes), dof, nlist, volume, omega, bondDamage, bondGeom, bondGeom, invShapeTensor, defGrad)
    for i in 1:nnodes
        @test defGrad[i, 1, 1] - 1 < 1e-7
        @test defGrad[i, 1, 2] < 1e-7
        @test defGrad[i, 2, 1] < 1e-7
        @test defGrad[i, 2, 2] - 1 < 1e-7
    end
    defcoor[1, 1] = -0.25
    defcoor[3, 1] = -0.25
    defcoor[2, 1] = 0.25
    defcoor[4, 1] = 0.25

    deformedBond = Geometry.bond_geometry(Vector(1:nnodes), dof, nlist, defcoor, deformedBond)
    defGrad = Geometry.deformation_gradient(Vector(1:nnodes), dof, nlist, volume, omega, bondDamage, deformedBond, bondGeom, invShapeTensor, defGrad)

    for i in 1:nnodes
        for j in 1:nn[i]
            testVec = defGrad[i, :, :] * bondGeom[i][j, 1:dof] - deformedBond[i][j, 1:dof]
            for k in 1:dof
                @test abs(testVec[k]) < 1e-7
            end
        end
    end

    defcoor = copy(coor)
    defcoor[3, 2] = 1.5
    defcoor[4, 2] = 1.5

    deformedBond = Geometry.bond_geometry(Vector(1:nnodes), dof, nlist, defcoor, deformedBond)
    defGrad = Geometry.deformation_gradient(Vector(1:nnodes), dof, nlist, volume, omega, bondDamage, deformedBond, bondGeom, invShapeTensor, defGrad)
    for i in 1:nnodes
        for j in nn[i]
            testVec = defGrad[i, :, :] * bondGeom[i][j, 1:dof] - deformedBond[i][j, 1:dof]
            for k in 1:dof
                @test abs(testVec[k]) < 1e-7
            end
        end
    end

    defcoor[1, 1] = 0.0
    defcoor[1, 2] = 0
    defcoor[2, 1] = 1
    defcoor[2, 2] = 0
    defcoor[3, 1] = 0.5
    defcoor[3, 2] = 0.5
    defcoor[4, 1] = 1.5
    defcoor[4, 2] = 0.5

    deformedBond = Geometry.bond_geometry(Vector(1:nnodes), dof, nlist, defcoor, deformedBond)
    defGrad = Geometry.deformation_gradient(Vector(1:nnodes), dof, nlist, volume, omega, bondDamage, deformedBond, bondGeom, invShapeTensor, defGrad)
    for i in 1:nnodes
        for j in nn[i]
            testVec = defGrad[i, :, :] * bondGeom[i][j, 1:dof] - deformedBond[i][j, 1:dof]
            for k in 1:dof
                @test abs(testVec[k]) < 1e-7
            end
        end
    end
end

@testset "ut_strain_increment" begin
    testDatamanager = Data_manager
    dof = testDatamanager.get_dof()
    nnodes = 4
    defGrad, defGradNP1 = testDatamanager.create_node_field("Deformation Gradient", Float32, "Matrix", dof)
    strainInc = testDatamanager.create_constant_node_field("Strain Increment", Float32, "Matrix", dof)
    nlist = testDatamanager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
    volume = testDatamanager.create_constant_node_field("Volume", Float32, 1)
    omega = testDatamanager.create_constant_bond_field("Influence Function", Float32, 1)
    bondDamage = testDatamanager.create_constant_bond_field("Bond Damage", Float32, 1)
    bondGeom = testDatamanager.create_constant_bond_field("Bond Geometry", Float32, dof + 1)
    shapeTensor = testDatamanager.create_constant_node_field("Shape Tensor", Float32, "Matrix", dof)
    invShapeTensor = testDatamanager.create_constant_node_field("Inverse Shape Tensor", Float32, "Matrix", dof)
    defGrad, defGradNP1 = testDatamanager.create_node_field("Deformation Gradient", Float32, "Matrix", dof)

    defGrad = Geometry.deformation_gradient(Vector(1:nnodes), dof, nlist, volume, omega, bondDamage, bondGeom, bondGeom, invShapeTensor, defGrad)
    defGradNP1 = Geometry.deformation_gradient(Vector(1:nnodes), dof, nlist, volume, omega, bondDamage, bondGeom, bondGeom, invShapeTensor, defGradNP1)
    strainInc = Geometry.strain_increment(Vector(1:nnodes), defGradNP1, defGrad, strainInc)

    for i in 1:nnodes
        @test strainInc[i, 1, 1] == 0
        @test strainInc[i, 2, 1] == 0
        @test strainInc[i, 1, 2] == 0
        @test strainInc[i, 2, 2] == 0
    end

    defGrad = zeros(4, 3, 3)
    defGradNP1 = zeros(4, 3, 3)
    strainInc = zeros(4, 3, 3)


    defGradNP1[1, 1, 1] = 2.0
    defGradNP1[1, 1, 2] = 1.0
    defGradNP1[1, 1, 3] = 2.0
    defGradNP1[1, 2, 1] = 2.0
    defGradNP1[1, 2, 2] = 1.0
    defGradNP1[1, 2, 3] = 2.3
    defGradNP1[1, 3, 1] = 2.0
    defGradNP1[1, 3, 2] = -1.0
    defGradNP1[1, 3, 3] = 3.0

    strainInc = Geometry.strain_increment(Vector(1:1), defGradNP1, defGrad, strainInc)
    for i in 1:dof
        for j in 1:dof
            @test strainInc[1, i, j] == defGradNP1[1, i, j]
        end
    end
end
