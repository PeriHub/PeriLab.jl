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
    coor = testDatamanager.create_constant_node_field("Coordinates", Float64, 2)
    bondGeom = testDatamanager.create_constant_bond_field("Bond Geometry", Float64, 3)
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
    nodes = Vector{Int64}(1:nnodes)
    testDatamanager = Data_manager
    testDatamanager.set_nmasters(nnodes)
    testDatamanager.set_dof(dof)
    nn = testDatamanager.create_constant_node_field("Number of Neighbors", Int32, 1)
    nn[:] = [3, 3, 3, 3]
    delete!(testDatamanager.fields[Int64], "Neighborhoodlist")
    delete!(testDatamanager.field_names, "Neighborhoodlist")
    delete!(testDatamanager.fields[Float64], "Bond Geometry")
    delete!(testDatamanager.field_names, "Bond Geometry")
    coor = testDatamanager.create_constant_node_field("Coordinates", Float64, 2)
    bondGeom = testDatamanager.create_constant_bond_field("Bond Geometry", Float64, dof + 1)
    # does not have to be NP1 for testing
    deformedBond = testDatamanager.create_constant_bond_field("Deformed Bond Geometry", Float64, dof + 1)

    nlist = testDatamanager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
    volume = testDatamanager.create_constant_node_field("Volume", Float64, 1)
    omega = testDatamanager.create_constant_bond_field("Influence Function", Float64, 1)

    bondDamage = testDatamanager.create_constant_bond_field("Bond Damage", Float64, 1)
    shapeTensor = testDatamanager.create_constant_node_field("Shape Tensor", Float64, "Matrix", dof)
    invShapeTensor = testDatamanager.create_constant_node_field("Inverse Shape Tensor", Float64, "Matrix", dof)
    defGrad = testDatamanager.create_constant_node_field("Deformation Gradient", Float64, "Matrix", dof)
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
    shapeTensor, invShapeTensor = Geometry.shape_tensor(view(nodes, eachindex(nodes)), dof, nlist, volume, omega, bondDamage, bondGeom, shapeTensor, invShapeTensor)

    defcoor = copy(coor)

    defGrad = Geometry.deformation_gradient(view(nodes, eachindex(nodes)), dof, nlist, volume, omega, bondDamage, bondGeom, bondGeom, invShapeTensor, defGrad)
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
    defGrad = Geometry.deformation_gradient(view(nodes, eachindex(nodes)), dof, nlist, volume, omega, bondDamage, deformedBond, bondGeom, invShapeTensor, defGrad)

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
    defGrad = Geometry.deformation_gradient(view(nodes, eachindex(nodes)), dof, nlist, volume, omega, bondDamage, deformedBond, bondGeom, invShapeTensor, defGrad)
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
    defGrad = Geometry.deformation_gradient(view(nodes, eachindex(nodes)), dof, nlist, volume, omega, bondDamage, deformedBond, bondGeom, invShapeTensor, defGrad)
    for i in 1:nnodes
        for j in nn[i]
            testVec = defGrad[i, :, :] * bondGeom[i][j, 1:dof] - deformedBond[i][j, 1:dof]
            for k in 1:dof
                @test abs(testVec[k]) < 1e-7
            end
        end
    end
end

@testset "ut_strain" begin
    testDatamanager = Data_manager
    dof = testDatamanager.get_dof()
    nnodes = 4
    nodes = Vector{Int64}(1:nnodes)
    defGrad = testDatamanager.create_constant_node_field("Deformation Gradient", Float64, "Matrix", dof)
    strainInc = testDatamanager.create_constant_node_field("Strain", Float64, "Matrix", dof)
    nlist = testDatamanager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
    volume = testDatamanager.create_constant_node_field("Volume", Float64, 1)
    omega = testDatamanager.create_constant_bond_field("Influence Function", Float64, 1)
    bondDamage = testDatamanager.create_constant_bond_field("Bond Damage", Float64, 1)
    bondGeom = testDatamanager.create_constant_bond_field("Bond Geometry", Float64, dof + 1)
    shapeTensor = testDatamanager.create_constant_node_field("Shape Tensor", Float64, "Matrix", dof)
    invShapeTensor = testDatamanager.create_constant_node_field("Inverse Shape Tensor", Float64, "Matrix", dof)


    defGrad = Geometry.deformation_gradient(view(nodes, eachindex(nodes)), dof, nlist, volume, omega, bondDamage, bondGeom, bondGeom, invShapeTensor, defGrad)
    strainInc = Geometry.strain(view(nodes, eachindex(nodes)), defGrad, strainInc)
    defGrad = Geometry.deformation_gradient(view(nodes, eachindex(nodes)), dof, nlist, volume, omega, bondDamage, bondGeom, bondGeom, invShapeTensor, defGrad)
    strainInc = Geometry.strain(view(nodes, eachindex(nodes)), defGrad, strainInc) - strainInc
    strainInc = Geometry.strain(view(nodes, eachindex(nodes)), defGrad, strainInc) - strainInc


    for i in 1:nnodes
        @test strainInc[i, 1, 1] == 0
        @test strainInc[i, 2, 1] == 0
        @test strainInc[i, 1, 2] == 0
        @test strainInc[i, 2, 2] == 0
    end

    defGrad = zeros(4, 3, 3)
    strainInc = zeros(4, 3, 3)


    defGrad[1, 1, 1] = 2.0
    defGrad[1, 1, 2] = 1.0
    defGrad[1, 1, 3] = 2.0
    defGrad[1, 2, 1] = 2.0
    defGrad[1, 2, 2] = 1.0
    defGrad[1, 2, 3] = 2.3
    defGrad[1, 3, 1] = 2.0
    defGrad[1, 3, 2] = -1.0
    defGrad[1, 3, 3] = 3.0

    strainInc = Geometry.strain(view(nodes, eachindex(nodes)), defGrad, strainInc)
    identity = zeros(3, 3)
    identity[1, 1] = 1
    identity[2, 2] = 1
    identity[3, 3] = 1
    test = 0.5 * (transpose(defGrad[1, :, :]) * defGrad[1, :, :] - identity)
    for i in 1:dof
        for j in 1:dof
            @test strainInc[1, i, j] == test[i, j]
        end
    end
end
@testset "ut_rotation_tensor" begin
    rot = Geometry.rotation_tensor(Float64(0))
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
    rot = Geometry.rotation_tensor(Float64(90))
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

