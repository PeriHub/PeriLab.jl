
using Test
include("../../../src/Support/data_manager.jl")
include("../../../src/Support/geometry.jl")
import .Geometry
import .Data_manager

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
    nlist = testDatamanager.create_constant_bond_field("Neighborhoodlist", Float32, 1)
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

    bondGeom = Geometry.bond_geometry(1:nnodes, dof, nlist, coor, bondGeom)

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
