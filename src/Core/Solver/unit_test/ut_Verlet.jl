using Test
include("../Verlet.jl")
include("../../../Support/geometry.jl")
include("../../../Support/data_manager.jl")
import .Data_manager
import .Geometry
@testset "ut_test_timestep" begin
    @test test_timestep(1, 2) == 1
    @test test_timestep(2, 1.1) == 1.1
    @test test_timestep(2, 2) == 2
end

@testset "ut_get_cs_denominator" begin
    volume = [1, 2, 3]
    bondgeometry = [1, 2, 3]
    @test get_cs_denominator(1, volume, bondgeometry) == 1
    @test get_cs_denominator(2, volume, bondgeometry) == 2
    @test get_cs_denominator(3, volume, bondgeometry) == 3
    bondgeometry = [2, 4, 6]
    @test get_cs_denominator(1, volume, bondgeometry) == 0.5
    @test get_cs_denominator(2, volume, bondgeometry) == 1
    @test get_cs_denominator(3, volume, bondgeometry) == 1.5
    bondgeometry = [1, 0.5, 2]
    @test get_cs_denominator(1, volume, bondgeometry) == 1
    @test get_cs_denominator(2, volume, bondgeometry) == 5
    @test get_cs_denominator(3, volume, bondgeometry) == 6.5
end

nnodes = 5
dof = 2
testDatamanager = Data_manager

testDatamanager.set_nnodes(5)
testDatamanager.set_dof(2)
blocks = testDatamanager.create_constant_node_field("Block_Id", Int64, 1)
horizon = testDatamanager.create_constant_node_field("Horizon", Float32, 1)
coor = testDatamanager.create_constant_node_field("Coordinates", Float32, 2)
density = testDatamanager.create_constant_node_field("Density", Float32, 1)
volume = testDatamanager.create_constant_node_field("Volume", Float32, 1)
lenNlist = testDatamanager.create_constant_node_field("Number of Neighbors", Int64, 1)
lenNlist[:] = [4, 4, 4, 4, 4]

nlist = testDatamanager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
bondGeom = testDatamanager.create_constant_bond_field("Bond Geometry", Float32, 3)
nlist[1] = [2, 3, 4, 5]
nlist[2] = [1, 3, 4, 5]
nlist[3] = [1, 2, 4, 5]
nlist[4] = [1, 2, 3, 5]
nlist[5] = [1, 2, 3, 4]

coor[1, 1] = 0;
coor[1, 2] = 0;
coor[2, 1] = 0.5;
coor[2, 2] = 0.5;
coor[3, 1] = 1;
coor[3, 2] = 0;
coor[4, 1] = 0;
coor[4, 2] = 1;
coor[5, 1] = 1;
coor[5, 2] = 1;

volume[:] = [0.5, 0.5, 0.5, 0.5, 0.5]
density[:] = [1e-6, 1e-6, 3e-6, 3e-6, 1e-6]
horizon[:] = [1.1, 1.1, 1.1, 1.1, 1.1]

bondGeom = Geometry.bond_geometry(nnodes, dof, nlist, coor, bondGeom)

blocks[:] = [1, 1, 2, 2, 1]
blocks = testDatamanager.set_block_list(blocks)
# from Peridigm
@testset "ut_crititical_time_step" begin
    testVal = 3.59255e-05
    t = compute_mechanical_crititical_time_step(testDatamanager, 140.0)
    @test testVal / t - 1 < 1e-6
    testVal = 1e50 # to take from Peridigm
    t = compute_thermodynamic_crititical_time_step(testDatamanager, 140.0, 5.1)
    #@test testVal / t - 1 < 1e-6

end