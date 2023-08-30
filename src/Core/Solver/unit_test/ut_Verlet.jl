using Test
using MPI
include("../Verlet.jl")
include("../Solver_control.jl")
include("../../../Support/geometry.jl")
include("../../../Support/data_manager.jl")
include("../../BC_manager.jl")
import .Data_manager
import .Geometry
import .Boundary_conditions
import .Solver
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
testValmech = 3.59255e-05
MPI.Init()
comm = MPI.COMM_WORLD
testDatamanager.set_comm(comm)
# from Peridigm
@testset "ut_crititical_time_step" begin

    t = compute_mechanical_crititical_time_step(testDatamanager, 140.0)
    @test testValmech / t - 1 < 1e-6
    testVal = 1e50 # to take from Peridigm
    t = compute_thermodynamic_crititical_time_step(testDatamanager, 140.0, 5.1)
    #@test testVal / t - 1 < 1e-6

end

testDatamanager.init_property()
testDatamanager.set_property(1, "Material Model", "Bulk Modulus", 140.0)
testDatamanager.set_property(2, "Material Model", "Bulk Modulus", 140.0)
@testset "ut_init_Verlet" begin
    params = Dict("Solver" => Dict("Initial Time" => 0.0, "Final Time" => 1.0, "Verlet" => Dict("Safety Factor" => 1.0)))
    start_time, dt, nsteps = init_Verlet(params, testDatamanager, true, false)

    @test start_time == params["Solver"]["Initial Time"]
    testStep = Int64(ceil((params["Solver"]["Final Time"] - params["Solver"]["Initial Time"]) / testValmech))
    @test nsteps == testStep
    testDt = (params["Solver"]["Final Time"] - params["Solver"]["Initial Time"]) / testStep

    @test testDt / dt - 1 < 1e-6
    params = Dict("Solver" => Dict("Initial Time" => 0.0, "Final Time" => 1.0, "Verlet" => Dict("Safety Factor" => 1.0, "Fixed dt" => 1e-5)))
    start_time, dt, nsteps = init_Verlet(params, testDatamanager, true, false)

    testStep = Int64(ceil((params["Solver"]["Final Time"] - params["Solver"]["Initial Time"]) / 1e-5))
    @test testStep == nsteps
    testFixdtVal = (params["Solver"]["Final Time"] - params["Solver"]["Initial Time"]) / testStep

    @test testFixdtVal / dt - 1 < 1e-6

end

nnodes = 5
dof = 2

testDatamanager = Data_manager
testDatamanager.set_comm(comm)
testDatamanager.set_nnodes(5)
testDatamanager.set_dof(2)

testDatamanager.set_glob_to_loc([1, 2, 3, 4, 5])
density = testDatamanager.create_constant_node_field("Density", Float32, 1)
force = testDatamanager.create_node_field("Forces", Float32, dof)
Y = testDatamanager.create_node_field("Deformed State", Float32, dof)
u = testDatamanager.create_node_field("Displacements", Float32, dof)
bu = testDatamanager.create_bond_field("Deformed Bond Geometry", Float32, dof + 1)
a = testDatamanager.create_constant_node_field("Acceleration", Float32, dof)
v = testDatamanager.create_node_field("Velocity", Float32, dof)

density[:] = [1e-6, 1e-6, 3e-6, 3e-6, 1e-6]
testDatamanager.set_nsets("Nset_1", [1, 2, 3])
testDatamanager.set_nsets("Nset_2", [3, 4, 7, 10])
blockNodes = [1, 1, 2, 2, 1]
params = Dict("Boundary Conditions" => Dict("BC_1" => Dict("Type" => "Force", "Node Set" => "Nset_1", "Coordinate" => "x", "Value" => "20*t"), "BC_2" => Dict("Type" => "Displacement", "Node Set" => "Nset_2", "Coordinate" => "y", "Value" => "5")))

bcs = Boundary_conditions.init_BCs(params, testDatamanager)
exos = []
outputs = Dict()
solver_options = Dict("Initial Time" => 0, "dt" => 3.59255e-05, "nsteps" => 2)
testDatamanager.set_rank(0)
exos = run_Verlet_solver(solver_options, Solver.get_blockNodes(blockNodes), bcs, testDatamanager, outputs, exos, Solver.write_results)
testDatamanager.set_rank(1)
# only if routine runs, if progress bar is not active
bcs = Boundary_conditions.init_BCs(params, testDatamanager)
exos = []
outputs = Dict()
solver_options = Dict("Initial Time" => 0, "dt" => 3.59255e-05, "nsteps" => 2)
exos = run_Verlet_solver(solver_options, Solver.get_blockNodes(blockNodes), bcs, testDatamanager, outputs, exos, Solver.write_results)

MPI.Finalize()