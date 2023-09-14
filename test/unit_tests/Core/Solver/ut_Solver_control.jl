using Test
using MPI
import TimerOutputs
include("../../../../src/Core/Solver/Solver_control.jl")
include("../../../../src/Support/geometry.jl")
if !isdefined(@__MODULE__, :Data_manager)
    include("../../../../src/Support/data_manager.jl")
end
include("../../../../src/Support/Parameters/parameter_handling.jl")
import .Solver

@testset "ut_get_blockNodes" begin
    blockIDs = [1, 1, 1, 2, 2, 3, 3, 3, 3, 1, 1, 2, 3, 3, 1, 1, 2]
    blockNodes = Solver.get_blockNodes(blockIDs, length(blockIDs))
    @test blockNodes[1] == [1, 2, 3, 10, 11, 15, 16]
    @test blockNodes[2] == [4, 5, 12, 17]
    @test blockNodes[3] == [6, 7, 8, 9, 13, 14]

end


@testset "ut_init" begin

    comm = MPI.COMM_WORLD
    params = Dict("Discretization" => Dict("Node Sets" => Dict("Nset_1" => [1 2 3])), "Boundary Conditions" => Dict("BC_1" => Dict("Type" => "Force", "Node Set" => "Nset_1", "Coordinate" => "x", "Value" => "20*t"), "BC_2" => Dict("Type" => "Displacement", "Node Set" => "Nset_2", "Coordinate" => "y", "Value" => "0")), "Blocks" => Dict("block_1" => Dict("Material Model" => "a", "Density" => 1e-6, "Horizon" => 2), "block_2" => Dict("Material Model" => "c", "Density" => 3e-6, "Horizon" => 2)), "Physics" => Dict("Material Models" => Dict("a" => Dict("Bulk Modulus" => 140.0), "c" => Dict("Bulk Modulus" => 140.0, "value2" => 1)), "Damage Models" => Dict("a" => Dict("value" => 3), "c" => Dict("value" => [1 2], "value2" => 1)), "Thermal Models" => Dict("therm" => Dict("Specific Heat Capacity" => 1e-9, "Lambda" => 1.1))), "Solver" => Dict("Material Models" => true, "Damage Models" => true, "Additive Models" => false, "Thermal Models" => true, "Initial Time" => 0.0, "Final Time" => 1.0, "Verlet" => Dict("Safety Factor" => 1.0)), "Outputs" => Dict("Output1" => Dict("Output Filename" => "test1.e", "Output Variables" => Dict("Velocity" => true, "Force" => false, "Displacements" => true)), "Output2" => Dict("Output Filename" => "test2.e", "Output Variables" => Dict("Displacements" => true))), "Boundary Conditions" => Dict("BC1" => Dict("Type" => "Force", "Node Set" => "Nset_1", "Coordinate" => "x", "Value" => "20*t"), "BC2" => Dict("Type" => "Force", "Node Set" => "Nset_2", "Coordinate" => "x", "Value" => "20*t")))
    nnodes = 5
    dof = 2
    testDatamanager = Data_manager
    testDatamanager.set_comm(comm)
    testDatamanager.set_dof(2)
    testDatamanager.set_nmasters(nnodes)
    testDatamanager.set_nnsets(2)
    nn = testDatamanager.create_constant_node_field("Number of Neighbors", Int64, 1)

    nn[1] = 4
    nn[2] = 4
    nn[3] = 4
    nn[4] = 4
    nn[5] = 4
    id = testDatamanager.create_constant_node_field("Block_Id", Int64, 1)

    id[:] = [1, 1, 2, 2, 1]
    testDatamanager.set_block_list(id)

    horizon = testDatamanager.create_constant_node_field("Horizon", Float32, 1)
    coor = testDatamanager.create_constant_node_field("Coordinates", Float32, 2)
    density = testDatamanager.create_constant_node_field("Density", Float32, 1)
    volume = testDatamanager.create_constant_node_field("Volume", Float32, 1)

    testDatamanager.set_nset("Nset_1", [2])
    testDatamanager.set_nset("Nset_2", [2, 3, 4])

    testDatamanager.set_glob_to_loc([1, 2, 3, 4, 5])

    nlist = testDatamanager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
    bondGeom = testDatamanager.create_constant_bond_field("Bond Geometry", Float32, 3)
    nlist[1] = [2, 3, 4, 5]
    nlist[2] = [1, 3, 4, 5]
    nlist[3] = [1, 2, 4, 5]
    nlist[4] = [1, 2, 3, 5]
    nlist[5] = [1, 2, 3, 4]

    coor[1, 1] = 0
    coor[1, 2] = 0
    coor[2, 1] = 0.5
    coor[2, 2] = 0.5
    coor[3, 1] = 1
    coor[3, 2] = 0
    coor[4, 1] = 0
    coor[4, 2] = 1
    coor[5, 1] = 1
    coor[5, 2] = 1

    volume[:] = [0.5, 0.5, 0.5, 0.5, 0.5]
    horizon[:] = [1.1, 1.1, 1.1, 1.1, 1.1]
    node_list = zeros(Int64, nnodes)
    node_list[:] = 1:nnodes
    bondGeom = Geometry.bond_geometry(node_list, dof, nlist, coor, bondGeom)

    blockNodes, bcs, datamanager, solver_options = Solver.init(params, testDatamanager)

    @test solver_options["Material Models"]
    @test solver_options["Damage Models"]
    @test solver_options["Thermal Models"]
    @test solver_options["Additive Models"] == false
    @test solver_options["Initial Time"] == 0.0
    @test solver_options["dt"] / 3.59247018e-5 - 1 < 1e-8
    @test solver_options["nsteps"] == 27836

    @test bcs["BC1"]["Node Set"] == [2]
    @test bcs["BC2"]["Node Set"] == [2, 3, 4]
    @test blockNodes[1] == [1, 2, 5]
    @test blockNodes[2] == [3, 4]
    @test "TemperatureNP1" in testDatamanager.get_all_field_keys()
    @test "DisplacementsNP1" in testDatamanager.get_all_field_keys()
    @test "DamageNP1" in testDatamanager.get_all_field_keys()
    @test ("ActivatedNP1" in testDatamanager.get_all_field_keys()) == false
    @test ("Activated" in testDatamanager.get_all_field_keys()) == false

end
