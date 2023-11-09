# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test


include("../../../../src/Core/BC_manager.jl")
include("../../../../src/Core/Solver/Verlet.jl")
include("../../../../src/Core/Solver/Solver_control.jl")
include("../../../../src/Support/geometry.jl")

if !isdefined(@__MODULE__, :Data_manager)
    include("../../../../src/Support/data_manager.jl")
end
include("../../../../src/Support/Parameters/parameter_handling.jl")
using .Data_manager
using .Geometry
using .Boundary_conditions
using .Solver
using .Verlet
using MPI
@testset "ut_test_timestep" begin
    @test Verlet.test_timestep(1.0, 2.0) == 1
    @test Verlet.test_timestep(2.0, 1.1) == 1.1
    @test Verlet.test_timestep(2.0, 2.0) == 2
end

@testset "ut_get_cs_denominator" begin
    volume = Float64[1, 2, 3]
    bond_geometryetry = [1, 2, 3]
    @test Verlet.get_cs_denominator(volume, bond_geometryetry) == 3
    bond_geometryetry = [2, 4, 6]
    @test Verlet.get_cs_denominator(volume, bond_geometryetry) == 1.5
    bond_geometryetry = [1, 0.5, 2]
    @test Verlet.get_cs_denominator(volume, bond_geometryetry) == 6.5
end

# nnodes = 5
# dof = 2
# test_Data_manager = Data_manager
# comm = MPI.COMM_WORLD
# test_Data_manager.set_comm(comm)
# test_Data_manager.set_nmasters(5)
# test_Data_manager.set_dof(2)
# blocks = test_Data_manager.create_constant_node_field("Block_Id", Int64, 1)
# horizon = test_Data_manager.create_constant_node_field("Horizon", Float64, 1)
# coor = test_Data_manager.create_constant_node_field("Coordinates", Float64, 2)
# density = test_Data_manager.create_constant_node_field("Density", Float64, 1)
# volume = test_Data_manager.create_constant_node_field("Volume", Float64, 1)
# lenNlist = test_Data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
# lenNlist .= 4

# nlist = test_Data_manager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
# bond_geometry = test_Data_manager.create_constant_bond_field("Bond Geometry", Float64, 3)
# nlist[1] = [2, 3, 4, 5]
# nlist[2] = [1, 3, 4, 5]
# nlist[3] = [1, 2, 4, 5]
# nlist[4] = [1, 2, 3, 5]
# nlist[5] = [1, 2, 3, 4]

# coor[1, 1] = 0;
# coor[1, 2] = 0;
# coor[2, 1] = 0.5;
# coor[2, 2] = 0.5;
# coor[3, 1] = 1;
# coor[3, 2] = 0;
# coor[4, 1] = 0;
# coor[4, 2] = 1;
# coor[5, 1] = 1;
# coor[5, 2] = 1;

# volume[:] = [0.5, 0.5, 0.5, 0.5, 0.5]
# density[:] = [1e-6, 1e-6, 3e-6, 3e-6, 1e-6]
# horizon[:] = [3.1, 3.1, 3.1, 3.1, 3.1]

# bond_geometry = Geometry.bond_geometry(Vector(1:nnodes), dof, nlist, coor, bond_geometry)

# blocks[:] = [1, 1, 2, 2, 1]
# blocks = test_Data_manager.set_block_list(blocks)
# # from Peridigm
# testValmech = 0.0002853254715348906
# testVal = 72.82376628733019

# # from Peridigm
# @testset "ut_mechanical_critical_time_step" begin

#     t = Verlet.compute_mechanical_critical_time_step(Vector{Int64}(1:nnodes), test_Data_manager, Float64(140.0))
#     @test abs(testValmech / t - 1) < 1e-6

# end
# # from Peridigm
# @testset "ut_thermodynamic_crititical_time_step" begin

#     t = Verlet.compute_thermodynamic_critical_time_step(Vector{Int64}(1:nnodes), test_Data_manager, Float64(0.12), Float64(1.8e9))
#     @test abs(testVal / t - 1) < 1e-6

# end

# test_Data_manager.init_property()
# test_Data_manager.set_property(1, "Material Model", "Bulk Modulus", Float64(140.0))
# test_Data_manager.set_property(2, "Material Model", "Bulk Modulus", Float64(140.0))
# @testset "ut_init_Verlet" begin
#     params = Dict("Solver" => Dict("Initial Time" => 0.0, "Final Time" => 1.0, "Verlet" => Dict("Safety Factor" => 1.0)))
#     start_time, dt, nsteps = Verlet.init_solver(params, test_Data_manager, Dict{Int64,Vector{Int64}}(1 => Vector{Int64}(1:nnodes)), true, false)

#     @test start_time == params["Solver"]["Initial Time"]
#     testStep = Int64(ceil((params["Solver"]["Final Time"] - params["Solver"]["Initial Time"]) / testValmech))
#     @test nsteps == testStep
#     testDt = (params["Solver"]["Final Time"] - params["Solver"]["Initial Time"]) / testStep

#     #nsteps = ceil((end_time - initial_time) / dt)
#     #dt = (end_time - initial_time) / nsteps

#     @test testDt / dt - 1 < 1e-6
#     params = Dict("Solver" => Dict("Initial Time" => 0.0, "Final Time" => 1.0, "Verlet" => Dict("Safety Factor" => 1.0, "Fixed dt" => 1e-5)))
#     start_time, dt, nsteps = Verlet.init_solver(params, test_Data_manager, Dict{Int64,Vector{Int64}}(1 => Vector{Int64}(1:nnodes)), true, false)

#     testStep = Int64(ceil((params["Solver"]["Final Time"] - params["Solver"]["Initial Time"]) / 1e-5))
#     @test testStep == nsteps
#     testFixdtVal = (params["Solver"]["Final Time"] - params["Solver"]["Initial Time"]) / testStep

#     @test testFixdtVal / dt - 1 < 1e-6

# end

# nnodes = 5
# dof = 2

# test_Data_manager = Data_manager
# test_Data_manager.set_comm(comm)
# test_Data_manager.set_nmasters(5)
# test_Data_manager.set_dof(2)

# test_Data_manager.set_glob_to_loc([1, 2, 3, 4, 5])
# density = test_Data_manager.create_constant_node_field("Density", Float64, 1)
# force = test_Data_manager.create_node_field("Forces", Float64, dof)
# Y = testDatama#nager.create_node_field("Deformed State", Float64, dof)
# u = test_Data_manager.create_node_field("Displacements", Float64, dof)
# bu = test_Data_manager.create_bond_field("Deformed Bond Geometry", Float64, dof + 1)
# a = test_Data_manager.create_constant_node_field("Acceleration", Float64, dof)
# v = test_Data_manager.create_node_field("Velocity", Float64, dof)

# density[:] = [1e-6, 1e-6, 3e-6, 3e-6, 1e-6]
# test_Data_manager.set_nset("Nset_1", [1, 2, 3])
# test_Data_manager.set_nset("Nset_2", [3, 4, 7, 10])
# blockNodes = [1, 1, 2, 2, 1]
# params = Dict("Boundary Conditions" => Dict("BC_1" => Dict("Type" => "Force", "Node Set" => "Nset_1", "Coordinate" => "x", "Value" => "20*t"), "BC_2" => Dict("Type" => "Displacement", "Node Set" => "Nset_2", "Coordinate" => "y", "Value" => "5")))

# bcs = Boundary_conditions.init_BCs(params, test_Data_manager)
# result_files = []
# outputs = Dict()
# solver_options = Dict("Initial Time" => 0, "dt" => 3.59255e-05, "nsteps" => 2)
# test_Data_manager.set_rank(0)
# result_files = run_Verlet_solver(solver_options, Solver.get_nodes(blockNodes), bcs, test_Data_manager, outputs, result_files, Solver.write_results)
# test_Data_manager.set_rank(1)
# # only if routine runs, if progress bar is not active
# bcs = Boundary_conditions.init_BCs(params, test_Data_manager)
# result_files = []
# outputs = Dict()
# solver_options = Dict("Initial Time" => 0, "dt" => 3.59255e-05, "nsteps" => 2)
# result_files = run_Verlet_solver(solver_options, Solver.get_nodes(blockNodes), bcs, test_Data_manager, outputs, result_files, Solver.write_results)

