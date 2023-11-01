# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../../src/IO/IO.jl")
include("../../../src/Support/data_manager.jl")
include("../../../src/Support/Parameters/parameter_handling.jl")
using Test
import .IO
using Exodus
using MPI

testDatamanager = Data_manager
filename1 = "test1"
filename2 = "test2"
dof = 2
nnodes = 5
comm = MPI.COMM_WORLD
testDatamanager.set_nmasters(nnodes)
testDatamanager.set_comm(comm)
testDatamanager.set_dof(dof)
testDatamanager.set_max_rank(1)
testDatamanager.create_constant_node_field("Coordinates", Float64, 2)
coordinates = testDatamanager.get_field("Coordinates")
testDatamanager.create_constant_node_field("Block_Id", Int64, 1)
block_Id = testDatamanager.get_field("Block_Id")
testDatamanager.create_node_field("Displacements", Float64, 2)
testDatamanager.create_node_field("Forces", Float64, 6)
params = Dict("Outputs" => Dict("Output1" => Dict("Output Filename" => filename1, "Output Variables" => Dict("Forces" => true)), "Output2" => Dict("Output Filename" => filename2, "Output Variables" => Dict("Displacements" => true, "Forces" => true))))
coordinates[1, 1] = 0
coordinates[1, 2] = 0
coordinates[2, 1] = 1
coordinates[2, 2] = 0
coordinates[3, 1] = 0
coordinates[3, 2] = 1
coordinates[4, 1] = 1
coordinates[4, 2] = 1
coordinates[5, 1] = 2
coordinates[5, 2] = 2
block_Id .+= 1
block_Id[end] = 2

@testset "ut_get_results_mapping" begin
    output = IO.get_results_mapping(params, testDatamanager)
    @test sort(collect(keys(output[1]["Fields"]))) == ["Forcesxx", "Forcesxy", "Forcesxz", "Forcesyx", "Forcesyy", "Forcesyz"]
    @test sort(collect(keys(output[2]["Fields"]))) == ["Displacementsx", "Displacementsy", "Forcesxx", "Forcesxy", "Forcesxz", "Forcesyx", "Forcesyy", "Forcesyz"]
    for i in 1:2
        dofForce = 0
        dofDisp = 0
        for entry in keys(sort(output[i]["Fields"]))
            if occursin("Forces", entry)
                dofForce += 1
                @test output[i]["Fields"][entry]["fieldname"] == "ForcesNP1"
                @test output[i]["Fields"][entry]["result_id"] == i
                @test output[i]["Fields"][entry]["dof"] == dofForce
                @test output[i]["Fields"][entry]["type"] == Float64
            else
                dofDisp += 1
                @test output[i]["Fields"][entry]["fieldname"] == "DisplacementsNP1"
                @test output[i]["Fields"][entry]["result_id"] == 1
                @test output[i]["Fields"][entry]["dof"] == dofDisp
                @test output[i]["Fields"][entry]["type"] == Float64
            end
        end
    end
end

@testset "ut_init_write_result_and_write_results" begin
    result_files, outputs = IO.init_write_results(params, "", testDatamanager, 2)
    @test length(result_files) == 2
    @test length(result_files[1].nodal_var_name_dict) == 6
    entries = collect(keys(result_files[1].nodal_var_name_dict))
    @test sort(entries) == ["Forcesxx", "Forcesxy", "Forcesxz", "Forcesyx", "Forcesyy", "Forcesyz"]
    @test length(result_files[2].nodal_var_name_dict) == 8
    entries = collect(keys(result_files[2].nodal_var_name_dict))
    @test sort(entries) == ["Displacementsx", "Displacementsy", "Forcesxx", "Forcesxy", "Forcesxz", "Forcesyx", "Forcesyy", "Forcesyz"]

    coords = vcat(transpose(coordinates))
    for exo in result_files
        @test exo.init.num_dim == dof
        exo_coords = read_coordinates(exo)
        exo_nsets = read_sets(exo, NodeSet)
        @test coords == exo_coords
        @test exo_nsets == []
    end
    for i in 1:2
        dofForce = 0
        dofDisp = 0
        for entry in keys(sort(outputs[i]["Fields"]))
            if occursin("Forces", entry)
                dofForce += 1
                @test outputs[i]["Fields"][entry]["fieldname"] == "ForcesNP1"
                @test outputs[i]["Fields"][entry]["result_id"] == i
                @test outputs[i]["Fields"][entry]["dof"] == dofForce
                @test outputs[i]["Fields"][entry]["type"] == Float64
            else
                dofDisp += 1
                @test outputs[i]["Fields"][entry]["fieldname"] == "DisplacementsNP1"
                @test outputs[i]["Fields"][entry]["result_id"] == 1
                @test outputs[i]["Fields"][entry]["dof"] == dofDisp
                @test outputs[i]["Fields"][entry]["type"] == Float64
            end
        end
    end
    IO.output_frequency = [Dict{String,Int64}("Counter" => 0, "Output Frequency" => 1, "Step" => 1), Dict{String,Int64}("Counter" => 0, "Output Frequency" => 1, "Step" => 1)]
    IO.write_results(result_files, 1.5, outputs, testDatamanager)

    @test read_time(result_files[1], 2) == 1.5
    @test read_time(result_files[2], 2) == 1.5
    IO.write_results(result_files, 1.6, outputs, testDatamanager)

    @test read_time(result_files[1], 3) == 1.6
    @test read_time(result_files[2], 3) == 1.6
    IO.write_results([], 1.6, outputs, testDatamanager)
    testBool = false
    try
        read_time(result_files[1], 4) == 1.6
    catch
        testBool = true
    end
    @test testBool
    testBool = false
    try
        read_time(result_files[2], 4) == 1.6
    catch
        testBool = true
    end
    @test testBool
    IO.close_result_files(result_files)

    rm(filename1 * ".e")
    rm(filename2 * ".e")
end

