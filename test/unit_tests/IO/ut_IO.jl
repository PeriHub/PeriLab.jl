# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# include("../../../src/IO/IO.jl")
# using PeriLab
using MPI
using TimerOutputs
using Test

test_data_manager = PeriLab.Data_Manager
test_data_manager.initialize_data()
filename1 = "test1"
filename2 = "test2"
filename3 = "test3"
dof = 2
nnodes = 5
comm = MPI.COMM_WORLD
test_data_manager.set_num_controller(nnodes)
test_data_manager.set_num_responder(0)
test_data_manager.set_comm(comm)
test_data_manager.set_dof(dof)
test_data_manager.set_max_rank(1)
test_data_manager.set_distribution([1, 2, 3, 4, 5])
test_data_manager.create_constant_node_vector_field("Coordinates", Float64, 2)
coordinates = test_data_manager.get_field("Coordinates")
test_data_manager.create_constant_node_scalar_field("Block_Id", Int64)
block_Id = test_data_manager.get_field("Block_Id")
test_data_manager.create_node_vector_field("Displacements", Float64, 2)
test_data_manager.create_node_vector_field("Forces", Float64, 6)

block_list = ["block_1", "block_2"]
test_data_manager.set_block_name_list(block_list)
test_data_manager.set_block_id_list([1, 2])

params = Dict("Outputs" => Dict("Output1" => Dict("Output Filename" => filename1,
                                                  "Flush File" => false,
                                                  "Output Variables" => Dict("Forces" => true)),
                                "Output2" => Dict("Output Filename" => filename2,
                                                  "Flush File" => false,
                                                  "Output Variables" => Dict("Displacements" => true,
                                                                             "Forces" => true)),
                                "Output3" => Dict("Output Filename" => filename3,
                                                  "Output File Type" => "CSV",
                                                  "Output Variables" => Dict("External_Displacement" => true))),
              "Compute Class Parameters" => Dict("External_Displacement" => Dict("Block" => "block_1",
                                                                                 "Calculation Type" => "Maximum",
                                                                                 "Compute Class" => "Block_Data",
                                                                                 "Variable" => "Displacements")))
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
    output = PeriLab.IO.get_results_mapping(params, "")
    @test sort(collect(keys(output[2]["Fields"]))) ==
          ["Forcesxx", "Forcesxy", "Forcesxz", "Forcesyx", "Forcesyy", "Forcesyz"]
    @test sort(collect(keys(output[3]["Fields"]))) == [
        "Displacementsx",
        "Displacementsy",
        "Forcesxx",
        "Forcesxy",
        "Forcesxz",
        "Forcesyx",
        "Forcesyy",
        "Forcesyz"
    ]
    for i in 2:3
        dof_force = 0
        dof_disp::Int64 = 0
        for entry in keys(sort(output[i]["Fields"]))
            if occursin("Forces", entry)
                dof_force += 1
                @test output[i]["Fields"][entry]["fieldname"] == "Forces"
                @test output[i]["Fields"][entry]["time"] == "NP1"
                @test output[i]["Fields"][entry]["dof"] == dof_force
                @test output[i]["Fields"][entry]["type"] == Float64
            else
                dof_disp += 1
                @test output[i]["Fields"][entry]["fieldname"] == "Displacements"
                @test output[i]["Fields"][entry]["time"] == "NP1"
                @test output[i]["Fields"][entry]["dof"] == dof_disp
                @test output[i]["Fields"][entry]["type"] == Float64
            end
        end
    end
end

@testset "ut_init_write_result_and_write_results" begin
    result_files,
    outputs = PeriLab.IO.init_write_results(params, "", "",
                                            "1.0.0")
    @test length(result_files) == 3
    @test length(result_files[2]["file"].nodal_var_name_dict) == 6
    entries = collect(keys(result_files[2]["file"].nodal_var_name_dict))
    @test sort(entries) ==
          ["Forcesxx", "Forcesxy", "Forcesxz", "Forcesyx", "Forcesyy", "Forcesyz"]
    @test length(result_files[3]["file"].nodal_var_name_dict) == 8
    entries = collect(keys(result_files[3]["file"].nodal_var_name_dict))
    @test sort(entries) == [
        "Displacementsx",
        "Displacementsy",
        "Forcesxx",
        "Forcesxy",
        "Forcesxz",
        "Forcesyx",
        "Forcesyy",
        "Forcesyz"
    ]

    coords = vcat(transpose(coordinates))
    for exo in result_files[2:3]
        exo_coords = read_coordinates(exo["file"])
        exo_nsets = read_sets(exo["file"], NodeSet)
        @test coords == exo_coords
        # @test exo_nsets == []
    end
    for i in 2:3
        dofForce = 0
        dofDisp = 0
        for entry in keys(sort(outputs[i]["Fields"]))
            if occursin("Forces", entry)
                dofForce += 1
                @test outputs[i]["Fields"][entry]["fieldname"] == "Forces"
                @test outputs[i]["Fields"][entry]["time"] == "NP1"
                @test outputs[i]["Fields"][entry]["dof"] == dofForce
                @test outputs[i]["Fields"][entry]["type"] == Float64
            else
                dofDisp += 1
                @test outputs[i]["Fields"][entry]["fieldname"] == "Displacements"
                @test outputs[i]["Fields"][entry]["time"] == "NP1"
                @test outputs[i]["Fields"][entry]["dof"] == dofDisp
                @test outputs[i]["Fields"][entry]["type"] == Float64
            end
        end
    end
    test_data_manager.data["output_frequency"] = [
        Dict{String,Int64}("Counter" => 0, "Output Frequency" => 1, "Step" => 1),
        Dict{String,Int64}("Counter" => 0, "Output Frequency" => 1, "Step" => 1),
        Dict{String,Int64}("Counter" => 0, "Output Frequency" => 1, "Step" => 1)
    ]
    PeriLab.IO.write_results(result_files, 1.5, 0.0, outputs)

    @test read_time(result_files[2]["file"], 2) == 1.5
    @test read_time(result_files[3]["file"], 2) == 1.5
    PeriLab.IO.write_results(result_files, 1.6, 0.0, outputs)

    @test read_time(result_files[2]["file"], 3) == 1.6
    @test read_time(result_files[3]["file"], 3) == 1.6
    testBool = false
    try
        read_time(result_files[2]["file"], 4) == 1.6
    catch
        testBool = true
    end
    @test testBool
    testBool = false
    try
        read_time(result_files[3]["file"], 4) == 1.6
    catch
        testBool = true
    end
    @test testBool
    @test PeriLab.IO.close_result_files(result_files[2:3])
    PeriLab.IO.merge_exodus_files(result_files[2:3], "")

    rm(filename1 * ".e")
    rm(filename2 * ".e")
    # rm(filename3 * ".csv")
end

# @testset "ut_show_block_summary" begin
#     test_data_manager.set_block_name_list([1])
#     solver_options = Dict("Material Models" => true, "Damage Models" => true, "Additive Models" => true, "Thermal Models" => true)
#     params = Dict("Blocks" => Dict("block_1" => Dict("Material Models" => true, "Damage Models" => true, "Additive Models" => true, "Thermal Models" => true)))
#     PeriLab.IO.show_block_summary(solver_options, params, "test.log", comm)
# end

# @testset "ut_show_mpi_summary" begin
#     test_data_manager.set_block_name_list([1])
#     PeriLab.IO.show_mpi_summary("test.log", comm)
# end

@testset "ut_init_orientations" begin
    angles = test_data_manager.create_constant_node_scalar_field("Angles", Float64;
                                                                 default_value = 90)
    test_data_manager.set_rotation(true)
    PeriLab.IO.init_orientations()
    orientations = test_data_manager.get_field("Orientations")
    @test isapprox(orientations[1, 1], 0; atol = 0.00001)
    @test isapprox(orientations[1, 2], 1; atol = 0.00001)
    @test isapprox(orientations[1, 3], 0; atol = 0.00001)
end

@testset "ut_init_orientations_3d" begin
    test_data_manager.initialize_data()
    dof = 3
    nnodes = 5
    comm = MPI.COMM_WORLD
    test_data_manager.set_num_controller(nnodes)
    test_data_manager.set_num_responder(0)
    test_data_manager.set_comm(comm)
    test_data_manager.set_dof(dof)
    # test_data_manager.set_max_rank(1)
    # test_data_manager.set_distribution([1, 2, 3, 4, 5])
    # test_data_manager.create_constant_node_field("Coordinates", Float64, 2)
    # coordinates = test_data_manager.get_field("Coordinates")
    # coordinates[1, 1] = 0
    # coordinates[1, 2] = 0
    # coordinates[1, 3] = 0
    # coordinates[2, 1] = 1
    # coordinates[2, 2] = 0
    # coordinates[2, 3] = 0
    # coordinates[3, 1] = 0
    # coordinates[3, 2] = 1
    # coordinates[3, 3] = 0
    # coordinates[4, 1] = 1
    # coordinates[4, 2] = 1
    # coordinates[4, 3] = 0
    # coordinates[5, 1] = 2
    # coordinates[5, 2] = 2
    # coordinates[5, 3] = 0
    angles = test_data_manager.create_constant_node_vector_field("Angles", Float64, 3;
                                                                 default_value = 90)
    test_data_manager.set_rotation(true)
    PeriLab.IO.init_orientations()
    orientations = test_data_manager.get_field("Orientations")
    @test isapprox(orientations[1, 1], 0; atol = 0.00001)
    @test isapprox(orientations[1, 2], 0; atol = 0.00001)
    @test isapprox(orientations[1, 3], 1.0; atol = 0.00001)
end

@testset "ut_get_mpi_rank_string" begin
    @test PeriLab.IO.get_mpi_rank_string(0, 2) == "0"
    @test PeriLab.IO.get_mpi_rank_string(1, 2) == "1"
    @test PeriLab.IO.get_mpi_rank_string(0, 11) == "00"
    @test PeriLab.IO.get_mpi_rank_string(10, 11) == "10"
    @test PeriLab.IO.get_mpi_rank_string(22, 100) == "022"
end

@testset "ut_show_block_summary" begin
    test_data_manager.initialize_data()
    dof = 2
    nnodes = 5
    comm = MPI.COMM_WORLD
    test_data_manager.set_num_controller(nnodes)
    test_data_manager.set_num_responder(0)
    test_data_manager.set_comm(comm)
    test_data_manager.set_dof(dof)
    test_data_manager.set_max_rank(1)
    test_data_manager.set_distribution([1, 2, 3, 4, 5])
    test_data_manager.create_constant_node_scalar_field("Block_Id", Int64)
    block_Id = test_data_manager.get_field("Block_Id")
    block_Id .+= 1
    block_Id[end] = 2
    test_data_manager.set_block_name_list(["block_1", "block_2"])
    solver_options = Dict("Models" => ("Material Models" => true,
                                       "Damage Models" => true,
                                       "Additive Models" => true,
                                       "Thermal Models" => true,
                                       "Degradation Models" => true))
    params = Dict("Blocks" => Dict("block_1" => Dict("Material Models" => true,
                                                     "Damage Models" => true,
                                                     "Additive Models" => true,
                                                     "Thermal Models" => true,
                                                     "Degradation Models" => true),
                                   "block_2" => Dict("Material Models" => true,
                                                     "Damage Models" => false,
                                                     "Additive Models" => false,
                                                     "Thermal Models" => false,
                                                     "Degradation Models" => false)))
    PeriLab.IO.show_block_summary(solver_options,
                                  params,
                                  "",
                                  false,
                                  comm)
end

@testset "ut_show_mpi_summary" begin
    # Not fully tested yet because MPI.size=1
    @test PeriLab.IO.show_mpi_summary("", false, comm) == 1
end

@testset "ut_clearNP1" begin
    @test PeriLab.IO.clearNP1("ForceNP1") == "Force"
    @test PeriLab.IO.clearNP1("Force") == "Force"
end
