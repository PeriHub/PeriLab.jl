# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
using TimerOutputs

using Exodus
test_data_manager = PeriLab.Data_Manager
test_data_manager.initialize_data()
@testset "ut_get_block_nodes" begin
    block_Id = [1, 1, 2, 2, 2, 3, 3, 3, 1, 3, 3, 4]
    test = PeriLab.IO.get_block_nodes(block_Id, 1)
    @test test == [1 2 9]
    test = PeriLab.IO.get_block_nodes(block_Id, 2)
    @test test == [3 4 5]
    test = PeriLab.IO.get_block_nodes(block_Id, 3)
    @test test == [6 7 8 10 11]
    test = PeriLab.IO.get_block_nodes(block_Id, 4)
    @test test[1] == 12
end

@testset "ut_paraview_specifics" begin
    @test PeriLab.IO.paraview_specifics(1) == "x"
    @test PeriLab.IO.paraview_specifics(2) == "y"
    @test PeriLab.IO.paraview_specifics(3) == "z"
end

@testset "ut_get_paraview_coordinates" begin
    for i in 1:3
        @test PeriLab.IO.get_paraview_coordinates(1, i) == "x"
        @test PeriLab.IO.get_paraview_coordinates(2, i) == "y"
        @test PeriLab.IO.get_paraview_coordinates(3, i) == "z"
    end
    @test isnothing(PeriLab.IO.get_paraview_coordinates(3, 10))

    for ref in 4:9
        for i in 1:3
            for j in 1:3
                @test PeriLab.IO.get_paraview_coordinates((i - 1) * 3 + j, ref) ==
                      PeriLab.IO.paraview_specifics(i) * PeriLab.IO.paraview_specifics(j)
            end
        end
    end
end

if !isdir("tmp")
    mkdir("tmp")
end

topology = test_data_manager.create_constant_free_size_field("FE Topology", Int64, (2, 4))
topology[1, 1] = 1
topology[1, 2] = 2
topology[1, 3] = 3
topology[1, 4] = 4
topology[2, 1] = 2
topology[2, 2] = 3
topology[2, 3] = 4
topology[2, 4] = 5
@testset "ut_create_result_file" begin
    filename1 = "./tmp/" * "test.e"
    filename2 = "./tmp/" * "test2.e"
    nnodes = 4
    dof = 3
    exo = PeriLab.IO.create_result_file(filename1, 5, dof, 1, 0, 2, topology)
    exo = PeriLab.IO.create_result_file(filename2, nnodes, dof, 1, 0)
    @test isfile(filename2)
    @test exo["file"].file_name == filename2
    # @test num_dim(exo["file"].init) == dof
    # @test num_nodes(exo["file"].init) == nnodes
    # @test num_node_sets(exo["file"].init) == 0
    # @test num_elem_blks(exo["file"].init) == 1

    close(exo["file"])
    rm(filename2)
    fid = open(filename2, "w")
    close(fid)
    nnodes = 300
    dof = 2
    @test isfile(filename2)
    exo = PeriLab.IO.create_result_file(filename2, nnodes, dof, 3, 2)
    @test isfile(filename2)
    @test exo["file"].file_name == filename2
    # @test num_dim(exo["file"].init) == dof
    # @test num_nodes(exo["file"].init) == nnodes
    # @test num_node_sets(exo["file"].init) == 2
    # @test num_elem_blks(exo["file"].init) == 3
    close(exo["file"])
    rm(filename1)
    rm(filename2)
end
filename = "./tmp/" * "test_2.e"
filename2 = "./tmp/" * "test_22.e"
nnodes = 5
dof = 2
test_data_manager.set_num_controller(nnodes)
test_data_manager.set_dof(dof)
coordinates = test_data_manager.create_constant_node_vector_field("Coordinates", Float64, 2)
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
test_data_manager.create_constant_node_scalar_field("Block_Id", Int64)
block_Id = test_data_manager.get_field("Block_Id")
block_Id .+= 1
block_Id[end] = 2
test_data_manager.create_constant_node_scalar_field("FEM", Bool; default_value = true)
fem_block = test_data_manager.get_field("FEM")
#outputs = ["Displacements", "Forces"]
test_data_manager.set_nset("Nset_1", [1, 2])
test_data_manager.set_nset("Nset_2", [5])

nsets = test_data_manager.get_nsets()
coords = vcat(transpose(coordinates))
outputs = Dict("Fields" => Dict("Forcesxx" => Dict("fieldname" => "Forces",
                                                   "time" => "NP1",
                                                   "global_var" => false,
                                                   "dof" => 1,
                                                   "type" => Float64),
                                "Forcesxy" => Dict("fieldname" => "Forces",
                                                   "time" => "NP1",
                                                   "global_var" => false,
                                                   "dof" => 1,
                                                   "type" => Float64),
                                "Forcesxz" => Dict("fieldname" => "Forces",
                                                   "time" => "NP1",
                                                   "global_var" => false,
                                                   "dof" => 1,
                                                   "type" => Float64),
                                "Forcesyx" => Dict("fieldname" => "Forces",
                                                   "time" => "NP1",
                                                   "global_var" => false,
                                                   "dof" => 1,
                                                   "type" => Float64),
                                "Forcesyy" => Dict("fieldname" => "Forces",
                                                   "time" => "NP1",
                                                   "global_var" => false,
                                                   "dof" => 1,
                                                   "type" => Float64),
                                "Forcesyz" => Dict("fieldname" => "Forces",
                                                   "time" => "NP1",
                                                   "global_var" => false,
                                                   "dof" => 1,
                                                   "type" => Float64),
                                "Displacements" => Dict("fieldname" => "Displacements",
                                                        "time" => "NP1",
                                                        "global_var" => false,
                                                        "dof" => 1,
                                                        "type" => Float64),
                                "External_Displacements" => Dict("fieldname" => "Displacements",
                                                                 "time" => "NP1",
                                                                 "global_var" => true,
                                                                 "dof" => 1,
                                                                 "type" => Float64,
                                                                 "compute_params" => Dict("Compute Class" => "Block_Data",
                                                                                          "Calculation Type" => "Maximum",
                                                                                          "Block" => "block_1",
                                                                                          "Variable" => "Displacements")),
                                "External_Forces" => Dict("fieldname" => "Forces",
                                                          "time" => "NP1",
                                                          "global_var" => true,
                                                          "dof" => 3,
                                                          "type" => Float64,
                                                          "compute_params" => Dict("Compute Class" => "Node_Set_Data",
                                                                                   "Calculation Type" => "Minimum",
                                                                                   "Node Set" => 1,
                                                                                   "Variable" => "DisplacementsNP1"))))
computes = Dict("Fields" => Dict("External_Displacements" => Dict("fieldname" => "Displacements",
                                                                  "time" => "NP1",
                                                                  "global_var" => true,
                                                                  "dof" => 1,
                                                                  "type" => Float64,
                                                                  "compute_params" => Dict("Compute Class" => "Block_Data",
                                                                                           "Calculation Type" => "Maximum",
                                                                                           "Block" => "block_1",
                                                                                           "Variable" => "DisplacementsNP1")),
                                 "External_Forces" => Dict("fieldname" => "Forces",
                                                           "time" => "NP1",
                                                           "global_var" => true,
                                                           "dof" => 3,
                                                           "type" => Float64,
                                                           "compute_params" => Dict("Compute Class" => "Node_Set_Data",
                                                                                    "Calculation Type" => "Minimum",
                                                                                    "Node Set" => 1,
                                                                                    "Variable" => "DisplacementsNP1"))))

exo1 = PeriLab.IO.create_result_file(filename2,
                                     nnodes,
                                     dof,
                                     maximum(block_Id),
                                     length(nsets),
                                     2,
                                     topology)
exo1["file"] = PeriLab.IO.init_results_in_exodus(exo1["file"],
                                                 outputs,
                                                 coords,
                                                 block_Id[1:nnodes],
                                                 ["Block_1", "Block_2"],
                                                 nsets,
                                                 [1, 2, 3, 4, 5],
                                                 "1.0.0",
                                                 fem_block,
                                                 topology,
                                                 [1, 2])
rm(filename2)
exo = PeriLab.IO.create_result_file(filename, nnodes, dof, maximum(block_Id), length(nsets))
exo["file"] = PeriLab.IO.init_results_in_exodus(exo["file"],
                                                outputs,
                                                coords,
                                                block_Id[1:nnodes],
                                                ["Block_1", "Block_2"],
                                                nsets,
                                                [1, 2, 3, 4, 5],
                                                "1.0.0")
result_files = []
push!(result_files, exo)
result_files[1]["file"] = PeriLab.IO.write_step_and_time(result_files[1]["file"], 2, 2.2)
result_files[1]["file"] = PeriLab.IO.write_step_and_time(result_files[1]["file"], 3, 3.7)
result_files[1]["file"] = PeriLab.IO.write_step_and_time(result_files[1]["file"], 4, 4.7)
result_files[1]["file"] = PeriLab.IO.write_step_and_time(result_files[1]["file"], 5, 5.7)
result_files[1]["file"] = PeriLab.IO.write_step_and_time(result_files[1]["file"], 6, 6.7)

@testset "ut_init_results_in_exodus" begin
    # @test exo["file"].init.num_dim == dof
    @test length(exo["file"].nodal_var_name_dict) == 7
    entries = collect(keys(exo["file"].nodal_var_name_dict))
    ref = collect(keys(outputs["Fields"]))
    @test sort(entries) == deleteat!(sort(ref), 2:3)
    exo_coords = read_coordinates(exo["file"])
    exo_nsets = read_sets(exo["file"], NodeSet)
    @test length(exo_nsets) == length(nsets)
    @test coords == exo_coords
    @warn "Info test deactivated"
    # @test ["PeriLab Version " * string(Pkg.project().version) * ", under BSD License", "Copyright (c) 2023, Christian Willberg, Jan-Timo Hesse", "compiled with Julia Version " * string(VERSION)] == read_info(exo["file"])
    @test read_number_of_time_steps(exo["file"]) == 6
    @test read_time(exo["file"], 2) == 2.2
    @test read_time(exo["file"], 3) == 3.7
    @test read_time(exo["file"], 4) == 4.7
    @test read_time(exo["file"], 5) == 5.7
    @test read_time(exo["file"], 6) == 6.7
    @test read_name(exo["file"], Block, 1) == "Block_1"
    @test read_name(exo["file"], Block, 2) == "Block_2"
end

test_data_manager.create_node_vector_field("Forces", Float64, 6)
test_data_manager.create_node_scalar_field("Displacements", Float64)
force = test_data_manager.get_field("Forces", "NP1")
disp = test_data_manager.get_field("Displacements", "NP1")
force[5, 1:6] .= 3.3
force[1:3, 6] .= 2.3
disp[1] = 3
disp[2] = 3.00001
disp[3] = 2.1
disp[4] = -1.8
disp[5] = 0

nodal_outputs = Dict(key => value
                     for (key, value) in outputs["Fields"] if (!value["global_var"]))
exo["file"] = PeriLab.IO.write_nodal_results_in_exodus(exo["file"], 2, nodal_outputs)

test_disp_step_zero = read_values(exo["file"], NodalVariable, 1, 1, "Displacements")

test_disp_step_one = read_values(exo["file"], NodalVariable, 2, 1, 1)
@testset "ut_write_results_in_exodus" begin
    @test isapprox(test_disp_step_zero .+ 1, ones(5))
    for id in eachindex(test_disp_step_one)
        if disp[id] != 0
            @test test_disp_step_one[id] / disp[id] - 1 < 1e-8
        else
            @test test_disp_step_one[id] == disp[id]
        end
    end

    ftest = read_values(exo["file"], NodalVariable, 2, 1, "Forcesxx")
    for id in 1:4
        @test ftest[id] == 0
    end
    @test ftest[5] / 3.3 - 1 < 1e-8
    ftest = read_values(exo["file"], NodalVariable, 2, 1, "Forcesxy")
    for id in 1:4
        @test ftest[id] == 0
    end
    @test ftest[5] / 3.3 - 1 < 1e-8
    ftest = read_values(exo["file"], NodalVariable, 2, 1, "Forcesxz")
    for id in 1:4
        @test ftest[id] == 0
    end
    @test ftest[5] / 3.3 - 1 < 1e-8
    ftest = read_values(exo["file"], NodalVariable, 2, 1, "Forcesyx")
    for id in 1:4
        @test ftest[id] == 0
    end
    @test ftest[5] / 3.3 - 1 < 1e-8

    ftest = read_values(exo["file"], NodalVariable, 2, 1, "Forcesyy")
    for id in 1:4
        @test ftest[id] == 0
    end
    @test ftest[5] / 3.3 - 1 < 1e-8

    ftest = read_values(exo["file"], NodalVariable, 2, 1, "Forcesyz")
    @test ftest[1] / 2.3 - 1 < 1e-8
    @test ftest[2] / 2.3 - 1 < 1e-8
    @test ftest[3] / 2.3 - 1 < 1e-8
    @test ftest[4] == 0
    @test ftest[5] / 3.3 - 1 < 1e-8
end

csvfilename = "./tmp/" * "test_2.csv"
csv_file = open(csvfilename, "w")
println(csv_file, "Test")
csv_file = PeriLab.IO.create_result_file(csvfilename, computes)
exo["file"] = PeriLab.IO.write_global_results_in_exodus(exo["file"], 2, [0.1, 0.2])

@testset "ut_write_global_results_in_exodus" begin
    global_vars = read_names(exo["file"], GlobalVariable)
    @test global_vars[1] == "External_Displacements"
    @test global_vars[2] == "External_Forces"

    ftest = read_values(exo["file"], GlobalVariable, 2)
    @test ftest[1] == 0.1
    @test ftest[2] == 0.2
end

@testset "ut_merge_exodus_file" begin
    merged = true
    try
        PeriLab.IO.merge_exodus_file(exo["filename"])
    catch
        merged = false
    end
    @test merged
end

PeriLab.IO.write_global_results_in_csv(csv_file["file"], 1.0, [0.1, 0.2])
#TODO: check if the csv file is correct
@testset "ut_write_global_results_in_csv" begin
    @test isfile(csv_file["filename"])
    @test csv_file["file"] isa IOStream
end

close(exo["file"])
close(csv_file["file"])
rm(filename)
rm(csvfilename)
