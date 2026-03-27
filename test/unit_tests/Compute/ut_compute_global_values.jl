# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

#using Test

#using PeriLab
@testset "ut_global_value_sum" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(4)
    PeriLab.Data_Manager.set_glob_to_loc(Dict(1 => 1, 2 => 2, 3 => 3, 4 => 4))
    nodes = Vector{Int64}(1:4)
    (forcesN,
     forcesNP1) = PeriLab.Data_Manager.create_node_vector_field("Forces", Float64,
                                                                3)
    forcesNP1[1, 1:3] .= 1:3
    forcesNP1[3, 1:3] .= 5.2
    forcesNP1[4, 1] = -20
    @test PeriLab.IO.global_value_sum(forcesNP1, 1, nodes) == -13.8
    @test PeriLab.IO.global_value_sum(forcesNP1, 2, nodes) == 7.2
    @test PeriLab.IO.global_value_sum(forcesNP1, 3, nodes) == 8.2
    @test PeriLab.IO.calculate_nodelist("Forces", 1, "Sum", nodes) == (-13.8, 4)
    @test PeriLab.IO.calculate_nodelist("Forces", 2, "Sum", nodes) == (7.2, 4)
    @test PeriLab.IO.calculate_nodelist("Forces", 3, "Sum", nodes) == (8.2, 4)

    nodes = Vector{Int64}(1:2)
    @test PeriLab.IO.global_value_sum(forcesNP1, 1, nodes) == 1
    @test PeriLab.IO.global_value_sum(forcesNP1, 2, nodes) == 2
    @test PeriLab.IO.global_value_sum(forcesNP1, 3, nodes) == 3
    @test PeriLab.IO.calculate_nodelist("Forces", 1, "Sum", nodes) == (1, 2)
    @test PeriLab.IO.calculate_nodelist("Forces", 2, "Sum", nodes) == (2, 2)
    @test PeriLab.IO.calculate_nodelist("Forces", 3, "Sum", nodes) == (3, 2)
    disp = PeriLab.Data_Manager.create_constant_node_scalar_field("Disp", Float64)
    disp[1:4] = 1:4
    nodes = Vector{Int64}(1:4)
    @test PeriLab.IO.global_value_sum(disp, 1, nodes) == 10
    @test PeriLab.IO.calculate_nodelist("Disp", 1, "Sum", nodes) == (10, 4)
    nodes = Vector{Int64}(2:3)
    @test PeriLab.IO.global_value_sum(disp, 1, nodes) == 5
    @test PeriLab.IO.calculate_nodelist("Disp", 1, "Sum", nodes) == (5, 2)
    @test_logs (:error, "Unknown calculation type ") PeriLab.IO.calculate_nodelist("Disp",
                                                                                   1, "",
                                                                                   nodes)
    @test_logs (:error, "Field not there does not exists for computation") PeriLab.IO.calculate_nodelist("not there",
                                                                                                         1,
                                                                                                         "Sum",
                                                                                                         nodes)

    PeriLab.Data_Manager.set_glob_to_loc(Dict(1 => 1, 2 => 2))
    nodes = Vector{Int64}(3:4)
    empty_nodes::Vector{Int64} = []

    @test PeriLab.IO.calculate_nodelist("Disp", 1, "Sum", nodes) == (7, 2)
    @test PeriLab.IO.calculate_nodelist("Disp", 1, "Sum", empty_nodes) == (0, 0)
    @test PeriLab.IO.calculate_nodelist("Disp", 1, "Average", nodes) == (3.5, 2)
    @test PeriLab.IO.calculate_nodelist("Disp", 1, "Average", empty_nodes) == (0, 0)
    @test PeriLab.IO.calculate_nodelist("Disp", 1, "Minimum", nodes) == (3, 2)
    @test PeriLab.IO.calculate_nodelist("Disp", 1, "Minimum", empty_nodes) ==
          (Inf, 0)
    @test PeriLab.IO.calculate_nodelist("Disp", 1, "Maximum", nodes) == (4, 2)
    @test PeriLab.IO.calculate_nodelist("Disp", 1, "Maximum", empty_nodes) ==
          (-Inf, 0)

    matrix = PeriLab.Data_Manager.create_constant_node_tensor_field("Matrix", Float64, 3)
    matrix[:, 1, 2] .= 4
    nodes = Vector{Int64}(1:2)
    @test PeriLab.IO.calculate_nodelist("Matrix", [1, 2], "Sum", nodes) == (8, 2)
end

@testset "ut_global_value_max" begin
    PeriLab.Data_Manager.set_num_controller(4)
    PeriLab.Data_Manager.set_glob_to_loc(Dict(1 => 1, 2 => 2, 3 => 3, 4 => 4))
    nodes = Vector{Int64}(1:4)
    (forcesN,
     forcesNP1) = PeriLab.Data_Manager.create_node_vector_field("Forces", Float64,
                                                                3)

    forcesNP1[1, 1:3] .= 1:3
    forcesNP1[3, 1:3] .= 5.2
    forcesNP1[4, 1] = -20

    @test PeriLab.IO.global_value_max(forcesNP1, 1, nodes) == 5.2
    @test PeriLab.IO.global_value_max(forcesNP1, 2, nodes) == 5.2
    @test PeriLab.IO.global_value_max(forcesNP1, 3, nodes) == 5.2
    @test PeriLab.IO.calculate_nodelist("Forces", 1, "Maximum", nodes) == (5.2, 4)
    @test PeriLab.IO.calculate_nodelist("Forces", 1, "Maximum", nodes) == (5.2, 4)
    @test PeriLab.IO.calculate_nodelist("Forces", 1, "Maximum", nodes) == (5.2, 4)
    nodes = Vector{Int64}(1:2)
    @test PeriLab.IO.global_value_max(forcesNP1, 1, nodes) == 1
    @test PeriLab.IO.global_value_max(forcesNP1, 2, nodes) == 2
    @test PeriLab.IO.global_value_max(forcesNP1, 3, nodes) == 3
    @test PeriLab.IO.calculate_nodelist("Forces", 1, "Maximum", nodes) == (1, 2)
    @test PeriLab.IO.calculate_nodelist("Forces", 2, "Maximum", nodes) == (2, 2)
    @test PeriLab.IO.calculate_nodelist("Forces", 3, "Maximum", nodes) == (3, 2)
    disp = PeriLab.Data_Manager.create_constant_node_scalar_field("Disp", Float64)
    disp[1:4] = 1:4
    nodes = Vector{Int64}(1:4)
    @test PeriLab.IO.global_value_max(disp, 1, nodes) == 4
    @test PeriLab.IO.calculate_nodelist("Disp", 1, "Maximum", nodes) == (4, 4)
    nodes = Vector{Int64}(2:3)
    @test PeriLab.IO.global_value_max(disp, 1, nodes) == 3
    @test PeriLab.IO.calculate_nodelist("Disp", 1, "Maximum", nodes) == (3, 2)
    @test_logs (:error, "Unknown calculation type ") PeriLab.IO.calculate_nodelist("Disp",
                                                                                   1, "",
                                                                                   nodes)
    @test_logs (:error, "Field not there does not exists for computation") PeriLab.IO.calculate_nodelist("not there",
                                                                                                         1,
                                                                                                         "Maximum",
                                                                                                         nodes)
end

@testset "ut_global_value_min" begin
    PeriLab.Data_Manager.set_num_controller(4)
    PeriLab.Data_Manager.set_glob_to_loc(Dict(1 => 1, 2 => 2, 3 => 3, 4 => 4))
    nodes = Vector{Int64}(1:4)
    (forcesN,
     forcesNP1) = PeriLab.Data_Manager.create_node_vector_field("Forces", Float64,
                                                                3)
    forcesNP1[1, 1:3] .= 1:3
    forcesNP1[3, 1:3] .= 5.2
    forcesNP1[4, 1] = -20
    @test PeriLab.IO.global_value_min(forcesNP1, 1, nodes) == -20
    @test PeriLab.IO.global_value_min(forcesNP1, 2, nodes) == 0
    @test PeriLab.IO.global_value_min(forcesNP1, 3, nodes) == 0
    @test PeriLab.IO.calculate_nodelist("Forces", 1, "Minimum", nodes) == (-20, 4)
    @test PeriLab.IO.calculate_nodelist("Forces", 2, "Minimum", nodes) == (0, 4)
    @test PeriLab.IO.calculate_nodelist("Forces", 3, "Minimum", nodes) == (0, 4)
    nodes = Vector{Int64}(1:2)
    @test PeriLab.IO.global_value_min(forcesNP1, 1, nodes) == 0
    @test PeriLab.IO.global_value_min(forcesNP1, 2, nodes) == 0
    @test PeriLab.IO.global_value_min(forcesNP1, 3, nodes) == 0
    @test PeriLab.IO.calculate_nodelist("Forces", 1, "Minimum", nodes) == (0, 2)
    @test PeriLab.IO.calculate_nodelist("Forces", 2, "Minimum", nodes) == (0, 2)
    @test PeriLab.IO.calculate_nodelist("Forces", 3, "Minimum", nodes) == (0, 2)
    disp = PeriLab.Data_Manager.create_constant_node_scalar_field("Disp", Float64)
    disp[1:4] = 1:4
    nodes = Vector{Int64}(1:4)
    @test PeriLab.IO.global_value_min(disp, 1, nodes) == 1
    @test PeriLab.IO.calculate_nodelist("Disp", 1, "Minimum", nodes) == (1, 4)
    nodes = Vector{Int64}(2:3)
    @test PeriLab.IO.global_value_min(disp, 1, nodes) == 2
    @test PeriLab.IO.calculate_nodelist("Disp", 1, "Minimum", nodes) == (2, 2)
    @test_logs (:error, "Unknown calculation type ") PeriLab.IO.calculate_nodelist("Disp",
                                                                                   1, "",
                                                                                   nodes)
    @test_logs (:error, "Field not there does not exists for computation") PeriLab.IO.calculate_nodelist("not there",
                                                                                                         1,
                                                                                                         "Minimum",
                                                                                                         nodes)
end
@testset "ut_global_value_avg" begin
    PeriLab.Data_Manager.set_num_controller(4)
    PeriLab.Data_Manager.set_glob_to_loc(Dict(1 => 1, 2 => 2, 3 => 3, 4 => 4))
    nodes = Vector{Int64}(1:4)
    (forcesN,
     forcesNP1) = PeriLab.Data_Manager.create_node_vector_field("Forces", Float64,
                                                                3)
    forcesNP1[1, 1:3] .= 1:3
    forcesNP1[3, 1:3] .= 5.2
    forcesNP1[4, 1] = -20
    @test PeriLab.IO.global_value_avg(forcesNP1, 1, nodes) == -13.8 / 4
    @test PeriLab.IO.global_value_avg(forcesNP1, 2, nodes) == 7.2 / 4
    @test PeriLab.IO.global_value_avg(forcesNP1, 3, nodes) == 8.2 / 4
    @test PeriLab.IO.calculate_nodelist("Forces", 1, "Average", nodes) ==
          (-13.8 / 4, 4)
    @test PeriLab.IO.calculate_nodelist("Forces", 2, "Average", nodes) ==
          (7.2 / 4, 4)
    @test PeriLab.IO.calculate_nodelist("Forces", 3, "Average", nodes) ==
          (8.2 / 4, 4)
    nodes = Vector{Int64}(1:2)
    @test PeriLab.IO.global_value_avg(forcesNP1, 1, nodes) == 1 / 2
    @test PeriLab.IO.global_value_avg(forcesNP1, 2, nodes) == 2 / 2
    @test PeriLab.IO.global_value_avg(forcesNP1, 3, nodes) == 3 / 2
    @test PeriLab.IO.calculate_nodelist("Forces", 1, "Average", nodes) == (1 / 2, 2)
    @test PeriLab.IO.calculate_nodelist("Forces", 2, "Average", nodes) == (2 / 2, 2)
    @test PeriLab.IO.calculate_nodelist("Forces", 3, "Average", nodes) == (3 / 2, 2)
    disp = PeriLab.Data_Manager.create_constant_node_scalar_field("Disp", Float64)
    disp[1:4] = 1:4
    nodes = Vector{Int64}(1:4)
    @test PeriLab.IO.global_value_avg(disp, 1, nodes) == 10 / 4
    @test PeriLab.IO.calculate_nodelist("Disp", 1, "Average", nodes) == (10 / 4, 4)
    nodes = Vector{Int64}(2:3)
    @test PeriLab.IO.global_value_avg(disp, 1, nodes) == 5 / 2
    @test PeriLab.IO.calculate_nodelist("Disp", 1, "Average", nodes) == (5 / 2, 2)
    @test_logs (:error, "Unknown calculation type ") PeriLab.IO.calculate_nodelist("Disp",
                                                                                   1, "",
                                                                                   nodes)
    @test_logs (:error,
                "Field not there does not exists for computation") PeriLab.IO.calculate_nodelist("not there",
                                                                                                 1,
                                                                                                 "Average",
                                                                                                 nodes)
end

@testset "ut_calculate_block" begin
    PeriLab.Data_Manager.create_constant_node_scalar_field("Block_Id", Int64;
                                                           default_value = 1)
    PeriLab.Data_Manager.set_glob_to_loc(Dict(1 => 1, 2 => 2, 3 => 3, 4 => 4))
    nodes = Vector{Int64}(1:4)
    (forcesN,
     forcesNP1) = PeriLab.Data_Manager.create_node_vector_field("Forces", Float64,
                                                                3)
    forcesNP1[1, 1:3] .= 1:3
    forcesNP1[3, 1:3] .= 5.2
    forcesNP1[4, 1] = -20
    @test_logs (:error,
                "Field no field does not exists for computation") PeriLab.IO.calculate_block("no field",
                                                                                             1,
                                                                                             "sum",
                                                                                             1)
    @test PeriLab.IO.calculate_block("Forces", 1, "Sum", 1) == (-13.8, 4)
end
