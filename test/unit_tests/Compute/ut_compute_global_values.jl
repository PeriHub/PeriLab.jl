# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
include("../../../src/Support/data_manager.jl")
include("../../../src/Compute/compute_global_values.jl")

@testset "ut_global_value_sum" begin
    test_Data_manager = Data_manager
    test_Data_manager.set_num_controller(4)
    nodes = Vector{Int64}(1:4)
    (forcesN, forcesNP1) = test_Data_manager.create_node_field("Forces", Float64, 3)
    forcesNP1[1, 1:3] .= 1:3
    forcesNP1[3, 1:3] .= 5.2
    forcesNP1[4, 1] = -20

    @test global_value_sum(forcesNP1, 1, nodes) == -13.8
    @test global_value_sum(forcesNP1, 2, nodes) == 7.2
    @test global_value_sum(forcesNP1, 3, nodes) == 8.2
    @test calculate_nodelist(test_Data_manager, "Forces", [1, 1], "Sum", nodes) == (-13.8, 4)
    @test calculate_nodelist(test_Data_manager, "Forces", [1, 2], "Sum", nodes) == (7.2, 4)
    @test calculate_nodelist(test_Data_manager, "Forces", [1, 3], "Sum", nodes) == (8.2, 4)
    nodes = Vector{Int64}(1:2)
    @test global_value_sum(forcesNP1, 1, nodes) == 1
    @test global_value_sum(forcesNP1, 2, nodes) == 2
    @test global_value_sum(forcesNP1, 3, nodes) == 3
    @test calculate_nodelist(test_Data_manager, "Forces", "Sum", 1, nodes) == (1, 2)
    @test calculate_nodelist(test_Data_manager, "Forces", "Sum", 2, nodes) == (2, 2)
    @test calculate_nodelist(test_Data_manager, "Forces", "Sum", 3, nodes) == (3, 2)
    disp = test_Data_manager.create_constant_node_field("Disp", Float64, 1)
    disp[1:4] = 1:4
    nodes = Vector{Int64}(1:4)
    @test global_value_sum(disp, 1, nodes) == 10
    @test calculate_nodelist(test_Data_manager, "Disp", "Sum", 1, nodes) == (10, 4)
    nodes = Vector{Int64}(2:3)
    @test global_value_sum(disp, 1, nodes) == 5
    @test calculate_nodelist(test_Data_manager, "Disp", "Sum", 1, nodes) == (5, 2)
    @test isnothing(calculate_nodelist(test_Data_manager, "Disp", "", 1, nodes))
    @test isnothing(calculate_nodelist(test_Data_manager, "not there", "Sum", 1, nodes))
end

@testset "ut_global_value_max" begin
    test_Data_manager = Data_manager
    test_Data_manager.set_num_controller(4)
    nodes = Vector{Int64}(1:4)
    (forcesN, forcesNP1) = test_Data_manager.create_node_field("Forces", Float64, 3)

    forcesNP1[1, 1:3] .= 1:3
    forcesNP1[3, 1:3] .= 5.2
    forcesNP1[4, 1] = -20

    @test global_value_max(forcesNP1, 1, nodes) == 5.2
    @test global_value_max(forcesNP1, 2, nodes) == 5.2
    @test global_value_max(forcesNP1, 3, nodes) == 5.2
    @test calculate_nodelist(test_Data_manager, "Forces", 1, "Maximum", nodes) == (5.2, 4)
    @test calculate_nodelist(test_Data_manager, "Forces", 1, "Maximum", nodes) == (5.2, 4)
    @test calculate_nodelist(test_Data_manager, "Forces", 1, "Maximum", nodes) == (5.2, 4)
    nodes = Vector{Int64}(1:2)
    @test global_value_max(forcesNP1, 1, nodes) == 1
    @test global_value_max(forcesNP1, 2, nodes) == 2
    @test global_value_max(forcesNP1, 3, nodes) == 3
    @test calculate_nodelist(test_Data_manager, "Forces", "Maximum", 1, nodes) == (1, 2)
    @test calculate_nodelist(test_Data_manager, "Forces", "Maximum", 2, nodes) == (2, 2)
    @test calculate_nodelist(test_Data_manager, "Forces", "Maximum", 3, nodes) == (3, 2)
    disp = test_Data_manager.create_constant_node_field("Disp", Float64, 1)
    disp[1:4] = 1:4
    nodes = Vector{Int64}(1:4)
    @test global_value_max(disp, 1, nodes) == 1
    @test calculate_nodelist(test_Data_manager, "Disp", "Maximum", 1, nodes) == (4, 4)
    nodes = Vector{Int64}(2:3)
    @test global_value_max(disp, 1, nodes) == 1
    @test calculate_nodelist(test_Data_manager, "Disp", "Maximum", 1, nodes) == (3, 2)
    @test isnothing(calculate_nodelist(test_Data_manager, "Disp", "", 1, nodes))
    @test isnothing(calculate_nodelist(test_Data_manager, "not there", "Maximum", 1, nodes))
end

@testset "ut_global_value_min" begin
    test_Data_manager = Data_manager
    test_Data_manager.set_num_controller(4)
    nodes = Vector{Int64}(1:4)
    (forcesN, forcesNP1) = test_Data_manager.create_node_field("Forces", Float64, 3)
    forcesNP1[1, 1:3] .= 1:3
    forcesNP1[3, 1:3] .= 5.2
    forcesNP1[4, 1] = -20
    @test global_value_min(forcesNP1, 1, nodes) == -20
    @test global_value_min(forcesNP1, 2, nodes) == 0
    @test global_value_min(forcesNP1, 3, nodes) == 0
    @test calculate_nodelist(test_Data_manager, "Forces", 1, "Minimum", nodes) == -20
    @test calculate_nodelist(test_Data_manager, "Forces", 1, "Minimum", nodes) == 0
    @test calculate_nodelist(test_Data_manager, "Forces", 1, "Minimum", nodes) == 0
    nodes = Vector{Int64}(1:2)
    @test global_value_min(forcesNP1, 1, nodes) == 3
    @test global_value_min(forcesNP1, 2, nodes) == 0
    @test global_value_min(forcesNP1, 3, nodes) == 0
    @test calculate_nodelist(test_Data_manager, "Forces", "Minimum", 1, nodes) == 3
    @test calculate_nodelist(test_Data_manager, "Forces", "Minimum", 2, nodes) == 0
    @test calculate_nodelist(test_Data_manager, "Forces", "Minimum", 3, nodes) == 0
    disp = test_Data_manager.create_constant_node_field("Disp", Float64, 1)
    disp[1:4] = 1:4
    nodes = Vector{Int64}(1:4)
    @test global_value_min(disp, 1, nodes) == 1
    @test calculate_nodelist(test_Data_manager, "Disp", "Minimum", 1, nodes) == 1
    nodes = Vector{Int64}(2:3)
    @test global_value_min(disp, 1, nodes) == 2
    @test calculate_nodelist(test_Data_manager, "Disp", "Minimum", 1, nodes) == 2
    @test isnothing(calculate_nodelist(test_Data_manager, "Disp", "", 1, nodes))
    @test isnothing(calculate_nodelist(test_Data_manager, "not there", "Minimum", 1, nodes))
end
@testset "ut_global_value_avg" begin
    test_Data_manager = Data_manager
    test_Data_manager.set_num_controller(4)
    nodes = Vector{Int64}(1:4)
    (forcesN, forcesNP1) = test_Data_manager.create_node_field("Forces", Float64, 3)
    forcesNP1[1, 1:3] .= 1:3
    forcesNP1[3, 1:3] .= 5.2
    forcesNP1[4, 1] = -20
    @test global_value_avg(forcesNP1, 1, nodes) == -13.8 / 4
    @test global_value_avg(forcesNP1, 2, nodes) == 7.2 / 4
    @test global_value_avg(forcesNP1, 3, nodes) == 8.2 / 4
    @test calculate_nodelist(test_Data_manager, "Forces", 1, "Average", nodes) == -13.8 / 4
    @test calculate_nodelist(test_Data_manager, "Forces", 1, "Average", nodes) == 7.2 / 4
    @test calculate_nodelist(test_Data_manager, "Forces", 1, "Average", nodes) == 8.2 / 4
    nodes = Vector{Int64}(1:2)
    @test global_value_avg(forcesNP1, 1, nodes) == 1 / 2
    @test global_value_avg(forcesNP1, 2, nodes) == 2 / 2
    @test global_value_avg(forcesNP1, 3, nodes) == 3 / 2
    @test calculate_nodelist(test_Data_manager, "Forces", "Average", 1, nodes) == 1 / 2
    @test calculate_nodelist(test_Data_manager, "Forces", "Average", 2, nodes) == 2 / 2
    @test calculate_nodelist(test_Data_manager, "Forces", "Average", 3, nodes) == 3 / 2
    disp = test_Data_manager.create_constant_node_field("Disp", Float64, 1)
    disp[1:4] = 1:4
    nodes = Vector{Int64}(1:4)
    @test global_value_avg(disp, 1, nodes) == 10 / 4
    @test calculate_nodelist(test_Data_manager, "Disp", "Average", 1, nodes) == 10 / 4
    nodes = Vector{Int64}(2:3)
    @test global_value_avg(disp, 1, nodes) == 5 / 2
    @test calculate_nodelist(test_Data_manager, "Disp", "Average", 1, nodes) == 5 / 2
    @test isnothing(calculate_nodelist(test_Data_manager, "Disp", "", 1, nodes))
    @test isnothing(calculate_nodelist(test_Data_manager, "not there", "Average", 1, nodes))
end

@testset "ut_calculate_block" begin
    test_Data_manager = Data_manager
    test_Data_manager.create_constant_node_field("Block_Id", Int64, 1)
    @test isnothing(calculate_block(test_Data_manager, "no field", 1, "sum", 1))
    @test isnothing(calculate_block(test_Data_manager, "Disp", 1, "no option", 1))
end
