# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test

include("../../../src/Compute/compute_global_values.jl")

include("../../../src/PeriLab.jl")
using .PeriLab
@testset "ut_global_value_sum" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.initialize_data()
    test_data_manager.set_num_controller(4)
    test_data_manager.set_glob_to_loc(Dict(1 => 1, 2 => 2, 3 => 3, 4 => 4))
    nodes = Vector{Int64}(1:4)
    (forcesN, forcesNP1) = test_data_manager.create_node_field("Forces", Float64, 3)
    forcesNP1[1, 1:3] .= 1:3
    forcesNP1[3, 1:3] .= 5.2
    forcesNP1[4, 1] = -20
    @test global_value_sum(forcesNP1, 1, nodes) == -13.8
    @test global_value_sum(forcesNP1, 2, nodes) == 7.2
    @test global_value_sum(forcesNP1, 3, nodes) == 8.2
    @test calculate_nodelist(test_data_manager, "Forces", 1, "Sum", nodes) == (-13.8, 4)
    @test calculate_nodelist(test_data_manager, "Forces", 2, "Sum", nodes) == (7.2, 4)
    @test calculate_nodelist(test_data_manager, "Forces", 3, "Sum", nodes) == (8.2, 4)

    nodes = Vector{Int64}(1:2)
    @test global_value_sum(forcesNP1, 1, nodes) == 1
    @test global_value_sum(forcesNP1, 2, nodes) == 2
    @test global_value_sum(forcesNP1, 3, nodes) == 3
    @test calculate_nodelist(test_data_manager, "Forces", 1, "Sum", nodes) == (1, 2)
    @test calculate_nodelist(test_data_manager, "Forces", 2, "Sum", nodes) == (2, 2)
    @test calculate_nodelist(test_data_manager, "Forces", 3, "Sum", nodes) == (3, 2)
    disp = test_data_manager.create_constant_node_field("Disp", Float64, 1)
    disp[1:4] = 1:4
    nodes = Vector{Int64}(1:4)
    @test global_value_sum(disp, 1, nodes) == 10
    @test calculate_nodelist(test_data_manager, "Disp", 1, "Sum", nodes) == (10, 4)
    nodes = Vector{Int64}(2:3)
    @test global_value_sum(disp, 1, nodes) == 5
    @test calculate_nodelist(test_data_manager, "Disp", 1, "Sum", nodes) == (5, 2)
    @test isnothing(calculate_nodelist(test_data_manager, "Disp", 1, "", nodes))
    @test isnothing(calculate_nodelist(test_data_manager, "not there", 1, "Sum", nodes))

    test_data_manager.set_glob_to_loc(Dict(1 => 1, 2 => 2))
    nodes = Vector{Int64}(3:4)

    @test calculate_nodelist(test_data_manager, "Disp", 1, "Sum", nodes) == (0, 0)
    @test calculate_nodelist(test_data_manager, "Disp", 1, "Average", nodes) == (0, 0)
    @test calculate_nodelist(test_data_manager, "Disp", 1, "Minimum", nodes) == (Inf, 0)
    @test calculate_nodelist(test_data_manager, "Disp", 1, "Maximum", nodes) == (-Inf, 0)

    matrix = test_data_manager.create_constant_node_field("Matrix", Float64, "Matrix", 3)
    matrix[:, 1, 2] .= 4
    nodes = Vector{Int64}(1:2)
    @test calculate_nodelist(test_data_manager, "Matrix", [1, 2], "Sum", nodes) == (8, 2)

end

@testset "ut_global_value_max" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.set_num_controller(4)
    test_data_manager.set_glob_to_loc(Dict(1 => 1, 2 => 2, 3 => 3, 4 => 4))
    nodes = Vector{Int64}(1:4)
    (forcesN, forcesNP1) = test_data_manager.create_node_field("Forces", Float64, 3)

    forcesNP1[1, 1:3] .= 1:3
    forcesNP1[3, 1:3] .= 5.2
    forcesNP1[4, 1] = -20

    @test global_value_max(forcesNP1, 1, nodes) == 5.2
    @test global_value_max(forcesNP1, 2, nodes) == 5.2
    @test global_value_max(forcesNP1, 3, nodes) == 5.2
    @test calculate_nodelist(test_data_manager, "Forces", 1, "Maximum", nodes) == (5.2, 4)
    @test calculate_nodelist(test_data_manager, "Forces", 1, "Maximum", nodes) == (5.2, 4)
    @test calculate_nodelist(test_data_manager, "Forces", 1, "Maximum", nodes) == (5.2, 4)
    nodes = Vector{Int64}(1:2)
    @test global_value_max(forcesNP1, 1, nodes) == 1
    @test global_value_max(forcesNP1, 2, nodes) == 2
    @test global_value_max(forcesNP1, 3, nodes) == 3
    @test calculate_nodelist(test_data_manager, "Forces", 1, "Maximum", nodes) == (1, 2)
    @test calculate_nodelist(test_data_manager, "Forces", 2, "Maximum", nodes) == (2, 2)
    @test calculate_nodelist(test_data_manager, "Forces", 3, "Maximum", nodes) == (3, 2)
    disp = test_data_manager.create_constant_node_field("Disp", Float64, 1)
    disp[1:4] = 1:4
    nodes = Vector{Int64}(1:4)
    @test global_value_max(disp, 1, nodes) == 4
    @test calculate_nodelist(test_data_manager, "Disp", 1, "Maximum", nodes) == (4, 4)
    nodes = Vector{Int64}(2:3)
    @test global_value_max(disp, 1, nodes) == 3
    @test calculate_nodelist(test_data_manager, "Disp", 1, "Maximum", nodes) == (3, 2)
    @test isnothing(calculate_nodelist(test_data_manager, "Disp", 1, "", nodes))
    @test isnothing(calculate_nodelist(test_data_manager, "not there", 1, "Maximum", nodes))
end

@testset "ut_global_value_min" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.set_num_controller(4)
    test_data_manager.set_glob_to_loc(Dict(1 => 1, 2 => 2, 3 => 3, 4 => 4))
    nodes = Vector{Int64}(1:4)
    (forcesN, forcesNP1) = test_data_manager.create_node_field("Forces", Float64, 3)
    forcesNP1[1, 1:3] .= 1:3
    forcesNP1[3, 1:3] .= 5.2
    forcesNP1[4, 1] = -20
    @test global_value_min(forcesNP1, 1, nodes) == -20
    @test global_value_min(forcesNP1, 2, nodes) == 0
    @test global_value_min(forcesNP1, 3, nodes) == 0
    @test calculate_nodelist(test_data_manager, "Forces", 1, "Minimum", nodes) == (-20, 4)
    @test calculate_nodelist(test_data_manager, "Forces", 2, "Minimum", nodes) == (0, 4)
    @test calculate_nodelist(test_data_manager, "Forces", 3, "Minimum", nodes) == (0, 4)
    nodes = Vector{Int64}(1:2)
    @test global_value_min(forcesNP1, 1, nodes) == 0
    @test global_value_min(forcesNP1, 2, nodes) == 0
    @test global_value_min(forcesNP1, 3, nodes) == 0
    @test calculate_nodelist(test_data_manager, "Forces", 1, "Minimum", nodes) == (0, 2)
    @test calculate_nodelist(test_data_manager, "Forces", 2, "Minimum", nodes) == (0, 2)
    @test calculate_nodelist(test_data_manager, "Forces", 3, "Minimum", nodes) == (0, 2)
    disp = test_data_manager.create_constant_node_field("Disp", Float64, 1)
    disp[1:4] = 1:4
    nodes = Vector{Int64}(1:4)
    @test global_value_min(disp, 1, nodes) == 1
    @test calculate_nodelist(test_data_manager, "Disp", 1, "Minimum", nodes) == (1, 4)
    nodes = Vector{Int64}(2:3)
    @test global_value_min(disp, 1, nodes) == 2
    @test calculate_nodelist(test_data_manager, "Disp", 1, "Minimum", nodes) == (2, 2)
    @test isnothing(calculate_nodelist(test_data_manager, "Disp", 1, "", nodes))
    @test isnothing(calculate_nodelist(test_data_manager, "not there", 1, "Minimum", nodes))
end
@testset "ut_global_value_avg" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.set_num_controller(4)
    test_data_manager.set_glob_to_loc(Dict(1 => 1, 2 => 2, 3 => 3, 4 => 4))
    nodes = Vector{Int64}(1:4)
    (forcesN, forcesNP1) = test_data_manager.create_node_field("Forces", Float64, 3)
    forcesNP1[1, 1:3] .= 1:3
    forcesNP1[3, 1:3] .= 5.2
    forcesNP1[4, 1] = -20
    @test global_value_avg(forcesNP1, 1, nodes) == -13.8 / 4
    @test global_value_avg(forcesNP1, 2, nodes) == 7.2 / 4
    @test global_value_avg(forcesNP1, 3, nodes) == 8.2 / 4
    @test calculate_nodelist(test_data_manager, "Forces", 1, "Average", nodes) == (-13.8 / 4, 4)
    @test calculate_nodelist(test_data_manager, "Forces", 2, "Average", nodes) == (7.2 / 4, 4)
    @test calculate_nodelist(test_data_manager, "Forces", 3, "Average", nodes) == (8.2 / 4, 4)
    nodes = Vector{Int64}(1:2)
    @test global_value_avg(forcesNP1, 1, nodes) == 1 / 2
    @test global_value_avg(forcesNP1, 2, nodes) == 2 / 2
    @test global_value_avg(forcesNP1, 3, nodes) == 3 / 2
    @test calculate_nodelist(test_data_manager, "Forces", 1, "Average", nodes) == (1 / 2, 2)
    @test calculate_nodelist(test_data_manager, "Forces", 2, "Average", nodes) == (2 / 2, 2)
    @test calculate_nodelist(test_data_manager, "Forces", 3, "Average", nodes) == (3 / 2, 2)
    disp = test_data_manager.create_constant_node_field("Disp", Float64, 1)
    disp[1:4] = 1:4
    nodes = Vector{Int64}(1:4)
    @test global_value_avg(disp, 1, nodes) == 10 / 4
    @test calculate_nodelist(test_data_manager, "Disp", 1, "Average", nodes) == (10 / 4, 4)
    nodes = Vector{Int64}(2:3)
    @test global_value_avg(disp, 1, nodes) == 5 / 2
    @test calculate_nodelist(test_data_manager, "Disp", 1, "Average", nodes) == (5 / 2, 2)
    @test isnothing(calculate_nodelist(test_data_manager, "Disp", 1, "", nodes))
    @test isnothing(calculate_nodelist(test_data_manager, "not there", 1, "Average", nodes))
end

@testset "ut_calculate_block" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.create_constant_node_field("Block_Id", Int64, 1, 1)
    test_data_manager.set_glob_to_loc(Dict(1 => 1, 2 => 2, 3 => 3, 4 => 4))
    nodes = Vector{Int64}(1:4)
    (forcesN, forcesNP1) = test_data_manager.create_node_field("Forces", Float64, 3)
    forcesNP1[1, 1:3] .= 1:3
    forcesNP1[3, 1:3] .= 5.2
    forcesNP1[4, 1] = -20
    @test isnothing(calculate_block(test_data_manager, "no field", 1, "sum", 1))
    @test calculate_block(test_data_manager, "Forces", 1, "Sum", 1) == (-13.8, 4)
end
