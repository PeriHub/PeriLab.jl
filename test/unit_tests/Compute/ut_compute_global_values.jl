# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
include("../../../src/Support/data_manager.jl")
include("../../../src/Compute/compute_global_values.jl")

@testset "ut_global_value_sum" begin
    test_Data_manager = Data_manager
    test_Data_manager.set_nmasters(4)
    nodes = Vector{Int64}(1:4)
    (forcesN, forcesNP1) = test_Data_manager.create_node_field("Forces", Float64, 3)
    forcesNP1[1, 1:3] .= 1:3
    forcesNP1[3, 1:3] .= 5.2
    forcesNP1[4, 1] = -20
    testValues = global_value_sum(forcesNP1, nodes)

    @test length(testValues) == 3
    @test testValues[1] == -13.8
    @test testValues[2] == 7.2
    @test testValues[3] == 8.2
    testValues = calculate_nodelist(test_Data_manager, "Forces", "Sum", nodes)
    @test length(testValues) == 3
    @test testValues[1] == -13.8
    @test testValues[2] == 7.2
    @test testValues[3] == 8.2

    nodes = Vector{Int64}(1:2)
    testValues = global_value_sum(forcesNP1, nodes)
    @test testValues[1] == 1
    @test testValues[2] == 2
    @test testValues[3] == 3
    testValues = calculate_nodelist(test_Data_manager, "Forces", "Sum", nodes)
    @test testValues[1] == 1
    @test testValues[2] == 2
    @test testValues[3] == 3
    disp = test_Data_manager.create_constant_node_field("Disp", Float64, 1)
    disp[1:4] = 1:4
    nodes = Vector{Int64}(1:4)
    testValues = global_value_sum(disp, nodes)
    @test length(testValues) == 1
    @test testValues[1] == 10
    testValues = calculate_nodelist(test_Data_manager, "Disp", "Sum", nodes)
    @test length(testValues) == 1
    @test testValues[1] == 10

    nodes = Vector{Int64}(2:3)
    testValues = global_value_sum(disp, nodes)
    @test length(testValues) == 1
    @test testValues[1] == 5

    testValues = calculate_nodelist(test_Data_manager, "Disp", "Sum", nodes)
    @test length(testValues) == 1
    @test testValues[1] == 5
    testValues = calculate_nodelist(test_Data_manager, "Disp", "", nodes)
    @test isnothing(testValues)
    testValues = calculate_nodelist(test_Data_manager, "not there", "Sum", nodes)
    @test isnothing(testValues)
end

@testset "ut_global_value_max" begin
    test_Data_manager = Data_manager
    test_Data_manager.set_nmasters(4)
    nodes = Vector{Int64}(1:4)
    (forcesN, forcesNP1) = test_Data_manager.create_node_field("Forces", Float64, 3)

    forcesNP1[1, 1:3] .= 1:3
    forcesNP1[3, 1:3] .= 5.2
    forcesNP1[4, 1] = -20
    testValues = global_value_max(forcesNP1, nodes)
    @test length(testValues) == 3
    @test testValues[1] == 5.2
    @test testValues[2] == 5.2
    @test testValues[3] == 5.2
    testValues = calculate_nodelist(test_Data_manager, "Forces", "Maximum", nodes)
    println()
    @test length(testValues) == 3
    @test testValues[1] == 5.2
    @test testValues[2] == 5.2
    @test testValues[3] == 5.2

    nodes = Vector{Int64}(1:2)
    testValues = global_value_max(forcesNP1, nodes)
    @test length(testValues) == 3
    @test testValues[1] == 1
    @test testValues[2] == 2
    @test testValues[3] == 3
    testValues = calculate_nodelist(test_Data_manager, "Forces", "Maximum", nodes)
    @test length(testValues) == 3
    @test testValues[1] == 1
    @test testValues[2] == 2
    @test testValues[3] == 3
    disp = test_Data_manager.create_constant_node_field("Disp", Float64, 1)
    disp[1:4] = 1:4
    nodes = Vector{Int64}(1:4)
    testValues = global_value_max(disp, nodes)
    @test length(testValues) == 1
    @test testValues[1] == 4
    testValues = calculate_nodelist(test_Data_manager, "Disp", "Maximum", nodes)
    testValues = global_value_max(disp, nodes)
    @test length(testValues) == 1
    @test testValues[1] == 4
    nodes = Vector{Int64}(2:3)
    testValues = global_value_max(disp, nodes)
    @test length(testValues) == 1
    @test testValues[1] == 3
    testValues = calculate_nodelist(test_Data_manager, "Disp", "Maximum", nodes)
    @test length(testValues) == 1
    @test testValues[1] == 3
    testValues = calculate_nodelist(test_Data_manager, "Disp", "", nodes)
    @test isnothing(testValues)
    testValues = calculate_nodelist(test_Data_manager, "not there", "Maximum", nodes)
    @test isnothing(testValues)
end

@testset "ut_global_value_min" begin
    test_Data_manager = Data_manager
    test_Data_manager.set_nmasters(4)
    nodes = Vector{Int64}(1:4)
    (forcesN, forcesNP1) = test_Data_manager.create_node_field("Forces", Float64, 3)
    forcesNP1[1, 1:3] .= 1:3
    forcesNP1[3, 1:3] .= 5.2
    forcesNP1[4, 1] = -20
    testValues = global_value_min(forcesNP1, nodes)
    @test length(testValues) == 3
    @test testValues[1] == -20
    @test testValues[2] == 0
    @test testValues[3] == 0
    testValues = calculate_nodelist(test_Data_manager, "Forces", "Minimum", nodes)
    @test length(testValues) == 3
    @test testValues[1] == -20
    @test testValues[2] == 0
    @test testValues[3] == 0
    nodes = Vector{Int64}(1:2)
    testValues = global_value_min(forcesNP1, nodes)
    @test length(testValues) == 3
    @test testValues[1] == 0
    @test testValues[2] == 0
    @test testValues[3] == 0
    testValues = calculate_nodelist(test_Data_manager, "Forces", "Minimum", nodes)
    @test length(testValues) == 3
    @test testValues[1] == 0
    @test testValues[2] == 0
    @test testValues[3] == 0
    disp = test_Data_manager.create_constant_node_field("Disp", Float64, 1)
    disp[1:4] = 1:4
    nodes = Vector{Int64}(1:4)
    testValues = global_value_min(disp, nodes)
    @test length(testValues) == 1
    @test testValues[1] == 1
    testValues = calculate_nodelist(test_Data_manager, "Disp", "Minimum", nodes)
    @test length(testValues) == 1
    @test testValues[1] == 1
    nodes = Vector{Int64}(2:3)
    testValues = global_value_min(disp, nodes)
    @test length(testValues) == 1
    @test testValues[1] == 2
    testValues = calculate_nodelist(test_Data_manager, "Disp", "Minimum", nodes)
    @test length(testValues) == 1
    @test testValues[1] == 2
    testValues = calculate_nodelist(test_Data_manager, "Disp", "", nodes)
    @test isnothing(testValues)
    testValues = calculate_nodelist(test_Data_manager, "not there", "Minimum", nodes)
    @test isnothing(testValues)
end
@testset "ut_global_value_avg" begin
    test_Data_manager = Data_manager
    test_Data_manager.set_nmasters(4)
    nodes = Vector{Int64}(1:4)
    (forcesN, forcesNP1) = test_Data_manager.create_node_field("Forces", Float64, 3)
    forcesNP1[1, 1:3] .= 1:3
    forcesNP1[3, 1:3] .= 5.2
    forcesNP1[4, 1] = -20
    testValues = global_value_avg(forcesNP1, nodes)
    @test length(testValues) == 3
    @test testValues[1] == -13.8 / 4
    @test testValues[2] == 7.2 / 4
    @test testValues[3] == 8.2 / 4
    testValues = calculate_nodelist(test_Data_manager, "Forces", "Average", nodes)
    @test length(testValues) == 3
    @test testValues[1] == -13.8 / 4
    @test testValues[2] == 7.2 / 4
    @test testValues[3] == 8.2 / 4
    nodes = Vector{Int64}(1:2)
    testValues = global_value_avg(forcesNP1, nodes)
    @test length(testValues) == 3
    @test testValues[1] == 1 / 2
    @test testValues[2] == 2 / 2
    @test testValues[3] == 3 / 2
    testValues = calculate_nodelist(test_Data_manager, "Forces", "Average", nodes)
    @test length(testValues) == 3
    @test testValues[1] == 1 / 2
    @test testValues[2] == 2 / 2
    @test testValues[3] == 3 / 2
    disp = test_Data_manager.create_constant_node_field("Disp", Float64, 1)
    disp[1:4] = 1:4
    nodes = Vector{Int64}(1:4)
    testValues = global_value_avg(disp, nodes)
    @test length(testValues) == 1
    @test testValues[1] == 10 / 4
    testValues = calculate_nodelist(test_Data_manager, "Disp", "Average", nodes)
    @test length(testValues) == 1
    @test testValues[1] == 10 / 4
    nodes = Vector{Int64}(2:3)
    testValues = global_value_avg(disp, nodes)
    @test length(testValues) == 1
    @test testValues[1] == 5 / 2
    testValues = calculate_nodelist(test_Data_manager, "Disp", "Average", nodes)
    @test length(testValues) == 1
    @test testValues[1] == 5 / 2
    testValues = calculate_nodelist(test_Data_manager, "Disp", "", nodes)
    @test isnothing(testValues)
    testValues = calculate_nodelist(test_Data_manager, "not there", "Average", nodes)
    @test isnothing(testValues)
end


