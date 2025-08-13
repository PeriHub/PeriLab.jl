# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test

include("../../../../src/Models/Thermal/thermal_flow.jl")
using .Thermal_Flow
#include("../../../../src/PeriLab.jl")
#using .PeriLab

@test Thermal_Flow.thermal_model_name() == "Thermal Flow"

@testset "ut_init_model" begin
    test_data_manager = PeriLab.Data_manager

    @test isnothing(Thermal_Flow.init_model(test_data_manager,
                                            Vector{Int64}(1:3),
                                            Dict("a" => 1, "Thermal Conductivity" => 100)))
    @test isnothing(Thermal_Flow.init_model(test_data_manager,
                                            Vector{Int64}(1:3),
                                            Dict("Type" => "a",
                                                 "Thermal Conductivity" => 100)))
    @test Thermal_Flow.init_model(test_data_manager,
                                  Vector{Int64}(1:3),
                                  Dict("Type" => "Bond based",
                                       "Thermal Conductivity" => 100)) == test_data_manager
    @test Thermal_Flow.init_model(test_data_manager,
                                  Vector{Int64}(1:3),
                                  Dict("Type" => "Correspondence",
                                       "Thermal Conductivity" => 100)) == test_data_manager
    @test isnothing(Thermal_Flow.init_model(test_data_manager,
                                            Vector{Int64}(1:3),
                                            Dict("Type" => "Bond based")))
    @test isnothing(Thermal_Flow.init_model(test_data_manager,
                                            Vector{Int64}(1:3),
                                            Dict("Type" => "Correspondence")))
    parameter = Dict("Type" => "Correspondence",
                     "Thermal Conductivity" => 100,
                     "Print Bed Temperature" => 2)

    coordinates = test_data_manager.create_constant_node_field("Coordinates", Float64, 3)
    coordinates[1, 1] = 0
    coordinates[1, 2] = 0
    coordinates[1, 3] = 1
    coordinates[2, 1] = 1
    coordinates[2, 2] = 0
    coordinates[2, 3] = 1
    coordinates[3, 1] = 1
    coordinates[3, 2] = 1
    coordinates[3, 3] = 1

    test_data_manager.set_dof(3)
    @test Thermal_Flow.init_model(test_data_manager, Vector{Int64}(1:3), parameter) ==
          test_data_manager

    @test haskey(parameter, "Print Bed Temperature")
    test_data_manager.set_dof(2)

    @test Thermal_Flow.init_model(test_data_manager, Vector{Int64}(1:3), parameter) ==
          test_data_manager

    @test !haskey(parameter, "Print Bed Temperature")
end
