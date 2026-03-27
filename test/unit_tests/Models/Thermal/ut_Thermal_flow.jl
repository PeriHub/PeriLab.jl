# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

#using Test

#include("../../../../src/PeriLab.jl")
#using .PeriLab

@test PeriLab.Solver_Manager.Model_Factory.Thermal.Thermal_Flow.thermal_model_name() ==
      "Thermal Flow"

@testset "ut_init_model" begin

    #@test_logs (:warn, "No model type has beed defined; ''Type'': ''Bond based'' or Type: ''Correspondence'; \n ''Bond based'' is set as default.'") PeriLab.Solver_Manager.Model_Factory.Thermal.Thermal_Flow.init_model(Vector{Int64}(1:3), Dict{String, Any}("a" => 1, "Thermal Conductivity" => 100))
    #@test_logs (:warn, "No model type has beed defined; ''Type'': ''Bond based'' or Type: ''Correspondence'; \n ''Bond based'' is set as default.'") PeriLab.Solver_Manager.Model_Factory.Thermal.Thermal_Flow.init_model(Vector{Int64}(1:3), Dict("Type" => "a", "Thermal Conductivity" => 100))
    PeriLab.Solver_Manager.Model_Factory.Thermal.Thermal_Flow.init_model(Vector{Int64}(1:3),
                                                                         Dict("Type" => "Bond based",
                                                                              "Thermal Conductivity" => 100))
    PeriLab.Solver_Manager.Model_Factory.Thermal.Thermal_Flow.init_model(Vector{Int64}(1:3),
                                                                         Dict("Type" => "Correspondence",
                                                                              "Thermal Conductivity" => 100))
    @test_logs (:error, "Thermal Conductivity not defined.") PeriLab.Solver_Manager.Model_Factory.Thermal.Thermal_Flow.init_model(Vector{Int64}(1:3),
                                                                                                                                  Dict("Type" => "Bond based"))
    @test_logs (:error, "Thermal Conductivity not defined.") PeriLab.Solver_Manager.Model_Factory.Thermal.Thermal_Flow.init_model(Vector{Int64}(1:3),
                                                                                                                                  Dict("Type" => "Correspondence"))
    parameter = Dict("Type" => "Correspondence",
                     "Thermal Conductivity" => 100,
                     "Print Bed Temperature" => 2, "Print Bed Z Coordinate" => -1)

    coordinates = PeriLab.Data_Manager.create_constant_node_vector_field("Coordinates",
                                                                         Float64, 3)
    coordinates[1, 1] = 0
    coordinates[1, 2] = 0
    coordinates[1, 3] = 1
    coordinates[2, 1] = 1
    coordinates[2, 2] = 0
    coordinates[2, 3] = 1
    coordinates[3, 1] = 1
    coordinates[3, 2] = 1
    coordinates[3, 3] = 1

    PeriLab.Data_Manager.set_dof(3)
    PeriLab.Solver_Manager.Model_Factory.Thermal.Thermal_Flow.init_model(Vector{Int64}(1:3),
                                                                         parameter)

    @test haskey(parameter, "Print Bed Temperature")
    PeriLab.Data_Manager.set_dof(2)

    PeriLab.Solver_Manager.Model_Factory.Thermal.Thermal_Flow.init_model(Vector{Int64}(1:3),
                                                                         parameter)
    @test !haskey(parameter, "Print Bed Temperature")
end
