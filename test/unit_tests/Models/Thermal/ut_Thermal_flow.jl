# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

#using Test

#include("../../../../src/PeriLab.jl")
#using .PeriLab

@test parentmodule(PeriLab.ParameterSpec.lookup_model(:thermal, "Thermal Flow")) ===
      PeriLab.Solver_Manager.Model_Factory.Thermal.Thermal_Flow

@testset "ut_init_model" begin

    TF = PeriLab.Solver_Manager.Model_Factory.Thermal.Thermal_Flow
    flow(raw) = typed_model(:thermal, merge(Dict{String,Any}("Thermal Model" => "Thermal Flow"), raw);
                            name_key = "Thermal Model")
    for type in ("Bond based", "Correspondence")
        th = flow(Dict("Type" => type, "Thermal Conductivity" => 100))
        TF.init_model(Vector{Int64}(1:3), th.model, th, 1)
        th = flow(Dict("Type" => type))
        @test_logs (:error,
                    "Thermal Conductivity not defined.") @test_throws PeriLab.PeriLabError TF.init_model(Vector{Int64}(1:3),
                                                                                                         th.model,
                                                                                                         th,
                                                                                                         1)
    end
    coordinates = PeriLab.Data_Manager.create_constant_node_vector_field("Coordinates",
                                                                         Float64, 3)
    coordinates[1:3, :] .= [0 0 1; 1 0 1; 1 1 1]
    th = flow(Dict("Type" => "Correspondence", "Thermal Conductivity" => 100,
                   "Print Bed Temperature" => 2, "Thermal Conductivity Print Bed" => 1,
                   "Print Bed Z Coordinate" => -1))
    PeriLab.Data_Manager.set_dof(3)
    TF.init_model(Vector{Int64}(1:3), th.model, th, 1)
    @test TF.print_bed_active(th.model, 3)
    PeriLab.Data_Manager.set_dof(2)
    TF.init_model(Vector{Int64}(1:3), th.model, th, 1)    # warns: warnings are off in the suite
    @test !TF.print_bed_active(th.model, 2)
end
