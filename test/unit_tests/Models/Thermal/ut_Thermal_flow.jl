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
    flow(raw) = typed_model(:thermal,
                            merge(Dict{String,Any}("Thermal Model" => "Thermal Flow"), raw);
                            name_key = "Thermal Model")
    DM = PeriLab.Data_Manager
    DM.initialize_data()
    DM.set_num_controller(3)
    DM.set_dof(2)
    nn = DM.create_constant_node_scalar_field("Number of Neighbors", Int64)
    nn .= 2
    nlist = DM.create_constant_bond_scalar_state("Neighborhoodlist", Int64)
    nlist[1] = [2, 3]
    nlist[2] = [1, 3]
    nlist[3] = [1, 2]
    bond_length = DM.create_constant_bond_scalar_state("Bond Length", Float64)
    bond_length[1] = [1.0, 2.0]
    bond_length[2] = [1.0, 1.0]
    bond_length[3] = [2.0, 1.0]
    DM.create_constant_node_scalar_field("Horizon", Float64; default_value = 2.5)
    volume = DM.create_constant_node_scalar_field("Volume", Float64)
    volume .= [1.0, 2.0, 3.0]
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
end

@testset "ut_bond_based_heat_flow" begin
    # the precomputed bond weights give the bond-based heat flow of OterkusS2014b
    TF = PeriLab.Solver_Manager.Model_Factory.Thermal.Thermal_Flow
    DM = PeriLab.Data_Manager
    th = typed_model(:thermal,
                     Dict("Thermal Model" => "Thermal Flow",
                          "Thermal Conductivity" => 100);
                     name_key = "Thermal Model")
    TF.init_model(Vector{Int64}(1:3), th.model, th, 1)
    weight = DM.get_field("Thermal Bond Weight")
    kernel = 6 / (pi * 2.5^3)
    @test weight[1] ≈ [100 * kernel * 2 / 1, 100 * kernel * 3 / 2]
    temperature = [1.0, 2.0, 4.0]
    damage = DM.create_constant_bond_scalar_state("Damage Test", Float64;
                                                  default_value = 1)
    damage[1][2] = 0
    heat_flow = zeros(3)
    TF.compute_heat_flow_state_bond_based!(1:3, DM.get_nlist(), weight, damage,
                                           [true, true, false], temperature, heat_flow)
    @test heat_flow[1] ≈ -100 * kernel * 2 * (2 - 1)
    @test heat_flow[2] ≈ -100 * kernel * (1 * (1 - 2) + 3 * (4 - 2))
    @test heat_flow[3] == 0
end
