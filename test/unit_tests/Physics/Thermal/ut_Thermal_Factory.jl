# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
using Test

include("../../../../src/Physics/Thermal/Thermal_Factory.jl")
using .Thermal
# include("../../../../src/Support/data_manager.jl")


@testset "init_thermal_model_fields" begin
    test_Data_manager = PeriLab.Data_manager
    test_Data_manager.clear_data_manager()
    test_Data_manager.set_dof(3)
    test_Data_manager.set_num_controller(4)
    nn = test_Data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1] = 2
    nn[2] = 3
    nn[3] = 1
    nn[4] = 2

    Thermal.init_thermal_model_fields(test_Data_manager)
    fieldkeys = test_Data_manager.get_all_field_keys()
    @test "TemperatureN" in fieldkeys
    @test "TemperatureNP1" in fieldkeys
    @test "Heat FlowN" in fieldkeys
    @test "Heat FlowNP1" in fieldkeys
    @test "Bond Heat Flow" in fieldkeys

end
@testset "init_thermal_model" begin
    test_Data_manager = PeriLab.Data_manager
    test_Data_manager.clear_data_manager()
    test_Data_manager.set_block_list([1, 2])
    test_Data_manager.init_property()
    test_Data_manager.set_properties(1, "Thermal Model", Dict("Thermal Model" => "Heat Transfer"))
    Thermal.init_thermal_model(test_Data_manager, [1], 1)
    test_Data_manager.set_properties(2, "Thermal Model", Dict("Thermal Model" => "Missing"))
    Thermal.init_thermal_model(test_Data_manager, [2], 1)
end

@testset "ut_distribute_heat_flows" begin

    nnodes = 2
    nodes = Vector{Int64}(1:nnodes)
    dof = 2
    test_Data_manager = PeriLab.Data_manager
    test_Data_manager.set_num_controller(2)
    nn = test_Data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    volume = test_Data_manager.create_constant_node_field("Volume", Float64, 1)

    nn[1] = 1
    nn[2] = 1
    volume[1:2] .= 1

    nlist = test_Data_manager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
    nlist[1][1] = 2
    nlist[2][1] = 1

    (heat_flowN, heat_flowNP1) = test_Data_manager.create_node_field("Heat Flow", Float64, 1)
    bond_heat_flow = test_Data_manager.create_constant_bond_field("Bond Heat Flow", Float64, 1)
    bond_heat_flow[1][1] = 1
    bond_heat_flow[2][1] = 1

    test_Data_manager = Thermal.distribute_heat_flows(test_Data_manager, nodes)

    @test heat_flowNP1[1] == -1.0
    @test heat_flowNP1[2] == -1.0
    heat_flowNP1[1:2] .= 0
    bond_heat_flow[1][1] = 1
    bond_heat_flow[2][1] = -1

    test_Data_manager = Thermal.distribute_heat_flows(test_Data_manager, nodes)

    @test heat_flowNP1[1] == -1
    @test heat_flowNP1[2] == 1

    volume[2] = 0

    test_Data_manager = Thermal.distribute_heat_flows(test_Data_manager, nodes)

    @test heat_flowNP1[1] == -1
    @test heat_flowNP1[2] == 2
end
