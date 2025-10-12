# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test

#include("../../../../src/PeriLab.jl")
#import .PeriLab

@test PeriLab.Solver_Manager.Model_Factory.Thermal.Heat_Transfer.thermal_model_name() == "Heat Transfer"

@testset "ut_calculate_specific_volume" begin
    nnodes = 10
    dof = 2
    nodes = Vector{Int64}(1:nnodes)
    test_data_manager = PeriLab.Data_Manager
    test_data_manager.initialize_data()
    test_data_manager.set_num_controller(nnodes)
    test_data_manager.set_dof(dof)
    nn = test_data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1] = 2
    nn[2] = 3
    nn[3] = 2
    nn[4] = 3
    nn[5] = 4
    nn[6] = 3
    nn[7] = 2
    nn[8] = 3
    nn[9] = 2
    nn[10] = 1
    nodes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    nlist = test_data_manager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
    nlist[1] = [4, 2]
    nlist[2] = [1, 3, 5]
    nlist[3] = [2, 6]
    nlist[4] = [1, 5, 7]
    nlist[5] = [2, 4, 6, 8]
    nlist[6] = [3, 5, 9]
    nlist[7] = [4, 8]
    nlist[8] = [5, 7, 9]
    nlist[9] = [6, 8]
    nlist[10] = [9]
    bond_norm = test_data_manager.create_constant_bond_field("Bond Norm", Float64, dof, 1)
    bond_norm[1] = [[1, 0], [-1, 0]]
    bond_norm[2] = [[1, 0], [-1, 0], [0, 1]]
    bond_norm[3] = [[1, 0], [0.5, 0.5]]
    bond_norm[4] = [[1, 0], [-1, 0], [0, -1]]
    bond_norm[5] = [[1, 0], [-1, 0], [0, 1], [0, -1]]
    bond_norm[6] = [[1, 0], [-1, 0], [0, 1]]
    bond_norm[7] = [[1, 0], [0, -1]]
    bond_norm[8] = [[1, 0], [-1, 0], [0, 1]]
    bond_norm[9] = [[1, 0], [0, 1]]
    bond_norm[10] = [[1, 0]]
    volume = test_data_manager.create_constant_node_field("Volume", Float64, 1, 0.25)
    specific_volume = test_data_manager.create_constant_node_field("specific_volume",
                                                                   Int64, 1)
    active = test_data_manager.create_constant_node_field("Active", Bool, 1, true)
    specific_volume_check = test_data_manager.create_constant_node_field("Specific Volume Check",
                                                                         Bool, 1, true)
    rotation_tensor = nothing
    PeriLab.Solver_Manager.Model_Factory.Thermal.Heat_Transfer.calculate_specific_volume!(specific_volume,
                                             nodes,
                                             nlist,
                                             active,
                                             bond_norm,
                                             rotation_tensor,
                                             specific_volume_check,
                                             dof)
    @test specific_volume == [
        2,
        1,
        3,
        1,
        0,
        1,
        2,
        1,
        2,
        3
    ]
end

@testset "ut_compute_thermal_model" begin
    test_data_manager = PeriLab.Data_Manager

    test_data_manager.create_node_field("Heat Flow", Float64, 1)
    test_data_manager.create_node_field("Temperature", Float64, 1)
    test_data_manager.create_constant_node_field("Specific Volume", Float64, 1)
    test_data_manager.create_constant_node_field("Surface_Nodes", Bool, 1)
    test_data_manager.create_bond_field("Bond Damage", Float64, 1)
    @test PeriLab.Solver_Manager.Model_Factory.Thermal.Heat_Transfer.compute_model(test_data_manager,
                                      Vector{Int64}(1:10),
                                      Dict("Heat Transfer Coefficient" => 1,
                                           "Environmental Temperature" => 1.2,
                                           "Allow Surface Change" => false),
                                      1,
                                      1.0,
                                      1.0) == test_data_manager

    dof = 3
    test_data_manager.set_dof(dof)

    @test PeriLab.Solver_Manager.Model_Factory.Thermal.Heat_Transfer.compute_model(test_data_manager,
                                      Vector{Int64}(1:10),
                                      Dict("Heat Transfer Coefficient" => 1,
                                           "Environmental Temperature" => 1.2,
                                           "Allow Surface Change" => false),
                                      1,
                                      1.0,
                                      1.0) == test_data_manager
end
