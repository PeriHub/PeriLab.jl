# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
#include("../../../../src/PeriLab.jl")
#import .PeriLab

@test PeriLab.Solver_Manager.Model_Factory.Thermal.Thermal_Expansion.thermal_model_name() ==
      "Thermal Expansion"

@testset "ut_thermal_deformation" begin
    nnodes = 2
    dof = 2
    test_data_manager = PeriLab.Data_Manager
    test_data_manager.initialize_data()
    test_data_manager.set_num_controller(2)
    nn = test_data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1] = 2
    nn[2] = 3
    temperature_N,
    temperature_NP1 = test_data_manager.create_node_field("Temperature",
                                                          Float64, 1)
    undeformed_bond = test_data_manager.create_constant_bond_field("Bond Geometry",
                                                                   Float64, dof)
    deformed_bond_N,
    deformed_bond_NP1 = test_data_manager.create_bond_field("Deformed Bond Geometry",
                                                            Float64,
                                                            dof)
    undeformed_bond_length = test_data_manager.create_constant_bond_field("Bond Length",
                                                                          Float64, 1)

    deformed_bond_length_N,
    deformed_bond_length_NP1 = test_data_manager.create_bond_field("Deformed Bond Length",
                                                                   Float64,
                                                                   1)

    undeformed_bond[1][1][1] = 0
    undeformed_bond[1][1][2] = 1
    undeformed_bond_length[1][1] = 1
    undeformed_bond[1][2][1] = 1
    undeformed_bond[1][2][2] = 1
    undeformed_bond_length[1][2] = sqrt(2)

    undeformed_bond[2][1][1] = -1
    undeformed_bond[2][1][2] = -1
    undeformed_bond_length[2][1] = sqrt(2)
    undeformed_bond[2][2][1] = 10
    undeformed_bond[2][2][2] = -10
    undeformed_bond_length[2][2] = sqrt(200)
    undeformed_bond[2][3][1] = 0.1
    undeformed_bond[2][3][2] = 0
    undeformed_bond_length[2][3] = 0.1

    deformed_bond_NP1[1][1][1] = 0
    deformed_bond_NP1[1][1][2] = 1
    deformed_bond_length_NP1[1][1] = 1
    deformed_bond_NP1[1][2][1] = 1
    deformed_bond_NP1[1][2][2] = 1
    deformed_bond_length_NP1[1][2] = sqrt(2)

    deformed_bond_NP1[2][1][1] = -1
    deformed_bond_NP1[2][1][2] = -1
    deformed_bond_length_NP1[2][1] = sqrt(2)
    deformed_bond_NP1[2][2][1] = 10
    deformed_bond_NP1[2][2][2] = -10
    deformed_bond_length_NP1[2][2] = sqrt(200)
    deformed_bond_NP1[2][3][1] = 0.1
    deformed_bond_NP1[2][3][2] = 0
    deformed_bond_length_NP1[2][3] = 0.1

    nodes = Vector{Int64}(1:nnodes)
    thermal_parameter = Dict("Thermal Expansion Coefficient" => 1.0,
                             "Reference Temperature" => 0.0)

    temperature_NP1 .= 0
    PeriLab.Solver_Manager.Model_Factory.Thermal.Thermal_Expansion.compute_model(nodes,
                                                                                 thermal_parameter,
                                                                                 1, 1.0,
                                                                                 1.0)

    for iID in nodes
        for jID in nn[iID]
            for i in 1:dof
                @test deformed_bond_NP1[iID][jID][i] == undeformed_bond[iID][jID][i]
            end
        end
    end

    # temperature_NP1 .= 1
    # PeriLab.Solver_Manager.Model_Factory.Thermal.Thermal_Expansion.compute_model(nodes, thermal_parameter, 1, 1.0, 1.0)

    # for iID in nodes
    #     for jID in nn[iID]
    #         @test deformed_bond_NP1[iID][jID] == .-undeformed_bond[iID][jID]
    #     end
    # end

    # temperature_NP1 .= 2
    # PeriLab.Solver_Manager.Model_Factory.Thermal.Thermal_Expansion.compute_model(nodes, thermal_parameter, 1, 1.0, 1.0)

    # for iID in nodes
    #     for jID in nn[iID]
    #         @test isapprox(deformed_bond_NP1[iID][jID] .+ 1,
    #                        .-undeformed_bond[iID][jID] .* 2 .+ 1)
    #     end
    # end

    # temperature_NP1[1] = 2
    # temperature_NP1[2] = -23

    # thermal_parameter = Dict("Thermal Expansion Coefficient" => [-1.1,2.1], "Reference Temperature" => 0.0)
    # PeriLab.Solver_Manager.Model_Factory.Thermal.Thermal_Expansion.compute_model(nodes, thermal_parameter, 1, 1.0, 1.0)

    # for iID in nodes
    #     for jID in nn[iID]
    #         @test isapprox(deformed_bond_NP1[iID][jID][1] + 1,
    #                        .-undeformed_bond[iID][jID][1] * -1.1 * temperature_NP1[iID] +
    #                        1)
    #         @test isapprox(deformed_bond_NP1[iID][jID][2] + 1,
    #                        .-undeformed_bond[iID][jID][2] * 2.1 * temperature_NP1[iID] +
    #                        1)
    #     end
    # end
end
