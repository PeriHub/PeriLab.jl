# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test

@testset "ut_get_forces_from_force_density" begin
    test_data_manager = PeriLab.Data_Manager
    test_data_manager.initialize_data()
    test_data_manager.set_num_controller(5)

    test_data_manager.create_node_field("Forces", Float64, 3)

    force_densityN,
    force_densityNP1 = test_data_manager.create_node_field("Force Densities",
                                                           Float64, 3)
    volume = test_data_manager.create_constant_node_field("Volume", Float64, 1)

    volume[1:5] = 1:5
    force_densityNP1 = rand(5, 3)

    test_data_manager = PeriLab.Solver_Manager.Verlet_Solver.get_forces_from_force_density(test_data_manager)
    forces = test_data_manager.get_field("Forces", "NP1")
    for i in 1:5
        for j in 1:3
            @test forces[i, j] / (force_densityNP1[i, j] * volume[i]) - 1 < 1e-8
        end
    end
end

@testset "ut_get_partial_stresses" begin
    test_data_manager = PeriLab.Data_Manager
    test_data_manager.set_dof(3)
    test_data_manager.set_num_controller(5)
    nn = test_data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1] = 1
    nn[2] = 2
    nn[3] = 1
    nn[4] = 2
    nn[5] = 3

    nlist = test_data_manager.create_constant_bond_field("Neighborhoodlist", Int64, 1)
    nlist[1] = [2]
    nlist[2] = [1, 3]
    nlist[3] = [1]
    nlist[4] = [1, 5]
    nlist[5] = [1, 2, 4]

    test_data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    bond_forces = test_data_manager.create_constant_bond_field("Bond Forces", Float64, 3, 1)
    bond_geometry = test_data_manager.create_constant_bond_field("Bond Geometry", Float64,
                                                                 3, 1)
    bond_length = test_data_manager.create_constant_bond_field("Bond Length", Float64, 1,
                                                               sqrt(3))
    test_data_manager.create_node_field("Cauchy Stress", Float64, 3,
                                        VectorOrMatrix = "Matrix")
    volume = test_data_manager.create_constant_node_field("Volume", Float64, 1)

    volume[1:5] = 1:5

    nodes = [1, 2, 3, 4, 5]

    test_data_manager = PeriLab.Solver_Manager.Verlet_Solver.get_partial_stresses(test_data_manager, nodes)
    stresses = test_data_manager.get_field("Cauchy Stress", "NP1")

    testval = zeros(5, 3, 3)
    for iID in 1:5
        for jID in eachindex(nlist[iID])
            for i in 1:3
                for j in 1:3
                    testval[iID, i,
                            j] += bond_forces[iID][jID][i] *
                                  bond_geometry[iID][jID][j] * volume[iID]
                end
            end
        end
    end

    for iID in 1:5
        for i in 1:3
            for j in 1:3
                @test testval[iID, i, j] == stresses[iID, i, j]
            end
        end
    end
end
