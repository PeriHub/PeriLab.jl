# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test

@testset "ut_get_forces_from_force_density" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(5)

    PeriLab.Data_Manager.create_node_vector_field("Forces", Float64, 3)

    force_densityN,
    force_densityNP1 = PeriLab.Data_Manager.create_node_vector_field("Force Densities",
                                                                     Float64, 3)
    volume = PeriLab.Data_Manager.create_constant_node_scalar_field("Volume", Float64)

    volume[1:5] = 1:5
    force_densityNP1 = rand(5, 3)

    PeriLab.Solver_Manager.Verlet_Solver.get_forces_from_force_density()
    forces = PeriLab.Data_Manager.get_field("Forces", "NP1")
    for i in 1:5
        for j in 1:3
            @test forces[i, j] / (force_densityNP1[i, j] * volume[i]) - 1 < 1e-8
        end
    end
end

@testset "ut_get_partial_stresses" begin
    PeriLab.Data_Manager.set_dof(3)
    PeriLab.Data_Manager.set_num_controller(5)
    nn = PeriLab.Data_Manager.create_constant_node_scalar_field("Number of Neighbors",
                                                                Int64)
    nn[1] = 1
    nn[2] = 2
    nn[3] = 1
    nn[4] = 2
    nn[5] = 3

    nlist = PeriLab.Data_Manager.create_constant_bond_scalar_state("Neighborhoodlist",
                                                                   Int64)
    nlist[1] = [2]
    nlist[2] = [1, 3]
    nlist[3] = [1]
    nlist[4] = [1, 5]
    nlist[5] = [1, 2, 4]

    PeriLab.Data_Manager.create_constant_node_scalar_field("Number of Neighbors", Int64)
    bond_forces = PeriLab.Data_Manager.create_constant_bond_vector_state("Bond Forces",
                                                                         Float64, 3)
    bond_geometry = PeriLab.Data_Manager.create_constant_bond_vector_state("Bond Geometry",
                                                                           Float64,
                                                                           3;
                                                                           default_value = 1)
    bond_length = PeriLab.Data_Manager.create_constant_bond_scalar_state("Bond Length",
                                                                         Float64;
                                                                         default_value = sqrt(3))
    PeriLab.Data_Manager.create_node_tensor_field("Cauchy Stress", Float64, 3)
    volume = PeriLab.Data_Manager.create_constant_node_scalar_field("Volume", Float64)

    volume[1:5] = 1:5

    nodes = [1, 2, 3, 4, 5]

    PeriLab.Solver_Manager.Verlet_Solver.get_partial_stresses(nodes)
    stresses = PeriLab.Data_Manager.get_field("Cauchy Stress", "NP1")

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
