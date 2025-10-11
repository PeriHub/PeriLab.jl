# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
#include("../../../src/PeriLab.jl")
@testset "ut_coupling_name" begin
    @test PeriLab.Solver_control.FEM.Coupling_PD_FEM.Arlequin_coupling.coupling_name() == "Arlequin"
end
@testset "ut_find_point_in_elements" begin
    dof = 2
    elements = 1:2  # FE        #@info elements
    nodes = 1:10   # beide
    nFEnodes = 6            # number of FE nodes
    nodesFE = nodes[1:nFEnodes, 1]
    nodesPD = nodes[(nFEnodes + 1):end, 1]
    topology = Array{Int64}(zeros(2, 4))
    topology[1, :] = [1, 2, 6, 4]
    topology[2, :] = [6, 4, 5, 3]

    coordinates = zeros(11, 2)
    coordinates[1, :] = [2, 0]
    coordinates[2, :] = [2, 1]      #FE
    coordinates[3, :] = [0, 1]
    coordinates[4, :] = [1, 1]      #FE
    coordinates[5, :] = [0, 0]
    coordinates[6, :] = [1, 0]      #FE
    coordinates[7, :] = [0.5, 0.5]
    coordinates[8, :] = [1.5, 0.5]     #PD
    coordinates[9, :] = [0.5, 1.5]
    coordinates[10, :] = [1.5, 1.5]    #PD

    test_dict = PeriLab.Solver_control.FEM.Coupling_PD_FEM.Arlequin_coupling.find_point_in_elements(coordinates,
                                                         topology,
                                                         nodesPD,
                                                         2)
    @test collect(keys(test_dict)) == [4, 2, 1]
    @test test_dict[4] == 1
    @test test_dict[2] == 1
    @test test_dict[1] == 1
end

@testset "ut_compute_coupling_matrix" begin
    dof = 2
    elements = 1:2  # FE        #@info elements
    nodes = 1:10   # beide
    nFEnodes = 6            # number of FE nodes
    nodesFE = nodes[1:nFEnodes, 1]
    nodesPD = nodes[(nFEnodes + 1):end, 1]
    topology = Array{Int64}(zeros(2, 4))
    topology[1, :] = [1, 2, 4, 6]
    topology[2, :] = [6, 4, 3, 5]

    coordinates = zeros(11, 2)
    coordinates[1, :] = [2, 0]
    coordinates[2, :] = [2, 1]      #FE
    coordinates[3, :] = [0, 1]
    coordinates[4, :] = [1, 1]      #FE
    coordinates[5, :] = [0, 0]
    coordinates[6, :] = [1, 0]      #FE
    coordinates[7, :] = [0.5, 0.5]
    coordinates[8, :] = [1.5, 0.5]     #PD
    coordinates[9, :] = [0.5, 1.5]
    coordinates[10, :] = [1.5, 1.5]    #PD
    kappa = 1

    p = [1, 1]
    dof = 2

    test_mat = PeriLab.Solver_control.FEM.Coupling_PD_FEM.Arlequin_coupling.compute_coupling_matrix(coordinates,
                                                         topology,
                                                         7,
                                                         2,
                                                         kappa,
                                                         p,
                                                         dof)

    @test test_mat == [1.0 -0.25 -0.25 -0.25 -0.25
           -0.25 0.0625 0.0625 0.0625 0.0625
           -0.25 0.0625 0.0625 0.0625 0.0625
           -0.25 0.0625 0.0625 0.0625 0.0625
           -0.25 0.0625 0.0625 0.0625 0.0625]
    test_mat = PeriLab.Solver_control.FEM.Coupling_PD_FEM.Arlequin_coupling.compute_coupling_matrix(coordinates,
                                                         topology,
                                                         8,
                                                         1,
                                                         kappa,
                                                         p,
                                                         dof)
    @test test_mat == [1.0 -0.25 -0.25 -0.25 -0.25
           -0.25 0.0625 0.0625 0.0625 0.0625
           -0.25 0.0625 0.0625 0.0625 0.0625
           -0.25 0.0625 0.0625 0.0625 0.0625
           -0.25 0.0625 0.0625 0.0625 0.0625]
end
