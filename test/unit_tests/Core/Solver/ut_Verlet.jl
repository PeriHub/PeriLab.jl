# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

#using Test
using MPI

@testset "ut_get_integration_steps" begin
    @test_logs (:error, "Time step -1.0 [s] is not valid") PeriLab.Solver_Manager.Verlet_Solver.get_integration_steps(0.0,
                                                                                                                      0.0,
                                                                                                                      -1.0)
    @test PeriLab.Solver_Manager.Verlet_Solver.get_integration_steps(0.0, 1.0, 1.0) ==
          (1, 1.0)
    @test PeriLab.Solver_Manager.Verlet_Solver.get_integration_steps(0.0, 2.0, 1.0) ==
          (2, 1.0)
    @test PeriLab.Solver_Manager.Verlet_Solver.get_integration_steps(0.0, 6.0, 2.0) ==
          (3, 2.0)
    @test PeriLab.Solver_Manager.Verlet_Solver.get_integration_steps(2.0, 6.0, 2.0) ==
          (2, 2.0)
end
