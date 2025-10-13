# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
@testset "get_quad_horizon" begin
    horizon::Float64 = 1.0
    thickness::Float64 = 2.0
    @test PeriLab.Solver_Manager.Model_Factory.Damage.Critical_Energy.get_quad_horizon(horizon, 3, thickness) ==
          Float64(4 / (pi * horizon^4))
    @test PeriLab.Solver_Manager.Model_Factory.Damage.Critical_Energy.get_quad_horizon(horizon, 2, thickness) ==
          Float64(3 / (pi * horizon^3 * thickness))

    horizon = 5.6
    thickness = 3.0
    @test PeriLab.Solver_Manager.Model_Factory.Damage.Critical_Energy.get_quad_horizon(horizon, 3, thickness) ==
          Float64(4 / (pi * horizon^4))
    @test PeriLab.Solver_Manager.Model_Factory.Damage.Critical_Energy.get_quad_horizon(horizon, 2, thickness) ==
          Float64(3 / (pi * horizon^3 * thickness))
end
