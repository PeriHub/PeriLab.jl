# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../../../src/Physics/Damage/Energy_release.jl")
include("../../../../src/Support/data_manager.jl")
using Test
using .Critical_Energy_Model
@testset "get_quad_horizon" begin
    horizon::Float64 = 1.0
    @test Critical_Energy_Model.get_quad_horizon(horizon, 3) == Float64(4 / (pi * horizon^4))
    @test Critical_Energy_Model.get_quad_horizon(horizon, 2) == Float64(3 / (pi * horizon^3))

    horizon = 5.6
    @test Critical_Energy_Model.get_quad_horizon(horizon, 3) == Float64(4 / (pi * horizon^4))
    @test Critical_Energy_Model.get_quad_horizon(horizon, 2) == Float64(3 / (pi * horizon^3))
end
