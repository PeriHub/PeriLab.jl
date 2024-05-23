# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test

include("../../../../src/Physics/Thermal/thermal_flow.jl")
using .Thermal_Flow
# include("../../../../src/Core/data_manager.jl")

@test Thermal_Flow.thermal_model_name() == "Thermal Flow"

@testset "ut_compute_thermal_model" begin
    test_Data_manager = PeriLab.Data_manager
    @test Thermal_Flow.compute_thermal_model(test_Data_manager, Vector{Int64}(1:3), Dict("a" => 1), 1.0, 1.0) == test_Data_manager
end
