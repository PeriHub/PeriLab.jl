# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test

include("../../../../src/Physics/Thermal/thermal_flow.jl")
using .Thermal_Flow
include("../../../../src/Support/data_manager.jl")

@test Thermal_Flow.thermal_model_name() == "Thermal Flow"

@testset "ut_compute_thermal_model" begin
    nnodes = 2
    dof = 2
    testDatamanager = Data_manager
    testDatamanager.set_nmasters(2)
    @test Thermal_Flow.compute_thermal_model(testDatamanager, Vector{Int64}(1:nnodes), Dict(), 1.0, 1.0) == testDatamanager
end
