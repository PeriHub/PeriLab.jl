# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
include("../../../src/Support/data_manager.jl")
include("../../../src/Compute/compute_forces.jl")
@testset "ut_get_forces_from_force_density" begin
    testDatamanager = Data_manager
    testDatamanager.set_nmasters(5)

    testDatamanager.create_node_field("Forces", Float32, 3)

    fdN, fdNP1 = testDatamanager.create_node_field("Force Densities", Float32, 3)
    volume = testDatamanager.create_constant_node_field("Volume", Float32, 1)

    volume[1:5] = 1:5
    fdNP1 = rand(5, 3)

    testDatamanager = get_forces_from_force_density(testDatamanager)
    forces = testDatamanager.get_field("Forces", "NP1")
    for i in 1:5
        for j in 1:3
            @test forces[i, j] / (fdNP1[i, j] * volume[i]) - 1 < 1e-8
        end
    end
end
