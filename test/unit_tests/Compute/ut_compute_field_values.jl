# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
include("../../../src/Support/data_manager.jl")
include("../../../src/Compute/compute_field_values.jl")
@testset "ut_get_forces_from_force_density" begin
    test_Data_manager = Data_manager
    test_Data_manager.set_num_controller(5)

    test_Data_manager.create_node_field("Forces", Float64, 3)

    fdN, fdNP1 = test_Data_manager.create_node_field("Force Densities", Float64, 3)
    volume = test_Data_manager.create_constant_node_field("Volume", Float64, 1)

    volume[1:5] = 1:5
    fdNP1 = rand(5, 3)

    test_Data_manager = get_forces_from_force_density(test_Data_manager)
    forces = test_Data_manager.get_field("Forces", "NP1")
    for i in 1:5
        for j in 1:3
            @test forces[i, j] / (fdNP1[i, j] * volume[i]) - 1 < 1e-8
        end
    end
end