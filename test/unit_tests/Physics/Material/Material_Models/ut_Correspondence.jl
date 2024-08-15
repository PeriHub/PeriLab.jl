# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

include("../../../../../src/Physics/Material/Material_Models/Correspondence.jl")
#include("../../../../../src/PeriLab.jl")
# include("../../../../../src/Core/data_manager.jl")
using Test
using .Correspondence
#using .PeriLab
@testset "zero_energy_mode_compensation_exception" begin

    test_data_manager = PeriLab.Data_manager
    nnodes = 2
    test_data_manager.set_num_controller(nnodes)
    nn = test_data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1] = 2
    nn[2] = 3
    nodes = Vector{Int64}(1:2)

    @test Correspondence.zero_energy_mode_compensation(
        test_data_manager,
        nodes,
        Dict(),
        0.0,
        0.0,
    ) == test_data_manager
end
