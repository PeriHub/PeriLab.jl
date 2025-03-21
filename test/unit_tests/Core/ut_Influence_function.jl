# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test

#include("../../../src/PeriLab.jl")
#using .PeriLab

@testset "init_influence_function" begin
    test_data_manager = PeriLab.Data_manager
    @test PeriLab.Solver_control.Influence_function.init_influence_function(Vector{Int64}(1:2),
                                                                            test_data_manager,
                                                                            Dict()) ==
          test_data_manager
end
