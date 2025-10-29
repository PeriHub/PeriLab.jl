# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test

#using PeriLab

@testset "init_influence_function" begin
    test_data_manager = PeriLab.Data_Manager
    PeriLab.Solver_Manager.Influence_Function.init_influence_function(Vector{Int64}(1:2),
                                                                      Dict())
end
