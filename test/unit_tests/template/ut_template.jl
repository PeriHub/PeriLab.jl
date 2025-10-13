# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
#include("../../../src/PeriLab.jl")
#using .PeriLab

"""
    Final test set should be included in runtest.jl
    Because PeriLab is initialized there it can be called inside the test. For local debugging, the upper comments have to be adapted. Before you commit, please comment these lines again
    #include("../../../src/PeriLab.jl")
    #using .PeriLab
"""
@testset "ut_compute_bond_level_deformation_gradient" begin
    # test_data_manager = PeriLab.Data_Manager
end
