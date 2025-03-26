# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
#include("../../../src/PeriLab.jl")

include("../../../../src/FEM/Coupling/Schwarz_coupling.jl")
@testset "ut_coupling_name" begin
    @test Schwarz_coupling.coupling_name() == "Schwarz"
end
@testset "ut_init_coupling_model" begin
    test_data_manager = PeriLab.Data_manager
    @test isnothing(Schwarz_coupling.init_coupling_model(test_data_manager, 1:2, Dict()))
end
