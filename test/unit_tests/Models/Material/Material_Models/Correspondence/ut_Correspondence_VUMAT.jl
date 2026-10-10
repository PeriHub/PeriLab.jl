# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

#using Test
using LinearAlgebra
# include("../../../../../../src/PeriLab.jl")
# using .PeriLab
@testset "get_name&fe_support" begin
    @test parentmodule(PeriLab.ParameterSpec.lookup_model(:material, "Correspondence VUMAT")) ===
      PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.Correspondence_VUMAT
    @test PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.Correspondence_VUMAT.fe_support()
end
