# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

#using Test
@testset "get_name&fe_support" begin
    @test parentmodule(PeriLab.ParameterSpec.lookup_model(:material, "PD Solid Elastic")) ===
      PeriLab.Solver_Manager.Model_Factory.Material.PD_Solid_Elastic
    @test !(PeriLab.Solver_Manager.Model_Factory.Material.PD_Solid_Elastic.fe_support())
end
