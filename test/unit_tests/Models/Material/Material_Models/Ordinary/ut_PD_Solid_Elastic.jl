# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
@testset "get_name&fe_support" begin
    @test PeriLab.Solver_control.Model_Factory.Material.PD_Solid_Elastic.material_name() == "PD Solid Elastic"
    @test !(PeriLab.Solver_control.Model_Factory.Material.PD_Solid_Elastic.fe_support())
end
