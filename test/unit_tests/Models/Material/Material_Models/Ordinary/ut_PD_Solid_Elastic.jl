# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
using Test
include("../../../../../../src/Models/Material/Material_Models/Ordinary/PD_Solid_Elastic.jl")

@testset "get_name&fe_support" begin
    @test PD_Solid_Elastic.material_name() == "PD Solid Elastic"
    @test !(PD_Solid_Elastic.fe_support())
end
