# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
using Test
include("../../../../../src/Physics/Material/Material_Models/Correspondence_Plastic.jl")


@testset "get_name&fe_support" begin
    @test Correspondence_Plastic.correspondence_name() == "Correspondence Plastic"
    @test !(Correspondence_Plastic.fe_support())
end

