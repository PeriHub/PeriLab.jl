# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
#include("../../../src/PeriLab.jl")

include("../../../../src/FEM/Coupling/Arlequin.jl")
@testset "ut_coupling_name" begin
    @test Arlequin.coupling_name() == "Arlequin"
end

Arlequin.find_point_in_elements(coordinates, topology, points_to_check)
