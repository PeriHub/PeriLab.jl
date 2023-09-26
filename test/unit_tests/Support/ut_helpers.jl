# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
include("../../../src/Support/helpers.jl")
@testset "ut_find_indices" begin
    @test find_indices([1, 1, 2, 3, 3, 4, 4, 4, 1], 1) == [1, 2, 9]
    @test find_indices([1, 1, 2, 3, 3, 4, 4, 4, 1], 2) == [3]
    @test find_indices([1, 1, 2, 3, 3, 4, 4, 4, 1], 3) == [4, 5]
    @test find_indices([1, 1, 2, 3, 3, 4, 4, 4, 1], 4) == [6, 7, 8]
    @test find_indices([1, 1, 2, 3, 3, 4, 4, 4, 1], 5) == []
end