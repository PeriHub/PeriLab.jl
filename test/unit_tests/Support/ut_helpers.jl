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

# Define a test case for the find_active function
@testset "ut_find_active" begin
    # Test case 1: Empty input
    @test isempty(find_active(Bool[]))
    # Test case 2: All elements are active
    @test find_active([true, true, true]) == [1, 2, 3]
    # Test case 3: No elements are active
    @test isempty(find_active([false, false, false]))
    # Test case 4: Mix of active and inactive elements
    @test find_active([false, true, false, true, true]) == [2, 4, 5]
end