using Test
include("../helpers.jl")
@testset "ut_find_indices" begin
    @test find_indices([1, 1, 2, 3, 3, 4, 4, 4, 1], 1) == [1, 2, 9]
    @test find_indices([1, 1, 2, 3, 3, 4, 4, 4, 1], 2) == [3]
    @test find_indices([1, 1, 2, 3, 3, 4, 4, 4, 1], 3) == [4, 5]
    @test find_indices([1, 1, 2, 3, 3, 4, 4, 4, 1], 4) == [6, 7, 8]
    @test find_indices([1, 1, 2, 3, 3, 4, 4, 4, 1], 5) == []
end