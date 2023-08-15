using Test
include("../Verlet.jl")
@testset "ut_test_timestep" begin
    @test test_timestep(1, 2) == 1
    @test test_timestep(2, 1.1) == 1.1
    @test test_timestep(2, 2) == 2
end

@testset "ut_get_cs_denominator" begin
    volume = [1, 2, 3]
    bondgeometry = [1, 2, 3]
    @test get_cs_denominator(1, volume, bondgeometry) == 1
    @test get_cs_denominator(2, volume, bondgeometry) == 2
    @test get_cs_denominator(3, volume, bondgeometry) == 3
    bondgeometry = [2, 4, 6]
    @test get_cs_denominator(1, volume, bondgeometry) == 0.5
    @test get_cs_denominator(2, volume, bondgeometry) == 1
    @test get_cs_denominator(3, volume, bondgeometry) == 1.5
    bondgeometry = [1, 0.5, 2]
    @test get_cs_denominator(1, volume, bondgeometry) == 1
    @test get_cs_denominator(2, volume, bondgeometry) == 5
    @test get_cs_denominator(3, volume, bondgeometry) == 6.5
end

