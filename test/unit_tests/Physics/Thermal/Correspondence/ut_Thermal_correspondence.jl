using Test
include("../../../../../src/Physics/Thermal/Correspondence/Thermal_correspondence.jl")

@testset "ut_thermal strain" begin
    strain = zeros(3, 3)
    alpha = ones(3, 3)
    temperature = 1.0
    @test -alpha == thermal_strain(alpha, temperature, strain)
    strain = ones(3, 3)
    @test zeros(3, 3) == thermal_strain(alpha, temperature, strain)
    @test ones(3, 3) == thermal_strain(alpha, 0.0, strain)

end
