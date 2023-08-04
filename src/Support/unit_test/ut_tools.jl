include("../tools.jl")

using Test

@testset "get_set_functions" begin
    a = ones(2, 2)
    @test check_inf_or_nan(a, "a") == false
    a[1, 1] = 1 / 0
    @test check_inf_or_nan(a, "a")
    a = 0
    @test check_inf_or_nan(a, "a") == false
end