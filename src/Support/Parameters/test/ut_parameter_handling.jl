include("../parameter_handling.jl")
using Test
@testset "ut_check_key_elements" begin
    params = Dict()
    @test check_key_elements(params) == false
    params = Dict("Damage" => Dict())
    @test check_key_elements(params) == false
    params = Dict("Damage" => Dict(), "Materials" => Dict())
    @test check_key_elements(params) == false
    params = Dict("Damage" => Dict(), "Materials" => Dict(), "Discretization" => Dict())
    @test check_key_elements(params) == false
    params = Dict("Damage" => Dict(), "Materials" => Dict(), "Discretization" => Dict(), "Blocks" => Dict())
    @test check_key_elements(params) == false
    params = Dict("Damage" => Dict(), "Materials" => Dict(), "Discretization" => Dict(), "Blocks" => Dict(), "Solver" => Dict())
    @test check_key_elements(params) == true
end