# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
include("../../../src/FEM/FEM_factory.jl")

using Test

@testset "ut_get_polynomial_degree" begin

    @test isnothing(FEM.get_polynomial_degree(Dict(), 1))
    @test isnothing(FEM.get_polynomial_degree(Dict(), 2))
    @test isnothing(FEM.get_polynomial_degree(Dict(), 3))

    params = Dict("Degree" => 1)
    println()
    @test FEM.get_polynomial_degree(params, 2) == [1, 1]
    @test FEM.get_polynomial_degree(params, 3) == [1, 1, 1]

    params = Dict("Degree" => 2)
    @test FEM.get_polynomial_degree(params, 2) == [2, 2]
    @test FEM.get_polynomial_degree(params, 3) == [2, 2, 2]

    params = Dict("Degree" => 2.1)
    println()
    @test FEM.get_polynomial_degree(params, 2) == [2, 2]
    @test FEM.get_polynomial_degree(params, 3) == [2, 2, 2]

    params = Dict("Degree" => [2 3 1])
    @test isnothing(FEM.get_polynomial_degree(params, 2))
    @test FEM.get_polynomial_degree(params, 3) == [2, 3, 1]

    params = Dict("Degree" => [2.1 2])
    @test FEM.get_polynomial_degree(params, 2) == [2, 2]
    @test isnothing(FEM.get_polynomial_degree(params, 3))

end
