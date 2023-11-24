# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
include("../../../src/FEM/FEM_routines.jl")

using Test
@testset "ut_get_weights_and_integration_points" begin
    @test get_weights_and_integration_points(2, [0, 0]) == ([0.0; 0.0;;], [2.0; 2.0;;])
    @test get_weights_and_integration_points(2, [1, 1]) == ([-0.5773502691896258 0.5773502691896258; -0.5773502691896258 0.5773502691896258], [1.0 1.0; 1.0 1.0])
    @test get_weights_and_integration_points(3, [1, 2, 3]) == ([-0.5773502691896258 0.5773502691896258 0.0 0.0; -0.7745966692414834 0.0 0.7745966692414834 0.0; -0.8611363115940526 -0.3399810435848563 0.3399810435848563 0.8611363115940526], [1.0 1.0 0.0 0.0; 0.5555555555555556 0.8888888888888888 0.5555555555555556 0.0; 0.34785484513745385 0.6521451548625462 0.6521451548625462 0.34785484513745385])
    @test get_weights_and_integration_points(3, [1, 0, 0]) == ([-0.5773502691896258 0.5773502691896258; 0.0 0.0; 0.0 0.0], [1.0 1.0; 2.0 0.0; 2.0 0.0])
end

@testset "ut_get_weights_and_integration_points" begin
    @test isnothing(define_lagarangian_grid_space(2, [0, 0]))
    @test isnothing(define_lagarangian_grid_space(3, [0, 0]))
    @test isnothing(define_lagarangian_grid_space(3, [0, 1]))
    @test isnothing(define_lagarangian_grid_space(3, [1, 0]))
    @test define_lagarangian_grid_space(2, [1, 1]) == [-1.0 1.0; -1.0 1.0]
    @test define_lagarangian_grid_space(3, [1, 1, 1]) == [-1.0 1.0; -1.0 1.0; -1.0 1.0]
    @test define_lagarangian_grid_space(2, [1, 2]) == [-1.0 1.0 0.0; -1.0 0.0 1.0]
    @test define_lagarangian_grid_space(2, [1, 4]) == [-1.0 1.0 0.0 0.0 0.0; -1.0 -0.5 0.0 0.5 1.0]
    @test define_lagarangian_grid_space(3, [2, 2, 1]) == [-1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 1.0 0.0]
    @test define_lagarangian_grid_space(3, [2, 1, 1]) == [-1.0 0.0 1.0; -1.0 1.0 0.0; -1.0 1.0 0.0]
    @test define_lagarangian_grid_space(3, [2, 1, 3]) == [-1.0 0.0 1.0 0.0; -1.0 1.0 0.0 0.0; -1.0 -0.33333333333333337 0.33333333333333326 1.0]
end