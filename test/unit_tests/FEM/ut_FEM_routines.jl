# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
include("../../../src/FEM/FEM_routines.jl")

using Test
@testset "ut_get_weights_and_integration_points" begin
    @test get_weights_and_integration_points(2, [0, 0]) == ([2.0; 2.0;;], [0.0; 0.0;;])
    @test get_weights_and_integration_points(2, [1, 1]) == ([1.0 1.0; 1.0 1.0], [-0.5773502691896258 0.5773502691896258; -0.5773502691896258 0.5773502691896258])
    @test get_weights_and_integration_points(3, [1, 2, 3]) == ([1.0 1.0 0.0 0.0; 0.5555555555555556 0.8888888888888888 0.5555555555555556 0.0; 0.34785484513745385 0.6521451548625462 0.6521451548625462 0.34785484513745385], [-0.5773502691896258 0.5773502691896258 0.0 0.0; -0.7745966692414834 0.0 0.7745966692414834 0.0; -0.8611363115940526 -0.3399810435848563 0.3399810435848563 0.8611363115940526])
    @test get_weights_and_integration_points(3, [1, 0, 0]) == ([1.0 1.0; 2.0 0.0; 2.0 0.0], [-0.5773502691896258 0.5773502691896258; 0.0 0.0; 0.0 0.0])
end

@testset "ut_get_weights_and_integration_points" begin
    @test isnothing(define_lagarangian_grid_space(2, [0, 0]))
    @test isnothing(define_lagarangian_grid_space(3, [0, 0, 0]))
    @test isnothing(define_lagarangian_grid_space(3, [0, 1, 0]))
    @test isnothing(define_lagarangian_grid_space(3, [1, 1, 0]))
    @test define_lagarangian_grid_space(2, [1, 1]) == [-1.0 1.0; -1.0 1.0]
    @test define_lagarangian_grid_space(3, [1, 1, 1]) == [-1.0 1.0; -1.0 1.0; -1.0 1.0]
    @test define_lagarangian_grid_space(2, [1, 2]) == [-1.0 1.0 0.0; -1.0 0.0 1.0]
    @test define_lagarangian_grid_space(2, [1, 4]) == [-1.0 1.0 0.0 0.0 0.0; -1.0 -0.5 0.0 0.5 1.0]
    @test define_lagarangian_grid_space(3, [2, 2, 1]) == [-1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 1.0 0.0]
    @test define_lagarangian_grid_space(3, [2, 1, 1]) == [-1.0 0.0 1.0; -1.0 1.0 0.0; -1.0 1.0 0.0]
    @test define_lagarangian_grid_space(3, [2, 1, 3]) == [-1.0 0.0 1.0 0.0; -1.0 1.0 0.0 0.0; -1.0 -0.33333333333333337 0.33333333333333326 1.0]
end
@testset "ut_get_recursive_lagrange_shape_functions" begin
    p = 1
    xi = define_lagarangian_grid_space(2, [p, p])
    @test get_recursive_lagrange_shape_functions(xi[1, :], 0.0, p) == [0.5, 0.5]
    weights, integration_points = get_weights_and_integration_points(2, [p, p])
    @test get_recursive_lagrange_shape_functions(xi[1, :], integration_points[1, 1], p) == [0.7886751345948129, 0.21132486540518708]
    @test get_recursive_lagrange_shape_functions(xi[1, :], integration_points[1, 2], p) == [0.21132486540518708, 0.7886751345948129]
    p = 2
    weights, integration_points = get_weights_and_integration_points(2, [p, p])
    xi = define_lagarangian_grid_space(2, [p, p])
    @test get_recursive_lagrange_shape_functions(xi[1, :], integration_points[1, 1], p) == [0.6872983346207417, 0.39999999999999997, -0.08729833462074169]
    @test get_recursive_lagrange_shape_functions(xi[1, :], integration_points[1, 2], p) == [0.0, 1.0, 0.0]
    @test get_recursive_lagrange_shape_functions(xi[1, :], integration_points[1, 3], p) == [-0.08729833462074169, 0.39999999999999997, 0.6872983346207417]
    p = 3
    weights, integration_points = get_weights_and_integration_points(2, [p, p])
    xi = define_lagarangian_grid_space(2, [p, p])
    @test get_recursive_lagrange_shape_functions(xi[1, :], integration_points[1, 1], p) == [0.6600056650728034, 0.5209376877117036, -0.23018790325073893, 0.04924455046623181]
    @test get_recursive_lagrange_shape_functions(xi[1, :], integration_points[1, 2], p) == [0.0033737364327725374, 1.0048858548256456, -0.009921353572324671, 0.0016617623139063968]
    @test get_recursive_lagrange_shape_functions(xi[1, :], integration_points[1, 3], p) == [0.0016617623139064258, -0.009921353572324838, 1.0048858548256463, 0.003373736432772594]
    @test get_recursive_lagrange_shape_functions(xi[1, :], integration_points[1, 4], p) == [0.04924455046623184, -0.23018790325073904, 0.5209376877117036, 0.6600056650728034]
end