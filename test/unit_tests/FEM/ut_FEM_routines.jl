# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
include("../../../src/FEM/FEM_routines.jl")

using Test
@testset "ut_get_weights_and_integration_points" begin
    @test get_weights_and_integration_points(2, [0, 0]) == ([2.0 2.0], [0.0 0.0])
    @test get_weights_and_integration_points(2, [1, 1]) == ([1.0 1.0; 1.0 1.0], [-0.5773502691896258 -0.5773502691896258; 0.5773502691896258 0.5773502691896258])
    @test get_weights_and_integration_points(3, [1, 2, 3]) == ([1.0 0.5555555555555556 0.34785484513745385; 1.0 0.8888888888888888 0.6521451548625462; 0.0 0.5555555555555556 0.6521451548625462; 0.0 0.0 0.34785484513745385], [-0.5773502691896258 -0.7745966692414834 -0.8611363115940526; 0.5773502691896258 0.0 -0.3399810435848563; 0.0 0.7745966692414834 0.3399810435848563; 0.0 0.0 0.8611363115940526])
    @test get_weights_and_integration_points(3, [1, 0, 0]) == ([1.0 2.0 2.0; 1.0 0.0 0.0], [-0.5773502691896258 0.0 0.0; 0.5773502691896258 0.0 0.0])
end
@testset "ut_get_multi_dimensional_integration_points" begin
    @test isnothing(get_multi_dimensional_integration_points(1, [1], zeros(2, 2)))
    @test isnothing(get_multi_dimensional_integration_points(4, [1], zeros(2, 2)))

    dof = 2
    p = 1
    weights, xi = get_weights_and_integration_points(dof, [p, p])
    integration_point_coordinates = get_multi_dimensional_integration_points(dof, [p, p], xi)
    @test length(integration_point_coordinates[:, 1]) == 4
    @test integration_point_coordinates[1, :] == [-0.5773502691896258, -0.5773502691896258]
    @test integration_point_coordinates[2, :] == [-0.5773502691896258, 0.5773502691896258]
    @test integration_point_coordinates[3, :] == [0.5773502691896258, -0.5773502691896258]
    @test integration_point_coordinates[4, :] == [0.5773502691896258, 0.5773502691896258]
    dof = 3
    weights, xi = get_weights_and_integration_points(dof, [p, p, p])
    integration_point_coordinates = get_multi_dimensional_integration_points(dof, [p, p, p], xi)
    @test length(integration_point_coordinates[:, 1]) == 8
    @test integration_point_coordinates[1, :] == [-0.5773502691896258, -0.5773502691896258, -0.5773502691896258]
    @test integration_point_coordinates[2, :] == [-0.5773502691896258, -0.5773502691896258, 0.5773502691896258]
    @test integration_point_coordinates[3, :] == [-0.5773502691896258, 0.5773502691896258, -0.5773502691896258]
    @test integration_point_coordinates[4, :] == [-0.5773502691896258, 0.5773502691896258, 0.5773502691896258]
    @test integration_point_coordinates[5, :] == [0.5773502691896258, -0.5773502691896258, -0.5773502691896258]
    @test integration_point_coordinates[6, :] == [0.5773502691896258, -0.5773502691896258, 0.5773502691896258]
    @test integration_point_coordinates[7, :] == [0.5773502691896258, 0.5773502691896258, -0.5773502691896258]
    @test integration_point_coordinates[8, :] == [0.5773502691896258, 0.5773502691896258, 0.5773502691896258]
    dof = 2
    weights, xi = get_weights_and_integration_points(dof, [p, p + 1])
    integration_point_coordinates = get_multi_dimensional_integration_points(dof, [p, p + 1], xi)
    @test length(integration_point_coordinates[:, 1]) == 6
    @test integration_point_coordinates[1, :] == [-0.5773502691896258, -0.7745966692414834]
    @test integration_point_coordinates[2, :] == [-0.5773502691896258, 0.0]
    @test integration_point_coordinates[3, :] == [-0.5773502691896258, 0.7745966692414834]
    @test integration_point_coordinates[4, :] == [0.5773502691896258, -0.7745966692414834]
    @test integration_point_coordinates[5, :] == [0.5773502691896258, 0.0]
    @test integration_point_coordinates[6, :] == [0.5773502691896258, 0.7745966692414834]
    weights, xi = get_weights_and_integration_points(dof, [p + 1, p])
    integration_point_coordinates = get_multi_dimensional_integration_points(dof, [p + 1, p], xi)

    @test length(integration_point_coordinates[:, 1]) == 6
    @test integration_point_coordinates[1, :] == [-0.7745966692414834, -0.5773502691896258]
    @test integration_point_coordinates[2, :] == [-0.7745966692414834, 0.5773502691896258]
    @test integration_point_coordinates[3, :] == [0.0, -0.5773502691896258]
    @test integration_point_coordinates[4, :] == [0.0, 0.5773502691896258]
    @test integration_point_coordinates[5, :] == [0.7745966692414834, -0.5773502691896258]
    @test integration_point_coordinates[6, :] == [0.7745966692414834, 0.5773502691896258]
end