# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
include("../../../src/FEM/FEM_routines.jl")
include("../../../src/FEM/Element_formulation/lagrange_element.jl")
using Test
@testset "ut_get_FE_material_model" begin
    params = Dict("FEM" => Dict("FE_1" => Dict("Degree" => 1, "Element Type" => "Lagrange", "Material Model" => "Elastic Model")),
        "Material Models" => Dict("Elastic Model 2" => Dict("Material Model" => "Correspondence Elastic", "Symmetry" => "isotropic plane strain", "Bulk Modulus" => 2.5e+3, "Shear Modulus" => 1.15e3)))

    @test isnothing(get_FE_material_model(params, "FE_1"))

    params = Dict("FEM" => Dict("FE_1" => Dict("Degree" => 1, "Element Type" => "Lagrange", "Material Model" => "Elastic Model")),
        "Material Models" => Dict("Elastic Model" => Dict("Material Model" => "Correspondence Elastic", "Symmetry" => "isotropic plane strain", "Bulk Modulus" => 2.5e+3, "Shear Modulus" => 1.15e3)))

    @test get_FE_material_model(params, "FE_1") == Dict("Material Model" => "Correspondence Elastic", "Symmetry" => "isotropic plane strain", "Bulk Modulus" => 2.5e+3, "Shear Modulus" => 1.15e3)
end


@testset "ut_get_number_of_integration_points" begin
    @test get_number_of_integration_points(Vector{Int64}([1, 1]), 2) == [2, 2]
    @test get_number_of_integration_points(Vector{Int64}([1, 1, 1]), 3) == [2, 2, 2]
    @test get_number_of_integration_points(Vector{Int64}([1, 2, 1]), 3) == [2, 3, 2]
    @test get_number_of_integration_points(Vector{Int64}([1, 3, 8]), 3) == [2, 5, 15]
    @test get_number_of_integration_points(Vector{Int64}([1, 3]), 2) == [2, 5]
end

@testset "ut_get_weights_and_integration_points" begin
    @test get_weights_and_integration_points(2, [1, 1]) == ([2.0 2.0], [0.0 0.0])
    @test get_weights_and_integration_points(2, [2, 2]) == ([1.0 1.0; 1.0 1.0], [-0.5773502691896258 -0.5773502691896258; 0.5773502691896258 0.5773502691896258])
    @test get_weights_and_integration_points(3, [2, 3, 4]) == ([1.0 0.5555555555555556 0.34785484513745385; 1.0 0.8888888888888888 0.6521451548625462; 0.0 0.5555555555555556 0.6521451548625462; 0.0 0.0 0.34785484513745385], [-0.5773502691896258 -0.7745966692414834 -0.8611363115940526; 0.5773502691896258 0.0 -0.3399810435848563; 0.0 0.7745966692414834 0.3399810435848563; 0.0 0.0 0.8611363115940526])
    @test get_weights_and_integration_points(3, [2, 1, 1]) == ([1.0 2.0 2.0; 1.0 0.0 0.0], [-0.5773502691896258 0.0 0.0; 0.5773502691896258 0.0 0.0])
end

@testset "ut_get_multi_dimensional_integration_points" begin
    @test isnothing(get_multi_dimensional_integration_point_data(1, [1], zeros(2, 2)))
    @test isnothing(get_multi_dimensional_integration_point_data(4, [1, 1, 1, 1], zeros(2, 2)))
    dof = 2
    num_int = 2
    weights, xi = get_weights_and_integration_points(dof, [num_int, num_int])
    integration_point_coordinates = get_multi_dimensional_integration_point_data(dof, [num_int, num_int], xi)
    @test length(integration_point_coordinates[:, 1]) == 4
    println()
    @test integration_point_coordinates[1, :] == [-0.5773502691896258, -0.5773502691896258]
    @test integration_point_coordinates[2, :] == [0.5773502691896258, -0.5773502691896258]
    @test integration_point_coordinates[3, :] == [-0.5773502691896258, 0.5773502691896258]
    @test integration_point_coordinates[4, :] == [0.5773502691896258, 0.5773502691896258]
    @test get_multi_dimensional_integration_point_data(dof, [num_int, num_int], weights) == [1.0 1.0; 1.0 1.0; 1.0 1.0; 1.0 1.0]
    dof = 3
    weights, xi = get_weights_and_integration_points(dof, [num_int, num_int, num_int])
    integration_point_coordinates = get_multi_dimensional_integration_point_data(dof, [num_int, num_int, num_int], xi)

    @test length(integration_point_coordinates[:, 1]) == 8
    @test integration_point_coordinates[1, :] == [-0.5773502691896258, -0.5773502691896258, -0.5773502691896258]
    @test integration_point_coordinates[2, :] == [0.5773502691896258, -0.5773502691896258, -0.5773502691896258]
    @test integration_point_coordinates[3, :] == [-0.5773502691896258, 0.5773502691896258, -0.5773502691896258]
    @test integration_point_coordinates[4, :] == [0.5773502691896258, 0.5773502691896258, -0.5773502691896258]
    @test integration_point_coordinates[5, :] == [-0.5773502691896258, -0.5773502691896258, 0.5773502691896258]
    @test integration_point_coordinates[6, :] == [0.5773502691896258, -0.5773502691896258, 0.5773502691896258]
    @test integration_point_coordinates[7, :] == [-0.5773502691896258, 0.5773502691896258, 0.5773502691896258]
    @test integration_point_coordinates[8, :] == [0.5773502691896258, 0.5773502691896258, 0.5773502691896258]
    integration_point_weights = get_multi_dimensional_integration_point_data(dof, [num_int, num_int, num_int], weights)

    for i in 1:8
        @test integration_point_weights[i, :] == [1.0, 1.0, 1.0]
    end

    dof = 2
    weights, xi = get_weights_and_integration_points(dof, [num_int, num_int + 1])
    integration_point_coordinates = get_multi_dimensional_integration_point_data(dof, [num_int, num_int + 1], xi)
    integration_point_weights = get_multi_dimensional_integration_point_data(dof, [num_int, num_int + 1], weights)

    @test length(integration_point_coordinates[:, 1]) == 6
    @test length(integration_point_weights[:, 1]) == 6

    @test integration_point_coordinates[1, :] == [-0.5773502691896258, -0.7745966692414834]
    @test integration_point_coordinates[2, :] == [0.5773502691896258, -0.7745966692414834]
    @test integration_point_coordinates[3, :] == [-0.5773502691896258, 0.0]
    @test integration_point_coordinates[4, :] == [0.5773502691896258, 0.0]
    @test integration_point_coordinates[5, :] == [-0.5773502691896258, 0.7745966692414834]
    @test integration_point_coordinates[6, :] == [0.5773502691896258, 0.7745966692414834]

    @test integration_point_weights[1, :] == [1.0, 0.5555555555555556]
    @test integration_point_weights[2, :] == [1.0, 0.5555555555555556]
    @test integration_point_weights[3, :] == [1.0, 0.8888888888888888]
    @test integration_point_weights[4, :] == [1.0, 0.8888888888888888]
    @test integration_point_weights[5, :] == [1.0, 0.5555555555555556]
    @test integration_point_weights[6, :] == [1.0, 0.5555555555555556]

    weights, xi = get_weights_and_integration_points(dof, [num_int + 1, num_int])
    integration_point_coordinates = get_multi_dimensional_integration_point_data(dof, [num_int + 1, num_int], xi)
    integration_point_weights = get_multi_dimensional_integration_point_data(dof, [num_int + 1, num_int], weights)

    @test length(integration_point_coordinates[:, 1]) == 6
    @test length(integration_point_weights[:, 1]) == 6
    @test integration_point_coordinates[1, :] == [-0.7745966692414834, -0.5773502691896258]
    @test integration_point_coordinates[2, :] == [0.0, -0.5773502691896258]
    @test integration_point_coordinates[3, :] == [0.7745966692414834, -0.5773502691896258]
    @test integration_point_coordinates[4, :] == [-0.7745966692414834, 0.5773502691896258]
    @test integration_point_coordinates[5, :] == [0.0, 0.5773502691896258]
    @test integration_point_coordinates[6, :] == [0.7745966692414834, 0.5773502691896258]

    @test integration_point_weights[1, :] == [0.5555555555555556, 1.0]
    @test integration_point_weights[2, :] == [0.8888888888888888, 1.0]
    @test integration_point_weights[3, :] == [0.5555555555555556, 1.0]
    @test integration_point_weights[4, :] == [0.5555555555555556, 1.0]
    @test integration_point_weights[5, :] == [0.8888888888888888, 1.0]
    @test integration_point_weights[6, :] == [0.5555555555555556, 1.0]
end

@testset "ut_create_element_matrices" begin
    dof::Int64 = 1
    p::Vector{Int64} = [1, 1]

    N, B = create_element_matrices(dof, p, Lagrange_element.create_element_matrices)

    @test isnothing(N)
    @test isnothing(B)
    dof = 4
    N, B = create_element_matrices(dof, Vector{Int64}([1, 1, 1, 1]), Lagrange_element.create_element_matrices)
    @test isnothing(N)
    @test isnothing(B)
    dof = 2
    N, B = create_element_matrices(dof, p, Lagrange_element.create_element_matrices)

    @test size(N) == (4, 8, 2)
    @test size(B) == (4, 8, 3)

    @test N[1, 1:4, :] == [0.6220084679281462 0.0; 0.0 0.6220084679281462; 0.16666666666666663 0.0; 0.0 0.16666666666666663]
    @test N[2, 1:4, :] == [0.16666666666666663 0.0; 0.0 0.16666666666666663; 0.6220084679281462 0.0; 0.0 0.6220084679281462]
    @test N[3, 1:4, :] == [0.16666666666666663 0.0; 0.0 0.16666666666666663; 0.044658198738520435 0.0; 0.0 0.044658198738520435]
    @test N[4, 1:4, :] == [0.044658198738520435 0.0; 0.0 0.044658198738520435; 0.16666666666666663 0.0; 0.0 0.16666666666666663]

    @test N[1, 5:8, :] == [0.16666666666666663 0.0; 0.0 0.16666666666666663; 0.044658198738520435 0.0; 0.0 0.044658198738520435]
    @test N[2, 5:8, :] == [0.044658198738520435 0.0; 0.0 0.044658198738520435; 0.16666666666666663 0.0; 0.0 0.16666666666666663]
    @test N[3, 5:8, :] == [0.6220084679281462 0.0; 0.0 0.6220084679281462; 0.16666666666666663 0.0; 0.0 0.16666666666666663]
    @test N[4, 5:8, :] == [0.16666666666666663 0.0; 0.0 0.16666666666666663; 0.6220084679281462 0.0; 0.0 0.6220084679281462]

    @test B[3, 1:6, 1] == [-0.10566243270259354, 0.0, 0.10566243270259354, 0.0, -0.39433756729740643, 0.0]
    @test B[3, 1:6, 2] == [0.0, -0.39433756729740643, 0.0, -0.10566243270259354, 0.0, 0.39433756729740643]
    @test B[3, 1:6, 3] == [-0.39433756729740643, -0.10566243270259354, -0.10566243270259354, 0.10566243270259354, 0.39433756729740643, -0.39433756729740643]

    dof = 3
    p = [1, 1, 1]
    N, B = create_element_matrices(dof, p, Lagrange_element.create_element_matrices)

    @test size(N) == (8, 24, 3)
    @test size(B) == (8, 24, 6)
    @test N[1, 1:4, :] == [0.4905626121623441 0.0 0.0; 0.0 0.4905626121623441 0.0; 0.0 0.0 0.4905626121623441; 0.13144585576580212 0.0 0.0]
    @test N[4, 8:10, :] == [0.0 0.13144585576580212 0.0; 0.0 0.0 0.13144585576580212; 0.4905626121623441 0.0 0.0]
    @test B[7, 1:5, 1] == [-0.022329099369260218, 0.0, 0.0, 0.022329099369260218, 0.0]
    @test B[7, 1:5, 2] == [0.0, -0.08333333333333331, 0.0, 0.0, -0.022329099369260218]
    @test B[7, 1:5, 3] == [0.0, 0.0, -0.08333333333333331, 0.0, 0.0]
    @test B[7, 1:5, 4] == [0.0, -0.08333333333333331, -0.08333333333333331, 0.0, -0.022329099369260218]
    @test B[7, 1:5, 5] == [-0.08333333333333331, 0.0, -0.022329099369260218, -0.022329099369260218, 0.0]
    @test B[7, 1:5, 6] == [-0.08333333333333331, -0.022329099369260218, 0.0, -0.022329099369260218, 0.022329099369260218]

    dof = 2
    p = [2, 1]
    N, B = create_element_matrices(dof, p, Lagrange_element.create_element_matrices)

    @test size(N) == (6, 12, 2)
    @test size(B) == (6, 12, 3)
    dof = 3
    p = [5, 3, 1]
    N, B = create_element_matrices(dof, p, Lagrange_element.create_element_matrices)

    @test size(N) == (90, 144, 3)
    @test size(B) == (90, 144, 6)

end