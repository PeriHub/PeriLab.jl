# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
include("../../../src/FEM/FEM_routines.jl")
include("../../../src/FEM/Element_formulation/lagrange_element.jl")
using .Lagrange_element
# include("../../../src/PeriLab.jl")
using .FEM_routines:
    create_element_matrices,
    get_Jacobian,
    get_polynomial_degree,
    get_number_of_integration_points,
    create_element_matrices,
    get_lumped_mass,
    get_weights_and_integration_points,
    get_multi_dimensional_integration_point_data,
    get_FE_material_model
# using .PeriLab
using Test

@testset "ut_jacobi" begin
    test_data_manager = PeriLab.Data_manager
    test_data_manager.initialize_data()
    dof = 2
    nelements = 1
    test_data_manager.set_dof(dof)
    test_data_manager.set_num_elements(nelements)
    test_data_manager.set_num_controller(4)
    test_data_manager.create_node_field("Force Densities", Float64, dof)

    coordinates = test_data_manager.create_constant_node_field("Coordinates", Float64, dof)
    params = Dict(
        "FEM" => Dict(
            "FE_1" => Dict(
                "Degree" => 1,
                "Element Type" => "Lagrange",
                "Material Model" => "Elastic Model",
            ),
        ),
        "Material Models" => Dict(
            "Elastic Model" => Dict(
                "Material Model" => "Correspondence Elastic",
                "Symmetry" => "isotropic plane strain",
                "Young's Modulus" => 2.5e+3,
                "Poisson's Ratio" => 0.33,
                "Shear Modulus" => 2.0e3,
            ),
        ),
    )

    topology =
        test_data_manager.create_constant_free_size_field("FE Topology", Int64, (2, 4))
    topology[1, 1] = 1
    topology[1, 2] = 2
    topology[1, 3] = 3
    topology[1, 4] = 4
    elements = Vector{Int64}(1:nelements)
    p = get_polynomial_degree(params["FEM"]["FE_1"], dof)
    num_int = get_number_of_integration_points(p, dof)

    N = test_data_manager.create_constant_free_size_field(
        "N Matrix",
        Float64,
        (prod(num_int), prod(p .+ 1) * dof, dof),
    )
    B = test_data_manager.create_constant_free_size_field(
        "B Matrix",
        Float64,
        (prod(num_int), prod(p .+ 1) * dof, 3 * dof - 3),
    )

    N, B = create_element_matrices(dof, p, Lagrange_element.create_element_matrices)

    jacobian = test_data_manager.create_constant_free_size_field(
        "Element Jacobi Matrix",
        Float64,
        (nelements, prod(num_int), dof, dof),
    )
    determinant_jacobian = test_data_manager.create_constant_free_size_field(
        "Element Jacobi Determinant",
        Float64,
        (nelements, prod(num_int)),
    )
    test_jacobian, test_determinant_jacobian = get_Jacobian(
        elements,
        dof,
        topology,
        coordinates,
        B,
        jacobian,
        determinant_jacobian,
    )
    @test isnothing(test_jacobian)
    @test isnothing(test_determinant_jacobian)
    coordinates[1, 1] = 0
    coordinates[1, 2] = 0
    coordinates[2, 1] = 1
    coordinates[2, 2] = 0
    coordinates[3, 1] = 0
    coordinates[3, 2] = 1
    coordinates[4, 1] = 1
    coordinates[4, 2] = 1

    jacobian, determinant_jacobian = get_Jacobian(
        elements,
        dof,
        topology,
        coordinates,
        B,
        jacobian,
        determinant_jacobian,
    )

    for i = 1:4
        @test determinant_jacobian[1, i] == 0.25
        @test jacobian[1, i, :, :] == [2.0 0.0; 0.0 2.0]
    end

    coordinates[1, 1] = 0
    coordinates[1, 2] = 0
    coordinates[2, 1] = 2
    coordinates[2, 2] = 0
    coordinates[3, 1] = 0
    coordinates[3, 2] = 0.5
    coordinates[4, 1] = 2
    coordinates[4, 2] = 0.5
    jacobian, determinant_jacobian = get_Jacobian(
        elements,
        dof,
        topology,
        coordinates,
        B,
        jacobian,
        determinant_jacobian,
    )

    for i = 1:4
        @test determinant_jacobian[1, i] == 0.25
        @test jacobian[1, i, :, :] == [1.0 0.0; 0.0 4.0]
    end

    coordinates[1, 1] = 0
    coordinates[1, 2] = 0
    coordinates[2, 1] = 1.0
    coordinates[2, 2] = 0.0
    coordinates[3, 1] = 0.5
    coordinates[3, 2] = 1
    coordinates[4, 1] = 1.5
    coordinates[4, 2] = 1.0


    jacobian, determinant_jacobian = get_Jacobian(
        elements,
        dof,
        topology,
        coordinates,
        B,
        jacobian,
        determinant_jacobian,
    )

    for i = 1:4
        @test isapprox(determinant_jacobian[1, i], 0.25)
        @test isapprox(jacobian[1, i, 1, 1], 2.0)
        @test isapprox(jacobian[1, i, 1, 2], -1.0)
        @test isapprox(jacobian[1, i, 2, 1], 0.0)
        @test isapprox(jacobian[1, i, 2, 2], 2)
    end

    coordinates[1, 1] = 0
    coordinates[1, 2] = 0
    coordinates[2, 1] = 1.0
    coordinates[2, 2] = 0.5
    coordinates[3, 1] = 0.0
    coordinates[3, 2] = 1
    coordinates[4, 1] = 1.0
    coordinates[4, 2] = 1.5

    jacobian, determinant_jacobian = get_Jacobian(
        elements,
        dof,
        topology,
        coordinates,
        B,
        jacobian,
        determinant_jacobian,
    )
    for i = 1:4
        @test isapprox(determinant_jacobian[1, i], 0.25)

        @test isapprox(jacobian[1, i, 1, 1], 2.0)
        @test isapprox(jacobian[1, i, 1, 2], 0.0)
        @test isapprox(jacobian[1, i, 2, 1], -1.0)
        @test isapprox(jacobian[1, i, 2, 2], 2)
    end
end

@testset "ut_lumped_mass" begin
    test_data_manager = PeriLab.Data_manager
    dof = 2
    nelements = 1
    test_data_manager.set_dof(dof)
    test_data_manager.set_num_elements(nelements)
    test_data_manager.set_num_controller(4)

    coordinates = test_data_manager.create_constant_node_field("Coordinates", Float64, dof)
    params = Dict(
        "FEM" => Dict(
            "FE_1" => Dict(
                "Degree" => 1,
                "Element Type" => "Lagrange",
                "Material Model" => "Elastic Model",
            ),
        ),
        "Material Models" => Dict(
            "Elastic Model" => Dict(
                "Material Model" => "Correspondence Elastic",
                "Symmetry" => "isotropic plane strain",
                "Young's Modulus" => 2.5e+3,
                "Poisson's Ratio" => 0.33,
                "Shear Modulus" => 2.0e3,
            ),
        ),
    )

    topology =
        test_data_manager.create_constant_free_size_field("FE Topology", Int64, (2, 4))
    topology[1, 1] = 1
    topology[1, 2] = 2
    topology[1, 3] = 3
    topology[1, 4] = 4
    elements = Vector{Int64}(1:nelements)
    p = get_polynomial_degree(params["FEM"]["FE_1"], dof)
    num_int = get_number_of_integration_points(p, dof)

    N = test_data_manager.create_constant_free_size_field(
        "N Matrix",
        Float64,
        (prod(num_int), prod(p .+ 1) * dof, dof),
    )
    B = test_data_manager.create_constant_free_size_field(
        "B Matrix",
        Float64,
        (prod(num_int), prod(p .+ 1) * dof, 3 * dof - 3),
    )

    N[:], B[:] = create_element_matrices(dof, p, Lagrange_element.create_element_matrices)
    lumped_mass =
        test_data_manager.create_constant_node_field("Lumped Mass Matrix", Float64, 1)

    matrix = zeros(3, 3)
    nu = 0.3
    E = 1000
    G = E / (2 * (1 + nu))
    temp = E / ((1 + nu) * (1 - 2 * nu))
    matrix[1, 1] = (1 - nu) * temp
    matrix[2, 2] = (1 - nu) * temp
    matrix[3, 3] = G
    matrix[1, 2] = nu * temp
    matrix[2, 1] = nu * temp

    matrix[1, 1] = E / (1 - nu * nu)
    matrix[1, 2] = E * nu / (1 - nu * nu)
    matrix[2, 1] = E * nu / (1 - nu * nu)
    matrix[2, 2] = E / (1 - nu * nu)
    matrix[3, 3] = G

    global K = zeros(8, 8)
    for i = 1:4
        K[:, :] += B[i, :, :] * matrix * B[i, :, :]'
    end

    Ktest = zeros(8, 8)
    # checked with separate code; without jacobian
    Ktest[1, :] = [
        494.50549450549454,
        178.57142857142853,
        -302.1978021978022,
        -13.73626373626374,
        54.94505494505495,
        13.736263736263744,
        -247.2527472527472,
        -178.57142857142853,
    ]
    Ktest[2, :] = [
        178.57142857142853,
        494.5054945054945,
        13.736263736263746,
        54.94505494505495,
        -13.736263736263744,
        -302.1978021978022,
        -178.57142857142853,
        -247.2527472527472,
    ]
    Ktest[3, :] = [
        -302.1978021978022,
        13.73626373626374,
        494.50549450549454,
        -178.57142857142853,
        -247.2527472527472,
        178.5714285714285,
        54.94505494505495,
        -13.736263736263744,
    ]
    Ktest[4, :] = [
        -13.736263736263746,
        54.94505494505495,
        -178.57142857142853,
        494.50549450549454,
        178.57142857142853,
        -247.2527472527472,
        13.736263736263744,
        -302.1978021978022,
    ]
    Ktest[5, :] = [
        54.94505494505495,
        -13.736263736263746,
        -247.2527472527472,
        178.5714285714285,
        494.50549450549454,
        -178.57142857142853,
        -302.1978021978022,
        13.736263736263737,
    ]
    Ktest[6, :] = [
        13.73626373626374,
        -302.1978021978023,
        178.57142857142853,
        -247.2527472527472,
        -178.57142857142853,
        494.5054945054945,
        -13.736263736263744,
        54.94505494505495,
    ]
    Ktest[7, :] = [
        -247.2527472527472,
        -178.57142857142853,
        54.94505494505495,
        13.736263736263744,
        -302.1978021978022,
        -13.736263736263744,
        494.50549450549454,
        178.57142857142853,
    ]
    Ktest[8, :] = [
        -178.57142857142853,
        -247.2527472527472,
        -13.736263736263744,
        -302.1978021978022,
        13.736263736263744,
        54.94505494505495,
        178.57142857142853,
        494.5054945054945,
    ]
    @test isapprox(Ktest, K)


    coordinates[1, 1] = 0
    coordinates[1, 2] = 0
    coordinates[2, 1] = 1
    coordinates[2, 2] = 0
    coordinates[3, 1] = 0
    coordinates[3, 2] = 1
    coordinates[4, 1] = 1
    coordinates[4, 2] = 1
    jacobian = test_data_manager.create_constant_free_size_field(
        "Element Jacobi Matrix",
        Float64,
        (nelements, prod(num_int), dof, dof),
    )
    determinant_jacobian = test_data_manager.create_constant_free_size_field(
        "Element Jacobi Determinant",
        Float64,
        (nelements, prod(num_int)),
    )
    jacobian, determinant_jacobian = get_Jacobian(
        elements,
        dof,
        topology,
        coordinates,
        B,
        jacobian,
        determinant_jacobian,
    )
    N = test_data_manager.get_field("N Matrix")
    rho = test_data_manager.create_constant_node_field("Density", Float64, 1)

    rho[:] .= 1.0
    lumped_mass =
        get_lumped_mass(elements, dof, topology, N, determinant_jacobian, rho, lumped_mass)
    for i = 1:4
        @test isapprox(lumped_mass[i], 0.25)
    end
    lumped_mass[:] .= 0
    rho[:] .= 2.0

    lumped_mass =
        get_lumped_mass(elements, dof, topology, N, determinant_jacobian, rho, lumped_mass)
    for i = 1:4
        @test isapprox(lumped_mass[i], 0.5)
    end
    lumped_mass[:] .= 0
    rho[:] .= 1.2

    lumped_mass =
        get_lumped_mass(elements, dof, topology, N, determinant_jacobian, rho, lumped_mass)
    @test lumped_mass == [0.3; 0.3; 0.3; 0.3]
end

@testset "ut_get_FE_material_model" begin
    params = Dict(
        "FEM" => Dict(
            "FE_1" => Dict(
                "Degree" => 1,
                "Element Type" => "Lagrange",
                "Material Model" => "Elastic Model",
            ),
        ),
        "Material Models" => Dict(
            "Elastic Model 2" => Dict(
                "Material Model" => "Correspondence Elastic",
                "Symmetry" => "isotropic plane strain",
                "Bulk Modulus" => 2.5e+3,
                "Shear Modulus" => 1.15e3,
            ),
        ),
    )

    @test isnothing(get_FE_material_model(params, "FE_1"))

    params = Dict(
        "FEM" => Dict(
            "FE_1" => Dict(
                "Degree" => 1,
                "Element Type" => "Lagrange",
                "Material Model" => "Elastic Model",
            ),
        ),
        "Material Models" => Dict(
            "Elastic Model" => Dict(
                "Material Model" => "Correspondence Elastic",
                "Symmetry" => "isotropic plane strain",
                "Bulk Modulus" => 2.5e+3,
                "Shear Modulus" => 1.15e3,
            ),
        ),
    )

    @test get_FE_material_model(params, "FE_1") == Dict(
        "Material Model" => "Correspondence Elastic",
        "Symmetry" => "isotropic plane strain",
        "Bulk Modulus" => 2.5e+3,
        "Shear Modulus" => 1.15e3,
    )
end

@testset "ut_get_polynomial_degree" begin

    @test isnothing(get_polynomial_degree(Dict(), 1))
    @test isnothing(get_polynomial_degree(Dict(), 2))
    @test isnothing(get_polynomial_degree(Dict(), 3))

    params = Dict("Degree" => 1)

    @test get_polynomial_degree(params, 2) == [1, 1]
    @test get_polynomial_degree(params, 3) == [1, 1, 1]

    params = Dict("Degree" => 2)
    @test get_polynomial_degree(params, 2) == [2, 2]
    @test get_polynomial_degree(params, 3) == [2, 2, 2]

    params = Dict("Degree" => 2.1)

    @test get_polynomial_degree(params, 2) == [2, 2]
    @test get_polynomial_degree(params, 3) == [2, 2, 2]

    params = Dict("Degree" => [2 3 1])
    @test isnothing(get_polynomial_degree(params, 2))
    @test get_polynomial_degree(params, 3) == [2, 3, 1]

    params = Dict("Degree" => [2.1 2])
    @test get_polynomial_degree(params, 2) == [2, 2]
    @test isnothing(get_polynomial_degree(params, 3))
    params = Dict("Degree" => "2")
    @test get_polynomial_degree(params, 3) == [2, 2, 2]
    params = Dict("Degree" => "2 2")
    @test get_polynomial_degree(params, 2) == [2, 2]
    params = Dict("Degree" => "2 1")
    @test get_polynomial_degree(params, 2) == [2, 1]
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
    @test get_weights_and_integration_points(2, [2, 2]) == (
        [1.0 1.0; 1.0 1.0],
        [-0.5773502691896258 -0.5773502691896258; 0.5773502691896258 0.5773502691896258],
    )
    @test get_weights_and_integration_points(3, [2, 3, 4]) == (
        [
            1.0 0.5555555555555556 0.34785484513745385
            1.0 0.8888888888888888 0.6521451548625462
            0.0 0.5555555555555556 0.6521451548625462
            0.0 0.0 0.34785484513745385
        ],
        [
            -0.5773502691896258 -0.7745966692414834 -0.8611363115940526
            0.5773502691896258 0.0 -0.3399810435848563
            0.0 0.7745966692414834 0.3399810435848563
            0.0 0.0 0.8611363115940526
        ],
    )
    @test get_weights_and_integration_points(3, [2, 1, 1]) == (
        [1.0 2.0 2.0; 1.0 0.0 0.0],
        [-0.5773502691896258 0.0 0.0; 0.5773502691896258 0.0 0.0],
    )
end

@testset "ut_get_multi_dimensional_integration_points" begin
    @test isnothing(get_multi_dimensional_integration_point_data(1, [1], zeros(2, 2)))
    @test isnothing(
        get_multi_dimensional_integration_point_data(4, [1, 1, 1, 1], zeros(2, 2)),
    )
    dof = 2
    num_int = 2
    weights, xi = get_weights_and_integration_points(dof, [num_int, num_int])
    integration_point_coordinates =
        get_multi_dimensional_integration_point_data(dof, [num_int, num_int], xi)
    @test length(integration_point_coordinates[:, 1]) == 4
    @test integration_point_coordinates[1, :] == [-0.5773502691896258, -0.5773502691896258]
    @test integration_point_coordinates[2, :] == [0.5773502691896258, -0.5773502691896258]
    @test integration_point_coordinates[3, :] == [-0.5773502691896258, 0.5773502691896258]
    @test integration_point_coordinates[4, :] == [0.5773502691896258, 0.5773502691896258]
    @test get_multi_dimensional_integration_point_data(dof, [num_int, num_int], weights) ==
          [1.0 1.0; 1.0 1.0; 1.0 1.0; 1.0 1.0]
    dof = 3
    weights, xi = get_weights_and_integration_points(dof, [num_int, num_int, num_int])
    integration_point_coordinates =
        get_multi_dimensional_integration_point_data(dof, [num_int, num_int, num_int], xi)

    @test length(integration_point_coordinates[:, 1]) == 8
    @test integration_point_coordinates[1, :] ==
          [-0.5773502691896258, -0.5773502691896258, -0.5773502691896258]
    @test integration_point_coordinates[2, :] ==
          [0.5773502691896258, -0.5773502691896258, -0.5773502691896258]
    @test integration_point_coordinates[3, :] ==
          [-0.5773502691896258, 0.5773502691896258, -0.5773502691896258]
    @test integration_point_coordinates[4, :] ==
          [0.5773502691896258, 0.5773502691896258, -0.5773502691896258]
    @test integration_point_coordinates[5, :] ==
          [-0.5773502691896258, -0.5773502691896258, 0.5773502691896258]
    @test integration_point_coordinates[6, :] ==
          [0.5773502691896258, -0.5773502691896258, 0.5773502691896258]
    @test integration_point_coordinates[7, :] ==
          [-0.5773502691896258, 0.5773502691896258, 0.5773502691896258]
    @test integration_point_coordinates[8, :] ==
          [0.5773502691896258, 0.5773502691896258, 0.5773502691896258]
    integration_point_weights = get_multi_dimensional_integration_point_data(
        dof,
        [num_int, num_int, num_int],
        weights,
    )

    for i = 1:8
        @test integration_point_weights[i, :] == [1.0, 1.0, 1.0]
    end

    dof = 2
    weights, xi = get_weights_and_integration_points(dof, [num_int, num_int + 1])
    integration_point_coordinates =
        get_multi_dimensional_integration_point_data(dof, [num_int, num_int + 1], xi)
    integration_point_weights =
        get_multi_dimensional_integration_point_data(dof, [num_int, num_int + 1], weights)

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
    integration_point_coordinates =
        get_multi_dimensional_integration_point_data(dof, [num_int + 1, num_int], xi)
    integration_point_weights =
        get_multi_dimensional_integration_point_data(dof, [num_int + 1, num_int], weights)

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
    N, B = create_element_matrices(
        dof,
        Vector{Int64}([1, 1, 1, 1]),
        Lagrange_element.create_element_matrices,
    )
    @test isnothing(N)
    @test isnothing(B)
    dof = 2
    N, B = create_element_matrices(dof, p, Lagrange_element.create_element_matrices)

    @test size(N) == (4, 8, 2)
    @test size(B) == (4, 8, 3)

    @test N[1, 1:4, :] == [
        0.6220084679281462 0.0
        0.0 0.6220084679281462
        0.16666666666666663 0.0
        0.0 0.16666666666666663
    ]
    @test N[2, 1:4, :] == [
        0.16666666666666663 0.0
        0.0 0.16666666666666663
        0.6220084679281462 0.0
        0.0 0.6220084679281462
    ]
    @test N[3, 1:4, :] == [
        0.16666666666666663 0.0
        0.0 0.16666666666666663
        0.044658198738520435 0.0
        0.0 0.044658198738520435
    ]
    @test N[4, 1:4, :] == [
        0.044658198738520435 0.0
        0.0 0.044658198738520435
        0.16666666666666663 0.0
        0.0 0.16666666666666663
    ]

    @test N[1, 5:8, :] == [
        0.16666666666666663 0.0
        0.0 0.16666666666666663
        0.044658198738520435 0.0
        0.0 0.044658198738520435
    ]
    @test N[2, 5:8, :] == [
        0.044658198738520435 0.0
        0.0 0.044658198738520435
        0.16666666666666663 0.0
        0.0 0.16666666666666663
    ]
    @test N[3, 5:8, :] == [
        0.6220084679281462 0.0
        0.0 0.6220084679281462
        0.16666666666666663 0.0
        0.0 0.16666666666666663
    ]
    @test N[4, 5:8, :] == [
        0.16666666666666663 0.0
        0.0 0.16666666666666663
        0.6220084679281462 0.0
        0.0 0.6220084679281462
    ]

    @test B[3, 1:6, 1] ==
          [-0.10566243270259354, 0.0, 0.10566243270259354, 0.0, -0.39433756729740643, 0.0]
    @test B[3, 1:6, 2] ==
          [0.0, -0.39433756729740643, 0.0, -0.10566243270259354, 0.0, 0.39433756729740643]
    @test B[3, 1:6, 3] == [
        -0.39433756729740643,
        -0.10566243270259354,
        -0.10566243270259354,
        0.10566243270259354,
        0.39433756729740643,
        -0.39433756729740643,
    ]

    dof = 3
    p = [1, 1, 1]
    N, B = create_element_matrices(dof, p, Lagrange_element.create_element_matrices)

    @test size(N) == (8, 24, 3)
    @test size(B) == (8, 24, 6)
    @test N[1, 1:4, :] == [
        0.4905626121623441 0.0 0.0
        0.0 0.4905626121623441 0.0
        0.0 0.0 0.4905626121623441
        0.13144585576580212 0.0 0.0
    ]
    @test N[4, 8:10, :] == [
        0.0 0.13144585576580212 0.0
        0.0 0.0 0.13144585576580212
        0.4905626121623441 0.0 0.0
    ]
    @test B[7, 1:5, 1] == [-0.022329099369260218, 0.0, 0.0, 0.022329099369260218, 0.0]
    @test B[7, 1:5, 2] == [0.0, -0.08333333333333331, 0.0, 0.0, -0.022329099369260218]
    @test B[7, 1:5, 3] == [0.0, 0.0, -0.08333333333333331, 0.0, 0.0]
    @test B[7, 1:5, 4] ==
          [0.0, -0.08333333333333331, -0.08333333333333331, 0.0, -0.022329099369260218]
    @test B[7, 1:5, 5] ==
          [-0.08333333333333331, 0.0, -0.022329099369260218, -0.022329099369260218, 0.0]
    @test B[7, 1:5, 6] == [
        -0.08333333333333331,
        -0.022329099369260218,
        0.0,
        -0.022329099369260218,
        0.022329099369260218,
    ]

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
