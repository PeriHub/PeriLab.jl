# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
include("../../../src/FEM/FEM_Factory.jl")

using Test
include("../../../src/Support/data_manager.jl")


@testset "ut_valid_models" begin
    @test isnothing(FEM.valid_models(Dict()))
    @test isnothing(FEM.valid_models(Dict("Additive Model" => "a")))
    @test isnothing(FEM.valid_models(Dict("Damage Model" => "a")))
    @test isnothing(FEM.valid_models(Dict("Damage Model" => "a", "Additive Model" => "a", "Material Model" => "a Correspondence")))
    @test isnothing(FEM.valid_models(Dict("Thermal Model" => "a")))
    @test isnothing(FEM.valid_models(Dict("Material Model" => "a")))
end
@testset "ut_init_FEM" begin
    test_Data_manager = Data_manager
    nelements = 2
    test_Data_manager.set_dof(2)
    test_Data_manager.set_num_elements(nelements)
    test_Data_manager.set_num_controller(6)
    rho = test_Data_manager.create_constant_node_field("Density", Float64, 1)
    rho .= 2
    test = FEM.init_FEM(Dict(), test_Data_manager)
    @test isnothing(test)
    test_Data_manager.set_dof(1)
    test = FEM.init_FEM(Dict("Degree" => 1), test_Data_manager)
    @test isnothing(test)
    test = FEM.init_FEM(Dict("Degree" => 4), test_Data_manager)
    @test isnothing(test)
    dof = 2
    test_Data_manager.set_dof(dof)
    test_Data_manager.create_node_field("Displacements", Float64, dof)
    coordinates = test_Data_manager.create_constant_node_field("Coordinates", Float64, dof)
    coordinates[1, 1] = 0
    coordinates[1, 2] = 0
    coordinates[2, 1] = 1
    coordinates[2, 2] = 0
    coordinates[3, 1] = 0
    coordinates[3, 2] = 1
    coordinates[4, 1] = 1
    coordinates[4, 2] = 1
    coordinates[5, 1] = 2
    coordinates[5, 2] = 0
    coordinates[6, 1] = 2
    coordinates[6, 2] = 1

    params = Dict("FEM" => Dict("FE_1" => Dict("Degree" => 1, "Element Type" => "Lagrange", "Material Model" => "Elastic Model")),
        "Material Models" => Dict("Elastic Model" => Dict("Material Model" => "Correspondence Elastic", "Symmetry" => "isotropic plane strain", "Young's Modulus" => 2.5e+3, "Poisson's Ratio" => 0.33, "Shear Modulus" => 2.0e3)))

    topology = test_Data_manager.create_constant_free_size_field("FE Topology", Int64, (2, 4))
    topology[1, 1] = 1
    topology[1, 2] = 2
    topology[1, 3] = 3
    topology[1, 4] = 4
    topology[2, 1] = 2
    topology[2, 2] = 5
    topology[2, 3] = 4
    topology[2, 4] = 6


    params = Dict("FEM" => Dict("FE_1" => Dict("Degree" => 1, "Element Type" => "Lagrange", "Material Model" => "Elastic Model")),
        "Material Models" => Dict("Elastic Model" => Dict("Material Model" => "Correspondence Elastic", "Symmetry" => "isotropic plane strain", "Young's Modulus" => 2.5e+3, "Poisson's Ratio" => 0.33, "Shear Modulus" => 2.0e3)))
    test_Data_manager = FEM.init_FEM(params["FEM"]["FE_1"], test_Data_manager)

    @test "N Matrix" in test_Data_manager.get_all_field_keys()
    @test "B Matrix" in test_Data_manager.get_all_field_keys()
    N = test_Data_manager.get_field("N Matrix")
    @test size(N) == (4, 8, 2)
    @test N[1, 1:4, :] == [0.6220084679281462 0.0; 0.0 0.6220084679281462; 0.16666666666666663 0.0; 0.0 0.16666666666666663]
    B = test_Data_manager.get_field("B Matrix")
    @test size(B) == (4, 8, 3)
    @test B[3, 1:6, 3] == [-0.39433756729740643, -0.10566243270259354, -0.10566243270259354, 0.10566243270259354, 0.39433756729740643, -0.39433756729740643]

    @test "Element StrainN" in test_Data_manager.get_all_field_keys()
    @test "Element StressN" in test_Data_manager.get_all_field_keys()
    @test "Element StrainNP1" in test_Data_manager.get_all_field_keys()
    @test "Element StressNP1" in test_Data_manager.get_all_field_keys()
    @test "Element Strain Increment" in test_Data_manager.get_all_field_keys()
    @test "Element Jacobi Matrix" in test_Data_manager.get_all_field_keys()
    @test "Element Jacobi Determinant" in test_Data_manager.get_all_field_keys()
    @test "Lumped Mass Matrix" in test_Data_manager.get_all_field_keys()
    lumped_mass = test_Data_manager.get_field("Lumped Mass Matrix")

    @test lumped_mass[:, 1] == [0.49999999999999994, 0.9999999999999998, 0.49999999999999994, 0.9999999999999998, 0.4999999999999999, 0.49999999999999994]
    @test lumped_mass[:, 2] == [0.49999999999999994, 0.9999999999999998, 0.49999999999999994, 0.9999999999999998, 0.4999999999999999, 0.49999999999999994]
end

@testset "ut_eval" begin

    test_Data_manager = Data_manager
    dof = 2
    test_Data_manager.set_dof(dof)
    test_Data_manager.set_num_elements(2)
    test_Data_manager.set_num_controller(6)
    test_Data_manager.create_node_field("Force Density", Float64, dof)
    test_Data_manager.create_node_field("Displacements", Float64, dof)
    coordinates = test_Data_manager.create_constant_node_field("Coordinates", Float64, dof)
    coordinates[1, 1] = 0
    coordinates[1, 2] = 0
    coordinates[2, 1] = 1
    coordinates[2, 2] = 0
    coordinates[3, 1] = 0
    coordinates[3, 2] = 1
    coordinates[4, 1] = 1
    coordinates[4, 2] = 1
    coordinates[5, 1] = 2
    coordinates[5, 2] = 0
    coordinates[6, 1] = 2
    coordinates[6, 2] = 1

    params = Dict("FEM" => Dict("Degree" => 1, "Element Type" => "Lagrange", "Material Model" => "Elastic Model"), "Physics" => Dict("Material Models" => Dict("Elastic Model" => Dict("Material Model" => "Correspondence Elastic", "Symmetry" => "isotropic plane strain", "Young's Modulus" => 1.5, "Poisson's Ratio" => 0.33, "Shear Modulus" => 0.5639))))

    topology = test_Data_manager.create_constant_free_size_field("FE Topology", Int64, (2, 4))
    topology[1, 1] = 1
    topology[1, 2] = 2
    topology[1, 3] = 3
    topology[1, 4] = 4
    topology[2, 1] = 2
    topology[2, 2] = 5
    topology[2, 3] = 4
    topology[2, 4] = 6

    test_Data_manager = FEM.init_FEM(params["FEM"], test_Data_manager)
    elements = Vector{Int64}([1, 2])
    test_Data_manager = FEM.eval(test_Data_manager, elements, params, 0.0, 1.0e-6)

    stress = test_Data_manager.get_field("Element Stress", "NP1")
    strain = test_Data_manager.get_field("Element Strain", "NP1")
    for iEl in 1:2
        for i_int in 1:4
            for i in 1:3
                @test stress[iEl, i_int, i] == 0
                @test strain[iEl, i_int, i] == 0
            end
        end
    end
    displacements = test_Data_manager.get_field("Displacements", "NP1")

    displacements[1, 1] = -1
    displacements[1, 2] = 0.5
    displacements[2, 1] = 0
    displacements[2, 2] = 0.5
    displacements[3, 1] = -1
    displacements[3, 2] = 0.5
    displacements[4, 1] = 0
    displacements[4, 2] = 0.5
    displacements[5, 1] = 1
    displacements[5, 2] = 0.5
    displacements[6, 1] = 1
    displacements[6, 2] = 0.5
    test_Data_manager = FEM.eval(test_Data_manager, elements, params, 0.0, 1.0e-6)
    stress = test_Data_manager.get_field("Element Stress", "NP1")
    strain = test_Data_manager.get_field("Element Strain", "NP1")

    for iEl in 1:2
        for i_int in 1:4
            @test isapprox(strain[iEl, i_int, 1], 1)
            @test isapprox(strain[iEl, i_int, 2] + 1, 1)
            @test isapprox(strain[iEl, i_int, 3] + 1, 1)
            @test isapprox(stress[iEl, i_int, 1], 2.222467934542238)
            @test isapprox(stress[iEl, i_int, 2], 1.0946483856700575)
            @test isapprox(stress[iEl, i_int, 3] + 1, 1)
        end

    end

end
@testset "ut_get_FEM_nodes" begin
    test_Data_manager = Data_manager
    topology = test_Data_manager.get_field("FE Topology")
    test_Data_manager = FEM.get_FEM_nodes(test_Data_manager, topology)
    @test "FE Nodes" in test_Data_manager.get_all_field_keys()
    fem_nodes = test_Data_manager.get_field("FE Nodes")
    for i in eachindex(fem_nodes)
        @test fem_nodes[i]
    end
    fem_nodes[:] .= false
    topology[1, 1] = 1
    topology[1, 2] = 5
    topology[1, 3] = 3
    topology[1, 4] = 4
    topology[2, 1] = 1
    topology[2, 2] = 5
    topology[2, 3] = 3
    topology[2, 4] = 4
    test_Data_manager = FEM.get_FEM_nodes(test_Data_manager, topology)
    @test fem_nodes[1]
    @test !fem_nodes[2]
    @test fem_nodes[3]
    @test fem_nodes[4]
    @test fem_nodes[5]
    @test !fem_nodes[6]
end