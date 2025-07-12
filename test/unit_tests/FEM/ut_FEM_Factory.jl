# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
#include("../../../src/PeriLab.jl")

#using .PeriLab
PeriLab.Data_manager.initialize_data()

@testset "ut_valid_models" begin
    @test isnothing(PeriLab.Solver_control.FEM.valid_models(Dict{String,Any}()))
    @test isnothing(PeriLab.Solver_control.FEM.valid_models(Dict{String,Any}("Additive Model" => "a")))
    @test isnothing(PeriLab.Solver_control.FEM.valid_models(Dict{String,Any}("Damage Model" => "a")))
    @test PeriLab.Solver_control.FEM.valid_models(Dict{String,Any}("Damage Model" => "a",
                                                                   "Additive Model" => "a",
                                                                   "Material Model" => "a Correspondence"))
    @test isnothing(PeriLab.Solver_control.FEM.valid_models(Dict{String,Any}("Thermal Model" => "a")))
    @test PeriLab.Solver_control.FEM.valid_models(Dict{String,Any}("Material Model" => "a"))
end
@testset "ut_init_FEM" begin
    test_data_manager = PeriLab.Data_manager
    nelements = 2
    test_data_manager.set_dof(2)
    test_data_manager.set_num_elements(nelements)
    test_data_manager.set_num_controller(6)
    rho = test_data_manager.create_constant_node_field("Density", Float64, 1)
    rho .= 2
    test = PeriLab.Solver_control.FEM.init_FEM(Dict{String,Any}(), test_data_manager)
    @test isnothing(test)
    test = PeriLab.Solver_control.FEM.init_FEM(Dict{String,Any}("FEM" => Dict()),
                                               test_data_manager)
    @test isnothing(test)
    test_data_manager.set_dof(1)
    test = PeriLab.Solver_control.FEM.init_FEM(Dict{String,Any}("FEM" => Dict("Degree" => 1)),
                                               test_data_manager)
    @test isnothing(test)
    test = PeriLab.Solver_control.FEM.init_FEM(Dict{String,Any}("FEM" => Dict("Degree" => 4)),
                                               test_data_manager)
    @test isnothing(test)
    dof = 2
    test_data_manager.set_dof(dof)
    test_data_manager.create_node_field("Displacements", Float64, dof)
    coordinates = test_data_manager.create_constant_node_field("Coordinates", Float64, dof)
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

    topology = test_data_manager.create_constant_free_size_field("FE Topology", Int64,
                                                                 (2, 4))
    topology[1, 1] = 1
    topology[1, 2] = 2
    topology[1, 3] = 3
    topology[1, 4] = 4
    topology[2, 1] = 2
    topology[2, 2] = 5
    topology[2, 3] = 4
    topology[2, 4] = 6
    test_data_manager.data["block_id_list"] = [1, 2]
    test_data_manager.init_properties()

    params = Dict{String,Any}("FEM" => Dict("Degree" => 1,
                                            "Element Type" => "Lagrange",
                                            "Material Model" => "Elastic Model"),
                              "Models" => Dict("Material Models" => Dict("Not FEM name" => Dict("Material Model" => "Correspondence Elastic",
                                                                                                "Symmetry" => "isotropic plane strain",
                                                                                                "Young's Modulus" => 2.5e+3,
                                                                                                "Poisson's Ratio" => 0.33,
                                                                                                "Shear Modulus" => 2.0e3))))
    @test isnothing(PeriLab.Solver_control.FEM.init_FEM(params, test_data_manager))
    params = Dict{String,Any}("FEM" => Dict("Degree" => 1,
                                            "Element Type" => "Lagrange",
                                            "Material Model" => "Elastic Model"),
                              "Models" => Dict{String,Any}("Material Models" => Dict("Elastic Model" => Dict("Material Model" => "Correspondence Elastic",
                                                                                                             "Symmetry" => "isotropic plane strain",
                                                                                                             "Young's Modulus" => 2.5e+3,
                                                                                                             "Poisson's Ratio" => 0.33,
                                                                                                             "Shear Modulus" => 2.0e3))))

    test_data_manager = PeriLab.Solver_control.FEM.init_FEM(params, test_data_manager)

    @test "N Matrix" in test_data_manager.get_all_field_keys()
    @test "B Matrix" in test_data_manager.get_all_field_keys()
    N = test_data_manager.get_field("N Matrix")
    @test size(N) == (4, 8, 2)
    @test N[1, 1:4, :] == [0.6220084679281462 0.0
           0.0 0.6220084679281462
           0.16666666666666663 0.0
           0.0 0.16666666666666663]
    B = test_data_manager.get_field("B Matrix")

    @test size(B) == (2, 4, 8, 3)
    @test B[1, 3, 1:6, 3] == [
        -0.7886751345948129,
        -0.21132486540518708,
        -0.21132486540518708,
        0.21132486540518708,
        0.7886751345948129,
        -0.7886751345948129
    ]
    @test "Element StrainN" in test_data_manager.get_all_field_keys()
    @test "Element StressN" in test_data_manager.get_all_field_keys()
    @test "Element StrainNP1" in test_data_manager.get_all_field_keys()
    @test "Element StressNP1" in test_data_manager.get_all_field_keys()
    @test "Element Strain Increment" in test_data_manager.get_all_field_keys()
    @test "Element Jacobi Matrix" in test_data_manager.get_all_field_keys()
    @test "Element Jacobi Determinant" in test_data_manager.get_all_field_keys()
    @test "Lumped Mass Matrix" in test_data_manager.get_all_field_keys()
    lumped_mass = test_data_manager.get_field("Lumped Mass Matrix")

    @test isapprox(lumped_mass[:], [2, 4, 2, 4, 2, 2])

    # only in tests for resize or redefinition reasons
    test_data_manager.fields[Int64]["FE Topology"] = zeros(Int64, 1, 6)
    @test isnothing(PeriLab.Solver_control.FEM.init_FEM(params, test_data_manager))
end

@testset "ut_eval" begin
    test_data_manager = PeriLab.Data_manager
    dof = 2
    test_data_manager.set_dof(dof)
    test_data_manager.set_num_elements(2)
    test_data_manager.set_num_controller(6)
    test_data_manager.create_node_field("Cauchy Stress", Float64, dof,
                                        VectorOrMatrix = "Matrix")
    test_data_manager.create_node_field("Force Densities", Float64, dof)
    test_data_manager.create_node_field("Displacements", Float64, dof)
    test_data_manager.create_node_field("Forces", Float64, dof)
    coordinates = test_data_manager.create_constant_node_field("Coordinates", Float64, dof)
    # only in tests for resize or redefinition reasons
    test_data_manager.fields[Int64]["FE Topology"] = zeros(Int64, 2, 4)
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

    topology = test_data_manager.create_constant_free_size_field("FE Topology", Int64,
                                                                 (2, 4))
    topology[1, 1] = 1
    topology[1, 2] = 2
    topology[1, 3] = 3
    topology[1, 4] = 4
    topology[2, 1] = 2
    topology[2, 2] = 5
    topology[2, 3] = 4
    topology[2, 4] = 6
    test_data_manager.data["block_id_list"] = [1, 2]
    test_data_manager.init_properties()

    params = Dict{String,Any}("FEM" => Dict("Degree" => 1,
                                            "Element Type" => "Lagrange",
                                            "Material Model" => "Elastic Model"),
                              "Models" => Dict("Material Models" => Dict("Elastic Model" => Dict("Material Model" => "Correspondence Elastic",
                                                                                                 "Symmetry" => "isotropic plane strain",
                                                                                                 "Young's Modulus" => 1.5,
                                                                                                 "Poisson's Ratio" => 0.33,
                                                                                                 "Shear Modulus" => 0.5639))))

    test_data_manager = PeriLab.Solver_control.FEM.init_FEM(params, test_data_manager)
    elements = Vector{Int64}([1, 2])
    test_data_manager = PeriLab.Solver_control.FEM.eval_FEM(test_data_manager,
                                                            elements,
                                                            test_data_manager.get_properties(1,
                                                                                             "FEM"),
                                                            0.0,
                                                            1.0e-6)

    stress = test_data_manager.get_field("Element Stress", "NP1")
    strain = test_data_manager.get_field("Element Strain", "NP1")
    for iEl in 1:2
        for i_int in 1:4
            for i in 1:3
                @test stress[iEl, i_int, i] == 0
                @test strain[iEl, i_int, i] == 0
            end
        end
    end
    displacements = test_data_manager.get_field("Displacements", "NP1")

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
    test_data_manager = PeriLab.Solver_control.FEM.eval_FEM(test_data_manager,
                                                            elements,
                                                            test_data_manager.get_properties(1,
                                                                                             "FEM"),
                                                            0.0,
                                                            1.0e-6)
    stress = test_data_manager.get_field("Element Stress", "NP1")
    strain = test_data_manager.get_field("Element Strain", "NP1")

    for iEl in 1:2
        for i_int in 1:4
            @test isapprox(strain[iEl, i_int, 1], 1)
            @test isapprox(strain[iEl, i_int, 2] + 1, 1)
            @test isapprox(strain[iEl, i_int, 3] + 1, 1)
            @test isapprox(stress[iEl, i_int, 1], 2.2224294117647054)
            @test isapprox(stress[iEl, i_int, 2], 1.0946294117647055)
            @test isapprox(stress[iEl, i_int, 3] + 1, 1)
        end
    end
end
@testset "ut_get_FEM_nodes" begin
    test_data_manager = PeriLab.Data_manager
    topology = test_data_manager.get_field("FE Topology")
    test_data_manager = PeriLab.Solver_control.FEM.get_FEM_nodes(test_data_manager,
                                                                 topology)
    @test "FE Nodes" in test_data_manager.get_all_field_keys()
    fem_nodes = test_data_manager.get_field("FE Nodes")
    for i in eachindex(fem_nodes)
        @test fem_nodes[i]
    end
    fem_nodes .= false
    topology[1, 1] = 1
    topology[1, 2] = 5
    topology[1, 3] = 3
    topology[1, 4] = 4
    topology[2, 1] = 1
    topology[2, 2] = 5
    topology[2, 3] = 3
    topology[2, 4] = 4
    test_data_manager = PeriLab.Solver_control.FEM.get_FEM_nodes(test_data_manager,
                                                                 topology)
    @test fem_nodes[1]
    @test !fem_nodes[2]
    @test fem_nodes[3]
    @test fem_nodes[4]
    @test fem_nodes[5]
    @test !fem_nodes[6]
end
