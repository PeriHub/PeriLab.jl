# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

@testset "init_influence_function" begin
    @testset "no 'Influence Function' key -> no-op" begin
        nodes = 2
        dof = 2
        Data_Manager.initialize_data()
        Data_Manager.set_num_controller(nodes)
        Data_Manager.set_dof(dof)
        nn = Data_Manager.create_constant_node_scalar_field("Number of Neighbors", Int64)
        nn[1] = 2
        nn[2] = 3
        bg = Data_Manager.create_constant_bond_vector_state("Bond Geometry", Float64,
                                                            dof + 1)
        omega = Data_Manager.create_constant_bond_scalar_state("Influence Function",
                                                               Float64)

        # should simply return without throwing, omega stays untouched
        PeriLab.Solver_Manager.Influence_Function.init_influence_function(Vector{Int64}(1:nodes),
                                                                          Dict())
        @test all(omega[1] .== 0.0)
        @test all(omega[2] .== 0.0)
    end
    @testset "predefined: 1/xi^2 (dof=2)" begin
        nodes = 2
        dof = 2
        Data_Manager.initialize_data()
        Data_Manager.set_num_controller(nodes)
        Data_Manager.set_dof(dof)
        nn = Data_Manager.create_constant_node_scalar_field("Number of Neighbors", Int64)
        nn[1] = 2
        nn[2] = 1
        bg = Data_Manager.create_constant_bond_vector_state("Bond Geometry", Float64,
                                                            dof + 1)
        omega = Data_Manager.create_constant_bond_scalar_state("Influence Function",
                                                               Float64)

        # node 1 has 2 bonds, node 2 has 1 bond
        # entries: [xiX, xiY, xi]
        bg[1][1] .= [1.0, 0.0, 1.0]
        bg[1][2] .= [0.0, 2.0, 2.0]
        bg[2][1] .= [3.0, 4.0, 5.0]

        params = Dict("Influence Function" => "1/xi^2")
        PeriLab.Solver_Manager.Influence_Function.init_influence_function(Vector{Int64}(1:nodes),
                                                                          params)

        @test omega[1][1] ≈ 1.0 / (1.0^2)
        @test omega[1][2] ≈ 1.0 / (2.0^2)
        @test omega[2][1] ≈ 1.0 / (5.0^2)
    end

    @testset "predefined: constant '1'" begin
        nodes = 1
        dof = 2
        Data_Manager.initialize_data()
        Data_Manager.set_num_controller(nodes)
        Data_Manager.set_dof(dof)
        nn = Data_Manager.create_constant_node_scalar_field("Number of Neighbors", Int64)
        nn[1] = 2
        bg = Data_Manager.create_constant_bond_vector_state("Bond Geometry", Float64,
                                                            dof + 1)
        omega = Data_Manager.create_constant_bond_scalar_state("Influence Function",
                                                               Float64)

        bg[1][1] .= [1.0, 1.0, 1.41421356]
        bg[1][2] .= [2.0, 0.0, 2.0]

        params = Dict("Influence Function" => "1")
        PeriLab.Solver_Manager.Influence_Function.init_influence_function(Vector{Int64}(1:nodes),
                                                                          params)

        @test all(omega[1] .≈ 1.0)
    end

    @testset "generic string expression: xi only (2D)" begin
        nodes = 1
        dof = 2
        Data_Manager.initialize_data()
        Data_Manager.set_num_controller(nodes)
        Data_Manager.set_dof(dof)
        nn = Data_Manager.create_constant_node_scalar_field("Number of Neighbors", Int64)
        nn[1] = 2
        bg = Data_Manager.create_constant_bond_vector_state("Bond Geometry", Float64,
                                                            dof + 1)
        omega = Data_Manager.create_constant_bond_scalar_state("Influence Function",
                                                               Float64)

        bg[1][1] .= [1.0, 0.0, 2.0]
        bg[1][2] .= [0.0, 1.0, 4.0]

        params = Dict("Influence Function" => "1/xi")
        PeriLab.Solver_Manager.Influence_Function.init_influence_function(Vector{Int64}(1:nodes),
                                                                          params)

        @test omega[1][1] ≈ 1.0 / 2.0
        @test omega[1][2] ≈ 1.0 / 4.0
    end

    @testset "generic string expression: uses xiX, xiY (2D)" begin
        nodes = 1
        dof = 2
        Data_Manager.initialize_data()
        Data_Manager.set_num_controller(nodes)
        Data_Manager.set_dof(dof)
        nn = Data_Manager.create_constant_node_scalar_field("Number of Neighbors", Int64)
        nn[1] = 1
        bg = Data_Manager.create_constant_bond_vector_state("Bond Geometry", Float64,
                                                            dof + 1)
        omega = Data_Manager.create_constant_bond_scalar_state("Influence Function",
                                                               Float64)

        bg[1][1] .= [3.0, 4.0, 5.0]

        params = Dict("Influence Function" => "xiX^2 + xiY^2")
        PeriLab.Solver_Manager.Influence_Function.init_influence_function(Vector{Int64}(1:nodes),
                                                                          params)

        @test omega[1][1] ≈ 3.0^2 + 4.0^2
    end

    @testset "generic string expression: exp() (2D)" begin
        nodes = 1
        dof = 2
        Data_Manager.initialize_data()
        Data_Manager.set_num_controller(nodes)
        Data_Manager.set_dof(dof)
        nn = Data_Manager.create_constant_node_scalar_field("Number of Neighbors", Int64)
        nn[1] = 1
        bg = Data_Manager.create_constant_bond_vector_state("Bond Geometry", Float64,
                                                            dof + 1)
        omega = Data_Manager.create_constant_bond_scalar_state("Influence Function",
                                                               Float64)

        bg[1][1] .= [1.0, 0.0, 1.0]

        params = Dict("Influence Function" => "exp(-xi)")
        PeriLab.Solver_Manager.Influence_Function.init_influence_function(Vector{Int64}(1:nodes),
                                                                          params)

        @test omega[1][1] ≈ exp(-1.0)
    end

    @testset "generic string expression: xiZ ignored/zero in 2D" begin
        nodes = 1
        dof = 2
        Data_Manager.initialize_data()
        Data_Manager.set_num_controller(nodes)
        Data_Manager.set_dof(dof)
        nn = Data_Manager.create_constant_node_scalar_field("Number of Neighbors", Int64)
        nn[1] = 1
        bg = Data_Manager.create_constant_bond_vector_state("Bond Geometry", Float64,
                                                            dof + 1)
        omega = Data_Manager.create_constant_bond_scalar_state("Influence Function",
                                                               Float64)

        bg[1][1] .= [1.0, 0.0, 1.0]

        params = Dict("Influence Function" => "xiZ + 1.0")
        PeriLab.Solver_Manager.Influence_Function.init_influence_function(Vector{Int64}(1:nodes),
                                                                          params)

        @test omega[1][1] ≈ 1.0  # xiZ should default to 0.0 in 2D
    end

    @testset "generic string expression: 3D with xiX, xiY, xiZ" begin
        nodes = 1
        dof = 3
        Data_Manager.initialize_data()
        Data_Manager.set_num_controller(nodes)
        Data_Manager.set_dof(dof)
        nn = Data_Manager.create_constant_node_scalar_field("Number of Neighbors", Int64)
        nn[1] = 1
        bg = Data_Manager.create_constant_bond_vector_state("Bond Geometry", Float64,
                                                            dof + 1)
        omega = Data_Manager.create_constant_bond_scalar_state("Influence Function",
                                                               Float64)

        # xiX=1, xiY=2, xiZ=2 -> xi = sqrt(1+4+4)=3
        bg[1][1] .= [1.0, 2.0, 2.0, 3.0]

        params = Dict("Influence Function" => "xiX + xiY + xiZ")
        PeriLab.Solver_Manager.Influence_Function.init_influence_function(Vector{Int64}(1:nodes),
                                                                          params)

        @test omega[1][1] ≈ 1.0 + 2.0 + 2.0
    end

    @testset "empty node list -> no error, omega unchanged" begin
        nodes = 1
        dof = 2
        Data_Manager.initialize_data()
        Data_Manager.set_num_controller(nodes)
        Data_Manager.set_dof(dof)
        nn = Data_Manager.create_constant_node_scalar_field("Number of Neighbors", Int64)
        nn[1] = 1
        bg = Data_Manager.create_constant_bond_vector_state("Bond Geometry", Float64,
                                                            dof + 1)
        omega = Data_Manager.create_constant_bond_scalar_state("Influence Function",
                                                               Float64)
        bg[1][1] .= [1.0, 0.0, 1.0]

        params = Dict("Influence Function" => "1/xi^2")
        PeriLab.Solver_Manager.Influence_Function.init_influence_function(Int64[],
                                                                          params)

        @test all(omega[1] .== 0.0)
    end
end
