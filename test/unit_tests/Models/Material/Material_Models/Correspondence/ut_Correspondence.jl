# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
#include("../../../../../../src/PeriLab.jl")
#using .PeriLab
@testset "zero_energy_mode_compensation_exception" begin
    test_data_manager = PeriLab.Data_Manager
    nnodes = 2
    test_data_manager.set_num_controller(nnodes)
    nn = test_data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn[1] = 2
    nn[2] = 3
    nodes = Vector{Int64}(1:2)

    PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.zero_energy_mode_compensation(nodes,
                                                                                               Dict{String,
                                                                                                    Any}(),
                                                                                               0.0,
                                                                                               0.0)
end

@testset "ut_init_model" begin
    test_data_manager = PeriLab.Data_Manager
    material_parameter = Dict{String,Any}()
    @test isnothing(PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.init_model([
                                                                                                1,
                                                                                                2
                                                                                            ],
                                                                                            1,
                                                                                            material_parameter))

    material_parameter = Dict{String,Any}("Material Model" => "Correspondence Non_Exist",
                                          "Symmetry" => "isotropic plane strain")
    @test isnothing(PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.init_model([
                                                                                                1,
                                                                                                2
                                                                                            ],
                                                                                            1,
                                                                                            material_parameter))
end
