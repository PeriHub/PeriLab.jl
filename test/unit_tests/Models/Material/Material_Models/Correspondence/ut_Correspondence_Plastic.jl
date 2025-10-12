# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
# include("../../../../../../src/PeriLab.jl")
# using .PeriLab

@testset "get_name&fe_support" begin
    @test PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.Correspondence_Plastic.correspondence_name() == "Correspondence Plastic"
    @test !(PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.Correspondence_Plastic.fe_support())
end

@testset "ut_init_model" begin
    nodes = 2
    test_data_manager = PeriLab.Data_Manager
    test_data_manager.initialize_data()
    material_parameter = Dict()
    test_data_manager.set_num_controller(nodes)
    dof = 3
    test_data_manager.set_dof(dof)
    nn = test_data_manager.create_constant_node_field("Number of Neighbors", Int64, 1)
    nn .= 2
    @test isnothing(PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.Correspondence_Plastic.init_model(test_data_manager,
                                                      Vector{Int64}(1:nodes),
                                                      material_parameter))

    material_parameter = Dict("Shear Modulus" => 10.5)
    @test isnothing(PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.Correspondence_Plastic.init_model(test_data_manager,
                                                      Vector{Int64}(1:nodes),
                                                      material_parameter))

    material_parameter = Dict("Shear Modulus" => 10.5,
                              "Yield Stress" => 3.4)
    test_data_manager = PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.Correspondence_Plastic.init_model(test_data_manager,
                                                          Vector{Int64}(1:nodes),
                                                          material_parameter)

    @test test_data_manager.has_key("von Mises Yield StressN")
    @test test_data_manager.has_key("Plastic StrainN")

    material_parameter = Dict("Shear Modulus" => 10.5,
                              "Yield Stress" => 3.4,
                              "Bond Associated" => true)
    test_data_manager = PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.Correspondence_Plastic.init_model(test_data_manager,
                                                          Vector{Int64}(1:nodes),
                                                          material_parameter)

    @test test_data_manager.has_key("von Mises Bond Yield StressN")
    @test test_data_manager.has_key("Plastic Bond StrainN")
end
