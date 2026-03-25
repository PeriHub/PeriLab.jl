# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
# include("../../../../../../src/PeriLab.jl")
# using .PeriLab

@testset "get_name&fe_support" begin
    @test PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.Correspondence_Plastic.correspondence_name() ==
          "Correspondence Plastic"
    @test !(PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.Correspondence_Plastic.fe_support())
end

@testset "ut_init_model" begin
    nodes = 2

    PeriLab.Data_Manager.initialize_data()
    material_parameter = Dict()
    PeriLab.Data_Manager.set_num_controller(nodes)
    dof = 3
    PeriLab.Data_Manager.set_dof(dof)
    nn = PeriLab.Data_Manager.create_constant_node_scalar_field("Number of Neighbors",
                                                                Int64)
    nn .= 2
    @test_logs (:error,
                "Shear Modulus must be defined to be able to run this plastic material") PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.Correspondence_Plastic.init_model(Vector{Int64}(1:nodes),
                                                                                                                                                                                        material_parameter)

    material_parameter = Dict("Shear Modulus" => 10.5)
    @test_logs (:error, "No ''Yield Stress'' is defined.") PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.Correspondence_Plastic.init_model(Vector{Int64}(1:nodes),
                                                                                                                                                          material_parameter)

    material_parameter = Dict("Shear Modulus" => 10.5,
                              "Yield Stress" => 3.4)
    PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.Correspondence_Plastic.init_model(Vector{Int64}(1:nodes),
                                                                                                   material_parameter)

    @test PeriLab.Data_Manager.has_key("von Mises Yield StressN")
    @test PeriLab.Data_Manager.has_key("Plastic StrainN")

    material_parameter = Dict("Shear Modulus" => 10.5,
                              "Yield Stress" => 3.4,
                              "Bond Associated" => true)
    PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.Correspondence_Plastic.init_model(Vector{Int64}(1:nodes),
                                                                                                   material_parameter)

    @test PeriLab.Data_Manager.has_key("von Mises Bond Yield StressN")
    @test PeriLab.Data_Manager.has_key("Plastic Bond StrainN")
end
