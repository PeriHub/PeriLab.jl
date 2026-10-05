# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

#using Test
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
    PeriLab.Data_Manager.set_num_controller(nodes)
    PeriLab.Data_Manager.set_dof(3)
    nn = PeriLab.Data_Manager.create_constant_node_scalar_field("Number of Neighbors",
                                                                Int64)
    nn .= 2
    PLASTIC = PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.Correspondence_Plastic
    m = typed_block_material(Dict("Material Model" => "Correspondence Elastic + Correspondence Plastic",
                                  "Symmetry" => "isotropic", "Bulk Modulus" => 10.0,
                                  "Shear Modulus" => 10.5, "Yield Stress" => 3.4))
    PLASTIC.init_model(Vector{Int64}(1:nodes), m.model.parts[2], m)
    @test PeriLab.Data_Manager.has_key("von Mises Yield StressN")
    @test PeriLab.Data_Manager.has_key("Plastic StrainN")
    ba = typed_block_material(Dict("Material Model" => "Correspondence Elastic + Correspondence Plastic",
                                   "Symmetry" => "isotropic", "Bulk Modulus" => 10.0,
                                   "Shear Modulus" => 10.5, "Yield Stress" => 3.4,
                                   "Bond Associated" => true))
    PLASTIC.init_model(Vector{Int64}(1:nodes), ba.model.parts[2], ba)
    @test PeriLab.Data_Manager.has_key("von Mises Bond Yield StressN")
    @test PeriLab.Data_Manager.has_key("Plastic Bond StrainN")
end
