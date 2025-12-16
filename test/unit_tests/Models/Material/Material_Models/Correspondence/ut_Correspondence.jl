# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
#include("../../../../../../src/PeriLab.jl")
#using .PeriLab

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
