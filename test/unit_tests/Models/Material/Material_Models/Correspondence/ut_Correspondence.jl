# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
#include("../../../../../../src/PeriLab.jl")
#using .PeriLab

@testset "ut_init_model" begin
    material_parameter = Dict{String,Any}()
    @test_logs (:error,
                "Symmetry for correspondence material is missing; options are 'isotropic plane strain', 'isotropic plane stress', 'anisotropic plane stress', 'anisotropic plane stress','isotropic' and 'anisotropic'. For 3D the plane stress or plane strain option is ignored.") PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.init_model([
                                                                                                                                                                                                                                                                                                                                                                 1,
                                                                                                                                                                                                                                                                                                                                                                 2
                                                                                                                                                                                                                                                                                                                                                             ],
                                                                                                                                                                                                                                                                                                                                                             1,
                                                                                                                                                                                                                                                                                                                                                             material_parameter)

    material_parameter = Dict{String,Any}("Material Model" => "Correspondence Non_Exist",
                                          "Symmetry" => "isotropic plane strain")
    @test_logs (:error,
                "No correspondence material of name Correspondence Non_Exist exists.") PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.init_model([
                                                                                                                                                                   1,
                                                                                                                                                                   2
                                                                                                                                                               ],
                                                                                                                                                               1,
                                                                                                                                                               material_parameter)
end
