# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
using LinearAlgebra
# include("../../../../../../src/PeriLab.jl")
# using .PeriLab
@testset "get_name&fe_support" begin
    @test PeriLab.Solver_control.Model_Factory.Material.Correspondence.Correspondence_VUMAT.correspondence_name() == "Correspondence VUMAT"
    @test PeriLab.Solver_control.Model_Factory.Material.Correspondence.Correspondence_VUMAT.fe_support()
end
@testset "init exceptions" begin
    nodes = 2
    test_data_manager = PeriLab.Data_manager
    test_data_manager.initialize_data()
    test_data_manager.set_num_controller(nodes)
    dof = 3
    test_data_manager.set_dof(dof)
    file = "./src/Models/Material/UMATs/libperuser.so"
    if !isfile(file)
        file = "../src/Models/Material/UMATs/libperuser.so"
    end

    @test !isnothing(PeriLab.Solver_control.Model_Factory.Material.Correspondence.Correspondence_VUMAT.init_model(test_data_manager,
                                                     Vector{Int64}(1:nodes),
                                                     Dict{String,Any}("File" => file,
                                                                      "Number of Properties" => 3,
                                                                      "Property_1" => 2,
                                                                      "Property_2" => 2,
                                                                      "Property_3" => 2.4,
                                                                      "Property_4" => 2)))
    @test isnothing(PeriLab.Solver_control.Model_Factory.Material.Correspondence.Correspondence_VUMAT.init_model(test_data_manager,
                                                    Vector{Int64}(1:nodes),
                                                    Dict{String,Any}("File" => file *
                                                                               "_not_there")))
    @test isnothing(PeriLab.Solver_control.Model_Factory.Material.Correspondence.Correspondence_VUMAT.init_model(test_data_manager,
                                                    Vector{Int64}(1:nodes),
                                                    Dict{String,Any}("File" => file)))

    @test isnothing(PeriLab.Solver_control.Model_Factory.Material.Correspondence.Correspondence_VUMAT.init_model(test_data_manager,
                                                    Vector{Int64}(1:nodes),
                                                    Dict{String,Any}()))

    @test isnothing(PeriLab.Solver_control.Model_Factory.Material.Correspondence.Correspondence_VUMAT.init_model(test_data_manager,
                                                    Vector{Int64}(1:nodes),
                                                    Dict{String,Any}("File" => file,
                                                                     "Number of Properties" => 3,
                                                                     "Property_1" => 2,
                                                                     "Property_2" => 2.4,
                                                                     "Property_3" => 2.4,
                                                                     "VUMAT Material Name" => "a"^81)))
    @test !isnothing(PeriLab.Solver_control.Model_Factory.Material.Correspondence.Correspondence_VUMAT.init_model(test_data_manager,
                                                     Vector{Int64}(1:nodes),
                                                     Dict{String,Any}("File" => file,
                                                                      "Number of Properties" => 3,
                                                                      "Property_1" => 2,
                                                                      "Property_2" => 2,
                                                                      "Property_3" => 2.4,
                                                                      "VUMAT Material Name" => "a"^80)))

    properties = test_data_manager.get_field("Properties")
    @test length(properties) == 3
    @test properties[1] == 2
    @test properties[2] == 2
    @test properties[3] == 2.4

    mat_dict = Dict{String,Any}("File" => file,
                                "VUMAT name" => "test_sub",
                                "Number of Properties" => 3,
                                "Property_1" => 2,
                                "Property_2" => 2,
                                "Property_3" => 2.4)
    PeriLab.Solver_control.Model_Factory.Material.Correspondence.Correspondence_VUMAT.init_model(test_data_manager, Vector{Int64}(1:nodes), mat_dict)
    @test mat_dict["VUMAT name"] == "test_sub"
    mat_dict = Dict{String,Any}("File" => file,
                                "Number of Properties" => 3,
                                "Property_1" => 2,
                                "Property_2" => 2,
                                "Property_3" => 2.4)
    @test !haskey(mat_dict, "VUMAT name")
    PeriLab.Solver_control.Model_Factory.Material.Correspondence.Correspondence_VUMAT.init_model(test_data_manager, Vector{Int64}(1:nodes), mat_dict)
    @test haskey(mat_dict, "VUMAT name")
    @test mat_dict["VUMAT name"] == "VUMAT"
end
