# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

#using Test
using LinearAlgebra
# include("../../../../../../src/PeriLab.jl")
# using .PeriLab
@testset "get_name&fe_support" begin
    @test PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.Correspondence_VUMAT.correspondence_name() ==
          "Correspondence VUMAT"
    @test PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.Correspondence_VUMAT.fe_support()
end
@testset "init exceptions" begin
    nodes = 2

    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(nodes)
    dof = 3
    PeriLab.Data_Manager.set_dof(dof)
    file = "./src/Models/Material/UMATs/libperuser.so"
    if !isfile(file)
        file = "../src/Models/Material/UMATs/libperuser.so"
    end
    directory = PeriLab.Data_Manager.get_directory()
    @test !isnothing(PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.Correspondence_VUMAT.init_model(Vector{Int64}(1:nodes),
                                                                                                                  Dict{String,
                                                                                                                       Any}("File" => file,
                                                                                                                            "Number of Properties" => 3,
                                                                                                                            "Property_1" => 2,
                                                                                                                            "Property_2" => 2,
                                                                                                                            "Property_3" => 2.4,
                                                                                                                            "Property_4" => 2)))
    @test_logs (:error,
                "File $(joinpath(pwd(), directory, file * "_not_there")) does not exist, please check name and directory.") PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.Correspondence_VUMAT.init_model(Vector{Int64}(1:nodes),
                                                                                                                                                                                                                         Dict{String,
                                                                                                                                                                                                                              Any}("File" => file *
                                                                                                                                                                                                                                             "_not_there"))
    @test_logs (:error, "Number of Properties must be at least equal 1") PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.Correspondence_VUMAT.init_model(Vector{Int64}(1:nodes),
                                                                                                                                                                      Dict{String,
                                                                                                                                                                           Any}("File" => file))

    @test_logs (:error, "VUMAT file is not defined.") PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.Correspondence_VUMAT.init_model(Vector{Int64}(1:nodes),
                                                                                                                                                   Dict{String,
                                                                                                                                                        Any}())

    @test_logs (:error,
                "Due to old Fortran standards only a name length of 80 is supported") PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.Correspondence_VUMAT.init_model(Vector{Int64}(1:nodes),
                                                                                                                                                                                   Dict{String,
                                                                                                                                                                                        Any}("File" => file,
                                                                                                                                                                                             "Number of Properties" => 3,
                                                                                                                                                                                             "Property_1" => 2,
                                                                                                                                                                                             "Property_2" => 2.4,
                                                                                                                                                                                             "Property_3" => 2.4,
                                                                                                                                                                                             "VUMAT Material Name" => "a"^81))
    @test !isnothing(PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.Correspondence_VUMAT.init_model(Vector{Int64}(1:nodes),
                                                                                                                  Dict{String,
                                                                                                                       Any}("File" => file,
                                                                                                                            "Number of Properties" => 3,
                                                                                                                            "Property_1" => 2,
                                                                                                                            "Property_2" => 2,
                                                                                                                            "Property_3" => 2.4,
                                                                                                                            "VUMAT Material Name" => "a"^80)))

    properties = PeriLab.Data_Manager.get_field("Properties")
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
    PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.Correspondence_VUMAT.init_model(Vector{Int64}(1:nodes),
                                                                                                 mat_dict)
    @test mat_dict["VUMAT name"] == "test_sub"
    mat_dict = Dict{String,Any}("File" => file,
                                "Number of Properties" => 3,
                                "Property_1" => 2,
                                "Property_2" => 2,
                                "Property_3" => 2.4)
    @test !haskey(mat_dict, "VUMAT name")
    PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.Correspondence_VUMAT.init_model(Vector{Int64}(1:nodes),
                                                                                                 mat_dict)
    @test haskey(mat_dict, "VUMAT name")
    @test mat_dict["VUMAT name"] == "VUMAT"
end
