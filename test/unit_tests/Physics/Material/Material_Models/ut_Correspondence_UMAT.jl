# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
using Test
include("../../../../../src/Physics/Material/Material_Models/Correspondence_UMAT.jl")
# include("../../../../../src/Support/data_manager.jl")
@testset "get_name&fe_support" begin
    @test Correspondence_UMAT.correspondence_name() == "Correspondence UMAT"
    @test Correspondence_UMAT.fe_support()
end
@testset "init exceptions" begin
    nodes = 2
    test_Data_manager = PeriLab.Data_manager
    test_Data_manager.set_num_controller(nodes)
    dof = 3
    test_Data_manager.set_dof(dof)
    file = "./src/Physics/Material/UMATs/libperuser.so"
    if !isfile(file)
        file = "../src/Physics/Material/UMATs/libperuser.so"
    end
    @test !isnothing(Correspondence_UMAT.init_material_model(test_Data_manager, Vector{Int64}(1:nodes), Dict{String,Any}("File" => file, "Number of Properties" => 3, "Property_1" => 2, "Property_2" => 2, "Property_3" => 2.4, "Property_4" => 2)))
    @test isnothing(Correspondence_UMAT.init_material_model(test_Data_manager, Vector{Int64}(1:nodes), Dict{String,Any}("File" => file * "_not_there")))
    @test isnothing(Correspondence_UMAT.init_material_model(test_Data_manager, Vector{Int64}(1:nodes), Dict{String,Any}("File" => file)))

    @test isnothing(Correspondence_UMAT.init_material_model(test_Data_manager, Vector{Int64}(1:nodes), Dict{String,Any}()))
    @test isnothing(Correspondence_UMAT.init_material_model(test_Data_manager, Vector{Int64}(1:nodes), Dict{String,Any}("File" => file, "Number of Properties" => 3, "Property_1" => 2, "Property_3" => 2.4)))
    @test isnothing(Correspondence_UMAT.init_material_model(test_Data_manager, Vector{Int64}(1:nodes), Dict{String,Any}("File" => file, "Number of Properties" => 3, "Property_1" => 2, "Property_2" => 2.4, "Property_4" => 2)))

    @test isnothing(Correspondence_UMAT.init_material_model(test_Data_manager, Vector{Int64}(1:nodes), Dict{String,Any}("File" => file, "Number of Properties" => 3, "Property_1" => 2, "Property_2" => 2.4, "Property_3" => 2.4, "UMAT Material Name" => "a"^81)))
    @test !isnothing(Correspondence_UMAT.init_material_model(test_Data_manager, Vector{Int64}(1:nodes), Dict{String,Any}("File" => file, "Number of Properties" => 3, "Property_1" => 2, "Property_2" => 2, "Property_3" => 2.4, "UMAT Material Name" => "a"^80)))

    properties = test_Data_manager.get_field("Properties")
    @test length(properties) == 3
    @test properties[1] == 2
    @test properties[2] == 2
    @test properties[3] == 2.4

    @test isnothing(Correspondence_UMAT.init_material_model(test_Data_manager, Vector{Int64}(1:nodes), Dict{String,Any}("File" => file, "Number of Properties" => 3, "Property_1" => 2, "Property_2" => 2, "Property_3" => 2.4, "Predefined Field Names" => "test_field_2 test_field_3")))

    test_1 = test_Data_manager.create_constant_node_field("test_field_2", Float64, 1)
    test_1[1] = 7.3
    test_2 = test_Data_manager.create_constant_node_field("test_field_3", Float64, 1)
    test_2 .= 3
    Correspondence_UMAT.init_material_model(test_Data_manager, Vector{Int64}(1:nodes), Dict{String,Any}("File" => file, "Number of Properties" => 3, "Property_1" => 2, "Property_2" => 2, "Property_3" => 2.4, "Predefined Field Names" => "test_field_2 test_field_3"))
    fields = test_Data_manager.get_field("Predefined Fields")
    inc = test_Data_manager.get_field("Predefined Fields Increment")
    @test size(fields) == (2, 2)
    @test size(inc) == (2, 2)
    @test fields[1, 1] == test_1[1]
    @test fields[2, 1] == test_1[2]
    @test fields[1, 2] == test_2[1]
    @test fields[2, 2] == test_2[2]
end