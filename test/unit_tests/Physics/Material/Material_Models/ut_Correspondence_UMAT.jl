# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
using Test
include("../../../../../src/Physics/Material/Material_Models/Correspondence_UMAT.jl")
include("../../../../../src/Support/data_manager.jl")
include("../../../../../src/Support/data_manager.jl")
@testset "get_name&fe_support" begin
    @test Correspondence_UMAT.correspondence_name() == "Correspondence UMAT"
    @test Correspondence_UMAT.fe_support()
end
@testset "init exceptions" begin
    nodes = 2
    test_Data_manager = Data_manager
    test_Data_manager.set_num_controller(nodes)
    dof = 3
    test_Data_manager.set_dof(dof)

    @test isnothing(Correspondence_UMAT.init_material_model(test_Data_manager, Vector{Int64}(1:nodes), Dict()))
    @test isnothing(Correspondence_UMAT.init_material_model(test_Data_manager, Vector{Int64}(1:nodes), Dict("Number of Properties" => 3, "Property_1" => 2, "Property_3" => 2)))
    @test isnothing(Correspondence_UMAT.init_material_model(test_Data_manager, Vector{Int64}(1:nodes), Dict("Number of Properties" => 3, "Property_1" => 2, "Property_3" => 2, "Property_4" => 2)))


    @test isnothing(Correspondence_UMAT.init_material_model(test_Data_manager, Vector{Int64}(1:nodes), Dict("Number of Properties" => 3, "Property_1" => 2, "Property_2" => 2, "Property_3" => 2.4, "UMAT Material Name" => "a"^81)))
    @test !isnothing(Correspondence_UMAT.init_material_model(test_Data_manager, Vector{Int64}(1:nodes), Dict("Number of Properties" => 3, "Property_1" => 2, "Property_2" => 2, "Property_3" => 2.4, "UMAT Material Name" => "a"^80)))

    @test isnothing(Correspondence_UMAT.init_material_model(test_Data_manager, Vector{Int64}(1:nodes), Dict("Number of Properties" => 3, "Property_1" => 2, "Property_2" => 2, "Property_3" => 2.4, "Property_4" => 2.4, "Predefined Field Names" => "s")))

    properties = test_Data_manager.get_field("Properties")
    @test length(properties) == 3
    @test properties[1] == 2
    @test properties[2] == 2
    @test properties[3] == 2.4

    test_Data_manager.create_constant_node_field("test_field", Float64, 2)
    Correspondence_UMAT.init_material_model(test_Data_manager, Vector{Int64}(1:nodes), Dict("Number of Properties" => 3, "Property_1" => 2, "Property_2" => 2, "Property_3" => 2.4, "Predefined Field Names" => "test_field"))
    @test ("test_field" in test_Data_manager.get_all_field_keys())
end