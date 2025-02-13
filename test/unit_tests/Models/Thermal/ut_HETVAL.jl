# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
using Test
using LinearAlgebra
include("../../../../src/Models/Thermal/hetval.jl")
# include("../../../../src/PeriLab.jl")
# using .PeriLab
@testset "get_name" begin
    @test HETVAL.thermal_model_name() == "HETVAL"
end
@testset "init exceptions" begin
    nodes = 2
    test_data_manager = PeriLab.Data_manager
    test_data_manager.initialize_data()
    test_data_manager.set_num_controller(nodes)
    dof = 3
    test_data_manager.set_dof(dof)
    file = "./src/Models/Thermal/HETVALs/hetval.so"
    if !isfile(file)
        file = "../src/Models/Thermal/HETVALs/hetval.so"
    end
    @test isnothing(
        HETVAL.init_model(
            test_data_manager,
            Vector{Int64}(1:nodes),
            Dict{String,Any}("File" => file * "_not_there"),
        ),
    )

    @test isnothing(
        HETVAL.init_model(test_data_manager, Vector{Int64}(1:nodes), Dict{String,Any}()),
    )

    @test isnothing(
        HETVAL.init_model(
            test_data_manager,
            Vector{Int64}(1:nodes),
            Dict{String,Any}(
                "File" => file,
                "Predefined Field Names" => "test_field_2 test_field_3",
            ),
        ),
    )

    test_1 = test_data_manager.create_constant_node_field("test_field_2", Float64, 1)
    test_1[1] = 7.3
    test_2 = test_data_manager.create_constant_node_field("test_field_3", Float64, 1)
    test_2 .= 3
    mat_dict = Dict{String,Any}(
        "File" => file,
        "HETVAL name" => "test_sub",
        "Number of Properties" => 3,
        "Property_1" => 2,
        "Property_2" => 2,
        "Property_3" => 2.4,
        "Predefined Field Names" => "test_field_2 test_field_3",
    )
    HETVAL.init_model(test_data_manager, Vector{Int64}(1:nodes), mat_dict)
    fields = test_data_manager.get_field("Predefined Fields")
    inc = test_data_manager.get_field("Predefined Fields Increment")
    @test size(fields) == (2, 2)
    @test size(inc) == (2, 2)
    @test fields[1, 1] == test_1[1]
    @test fields[2, 1] == test_1[2]
    @test fields[1, 2] == test_2[1]
    @test fields[2, 2] == test_2[2]
    @test mat_dict["HETVAL name"] == "test_sub"
    mat_dict = Dict{String,Any}(
        "File" => file,
        "Number of Properties" => 3,
        "Property_1" => 2,
        "Property_2" => 2,
        "Property_3" => 2.4,
        "Predefined Field Names" => "test_field_2 test_field_3",
    )
    @test !haskey(mat_dict, "HETVAL name")
    HETVAL.init_model(test_data_manager, Vector{Int64}(1:nodes), mat_dict)
    @test haskey(mat_dict, "HETVAL name")
    @test mat_dict["HETVAL name"] == "HETVAL"
end
@testset "ut_malloc_cstring" begin
    CMNAME = HETVAL.malloc_cstring("HETVAL_TEST")
    @test typeof(CMNAME) == Cstring
end

# Test wrapper function for HETVAL_interface
@testset "HETVAL_interface Tests" begin
    # Example test case (you should define your own)
    test_data_manager = PeriLab.Data_manager
    file = "./src/Models/Thermal/HETVALs/hetval.so"
    if !isfile(file)
        file = "../src/Models/Thermal/HETVALs/hetval.so"
    end
    CMNAME::Cstring = HETVAL.malloc_cstring("HETVAL_TEST")  # Adjust as needed
    TEMP::Vector{Float64} = [200.0, 50]
    TIME::Vector{Float64} = [0.0, 0.1]
    DTIME::Float64 = 0.1
    STATEV::Vector{Float64} = zeros(Float64, 2)  # Adjust the size and values as needed
    FLUX::Vector{Float64} = [0.0, 0.0]
    PREDEF::Vector{Float64} = zeros(Float64, 1)  # Adjust as needed
    DPRED::Vector{Float64} = zeros(Float64, 1)  # Adjust as needed
    HETVAL.hetval_file_path = file
    HETVAL.HETVAL_interface(CMNAME, TEMP, TIME, DTIME, STATEV, FLUX, PREDEF, DPRED)
    @test FLUX[1] == 100
    @test FLUX[2] == 100
    @test TEMP[1] == 2.5
    @test STATEV[2] == 200
end
