# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
using Test
using LinearAlgebra
include("../../../../../src/Physics/Material/Material_Models/Correspondence_UMAT.jl")
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
@testset "ut_malloc_cstring" begin
    CMNAME = Correspondence_UMAT.malloc_cstring("UMAT_TEST")
    @test typeof(CMNAME) == Cstring
end

# Test wrapper function for UMAT_interface
#@testset "UMAT_interface Tests" begin
# Example test case (you should define your own)
file = "./src/Physics/Material/UMATs/libusertest.so"
if !isfile(file)
    file = "../src/Physics/Material/UMATs/libusertest.so"
end
STRESS::Vector{Float64} = zeros(Float64, 6)  # Example initialization, adjust the size as needed
STATEV::Vector{Float64} = zeros(Float64, 10)  # Adjust the size and values as needed
DDSDDE::Matrix{Float64} = Matrix{Float64}(I, 6, 6)  # Identity matrix for simplicity, adjust as needed
SSE::Float64 = 0.0
SPD::Float64 = 0.0
SCD::Float64 = 0.0
RPL::Float64 = 0.0
DDSDDT::Vector{Float64} = zeros(Float64, 6)  # Adjust as needed
DRPLDE::Vector{Float64} = zeros(Float64, 6)  # Adjust as needed
DRPLDT::Float64 = 0.0  # Adjust as needed
STRAN::Vector{Float64} = zeros(Float64, 6)  # Adjust as needed
DSTRAN::Vector{Float64} = zeros(Float64, 6)  # Adjust as needed
TIME::Float64 = 0.0
DTIME::Float64 = 0.1
TEMP::Float64 = 25.0
DTEMP::Float64 = 0.1
PREDEF::Vector{Float64} = zeros(Float64, 1)  # Adjust as needed
DPRED::Vector{Float64} = zeros(Float64, 1)  # Adjust as needed
CMNAME::Cstring = Correspondence_UMAT.malloc_cstring("UMAT_TEST")  # Adjust as needed
NDI::Int64 = 3
NSHR::Int64 = 3
NTENS::Int64 = 6
NSTATEV::Int64 = length(STATEV)
PROPS::Vector{Float64} = zeros(Float64, 2)  # Adjust as needed
NPROPS::Int64 = 2
COORDS::Vector{Float64} = zeros(Float64, 3)  # Adjust as needed
DROT::Matrix{Float64} = Matrix{Float64}(I, 3, 3)  # Adjust as needed
PNEWDT::Float64 = 0.1
CELENT::Float64 = 1.0
DFGRD0::Matrix{Float64} = Matrix{Float64}(I, 3, 3)  # Adjust as needed
DFGRD1::Matrix{Float64} = Matrix{Float64}(I, 3, 3)  # Adjust as needed
NOEL::Int64 = 1
NPT::Int64 = 1
LAYER::Int64 = 1
KSPT::Int64 = 1
JSTEP::Int64 = 1
KINC::Int64 = 1

# Call the UMAT_interface
Correspondence_UMAT.UMAT_interface(file, STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATEV, PROPS, NPROPS, COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER, KSPT, JSTEP, KINC)
println()
# Example assertion (you should define appropriate checks based on your UMAT's expected behavior)
# @test STRESS != zeros(Float64, 6)  # Check if STRESS was updated (or any other appropriate check)
#end