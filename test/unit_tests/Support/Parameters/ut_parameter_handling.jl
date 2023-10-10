# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

if !isdefined(@__MODULE__, :Data_manager)
    include("../../../../src/Support/data_manager.jl")
end
include("../../../../src/Support/Parameters/parameter_handling.jl")

using Test
using Random

@testset "ut_check_key_elements" begin
    params = Dict()
    @info "Error messages are tested and therefore okay."
    @test check_key_elements(params) == false
    params = Dict("Physics" => Dict())
    @test check_key_elements(params) == false
    params = Dict("Physics" => Dict("Material Models" => Dict()), "Discretization" => Dict())
    @test check_key_elements(params) == false
    params = Dict("Physics" => Dict("Material Models" => Dict()), "Discretization" => Dict(), "Blocks" => Dict())
    @test check_key_elements(params) == false
    params = Dict("Physics" => Dict("Material Models" => Dict()), "Blocks" => Dict())
    @test check_key_elements(params) == false
    params = Dict("Blocks" => Dict())
    @test check_key_elements(params) == false
    params = Dict("Physics" => Dict(), "Discretization" => Dict(), "Blocks" => Dict(), "Solver" => Dict())
    @test check_key_elements(params) == false
    @info "No error messages are okay for this test until now."
    params = Dict("Physics" => Dict("Material Models" => Dict()), "Discretization" => Dict(), "Blocks" => Dict(), "Solver" => Dict())
    @test check_key_elements(params) == true
end

@testset "ut_get_mesh_name" begin
    params = Dict("Discretization" => Dict())
    @test get_mesh_name(params) === nothing
    name = randstring(12)
    params = Dict("Discretization" => Dict("Input Mesh File" => name))
    @test get_mesh_name(params) == name
end

@testset "ut_search_for_duplicates" begin
    filenames = ["a", "b", "c"]
    newfilenames = search_for_duplicates(filenames)
    for i in 1:3
        @test filenames[i] == newfilenames[i]
    end
    filenames = ["a", "b", "b", "b", "a", "c"]
    newfilenames = search_for_duplicates(filenames)
    @test newfilenames[1] == "a_1"
    @test newfilenames[2] == "a_2"
    @test newfilenames[3] == "b_1"
    @test newfilenames[4] == "b_2"
    @test newfilenames[5] == "b_3"
    @test newfilenames[6] == "c"
end

@testset "ut_get_output_filenames" begin
    params = Dict()
    filenames = get_output_filenames(params)
    @test filenames == []
    params = Dict("Outputs" => Dict("Output1" => Dict("Output Filename" => "1.e"), "Output2" => Dict("Output Filename" => "2.e")))
    filenames = get_output_filenames(params)
    @test filenames[1] == "1.e"
    @test filenames[2] == "2.e"
end

@testset "get_output_frequency" begin
    nsteps = 40

    params = Dict()
    params = Dict("Outputs" => Dict("Output1" => Dict("Output Frequency" => 2), "Output2" => Dict("Number of Outputs" => 1, "Output Frequency" => 1)))
    freq = get_output_frequency(params, nsteps)
    @test freq[1] == 2
    @test freq[2] == 40
    params = Dict("Outputs" => Dict("Output1" => Dict("Output Frequency" => 20), "Output2" => Dict("Number of Outputs" => 10)))
    freq = get_output_frequency(params, nsteps)
    @test freq[1] == 20
    @test freq[2] == 4
    nsteps = 1000
    freq = get_output_frequency(params, nsteps)
    @test freq[1] == 20
    @test freq[2] == 100
end

@testset "ut_get_outputs" begin
    testDatamanager = Data_manager
    testDatamanager.set_nmasters(5)
    testDatamanager.create_constant_node_field("A", Float32, 1)
    testDatamanager.create_node_field("B", Bool, 1)
    testDatamanager.create_constant_node_field("C", Float32, 4)
    testDatamanager.create_node_field("D", Int64, 7)
    testDatamanager.create_node_field("F", Float32, 1)
    testDatamanager.create_constant_node_field("E", Float32, 4)
    testfield_keys = testDatamanager.get_all_field_keys()

    params = Dict("Outputs" => Dict("Output1" => Dict("Output Variables" => Dict("A" => true, "B" => false, "C" => true)), "Output2" => Dict("Output Variables" => Dict("A" => true, "B" => true, "D" => false, "E" => true))))

    outputs = get_outputs(params, testfield_keys)

    @test "A" in outputs[1]
    @test ("BNP1" in outputs[1]) == false
    @test "C" in outputs[1]
    @test "A" in outputs[2]
    @test "BNP1" in outputs[2]
    @test ("D" in outputs[2]) == false
    @test "E" in outputs[2]
end
@testset "ut_node_sets" begin
    numbers = [11, 12, 13, 44, 125]
    lenNumbers = length(numbers)
    filename = "test.txt"
    file = open(filename, "w")
    println(file, "header: global_id")
    for number in numbers
        println(file, number)
    end
    close(file)

    params = Dict("Discretization" => Dict("Node Sets" => Dict("Nset_1" => "1 2 3 4 5 6 7", "Nset_2" => filename)))
    nsets = get_node_sets(params)
    @test "Nset_1" in keys(nsets)
    @test "Nset_2" in keys(nsets)
    @test length(nsets["Nset_1"]) == 7
    for i in 1:7
        @test nsets["Nset_1"][i] == i
    end
    @test length(nsets["Nset_2"]) == lenNumbers
    for i in 1:lenNumbers
        @test nsets["Nset_2"][i] == numbers[i]
    end
    rm(filename)
end
@testset "ut_get_bc_definitions" begin
    params = Dict()
    bcs = get_bc_definitions(params)
    @test length(bcs) == 0
    params = Dict("Boundary Conditions" => Dict())
    bcs = get_bc_definitions(params)
    @test length(bcs) == 0
    params = Dict("Boundary Conditions" => Dict("BC_1" => Dict("Type" => "Force", "Node Set" => "Nset_1", "Coordinate" => "x", "Value" => "20*t"), "BC_2" => Dict("Type" => "Displacement", "Node Set" => "Nset_2", "Coordinate" => "y", "Value" => "0")))
    bcs = get_bc_definitions(params)
    @test length(bcs) == 2
    @test bcs["BC_1"] == Dict("Type" => "Force", "Node Set" => "Nset_1", "Coordinate" => "x", "Value" => "20*t")
    @test bcs["BC_2"] == Dict("Type" => "Displacement", "Node Set" => "Nset_2", "Coordinate" => "y", "Value" => "0")
end
@testset "ut_get_solver_options" begin
    params = Dict("Solver" => Dict("Material Models" => true, "Damage Models" => true, "Additive Models" => true, "Thermal Models" => true))
    solver_options = get_solver_options(params)
    @test solver_options["Additive Models"]
    @test solver_options["Damage Models"]
    @test solver_options["Material Models"]
    @test solver_options["Thermal Models"]
    params = Dict("Solver" => Dict())
    solver_options = get_solver_options(params)
    @test solver_options["Additive Models"] == false
    @test solver_options["Damage Models"] == false
    @test solver_options["Material Models"]
    @test solver_options["Thermal Models"] == false

    params = Dict("Solver" => Dict("Material Models" => false, "Damage Models" => true, "Thermal Models" => true))
    solver_options = get_solver_options(params)
    @test solver_options["Additive Models"] == false
    @test solver_options["Damage Models"]
    @test solver_options["Material Models"] == false
    @test solver_options["Thermal Models"]
end

@testset "ut_get_number_of_blocks" begin
    params = Dict("Blocks" => Dict())
    @test get_number_of_blocks(params) == 0
    params = Dict("Blocks" => Dict("block_1" => Dict(), "block_2" => Dict()))
    @test get_number_of_blocks(params) == 2
    params = Dict("Blocks" => Dict("block_1" => Dict(), "block_2" => Dict(), "block_3" => Dict()))
    @test get_number_of_blocks(params) == 3
    params = Dict("Blocks" => Dict("block_1" => Dict(), "block_2" => Dict(), "block_3" => Dict(), "block_4" => Dict()))
    @test get_number_of_blocks(params) == 4
end

@testset "ut_solver" begin
    params = Dict("Solver" => Dict("Initial Time" => 0.0, "Final Time" => 1.0, "Verlet" => Dict("Safety Factor" => 0.95, "Fixed dt" => 1e-3), "Numerical Damping" => 5e-6))
    @test get_solver_name(params) == "Verlet"
    @test get_final_time(params) == params["Solver"]["Final Time"]
    @test get_initial_time(params) == params["Solver"]["Initial Time"]
    @test get_safety_factor(params) == params["Solver"]["Verlet"]["Safety Factor"]
    @test get_fixed_dt(params) == params["Solver"]["Verlet"]["Fixed dt"]
    @test get_numerical_damping(params) == params["Solver"]["Numerical Damping"]
    params = Dict("Solver" => Dict("Verlet" => Dict()))
    @test get_safety_factor(params) == 1
    @test get_fixed_dt(params) == true
    @test get_numerical_damping(params) == 0.0
end


params = Dict("Physics" => Dict("Material Models" => Dict("A" => Dict("s" => 0, "d" => true), "B" => Dict("sa" => [3.2, 2, 3], "d" => "true")), "Damage Models" => Dict("E" => Dict("ss" => 0, "d" => 1.1))), "Blocks" => Dict("block_1" => Dict("Material Model" => "A", "Damage Model" => "E"), "block_2" => Dict("Material Model" => "B")))

@testset "ut_get_model_parameter" begin
    blockModels = Dict{Int32,Dict{String,String}}()
    for id in 1:2
        blockModels[id] = get_block_models(params, id)
    end
    @test blockModels[1]["Material Model"] == "A"
    @test blockModels[1]["Damage Model"] == "E"
    @test blockModels[2]["Material Model"] == "B"
    testData = Dict("Material Model" => Dict(), "Damage Model" => Dict())

    testData["Material Model"] = get_model_parameter(params, "Material Model", blockModels[1]["Material Model"])
    @test testData["Material Model"]["s"] == 0
    @test testData["Material Model"]["d"] == true
    testData["Material Model"] = get_model_parameter(params, "Material Model", blockModels[2]["Material Model"])
    @test testData["Material Model"]["sa"] == [3.2, 2, 3]
    @test testData["Material Model"]["d"] == "true"
    testData["Damage Model"] = get_model_parameter(params, "Damage Model", blockModels[1]["Damage Model"])
    @test testData["Damage Model"]["ss"] == 0
    @test testData["Damage Model"]["d"] == 1.1
end
@testset "get_physics_option" begin
    params = Dict("Physics" => Dict("Pre Calculation" => Dict{String,Bool}("Deformed Bond Geometry" => false,
            "Deformation Gradient" => false,
            "Shape Tensor" => true,
            "Bond Associated Shape Tensor" => false,
            "Bond Associated Deformation Gradient" => false),
        "Material Models" => Dict("a" => Dict("value" => 1, "name" => "t4"), "c" => Dict("value" => [1 2], "value2" => 1, "name" => "t4"))))

    options = Dict{String,Bool}("Deformed Bond Geometry" => true,
        "Deformation Gradient" => false,
        "Shape Tensor" => false,
        "Bond Associated Shape Tensor" => false,
        "Bond Associated Deformation Gradient" => false)

    optionTest = get_physics_options(params, options)

    @test optionTest == params["Physics"]["Pre Calculation"]

    params = Dict("Physics" => Dict("Pre Calculation" => Dict{String,Bool}("Deformed Bond Geometry" => true,
            "Deformation Gradient" => true,
            "Shape Tensor" => false,
            "Bond Associated Shape Tensor" => false,
            "Bond Associated Deformation Gradient" => true),
        "Material Models" => Dict("a" => Dict("value" => 1, "name" => "t4"), "c" => Dict("value" => [1 2], "value2" => 1, "name" => "t4"))))
    optionTest = get_physics_options(params, options)

    @test optionTest == params["Physics"]["Pre Calculation"]
end