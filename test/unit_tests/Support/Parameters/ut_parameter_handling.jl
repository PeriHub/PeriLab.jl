# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

if !isdefined(@__MODULE__, :Data_manager)
    include("../../../../src/Support/data_manager.jl")
end
include("../../../../src/Support/Parameters/parameter_handling.jl")
include("../../../../src/Support/helpers.jl")

using Test
using Random



@testset "ut_get_output_type" begin
    @test get_output_type(Dict("Output1" => Dict()), "Output1") == "Exodus"
    @test get_output_type(Dict("Output1" => Dict("Output File Type" => "CSV")), "Output1") == "CSV"
    @test get_output_type(Dict("Output1" => Dict("Output File Type" => "Exodus"), "Output2" => Dict("Output File Type" => "CSV")), "Output1") == "Exodus"
    @test get_output_type(Dict("Output1" => Dict("Output File Type" => "Exodus"), "Output2" => Dict("Output File Type" => "CSV")), "Output2") == "CSV"
    @test get_output_type(Dict("Output1" => Dict("Output File Type" => "CSV"), "Output2" => Dict()), "Output2") == "Exodus"
    @test get_output_type(Dict("Output1" => Dict("Output File Type" => "Exodus")), "Output1") == "Exodus"
end
@testset "ut_get_bond_filters" begin
    params = Dict("Discretization" => Dict())
    check, bfList = get_bond_filters(params)
    @test !check
    @test bfList == Dict{String,Dict{String,Any}}()
    params = Dict("Discretization" => Dict("Bond Filters" => Dict()))
    check, bfList = get_bond_filters(params)
    @test check
    @test bfList == Dict()
    params = Dict("Discretization" => Dict("Bond Filters" => Dict("a" => Dict("a" => 1))))
    check, bfList = get_bond_filters(params)
    @test check
    @test bfList == Dict("a" => Dict("a" => 1))
    params = Dict("Discretization" => Dict("Bond Filters" => Dict("a" => Dict("a" => 1), "g" => Dict("a" => 1), "adas" => Dict("a" => 1))))
    check, bfList = get_bond_filters(params)
    @test check
    @test bfList == Dict("a" => Dict("a" => 1), "g" => Dict("a" => 1), "adas" => Dict("a" => 1))
end
@testset "ut_node_sets" begin

    filename = "test.txt"
    numbers = [11, 12, 13, 44, 125]
    lenNumbers = length(numbers)
    params = Dict("Discretization" => Dict())

    @test get_node_sets(params, "") == Dict{String,Any}()
    params = Dict("Discretization" => Dict("Node Sets" => Dict("Nset_1" => "1 2 3 4 5 6 7", "Nset_2" => filename)))

    file = open(filename, "w")
    println(file, "header: global_id")
    for number in numbers
        println(file, number)
    end
    close(file)

    nsets = get_node_sets(params, "")
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

    filename = "test.txt"
    file = open(filename, "w")
    println(file, "header: global_id")
    close(file)
    nsets = get_node_sets(params, "")
    @test haskey(nsets, "Nset_1")
    @test !haskey(nsets, "Nset_2")
    rm(filename)
    filename = "test.txt"
    file = open(filename, "w")
    close(file)
    nsets = get_node_sets(params, "")
    @test haskey(nsets, "Nset_1")
    @test !haskey(nsets, "Nset_2")
    rm(filename)
end
@testset "ut_node_set" begin

    filename = "test.txt"
    numbers = [11, 12, 13, 44, 125]
    file = open(filename, "w")
    println(file, "header: global_id")
    for number in numbers
        println(file, number)
    end
    close(file)

    params = Dict("Discretization" => Dict("Type" => "Text File"))
    computes = Dict("Node Setas" => 12)
    @test [] == get_node_set(computes, "", params)
    computes = Dict("Node Set" => 12)
    @test [12] == get_node_set(computes, "", params)
    computes = Dict("Node Set" => filename)
    @test [11, 12, 13, 44, 125] == get_node_set(computes, "", params)
    computes = Dict("Node Set" => "13 44 125")
    @test [13, 44, 125] == get_node_set(computes, "", params)
    rm(filename)
end
@testset "ut_validate_yaml" begin
    params = Dict{Any,Any}()
    @info "Error messages are tested and therefore okay."
    @test isnothing(validate_yaml(params))
    params = Dict{Any,Any}("PeriLab" => Dict{Any,Any}("Physics" => Dict{Any,Any}()))
    @test isnothing(validate_yaml(params))
    params = Dict{Any,Any}("PeriLab" => Dict{Any,Any}("Physics" => Dict{Any,Any}("Material Models" => Dict{Any,Any}()), "Discretization" => Dict{Any,Any}()))
    @test isnothing(validate_yaml(params))
    params = Dict{Any,Any}("PeriLab" => Dict{Any,Any}("Physics" => Dict{Any,Any}("Material Models" => Dict{Any,Any}()), "Discretization" => Dict{Any,Any}(), "Blocks" => Dict{Any,Any}()))
    @test isnothing(validate_yaml(params))
    params = Dict{Any,Any}("PeriLab" => Dict{Any,Any}("Physics" => Dict{Any,Any}("Material Models" => Dict{Any,Any}()), "Blocks" => Dict{Any,Any}()))
    @test isnothing(validate_yaml(params))
    params = Dict{Any,Any}("PeriLab" => Dict{Any,Any}("Blocks" => Dict{Any,Any}()))
    @test isnothing(validate_yaml(params))
    params = Dict{Any,Any}("PeriLab" => Dict{Any,Any}("Physics" => Dict{Any,Any}(), "Discretization" => Dict{Any,Any}(), "Blocks" => Dict{Any,Any}(), "Solver" => Dict{Any,Any}()))
    @test isnothing(validate_yaml(params))
    params = Dict{Any,Any}("PeriLab" => Dict{Any,Any}("Physics" => Dict{Any,Any}("Material Models" => Dict{Any,Any}("mat_1" => Dict{Any,Any}("Material Model" => "a"))), "Discretization" => Dict{Any,Any}("Input Mesh File" => "test", "Type" => "test"), "Blocks" => Dict{Any,Any}("Block_1" => Dict{Any,Any}("Block Names" => "Block_1", "Density" => 1.0, "Horizon" => "1.0")), "Solver" => Dict{Any,Any}("Final Time" => 1.0, "Initial Time" => 0.0)))
    @test isnothing(validate_yaml(params))
    params = Dict{Any,Any}("PeriLab" => Dict{Any,Any}("Physics" => Dict{Any,Any}("Material Models" => Dict{Any,Any}("mat_1" => Dict{Any,Any}("Material Model" => "a"))), "Discretization" => Dict{Any,Any}("Input Mesh File" => "test", "Type" => "test"), "Blocks" => Dict{Any,Any}("Block_1" => Dict{Any,Any}("Block Names" => "Block_1", "Density" => 1.0, "Horizon" => 1.0)), "Solver" => Dict{Any,Any}("Final Time" => 1.0, "Initial Time" => 0.0)))
    @test validate_yaml(params) == params["PeriLab"]
end

@testset "ut_get_mesh_name" begin
    params = Dict("Discretization" => Dict())
    @test get_mesh_name(params) === nothing
    name = randstring(12)
    params = Dict("Discretization" => Dict("Input Mesh File" => name))
    @test get_mesh_name(params) == name
end

@testset "ut_get_output_filenames" begin
    params = Dict()
    filenames = get_output_filenames(params, "")
    @test filenames == []
    params = Dict("Outputs" => Dict("Output1" => Dict("Output Filename" => "1"), "Output2" => Dict("Output Filename" => "2")))
    filenames = get_output_filenames(params, "")
    @test filenames[1] == "1.e"
    @test filenames[2] == "2.e"
    params = Dict("Outputs" => Dict("Output1" => Dict("Output Filename" => "3", "Output File Type" => "CSV"), "Output2" => Dict("Output Filename" => "4", "Output File Type" => "Exodus")))
    filenames = get_output_filenames(params, "test")
    @test filenames[1] == "test/3.csv"
    @test filenames[2] == "test/4.e"
end

@testset "get_output_frequency" begin
    nsteps = 40
    params = Dict()
    params = Dict("Outputs" => Dict("Output1" => Dict("Output Frequency" => 2), "Output2" => Dict("Number of Output Steps" => 1, "Output Frequency" => 1)))
    freq = get_output_frequency(params, nsteps)
    @test freq[1] == 2
    @test freq[2] == 40
    params = Dict("Outputs" => Dict("Output1" => Dict("Output Frequency" => 20), "Output2" => Dict("Number of Output Steps" => 10)))
    freq = get_output_frequency(params, nsteps)
    @test freq[1] == 20
    @test freq[2] == 4
    nsteps = 1000
    freq = get_output_frequency(params, nsteps)
    @test freq[1] == 20
    @test freq[2] == 100
    nsteps = 2
    freq = get_output_frequency(params, nsteps)
    @test freq[1] == 2
    @test freq[2] == 1
    params = Dict("Outputs" => Dict("Output1" => Dict("Output Frequency" => 20, "Number of Output Steps" => 10), "Output2" => Dict("Number of Output Steps" => 10, "Output Frequency" => 20)))
    nsteps = 1000
    freq = get_output_frequency(params, nsteps)
    @test (freq[1] == 100) || (freq[1] == 20)
    @test (freq[2] == 100) || (freq[2] == 20)
end

test_Data_manager = Data_manager
@testset "ut_get_outputs" begin
    test_Data_manager.set_num_controller(5)
    test_Data_manager.create_constant_node_field("A", Float64, 1)
    test_Data_manager.create_node_field("B", Bool, 1)
    test_Data_manager.create_constant_node_field("C", Float64, 4)
    test_Data_manager.create_node_field("D", Int64, 7)
    test_Data_manager.create_node_field("F", Float64, 1)
    test_Data_manager.create_constant_node_field("E", Float64, 4)
    testfield_keys = test_Data_manager.get_all_field_keys()

    params = Dict("Outputs" => Dict("Output1" => Dict("fieldnames" => [], "Output Variables" => Dict("A" => true, "B" => false, "C" => true)), "Output2" => Dict("fieldnames" => [], "Output Variables" => Dict("A" => true, "B" => true, "D" => false, "E" => true, "M" => true))))

    outputs = get_outputs(params, testfield_keys, String[])

    @test "A" in outputs["Output1"]["fieldnames"]
    @test ("BNP1" in outputs["Output1"]["fieldnames"]) == false
    @test "C" in outputs["Output1"]["fieldnames"]
    @test "A" in outputs["Output2"]["fieldnames"]
    @test "BNP1" in outputs["Output2"]["fieldnames"]
    @test ("D" in outputs["Output2"]["fieldnames"]) == false
    @test "E" in outputs["Output2"]["fieldnames"]
    @test !("M" in outputs["Output2"]["fieldnames"])
    params = Dict("Outputs" => Dict("Output1" => Dict("fieldnames" => [], "Output File Type" => "CSV", "Output Variables" => Dict("E" => true, "B" => false, "C" => true)), "Output2" => Dict("fieldnames" => [], "Output Variables" => Dict("A" => true, "B" => true, "D" => false, "E" => true, "M" => true))))
    outputs = get_outputs(params, testfield_keys, String["M"])
    @test !("A" in outputs["Output1"]["fieldnames"])
    @test !("BNP1" in outputs["Output1"]["fieldnames"])
    @test !("C" in outputs["Output1"]["fieldnames"])
    @test ("BNP1" in outputs["Output2"]["fieldnames"])
    @test !("D" in outputs["Output2"]["fieldnames"])
    @test ("E" in outputs["Output2"]["fieldnames"])
    @test ("M" in outputs["Output2"]["fieldnames"])
    params = Dict("Outputs" => Dict("Output1" => Dict("fieldnames" => [], "Output File Type" => "CSV", "Output Variables" => Dict("M" => true, "A" => true))))
    outputs = get_outputs(params, testfield_keys, String["M"])
    @test !("A" in outputs["Output1"]["fieldnames"])
    @test "M" in outputs["Output1"]["fieldnames"]

    params = Dict("Outputs" => Dict("Output1" => Dict("fieldnames" => [], "Output Variables" => Dict())))
    outputs = get_outputs(params, testfield_keys, String[])
    @test outputs["Output1"]["fieldnames"] == []
    params = Dict("Outputs" => Dict("Output1" => Dict("fieldnames" => [])))
    outputs = get_outputs(params, testfield_keys, String[])
    @test outputs["Output1"]["fieldnames"] == []
end
@testset "ut_get_computes" begin
    params = Dict()
    testfield_keys = test_Data_manager.get_all_field_keys()
    @test get_computes(params, testfield_keys) == Dict()

    params = Dict("Compute Class Parameters" => Dict("External_Forces" => Dict("Compute Class" => "Block_Data", "Calculation Type" => "Sum", "Block" => "block_2", "Variable" => "A"), "External_Displacements" => Dict("Compute Class" => "Block_Data", "Calculation Type" => "Maximum", "Block" => "block_1", "Variable" => "B"), "warn_test" => Dict("Compute Class" => "Block_Data", "Calculation Type" => "Maximum", "Block" => "block_1")))

    computes = get_computes(params, testfield_keys)

    @test haskey(computes, "External_Forces")
    @test haskey(computes, "External_Displacements")
    @test !haskey(computes, "warn_test")
    @test computes["External_Forces"]["Variable"] == "A"
    @test computes["External_Displacements"]["Variable"] == "BNP1"
end
@testset "ut_get_computes_names" begin
    testfield_keys = test_Data_manager.get_all_field_keys()

    params = Dict("Compute Class Parameters" => Dict("External_Forces" => Dict("Compute Class" => "Block_Data", "Calculation Type" => "Sum", "Block" => "block_2", "Variable" => "A"), "External_Displacements" => Dict("Compute Class" => "Block_Data", "Calculation Type" => "Maximum", "Block" => "block_1", "Variable" => "B")))

    computes_names = get_computes_names(params)

    @test "External_Forces" in computes_names
    @test "External_Displacements" in computes_names
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
    @test isnothing(get_number_of_blocks(Dict()))
    params = Dict("Blocks" => Dict())
    @test isnothing(get_number_of_blocks(params))
    params = Dict("Blocks" => Dict("block_1" => Dict(), "block_2" => Dict()))
    @test get_number_of_blocks(params) == 2
    params = Dict("Blocks" => Dict("block_1" => Dict(), "block_2" => Dict(), "block_3" => Dict()))
    @test get_number_of_blocks(params) == 3
    params = Dict("Blocks" => Dict("block_1" => Dict(), "block_2" => Dict(), "block_3" => Dict(), "block_4" => Dict()))
    @test get_number_of_blocks(params) == 4
end



function get_density(params::Dict, block_id::Int64)
    return get_values(params, block_id, "Density")
end

function get_heatcapacity(params::Dict, block_id::Int64)
    return get_values(params, block_id, "Specific Heat Capacity")
end

function get_horizon(params::Dict, block_id::Int64)
    return get_values(params, block_id, "Horizon")
end

function get_values(params::Dict, block_id::Int64, valueName::String)
    if haskey(params["Blocks"], "block_" * string(block_id))
        if haskey(params["Blocks"]["block_"*string(block_id)], valueName)
            return params["Blocks"]["block_"*string(block_id)][valueName]
        end
        @error "$valueName of Block $block_id is not defined"
        return
    end
    @error "Block $block_id is not defined"
    return
end

@testset "ut_block_values" begin
    params = Dict("Blocks" => Dict())
    @test isnothing(get_horizon(params, 1))
    @test isnothing(get_density(params, 1))
    @test isnothing(get_heatcapacity(params, 1))
    @test isnothing(get_values(params, 1, "Density"))
    @test isnothing(get_values(params, 1, "not there"))
    params = Dict("Blocks" => Dict("block_1" => Dict(), "block_2" => Dict()))
    @test isnothing(get_horizon(params, 1))
    @test isnothing(get_density(params, 1))
    @test isnothing(get_heatcapacity(params, 1))
    @test isnothing(get_values(params, 1, "Density"))
    @test isnothing(get_horizon(params, 2))
    @test isnothing(get_density(params, 2))
    @test isnothing(get_heatcapacity(params, 2))
    @test isnothing(get_values(params, 2, "Density"))
    params = Dict("Blocks" => Dict("block_1" => Dict("Density" => 1, "Specific Heat Capacity" => 3), "block_2" => Dict("Density" => 12.3, "Horizon" => 2)))
    @test get_values(params, 1, "Density") == 1
    @test get_values(params, 2, "Density") == 12.3
    @test isnothing(get_values(params, 3, "Density"))
    @test get_values(params, 1, "Specific Heat Capacity") == 3
    @test isnothing(get_values(params, 2, "Specific Heat Capacity"))
    @test isnothing(get_values(params, 3, "Specific Heat Capacity"))
    @test isnothing(get_values(params, 1, "Horizon"))
    @test get_values(params, 2, "Horizon") == 2
    @test isnothing(get_values(params, 3, "Horizon"))
    @test get_values(params, 1, "Density") == get_density(params, 1)
    @test get_values(params, 2, "Density") == get_density(params, 2)
    @test get_values(params, 3, "Density") == get_density(params, 3)

    @test get_values(params, 1, "Horizon") == get_horizon(params, 1)
    @test get_values(params, 2, "Horizon") == get_horizon(params, 2)
    @test get_values(params, 3, "Horizon") == get_horizon(params, 3)

    @test get_values(params, 1, "Specific Heat Capacity") == get_heatcapacity(params, 1)
    @test get_values(params, 2, "Specific Heat Capacity") == get_heatcapacity(params, 2)
    @test get_values(params, 3, "Specific Heat Capacity") == get_heatcapacity(params, 3)
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
    @test isnothing(get_initial_time(Dict("Solver" => Dict())))
    @test isnothing(get_final_time(Dict("Solver" => Dict())))
    @test isnothing(get_final_time(Dict("Solver" => Dict())))
    @test isnothing(get_solver_name(Dict("Solver" => Dict("Solvername" => Dict()))))
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
    @test isnothing(get_model_parameter(params, "Does not exist Model", blockModels[1]["Material Model"]))
    @test isnothing(get_model_parameter(params, "Does not exist Model", "s"))
    @test isnothing(get_model_parameter(params, "Material Model", "s"))
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
        "Bond Associated Deformation Gradient" => false)))
    options = Dict{String,Bool}("Deformed Bond Geometry" => true,
        "Deformation Gradient" => false,
        "Shape Tensor" => false,
        "Bond Associated Shape Tensor" => false,
        "Bond Associated Deformation Gradient" => false)
    optionTest = get_physics_option(params, options)
    # no material models included
    @test optionTest == params["Physics"]["Pre Calculation"]

    params = Dict("Physics" => Dict("Pre Calculation" => Dict{String,Bool}("Deformed Bond Geometry" => true,
            "Deformation Gradient" => true,
            "Shape Tensor" => false,
            "Bond Associated Shape Tensor" => false,
            "Bond Associated Deformation Gradient" => true),
        "Material Models" => Dict("a" => Dict("value" => 1), "c" => Dict("value" => [1 2], "value2" => 1))))
    optionTest = get_physics_option(params, options)

    @test isnothing(optionTest)

    params = Dict("Physics" => Dict(
        "Material Models" => Dict("a" => Dict("Material Model" => "adaCoB", "value" => 1), "c" => Dict("value" => [1 2], "value2" => 1))))
    optionTest = get_physics_option(params, options)
    @test isnothing(optionTest)
    params = Dict("Physics" => Dict(
        "Material Models" => Dict("a" => Dict("value" => 1), "c" => Dict("value" => [1 2], "value2" => 1, "Material Model" => "adaCoB"))))
    optionTest = get_physics_option(params, options)
    @test isnothing(optionTest)

    params = Dict("Physics" => Dict(
        "Material Models" => Dict("a" => Dict("value" => 1, "Material Model" => "adaCoB"), "c" => Dict("value" => [1 2], "value2" => 1, "Material Model" => "adaCoB"))))
    optionTest = get_physics_option(params, options)
    @test optionTest == options

    params = Dict("Physics" => Dict(
        "Material Models" => Dict("a" => Dict("Material Model" => "adaCorrespondence", "value" => 1), "c" => Dict("value" => [1 2], "value2" => 1, "Material Model" => "adaCorresponde"))))
    options = Dict{String,Bool}("Deformed Bond Geometry" => false,
        "Deformation Gradient" => false,
        "Shape Tensor" => false,
        "Bond Associated Shape Tensor" => false,
        "Bond Associated Deformation Gradient" => false)
    optionTest = get_physics_option(params, options)
    @test optionTest["Shape Tensor"]
    @test !optionTest["Bond Associated Shape Tensor"]
    @test !optionTest["Bond Associated Deformation Gradient"]
    @test optionTest["Deformation Gradient"]
    @test optionTest["Deformed Bond Geometry"]

    params = Dict("Physics" => Dict(
        "Material Models" => Dict("aBond Associated" => Dict("Material Model" => "adaCoBond Associated", "value" => 1), "c" => Dict("value" => [1 2], "value2" => 1, "Material Model" => "adaCoB"))))
    optionTest = get_physics_option(params, options)

    @test optionTest["Shape Tensor"]
    @test optionTest["Bond Associated Shape Tensor"]
    @test optionTest["Bond Associated Deformation Gradient"]
    @test optionTest["Deformation Gradient"]
    @test optionTest["Deformed Bond Geometry"]
    params = Dict("Physics" => Dict(
        "Material Models" => Dict("aCorrespondence" => Dict("Material Model" => "adaCorrespondence", "value" => 1), "c Bond Associated" => Dict("Material Model" => "adaCoBond Associated", "value" => [1 2], "value2" => 1))))
    options = Dict{String,Bool}("Deformed Bond Geometry" => false,
        "Deformation Gradient" => false,
        "Shape Tensor" => false,
        "Bond Associated Shape Tensor" => false,
        "Bond Associated Deformation Gradient" => false)
end

@testset "ut_check_for_duplicates" begin
    check_for_duplicates(["a", "b", "c"])
    check_for_duplicates(["a", "b", "c", "a"])
end