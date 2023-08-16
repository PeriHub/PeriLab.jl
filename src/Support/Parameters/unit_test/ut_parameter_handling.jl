include("../parameter_handling.jl")
include("../../data_manager.jl")
using Test
import .Data_manager
using Random

@testset "ut_check_key_elements" begin
    params = Dict()
    @info "Error messages are tested and therefore okay."
    @test check_key_elements(params) == false
    params = Dict("Damage Models" => Dict())
    @test check_key_elements(params) == false
    params = Dict("Damage Models" => Dict(), "Materials Models" => Dict())
    @test check_key_elements(params) == false
    params = Dict("Damage Models" => Dict(), "Materials Models" => Dict(), "Discretization" => Dict())
    @test check_key_elements(params) == false
    params = Dict("Damage Models" => Dict(), "Materials" => Dict(), "Discretization" => Dict(), "Blocks" => Dict())
    @test check_key_elements(params) == false
    @info "No error messages are okay for this test until now."
    params = Dict("Damage Models" => Dict(), "Materials Models" => Dict(), "Discretization" => Dict(), "Blocks" => Dict(), "Solver" => Dict())
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
    params = Dict("Output" => Dict("Output1" => Dict("Output Filename" => "1.e"), "Output2" => Dict("Output Filename" => "2.e")))
    filenames = get_output_filenames(params)
    @test filenames[1] == "1.e"
    @test filenames[2] == "2.e"
end
@testset "ut_get_outputs" begin
    testDatamanager = Data_manager
    testDatamanager.set_nnodes(5)
    testDatamanager.create_constant_node_field("A", Float32, 1)
    testDatamanager.create_node_field("B", Bool, 1)
    testDatamanager.create_constant_node_field("C", Float32, 4)
    testDatamanager.create_node_field("D", Int64, 7)
    testDatamanager.create_node_field("F", Float32, 1)
    testDatamanager.create_constant_node_field("E", Float32, 4)
    testfield_keys = testDatamanager.get_all_field_keys()

    params = Dict("Output" => Dict("Output1" => Dict("Output Variables" => Dict("A" => true, "B" => false, "C" => true)), "Output2" => Dict("Output Variables" => Dict("A" => true, "B" => true, "D" => false, "E" => true))))

    outputs = get_outputs(params, testfield_keys)

    @test "A" in outputs["Output1"]
    @test ("B" in outputs["Output1"]) == false
    @test "C" in outputs["Output1"]
    @test "A" in outputs["Output2"]
    @test "B" in outputs["Output2"]
    @test ("D" in outputs["Output2"]) == false
    @test "E" in outputs["Output2"]
end
@testset "ut_node_sets" begin
    numbers = [11, 12, 13, 44, 125]
    lenNumbers = length(numbers)
    filename = "test.txt"
    file = open(filename, "w")
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
testDatamanager = Data_manager
block_list = [1, 2, 3]
testDatamanager.set_block_list(block_list)
prop_keys = testDatamanager.init_property()
params = Dict("Blocks" => Dict("block_1" => Dict("Material Models" => "a"), "block_2" => Dict("Material  Models" => "c"), "block_3" => Dict("Material Models" => "a", "Damage Models" => "a", "Thermal Models" => "therm"), "Material  Models" => Dict("a" => Dict("value" => 1), "c" => Dict("value" => [1 2], "value2" => 1)), "Damage Models" => Dict("a" => Dict("value" => 3), "c" => Dict("value" => [1 2], "value2" => 1)), "Thermal Models" => Dict("therm" => Dict("value" => "hot", "bool" => true))))

for block in block_list
    get_block_model_definition(params, block, prop_keys, testDatamanager.set_properties)
end
"""
function get_block_model_definition(params, blockID, prop_keys, properties)
    if check_element(params, "block_" * string(blockID))
        block = params["block_"*string(blockID)]
        for model in prop_keys
            if check_element(block, model)
                properties(blockID, model) = get_model_parameter(params, model, block[model])
            end
        end
    end
    return property
end
"""