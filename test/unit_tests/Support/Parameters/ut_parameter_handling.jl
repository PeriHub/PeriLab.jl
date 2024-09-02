# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause


using Test
using Random
using Dierckx

#include("../../../../src/PeriLab.jl")
#using .PeriLab
@testset "ut_get_element_degree" begin
    @test isnothing(PeriLab.Solver.Parameter_Handling.get_element_degree(Dict()))
    @test isnothing(
        PeriLab.Solver.Parameter_Handling.get_element_degree(Dict("Degree" => "ABC")),
    )
    @test isnothing(
        PeriLab.Solver.Parameter_Handling.get_element_degree(Dict("Degree" => "1")),
    )
    @test PeriLab.Solver.Parameter_Handling.get_element_degree(Dict("Degree" => 1)) == 1
    @test PeriLab.Solver.Parameter_Handling.get_element_degree(
        Dict("Degree" => [1, 2, 3]),
    ) == [1, 2, 3]
    @test PeriLab.Solver.Parameter_Handling.get_element_degree(
        Dict("Degree" => [1, 2, 3, 5]),
    ) == [1, 2, 3, 5]
end

@testset "ut_get_element_type" begin
    @test isnothing(PeriLab.Solver.Parameter_Handling.get_element_type(Dict()))
    @test PeriLab.Solver.Parameter_Handling.get_element_type(
        Dict("Element Type" => "ABC"),
    ) == "ABC"
    @test PeriLab.Solver.Parameter_Handling.get_element_type(Dict("Element Type" => 12)) ==
          "12"
end

@testset "ut_get_output_type" begin
    @test PeriLab.Solver.Parameter_Handling.get_output_type(
        Dict("Output1" => Dict()),
        "Output1",
    ) == "Exodus"
    @test PeriLab.Solver.Parameter_Handling.get_output_type(
        Dict("Output1" => Dict("Output File Type" => "CSV")),
        "Output1",
    ) == "CSV"
    @test PeriLab.Solver.Parameter_Handling.get_output_type(
        Dict(
            "Output1" => Dict("Output File Type" => "Exodus"),
            "Output2" => Dict("Output File Type" => "CSV"),
        ),
        "Output1",
    ) == "Exodus"
    @test PeriLab.Solver.Parameter_Handling.get_output_type(
        Dict(
            "Output1" => Dict("Output File Type" => "Exodus"),
            "Output2" => Dict("Output File Type" => "CSV"),
        ),
        "Output2",
    ) == "CSV"
    @test PeriLab.Solver.Parameter_Handling.get_output_type(
        Dict("Output1" => Dict("Output File Type" => "CSV"), "Output2" => Dict()),
        "Output2",
    ) == "Exodus"
    @test PeriLab.Solver.Parameter_Handling.get_output_type(
        Dict("Output1" => Dict("Output File Type" => "Exodus")),
        "Output1",
    ) == "Exodus"
end
@testset "ut_get_bond_filters" begin
    params = Dict("Discretization" => Dict())
    check, bfList = PeriLab.Solver.Parameter_Handling.get_bond_filters(params)
    @test !check
    @test bfList == Dict{String,Dict{String,Any}}()
    params = Dict("Discretization" => Dict("Bond Filters" => Dict()))
    check, bfList = PeriLab.Solver.Parameter_Handling.get_bond_filters(params)
    @test check
    @test bfList == Dict()
    params = Dict("Discretization" => Dict("Bond Filters" => Dict("a" => Dict("a" => 1))))
    check, bfList = PeriLab.Solver.Parameter_Handling.get_bond_filters(params)
    @test check
    @test bfList == Dict("a" => Dict("a" => 1))
    params = Dict(
        "Discretization" => Dict(
            "Bond Filters" => Dict(
                "a" => Dict("a" => 1),
                "g" => Dict("a" => 1),
                "adas" => Dict("a" => 1),
            ),
        ),
    )
    check, bfList = PeriLab.Solver.Parameter_Handling.get_bond_filters(params)
    @test check
    @test bfList ==
          Dict("a" => Dict("a" => 1), "g" => Dict("a" => 1), "adas" => Dict("a" => 1))
end

@testset "ut_validate_yaml" begin
    params = Dict{Any,Any}()
    @info "Error messages are tested and therefore okay."

    @test isnothing(PeriLab.Solver.Parameter_Handling.validate_yaml(params))
    params = Dict{Any,Any}("PeriLab" => Dict{Any,Any}("Models" => Dict{Any,Any}()))
    @test isnothing(PeriLab.Solver.Parameter_Handling.validate_yaml(params))
    params = Dict{Any,Any}(
        "PeriLab" => Dict{Any,Any}(
            "Models" => Dict{Any,Any}("Material Models" => Dict{Any,Any}()),
            "Discretization" => Dict{Any,Any}(),
        ),
    )
    @test isnothing(PeriLab.Solver.Parameter_Handling.validate_yaml(params))
    params = Dict{Any,Any}(
        "PeriLab" => Dict{Any,Any}(
            "Models" => Dict{Any,Any}("Material Models" => Dict{Any,Any}()),
            "Discretization" => Dict{Any,Any}(),
            "Blocks" => Dict{Any,Any}(),
        ),
    )
    @test isnothing(PeriLab.Solver.Parameter_Handling.validate_yaml(params))
    params = Dict{Any,Any}(
        "PeriLab" => Dict{Any,Any}(
            "Models" => Dict{Any,Any}("Material Models" => Dict{Any,Any}()),
            "Blocks" => Dict{Any,Any}(),
        ),
    )
    @test isnothing(PeriLab.Solver.Parameter_Handling.validate_yaml(params))

    params = Dict{Any,Any}(
        "PeriLab" => Dict{Any,Any}(
            "Models" => Dict{Any,Any}(),
            "Discretization" => Dict{Any,Any}(),
            "Blocks" => Dict{Any,Any}(),
            "Solver" => Dict{Any,Any}(),
        ),
    )

    @test isnothing(PeriLab.Solver.Parameter_Handling.validate_yaml(params))
    params = Dict{Any,Any}(
        "PeriLab" => Dict{Any,Any}(
            "Models" => Dict{Any,Any}(
                "Material Models" => Dict{Any,Any}(
                    "mat_1" => Dict{Any,Any}("Material Model" => "a"),
                ),
            ),
            "Discretization" =>
                Dict{Any,Any}("Input Mesh File" => "test", "Type" => "test"),
            "Blocks" => Dict{Any,Any}(
                "Block_1" => Dict{Any,Any}(
                    "Block Names" => "Block_1",
                    "Density" => 1.0,
                    "Horizon" => "1.0",
                ),
            ),
            "Solver" => Dict{Any,Any}("Final Time" => 1.0, "Initial Time" => 0.0),
        ),
    )
    @test isnothing(PeriLab.Solver.Parameter_Handling.validate_yaml(params))
    params = Dict{Any,Any}(
        "PeriLab" => Dict{Any,Any}(
            "Models" => Dict{Any,Any}(
                "Material Models" => Dict{Any,Any}(
                    "mat_1" => Dict{Any,Any}("Material Model" => "a"),
                ),
            ),
            "Discretization" =>
                Dict{Any,Any}("Input Mesh File" => "test", "Type" => "test"),
            "Blocks" => Dict{Any,Any}(
                "Block_1" => Dict{Any,Any}(
                    "Block Names" => "Block_1",
                    "Density" => 1.0,
                    "Horizon" => 1.0,
                ),
            ),
            "Solver" => Dict{Any,Any}("Final Time" => 1.0, "Initial Time" => 0.0),
        ),
    )
    @test PeriLab.Solver.Parameter_Handling.validate_yaml(params) == params["PeriLab"]
end

@testset "ut_get_external_topology_name" begin
    params = Dict("Discretization" => Dict())
    @test isnothing(PeriLab.Solver.Parameter_Handling.get_external_topology_name(params))
    params = Dict("Discretization" => Dict("Input External Topology" => Dict()))
    @test isnothing(PeriLab.Solver.Parameter_Handling.get_external_topology_name(params))
    name = randstring(12)
    params =
        Dict("Discretization" => Dict("Input External Topology" => Dict("File" => name)))
    @test PeriLab.Solver.Parameter_Handling.get_external_topology_name(params) == name
end
@testset "ut_get_mesh_name" begin
    params = Dict("Discretization" => Dict())
    @test PeriLab.Solver.Parameter_Handling.get_mesh_name(params) === nothing
    name = randstring(12)
    params = Dict("Discretization" => Dict("Input Mesh File" => name))
    @test PeriLab.Solver.Parameter_Handling.get_mesh_name(params) == name
end

@testset "ut_get_output_filenames" begin
    params = Dict()
    filenames = PeriLab.Solver.Parameter_Handling.get_output_filenames(params, "")
    @test filenames == []
    params = Dict(
        "Outputs" => Dict(
            "Output1" => Dict("Output Filename" => "1"),
            "Output2" => Dict("Output Filename" => "2"),
        ),
    )
    filenames = PeriLab.Solver.Parameter_Handling.get_output_filenames(params, "")
    @test filenames[1] == "1.e"
    @test filenames[2] == "2.e"
    params = Dict(
        "Outputs" => Dict(
            "Output1" => Dict("Output Filename" => "3", "Output File Type" => "CSV"),
            "Output2" => Dict("Output Filename" => "4", "Output File Type" => "Exodus"),
        ),
    )
    filenames = PeriLab.Solver.Parameter_Handling.get_output_filenames(params, "test")
    @test filenames[1] == "test/3.csv"
    @test filenames[2] == "test/4.e"
end

@testset "ut_get_output_fieldnames" begin
    outputs =
        Dict("Displacements" => true, "Forces" => true, "External_Displacements" => true)
    variables = ["DisplacementsNP1", "Forces"]
    computes = ["External_Displacements"]
    output_type = "Exodus"
    fieldnames = PeriLab.Solver.Parameter_Handling.get_output_fieldnames(
        outputs,
        variables,
        computes,
        output_type,
    )
    @test fieldnames == ["DisplacementsNP1", "External_Displacements", "Forces"]

    outputs = Dict("Displacements" => "true")
    @test isnothing(
        PeriLab.Solver.Parameter_Handling.get_output_fieldnames(
            outputs,
            variables,
            computes,
            output_type,
        ),
    )

    outputs = Dict("External_Displacements" => true)
    output_type = "CSV"
    fieldnames = PeriLab.Solver.Parameter_Handling.get_output_fieldnames(
        outputs,
        variables,
        computes,
        output_type,
    )
    @test fieldnames == ["External_Displacements"]

    outputs = Dict("External_Forces" => true)
    fieldnames = PeriLab.Solver.Parameter_Handling.get_output_fieldnames(
        outputs,
        variables,
        computes,
        output_type,
    )
    @test fieldnames == []
end

@testset "get_output_frequency" begin
    nsteps = 40
    params = Dict()
    params = Dict(
        "Outputs" => Dict(
            "Output1" => Dict("Output Frequency" => 2),
            "Output2" => Dict("Number of Output Steps" => 1, "Output Frequency" => 1),
        ),
    )
    freq = PeriLab.Solver.Parameter_Handling.get_output_frequency(params, nsteps)
    @test freq[1] == 2
    @test freq[2] == 40

    params = Dict(
        "Outputs" => Dict(
            "Output1" => Dict("Output Frequency" => 20),
            "Output2" => Dict("Number of Output Steps" => 10),
        ),
    )
    freq = PeriLab.Solver.Parameter_Handling.get_output_frequency(params, nsteps)
    @test freq[1] == 20
    @test freq[2] == 4

    nsteps = 1000
    freq = PeriLab.Solver.Parameter_Handling.get_output_frequency(params, nsteps)
    @test freq[1] == 20
    @test freq[2] == 100

    nsteps = 2
    freq = PeriLab.Solver.Parameter_Handling.get_output_frequency(params, nsteps)
    @test freq[1] == 2
    @test freq[2] == 1

    params = Dict(
        "Outputs" => Dict(
            "Output1" => Dict("Output Frequency" => 20, "Number of Output Steps" => 10),
            "Output2" => Dict("Number of Output Steps" => 10, "Output Frequency" => 20),
        ),
    )
    nsteps = 1000
    freq = PeriLab.Solver.Parameter_Handling.get_output_frequency(params, nsteps)
    @test (freq[1] == 100) || (freq[1] == 20)
    @test (freq[2] == 100) || (freq[2] == 20)
end

test_data_manager = PeriLab.Data_manager
@testset "ut_get_outputs" begin
    test_data_manager.set_num_controller(5)
    test_data_manager.create_constant_node_field("A", Float64, 1)
    test_data_manager.create_node_field("B", Bool, 1)
    test_data_manager.create_constant_node_field("C", Float64, 4)
    test_data_manager.create_node_field("D", Int64, 7)
    test_data_manager.create_node_field("F", Float64, 1)
    test_data_manager.create_constant_node_field("E", Float64, 4)
    testfield_keys = test_data_manager.get_all_field_keys()

    params = Dict(
        "Outputs" => Dict(
            "Output1" => Dict(
                "fieldnames" => [],
                "Output Variables" => Dict("A" => true, "B" => false, "C" => true),
            ),
            "Output2" => Dict(
                "fieldnames" => [],
                "Output Variables" => Dict(
                    "A" => true,
                    "B" => true,
                    "D" => false,
                    "E" => true,
                    "M" => true,
                ),
            ),
        ),
    )

    outputs =
        PeriLab.Solver.Parameter_Handling.get_outputs(params, testfield_keys, String[])

    @test "A" in outputs["Output1"]["fieldnames"]
    @test ("BNP1" in outputs["Output1"]["fieldnames"]) == false
    @test "C" in outputs["Output1"]["fieldnames"]
    @test "A" in outputs["Output2"]["fieldnames"]
    @test "BNP1" in outputs["Output2"]["fieldnames"]
    @test ("D" in outputs["Output2"]["fieldnames"]) == false
    @test "E" in outputs["Output2"]["fieldnames"]
    @test !("M" in outputs["Output2"]["fieldnames"])
    params = Dict(
        "Outputs" => Dict(
            "Output1" => Dict(
                "fieldnames" => [],
                "Output File Type" => "CSV",
                "Output Variables" => Dict("E" => true, "B" => false, "C" => true),
            ),
            "Output2" => Dict(
                "fieldnames" => [],
                "Output Variables" => Dict(
                    "A" => true,
                    "B" => true,
                    "D" => false,
                    "E" => true,
                    "M" => true,
                ),
            ),
        ),
    )
    outputs =
        PeriLab.Solver.Parameter_Handling.get_outputs(params, testfield_keys, String["M"])
    @test !("A" in outputs["Output1"]["fieldnames"])
    @test !("BNP1" in outputs["Output1"]["fieldnames"])
    @test !("C" in outputs["Output1"]["fieldnames"])
    @test ("BNP1" in outputs["Output2"]["fieldnames"])
    @test !("D" in outputs["Output2"]["fieldnames"])
    @test ("E" in outputs["Output2"]["fieldnames"])
    @test ("M" in outputs["Output2"]["fieldnames"])
    params = Dict(
        "Outputs" => Dict(
            "Output1" => Dict(
                "fieldnames" => [],
                "Output File Type" => "CSV",
                "Output Variables" => Dict("M" => true, "A" => true),
            ),
        ),
    )
    outputs =
        PeriLab.Solver.Parameter_Handling.get_outputs(params, testfield_keys, String["M"])
    @test !("A" in outputs["Output1"]["fieldnames"])
    @test "M" in outputs["Output1"]["fieldnames"]

    params = Dict(
        "Outputs" =>
            Dict("Output1" => Dict("fieldnames" => [], "Output Variables" => Dict())),
    )
    outputs =
        PeriLab.Solver.Parameter_Handling.get_outputs(params, testfield_keys, String[])
    @test outputs["Output1"]["fieldnames"] == []
    params = Dict("Outputs" => Dict("Output1" => Dict("fieldnames" => [])))
    outputs =
        PeriLab.Solver.Parameter_Handling.get_outputs(params, testfield_keys, String[])
    @test outputs["Output1"]["fieldnames"] == []
end
@testset "ut_get_computes" begin
    params = Dict()
    testfield_keys = test_data_manager.get_all_field_keys()
    @test PeriLab.Solver.Parameter_Handling.get_computes(params, testfield_keys) == Dict()

    params = Dict(
        "Compute Class Parameters" => Dict(
            "External_Forces" => Dict(
                "Compute Class" => "Block_Data",
                "Calculation Type" => "Sum",
                "Block" => "block_2",
                "Variable" => "A",
            ),
            "External_Displacements" => Dict(
                "Compute Class" => "Block_Data",
                "Calculation Type" => "Maximum",
                "Block" => "block_1",
                "Variable" => "B",
            ),
            "warn_test" => Dict(
                "Compute Class" => "Block_Data",
                "Calculation Type" => "Maximum",
                "Block" => "block_1",
            ),
        ),
    )

    computes = PeriLab.Solver.Parameter_Handling.get_computes(params, testfield_keys)

    @test haskey(computes, "External_Forces")
    @test haskey(computes, "External_Displacements")
    @test !haskey(computes, "warn_test")
    @test computes["External_Forces"]["Variable"] == "A"
    @test computes["External_Displacements"]["Variable"] == "BNP1"
end
@testset "ut_get_computes_names" begin
    testfield_keys = test_data_manager.get_all_field_keys()

    params = Dict(
        "Compute Class Parameters" => Dict(
            "External_Forces" => Dict(
                "Compute Class" => "Block_Data",
                "Calculation Type" => "Sum",
                "Block" => "block_2",
                "Variable" => "A",
            ),
            "External_Displacements" => Dict(
                "Compute Class" => "Block_Data",
                "Calculation Type" => "Maximum",
                "Block" => "block_1",
                "Variable" => "B",
            ),
        ),
    )

    computes_names = PeriLab.Solver.Parameter_Handling.get_computes_names(params)

    @test "External_Forces" in computes_names
    @test "External_Displacements" in computes_names
end

@testset "ut_get_output_variables" begin
    @test PeriLab.Solver.Parameter_Handling.get_output_variables(
        "Force",
        ["Force", "DisplacementsNP1"],
    ) == "Force"
    @test PeriLab.Solver.Parameter_Handling.get_output_variables(
        "Displacements",
        ["Force", "DisplacementsNP1"],
    ) == "DisplacementsNP1"
    @test isnothing(
        PeriLab.Solver.Parameter_Handling.get_output_variables(
            "Stress",
            ["Force", "DisplacementsNP1"],
        ),
    )
end

@testset "ut_get_bc_definitions" begin
    params = Dict()
    bcs = PeriLab.Solver.Parameter_Handling.get_bc_definitions(params)
    @test length(bcs) == 0
    params = Dict("Boundary Conditions" => Dict())
    bcs = PeriLab.Solver.Parameter_Handling.get_bc_definitions(params)
    @test length(bcs) == 0
    params = Dict(
        "Boundary Conditions" => Dict(
            "BC_1" => Dict(
                "Variable" => "Force",
                "Node Set" => "Nset_1",
                "Coordinate" => "x",
                "Value" => "20*t",
            ),
            "BC_2" => Dict(
                "Variable" => "Displacement",
                "Node Set" => "Nset_2",
                "Coordinate" => "y",
                "Value" => "0",
            ),
        ),
    )
    bcs = PeriLab.Solver.Parameter_Handling.get_bc_definitions(params)
    @test length(bcs) == 2
    @test bcs["BC_1"] == Dict(
        "Variable" => "Force",
        "Node Set" => "Nset_1",
        "Coordinate" => "x",
        "Value" => "20*t",
    )
    @test bcs["BC_2"] == Dict(
        "Variable" => "Displacement",
        "Node Set" => "Nset_2",
        "Coordinate" => "y",
        "Value" => "0",
    )
end
@testset "ut_get_solver_options" begin
    params = Dict(
        "Solver" => Dict(
            "Material Models" => true,
            "Damage Models" => true,
            "Additive Models" => true,
            "Thermal Models" => true,
            "Corrosion Models" => true,
        ),
    )
    solver_options = PeriLab.Solver.Parameter_Handling.get_model_options(params)
    println()
    @test solver_options ==
          ["Additive", "Pre_Calculation", "Damage", "Thermal", "Corrosion", "Material"]
    params = Dict("Solver" => Dict())
    solver_options = PeriLab.Solver.Parameter_Handling.get_model_options(params)
    @test solver_options == ["Pre_Calculation", "Material"]
    params = Dict(
        "Solver" => Dict(
            "Material Models" => false,
            "Damage Models" => true,
            "Thermal Models" => true,
        ),
    )

    solver_options = PeriLab.Solver.Parameter_Handling.get_model_options(params)
    @test solver_options == ["Pre_Calculation", "Damage", "Thermal"]
end
@testset "ut_get_number_of_blocks" begin
    @test isnothing(PeriLab.Solver.Parameter_Handling.get_number_of_blocks(Dict()))
    params = Dict("Blocks" => Dict())
    @test isnothing(PeriLab.Solver.Parameter_Handling.get_number_of_blocks(params))
    params = Dict("Blocks" => Dict("block_1" => Dict(), "block_2" => Dict()))
    @test PeriLab.Solver.Parameter_Handling.get_number_of_blocks(params) == 2
    params = Dict(
        "Blocks" => Dict("block_1" => Dict(), "block_2" => Dict(), "block_3" => Dict()),
    )
    @test PeriLab.Solver.Parameter_Handling.get_number_of_blocks(params) == 3
    params = Dict(
        "Blocks" => Dict(
            "block_1" => Dict(),
            "block_2" => Dict(),
            "block_3" => Dict(),
            "block_4" => Dict(),
        ),
    )
    @test PeriLab.Solver.Parameter_Handling.get_number_of_blocks(params) == 4
end

@testset "ut_block_values" begin
    params = Dict("Blocks" => Dict())
    @test isnothing(PeriLab.Solver.Parameter_Handling.get_horizon(params, 1))
    @test isnothing(PeriLab.Solver.Parameter_Handling.get_density(params, 1))
    @test isnothing(PeriLab.Solver.Parameter_Handling.get_heat_capacity(params, 1))
    @test isnothing(PeriLab.Solver.Parameter_Handling.get_values(params, 1, "Density"))
    @test isnothing(PeriLab.Solver.Parameter_Handling.get_values(params, 1, "not there"))
    params = Dict("Blocks" => Dict("block_1" => Dict(), "block_2" => Dict()))
    @test isnothing(PeriLab.Solver.Parameter_Handling.get_horizon(params, 1))
    @test isnothing(PeriLab.Solver.Parameter_Handling.get_density(params, 1))
    @test isnothing(PeriLab.Solver.Parameter_Handling.get_heat_capacity(params, 1))
    @test isnothing(PeriLab.Solver.Parameter_Handling.get_values(params, 1, "Density"))
    @test isnothing(PeriLab.Solver.Parameter_Handling.get_horizon(params, 2))
    @test isnothing(PeriLab.Solver.Parameter_Handling.get_density(params, 2))
    @test isnothing(PeriLab.Solver.Parameter_Handling.get_heat_capacity(params, 2))
    @test isnothing(PeriLab.Solver.Parameter_Handling.get_values(params, 2, "Density"))
    params = Dict(
        "Blocks" => Dict(
            "block_1" => Dict("Density" => 1, "Specific Heat Capacity" => 3),
            "block_2" => Dict("Density" => 12.3, "Horizon" => 2),
        ),
    )
    @test PeriLab.Solver.Parameter_Handling.get_values(params, 1, "Density") == 1
    @test PeriLab.Solver.Parameter_Handling.get_values(params, 2, "Density") == 12.3
    @test isnothing(PeriLab.Solver.Parameter_Handling.get_values(params, 3, "Density"))
    @test PeriLab.Solver.Parameter_Handling.get_values(
        params,
        1,
        "Specific Heat Capacity",
    ) == 3
    @test isnothing(
        PeriLab.Solver.Parameter_Handling.get_values(params, 2, "Specific Heat Capacity"),
    )
    @test isnothing(
        PeriLab.Solver.Parameter_Handling.get_values(params, 3, "Specific Heat Capacity"),
    )
    @test isnothing(PeriLab.Solver.Parameter_Handling.get_values(params, 1, "Horizon"))
    @test PeriLab.Solver.Parameter_Handling.get_values(params, 2, "Horizon") == 2
    @test isnothing(PeriLab.Solver.Parameter_Handling.get_values(params, 3, "Horizon"))
    @test PeriLab.Solver.Parameter_Handling.get_values(params, 1, "Density") ==
          PeriLab.Solver.Parameter_Handling.get_density(params, 1)
    @test PeriLab.Solver.Parameter_Handling.get_values(params, 2, "Density") ==
          PeriLab.Solver.Parameter_Handling.get_density(params, 2)
    @test PeriLab.Solver.Parameter_Handling.get_values(params, 3, "Density") ==
          PeriLab.Solver.Parameter_Handling.get_density(params, 3)

    @test PeriLab.Solver.Parameter_Handling.get_values(params, 1, "Horizon") ==
          PeriLab.Solver.Parameter_Handling.get_horizon(params, 1)
    @test PeriLab.Solver.Parameter_Handling.get_values(params, 2, "Horizon") ==
          PeriLab.Solver.Parameter_Handling.get_horizon(params, 2)
    @test PeriLab.Solver.Parameter_Handling.get_values(params, 3, "Horizon") ==
          PeriLab.Solver.Parameter_Handling.get_horizon(params, 3)

    @test PeriLab.Solver.Parameter_Handling.get_values(
        params,
        1,
        "Specific Heat Capacity",
    ) == PeriLab.Solver.Parameter_Handling.get_heat_capacity(params, 1)
    @test PeriLab.Solver.Parameter_Handling.get_values(
        params,
        2,
        "Specific Heat Capacity",
    ) == PeriLab.Solver.Parameter_Handling.get_heat_capacity(params, 2)
    @test PeriLab.Solver.Parameter_Handling.get_values(
        params,
        3,
        "Specific Heat Capacity",
    ) == PeriLab.Solver.Parameter_Handling.get_heat_capacity(params, 3)
end

@testset "ut_solver" begin
    params = Dict(
        "Solver" => Dict(
            "Initial Time" => 0.0,
            "Final Time" => 1.0,
            "Verlet" => Dict("Safety Factor" => 0.95, "Fixed dt" => 1e-3),
            "Numerical Damping" => 5e-6,
        ),
    )
    @test PeriLab.Solver.Parameter_Handling.get_solver_name(params) == "Verlet"
    @test PeriLab.Solver.Parameter_Handling.get_final_time(params) ==
          params["Solver"]["Final Time"]
    @test PeriLab.Solver.Parameter_Handling.get_initial_time(params) ==
          params["Solver"]["Initial Time"]
    @test PeriLab.Solver.Parameter_Handling.get_safety_factor(params) ==
          params["Solver"]["Verlet"]["Safety Factor"]
    @test PeriLab.Solver.Parameter_Handling.get_fixed_dt(params) ==
          params["Solver"]["Verlet"]["Fixed dt"]
    @test PeriLab.Solver.Parameter_Handling.get_numerical_damping(params) ==
          params["Solver"]["Numerical Damping"]
    params = Dict("Solver" => Dict("Verlet" => Dict()))
    @test PeriLab.Solver.Parameter_Handling.get_safety_factor(params) == 1
    @test PeriLab.Solver.Parameter_Handling.get_fixed_dt(params) == -1.0
    @test PeriLab.Solver.Parameter_Handling.get_nsteps(params) == 1
    @test PeriLab.Solver.Parameter_Handling.get_nsteps(
        Dict(
            "Solver" => Dict(
                "Verlet" => Dict("Safety Factor" => 0.95, "Number of Steps" => 6),
                "Numerical Damping" => 5e-6,
            ),
        ),
    ) == 6
    @test PeriLab.Solver.Parameter_Handling.get_numerical_damping(params) == 0.0
    @test isnothing(
        PeriLab.Solver.Parameter_Handling.get_initial_time(Dict("Solver" => Dict())),
    )
    @test isnothing(
        PeriLab.Solver.Parameter_Handling.get_final_time(Dict("Solver" => Dict())),
    )
    @test isnothing(
        PeriLab.Solver.Parameter_Handling.get_final_time(Dict("Solver" => Dict())),
    )
    @test isnothing(
        PeriLab.Solver.Parameter_Handling.get_solver_name(
            Dict("Solver" => Dict("Solvername" => Dict())),
        ),
    )
    params = Dict(
        "Solver" => Dict(
            "Initial Time" => 0.0,
            "Final Time" => 1.0,
            "External" => Dict("Safety Factor" => 0.95, "Fixed dt" => 1e-3),
            "Numerical Damping" => 5e-6,
        ),
    )
    @test PeriLab.Solver.Parameter_Handling.get_solver_name(params) == "External"
end

path = "./test/unit_tests/Support/Parameters/"
if !isfile(path * "test_data_file.txt")
    path = "./unit_tests/Support/Parameters/"
end
params = Dict(
    "Models" => Dict(
        "Material Models" => Dict(
            "A" => Dict(
                "s" => 0,
                "d" => true,
                "test_file" => path * "test_data_file.txt",
                "test_file_2" => path * "test_data_file.txt",
            ),
            "B" => Dict(
                "sa" => [3.2, 2, 3],
                "d" => "true",
                "test_file_B" => path * "test_data_file.txt",
            ),
        ),
        "Damage Models" => Dict("E" => Dict("ss" => 0, "d" => 1.1)),
    ),
    "Blocks" => Dict(
        "block_1" => Dict("Material Model" => "A", "Damage Model" => "E"),
        "block_2" => Dict("Material Model" => "B"),
    ),
)

@testset "ut_find_data_files" begin
    @test sort(
        PeriLab.Solver.Parameter_Handling.find_data_files(
            params["Models"]["Material Models"]["A"],
        ),
    ) == ["test_file", "test_file_2"]
    @test PeriLab.Solver.Parameter_Handling.find_data_files(
        params["Models"]["Material Models"]["B"],
    ) == ["test_file_B"]
    @test PeriLab.Solver.Parameter_Handling.find_data_files(
        params["Models"]["Damage Models"]["E"],
    ) == []
end
@testset "ut_get_model_parameter" begin
    block_models = Dict{Int32,Dict{String,String}}()
    for id = 1:2
        block_models[id] = PeriLab.Solver.Parameter_Handling.get_block_models(params, id)
    end
    @test block_models[1]["Material Model"] == "A"
    @test block_models[1]["Damage Model"] == "E"
    @test block_models[2]["Material Model"] == "B"
    testData = Dict("Material Model" => Dict(), "Damage Model" => Dict())
    @test isnothing(
        PeriLab.Solver.Parameter_Handling.get_model_parameter(
            params,
            "Does not exist Model",
            block_models[1]["Material Model"],
        ),
    )
    @test isnothing(
        PeriLab.Solver.Parameter_Handling.get_model_parameter(
            params,
            "Does not exist Model",
            "s",
        ),
    )
    @test isnothing(
        PeriLab.Solver.Parameter_Handling.get_model_parameter(
            params,
            "Material Model",
            "s",
        ),
    )
    testData["Material Model"] = PeriLab.Solver.Parameter_Handling.get_model_parameter(
        params,
        "Material Model",
        block_models[1]["Material Model"],
    )
    @test testData["Material Model"]["s"] == 0
    @test testData["Material Model"]["d"] == true
    testData["Material Model"] = PeriLab.Solver.Parameter_Handling.get_model_parameter(
        params,
        "Material Model",
        block_models[2]["Material Model"],
    )
    @test testData["Material Model"]["sa"] == [3.2, 2, 3]
    @test testData["Material Model"]["d"] == "true"
    test_dict = testData["Material Model"]["test_file_B"]
    @test test_dict["max"] == 0.5
    @test test_dict["min"] == -2.5

    @test typeof(test_dict["spl"]) == Dierckx.Spline1D
    testData["Damage Model"] = PeriLab.Solver.Parameter_Handling.get_model_parameter(
        params,
        "Damage Model",
        block_models[1]["Damage Model"],
    )
    @test testData["Damage Model"]["ss"] == 0
    @test testData["Damage Model"]["d"] == 1.1
end

@testset "ut_check_for_duplicates" begin
    @test !(PeriLab.Solver.Parameter_Handling.check_for_duplicates(["a", "b", "c"]))
    @test isnothing(
        PeriLab.Solver.Parameter_Handling.check_for_duplicates(["a", "b", "c", "a"]),
    )
end

@testset "ut_validate_structure_recursive" begin
    expected_structure = Dict(
        "PeriLab" => [
            Dict{Any,Any}(
                "Blocks" => [
                    Dict{Any,Any}(
                        "Any" => [
                            Dict{Any,Any}(
                                "Density" => [Union{Float64,Int64}, true],
                                "Material Model" => [String, false],
                            ),
                            true,
                        ],
                    ),
                    true,
                ],
            ),
            true,
        ],
    )
    params = Dict(
        "PeriLab" => Dict{Any,Any}(
            "Blocks" => Dict{Any,Any}(
                "Any" => Dict{Any,Any}("Density" => 2.1, "Material Model" => "Test"),
            ),
        ),
    )
    validate = true
    checked_keys = []
    validate, checked_keys = PeriLab.Solver.Parameter_Handling.validate_structure_recursive(
        expected_structure,
        params,
        validate,
        checked_keys,
    )
    @test validate
    @test checked_keys == ["PeriLab", "Blocks", "Material Model", "Density", "Any"]

    params = Dict(
        "PeriLab" => Dict{Any,Any}(
            "Blocks" => Dict{Any,Any}(
                "Any" => Dict{Any,Any}("Density" => "Test", "Material Model" => "Test"),
            ),
        ),
    )
    validate, checked_keys = PeriLab.Solver.Parameter_Handling.validate_structure_recursive(
        expected_structure,
        params,
        validate,
        checked_keys,
    )
    @test !validate

    params = Dict(
        "PeriLab" => Dict{Any,Any}(
            "Blocks" =>
                Dict{Any,Any}("Any" => Dict{Any,Any}("Material Model" => "Test")),
        ),
    )
    validate, checked_keys = PeriLab.Solver.Parameter_Handling.validate_structure_recursive(
        expected_structure,
        params,
        validate,
        checked_keys,
    )
    @test !validate
end
@testset "ut_node_sets" begin

    filename = "test.txt"
    numbers = [11, 12, 13, 44, 125]
    lenNumbers = length(numbers)
    params = Dict("Discretization" => Dict())

    @test PeriLab.Solver.Parameter_Handling.get_node_sets(params, "") == Dict{String,Any}()
    params = Dict(
        "Discretization" => Dict(
            "Node Sets" => Dict("Nset_1" => "1 2 3 4 5 6 7", "Nset_2" => filename),
        ),
    )

    file = open(filename, "w")
    println(file, "header: global_id")
    for number in numbers
        println(file, number)
    end
    close(file)

    nsets = PeriLab.Solver.Parameter_Handling.get_node_sets(params, "")
    @test "Nset_1" in keys(nsets)
    @test "Nset_2" in keys(nsets)
    @test length(nsets["Nset_1"]) == 7
    for i = 1:7
        @test nsets["Nset_1"][i] == i
    end
    @test length(nsets["Nset_2"]) == lenNumbers
    for i = 1:lenNumbers
        @test nsets["Nset_2"][i] == numbers[i]
    end
    rm(filename)

    filename = "test.txt"
    file = open(filename, "w")
    println(file, "header: global_id")
    close(file)
    nsets = PeriLab.Solver.Parameter_Handling.get_node_sets(params, "")
    @test haskey(nsets, "Nset_1")
    @test !haskey(nsets, "Nset_2")
    rm(filename)
    filename = "test.txt"
    file = open(filename, "w")
    close(file)
    nsets = PeriLab.Solver.Parameter_Handling.get_node_sets(params, "")
    @test haskey(nsets, "Nset_1")
    @test !haskey(nsets, "Nset_2")
    rm(filename)

    filename = "example_mesh.g"
    params = Dict(
        "Discretization" => Dict(
            "Type" => "Exodus",
            "Input Mesh File" => filename,
            "Node Sets" => Dict("Nset_1" => "1 2 3 4 5 6 7", "Nset_2" => filename),
        ),
    )
    nsets = PeriLab.Solver.Parameter_Handling.get_node_sets(
        params,
        "unit_tests/Support/Parameters",
    )
    @test "Set-1" in keys(nsets)
    @test "Set-2" in keys(nsets)
    @test "Set-3" in keys(nsets)
    @test length(nsets["Set-1"]) == 297
    @test length(nsets["Set-2"]) == 27
    @test length(nsets["Set-3"]) == 3
end
