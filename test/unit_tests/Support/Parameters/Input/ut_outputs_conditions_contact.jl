# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const PS = PeriLab.ParameterSpec
const ID = PeriLab.InputDeck

function ut_occ(T, dict)
    ctx = PS.ParseContext()
    return PS.parse_section(T, dict, "X", ctx), ctx
end

const UT_VARS = Dict{String,Any}("Displacements" => true, "Forces" => false)

@testset "Outputs" begin
    o, ctx = ut_occ(ID.OutputParams,
                    Dict{String,Any}("Output Filename" => "out", "Output Frequency" => 10,
                                     "Output Variables" => UT_VARS))
    @test isempty(ctx.errors)
    @test o.output_file_type == "Exodus" && o.flush_file && !o.write_after_damage
    @test o.start_time === 0.0 && o.end_time === Inf && !o.bond_export
    @test o.output_frequency === 10 && o.number_of_output_steps === nothing
    @test o.output_variables == Dict("Displacements" => true, "Forces" => false)
    o, ctx = ut_occ(ID.OutputParams,
                    Dict{String,Any}("Output Filename" => "out", "Output File Type" => "VTK",
                                     "Output Frequency" => 1, "Output Variables" => UT_VARS))
    @test ctx.errors[1].message == "\"VTK\" is not one of: \"Exodus\", \"CSV\""
end

@testset "output frequency rules" begin
    o, ctx = ut_occ(ID.OutputParams,
                    Dict{String,Any}("Output Filename" => "out",
                                     "Output Variables" => UT_VARS))
    @test ctx.errors[1].path == "X"
    @test ctx.errors[1].message ==
          "\"Output Frequency\" or \"Number of Output Steps\" is required"
    o, ctx = ut_occ(ID.OutputParams,
                    Dict{String,Any}("Output Filename" => "out", "Output Frequency" => 1,
                                     "Number of Output Steps" => "10 100",
                                     "Output Variables" => UT_VARS))
    @test isempty(ctx.errors)
    @test o.number_of_output_steps == "10 100"
end

@testset "boundary condition values" begin
    for (raw, expected) in [(0.0, 0.0), (-25, -25.0), ("0.0", "0.0"), ("0*t", "0*t")]
        bc, ctx = ut_occ(ID.BoundaryConditionParams,
                         Dict{String,Any}("Variable" => "Displacements",
                                          "Node Set" => "Set-1", "Value" => raw,
                                          "Coordinate" => "x", "Type" => "Dirichlet"))
        @test isempty(ctx.errors)
        @test bc.value == expected && typeof(bc.value) == typeof(expected)
    end
    bc, ctx = ut_occ(ID.BoundaryConditionParams,
                     Dict{String,Any}("Variable" => "Displacements", "Node Set" => "Set-1",
                                      "Value" => 1, "Type" => "dirichlet"))
    @test ctx.errors[1].message ==
          "\"dirichlet\" is not one of: \"Initial\", \"Dirichlet\", \"Neumann\""
    bc, ctx = ut_occ(ID.BoundaryConditionParams,
                     Dict{String,Any}("Variable" => "Temperature", "Node Set" => "a+b",
                                      "Value" => 1, "Step ID" => 2))
    @test isempty(ctx.errors) && bc.type === nothing && bc.step_id === 2
end

@testset "Compute Class Parameters" begin
    c, ctx = ut_occ(ID.ComputeClassParams,
                    Dict{String,Any}("Compute Class" => "Block_Data", "Variable" => "Forces",
                                     "Calculation Type" => "Sum", "Block" => "block_1",
                                     "X" => 1.0))
    @test isempty(ctx.errors)
    @test c.compute_class == "Block_Data" && c.x === 1.0 && c.node_set === nothing
end

function ut_contact(raw)
    ctx = PS.ParseContext()
    return ID.parse_contact(raw, "Contact", ctx), ctx
end

const UT_CONTACT_MODEL = Dict{String,Any}("Type" => "Penalty Contact",
                                          "Contact Radius" => 0.005,
                                          "Contact Stiffness" => 1e8,
                                          "Contact Groups" => Dict{String,Any}("Group 1" => Dict{String,Any}("Master Block ID" => 2,
                                                                                                             "Slave Block ID" => 1,
                                                                                                             "Search Radius" => 0.005)))

@testset "Contact" begin
    c, ctx = ut_contact(Dict{String,Any}("Contact_1" => UT_CONTACT_MODEL))
    @test isempty(ctx.errors)
    @test c.globals.global_search_frequency === 1 && c.globals.only_surface_contact_nodes
    m = c.models["Contact_1"]
    @test m.contact_stiffness === 1e8 && m.friction_coefficient === nothing
    @test m.contact_groups["Group 1"].master_block_id === 2
    c, ctx = ut_contact(Dict{String,Any}("Globals" => Dict{String,Any}("Global Search Frequency" => 2,
                                                                        "Only Surface Contact Nodes" => false),
                                         "Contact_1" => UT_CONTACT_MODEL))
    @test isempty(ctx.errors)
    @test c.globals.global_search_frequency === 2 && !c.globals.only_surface_contact_nodes
    @test collect(keys(c.models)) == ["Contact_1"]
    c, ctx = ut_contact(Dict{String,Any}("Globals" => Dict{String,Any}("Global Search Frequncy" => 2),
                                         "Contact_1" => 5))
    @test c === nothing
    msgs = Dict(e.path => e.message for e in ctx.errors)
    @test msgs["Contact.Globals.\"Global Search Frequncy\""] ==
          "unknown key — did you mean \"Global Search Frequency\"?"
    @test msgs["Contact.Contact_1"] == "expected a section of `key: value` entries, got 5"
end
