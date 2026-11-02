# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const PS = PeriLab.ParameterSpec
const ID = PeriLab.InputDeck

function ut_deck()
    return Dict{String,Any}("Discretization" => Dict{String,Any}("Type" => "Text File",
                                                                 "Input Mesh File" => "m.txt"),
                            "Blocks" => Dict{String,Any}("block_1" => Dict{String,Any}("Block ID" => 1,
                                                                                       "Density" => 1.0,
                                                                                       "Horizon" => 1.0,
                                                                                       "Material Model" => "Mat")),
                            "Models" => Dict{String,Any}("Material Models" => Dict{String,Any}("Mat" => Dict{String,Any}("Material Model" => "PD Solid Elastic",
                                                                                                                         "Bulk Modulus" => 1.0,
                                                                                                                         "Shear Modulus" => 1.0))),
                            "Solver" => Dict{String,Any}("Initial Time" => 0.0,
                                                         "Final Time" => 1.0,
                                                         "Verlet" => Dict{String,Any}()))
end

@testset "minimal deck" begin
    input, ctx = ID.read_input(ut_deck())
    @test isempty(ctx.errors)
    @test input isa ID.PeriLabInput
    @test input.sections.blocks["block_1"].block_id === 1
    @test input.sections.solver.options.safety_factor === 1.0
    @test input.contact === nothing
    @test input.materials["Mat"] isa PS.WithBase
    @test isempty(input.globals)
    @test input.sections.strict_validation
end

@testset "unknown top-level section" begin
    deck = ut_deck()
    deck["Boundary Condition"] = Dict{String,Any}()
    input, ctx = ID.read_input(deck)
    @test input === nothing
    @test ctx.errors[1].path == "\"Boundary Condition\""
    @test ctx.errors[1].message == "unknown key — did you mean \"Boundary Conditions\"?"
    input, ctx = ID.read_input(deck; strict = false)
    @test input isa ID.PeriLabInput
    @test ctx.errors[1].severity == :warning
end

@testset "top-level rules" begin
    deck = ut_deck()
    deck["Blocks"] = Dict{String,Any}()
    delete!(deck, "Solver")
    delete!(deck, "Models")
    input, ctx = ID.read_input(deck)
    msgs = Dict(e.path => e.message for e in ctx.errors)
    @test msgs["Blocks"] == "at least one block is required"
    @test msgs[""] == "\"Solver\" or \"Multistep Solver\" is required"
    @test msgs["Models"] == "missing (required)"
    deck = ut_deck()
    delete!(deck["Solver"], "Initial Time")
    input, ctx = ID.read_input(deck)
    @test ctx.errors[1].path == "Solver.\"Initial Time\""
    @test ctx.errors[1].message == "missing (required for \"Solver\")"
    deck = ut_deck()
    delete!(deck, "Solver")
    deck["Multistep Solver"] = Dict{String,Any}("Step_1" => Dict{String,Any}("Final Time" => 1.0,
                                                                             "Verlet" => Dict{String,Any}()))
    input, ctx = ID.read_input(deck)
    @test ctx.errors[1].path == "\"Multistep Solver\".Step_1.\"Step ID\""
    @test ctx.errors[1].message == "missing (required in a multistep solver step)"
end

@testset "required groups" begin
    deck = ut_deck()
    deck["Multistep Solver"] = Dict{String,Any}("S" => Dict{String,Any}("Step ID" => 1,
                                                                        "Final Time" => 1.0,
                                                                        "Verlet" => Dict{String,Any}()))
    _, ctx = ID.read_input(deck)
    @test only(ctx.errors).path == ""
    @test only(ctx.errors).message == "give only one of \"Solver\" or \"Multistep Solver\""
    delete!(deck, "Solver")
    @test isempty(last(ID.read_input(deck)).errors)
    delete!(deck["Multistep Solver"]["S"], "Final Time")
    _, ctx = ID.read_input(deck)
    @test only(ctx.errors).path == "\"Multistep Solver\".S"
    @test only(ctx.errors).message == "\"Final Time\" or \"Additional Time\" is required"
    deck["Multistep Solver"]["S"]["Additional Time"] = 1.0
    @test isempty(last(ID.read_input(deck)).errors)
    deck = ut_deck()
    deck["Outputs"] = Dict{String,Any}("O" => Dict{String,Any}("Output Filename" => "o",
                                                               "Output Variables" => Dict{String,Any}("Displacements" => true)))
    _, ctx = ID.read_input(deck)
    @test only(ctx.errors).path == "Outputs.O"
    @test only(ctx.errors).message ==
          "\"Output Frequency\" or \"Number of Output Steps\" is required"
    deck["Outputs"]["O"]["Number of Output Steps"] = 10
    @test isempty(last(ID.read_input(deck)).errors)
end

@testset "a text mesh header gives node fields" begin
    mktempdir() do dir
        write(joinpath(dir, "m.txt"), "# comment\nheader: x y block_id volume Bulk_Modulus\n0 0 1 1 1\n")
        deck = ut_deck()
        delete!(deck["Models"]["Material Models"]["Mat"], "Bulk Modulus")
        _, ctx = ID.read_input(deck, dir)
        @test isempty(ctx.errors)
        @test ctx.node_fields == Set(["x", "y", "block_id", "volume", "Bulk_Modulus"])
        _, ctx = ID.read_input(deck)            # no mesh file: no node fields
        @test only(ctx.errors).path == "Models.\"Material Models\".Mat"
        deck["Discretization"]["Type"] = "Exodus"
        _, ctx = ID.read_input(deck)            # nodal variables unknown: not checked
        @test isempty(ctx.errors)
        @test ctx.node_fields === nothing
    end
end

@testset "Contact, Globals and Strict Validation keys" begin
    deck = ut_deck()
    deck["Globals"] = Dict{String,Any}("anything" => 1)
    deck["Strict Validation"] = false
    deck["Contact"] = Dict{String,Any}("C" => Dict{String,Any}("Type" => "Penalty Contact",
                                                               "Contact Radius" => 0.1,
                                                               "Contact Stiffness" => 1.0,
                                                               "Contact Groups" => Dict{String,Any}()))
    input, ctx = ID.read_input(deck)
    @test isempty(ctx.errors)
    @test input.globals == Dict("anything" => 1)
    @test !input.sections.strict_validation
    @test input.contact.models["C"].name == "Penalty Contact"
end

@testset "case variants of Contact and Models are caught" begin
    deck = ut_deck()
    deck["contact"] = Dict{String,Any}()
    input, ctx = ID.read_input(deck; strict = false)
    @test PS.has_errors(ctx)
    @test ctx.errors[1].path == "contact"
    @test ctx.errors[1].message == "unknown key — did you mean \"Contact\"?"
    deck = ut_deck()
    deck["Contacts"] = Dict{String,Any}()
    input, ctx = ID.read_input(deck)
    @test ctx.errors[1].message == "unknown key — did you mean \"Contact\"?"
end
