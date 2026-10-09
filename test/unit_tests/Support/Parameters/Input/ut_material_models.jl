# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

const MMID = PeriLab.InputDeck

function ut_material_deck(materials)
    return Dict{String,Any}("Discretization" => Dict{String,Any}("Type" => "Text File",
                                                                 "Input Mesh File" => "m.txt"),
                            "Blocks" => Dict{String,Any}("block_1" => Dict{String,Any}("Block ID" => 1,
                                                                                       "Density" => 1.0,
                                                                                       "Horizon" => 1.0,
                                                                                       "Material Model" => "Mat")),
                            "Models" => Dict{String,Any}("Material Models" => materials),
                            "Solver" => Dict{String,Any}("Initial Time" => 0.0,
                                                         "Final Time" => 1.0,
                                                         "Verlet" => Dict{String,Any}()))
end

const MM_PATH = "Models.\"Material Models\".Mat"

@testset "material models are parsed into typed structs" begin
    input, ctx = MMID.read_input(ut_material_deck(Dict{String,Any}("Mat" => Dict{String,Any}("Material Model" => "Bond-based Elastic",
                                                                                             "Bulk Modulus" => 2.0e5))))
    @test isempty(ctx.errors)
    @test input.materials["Mat"] isa PeriLab.ParameterSpec.WithBase
    @test input.materials["Mat"].base.bulk_modulus == 2.0e5
    @test !hasfield(PeriLab.InputDeck.PeriLabInput, :models)   # no raw model dict
end

@testset "misspelled key in a composite" begin
    input, ctx = MMID.read_input(ut_material_deck(Dict{String,Any}("Mat" => Dict{String,Any}("Material Model" => "Correspondence Elastic + Correspondence Plastic",
                                                                                             "Bulk Modulus" => 1.0,
                                                                                             "Shear Modulus" => 1.0,
                                                                                             "Yeild Stress" => 5.0))))
    @test input === nothing
    msgs = Dict(e.path => e.message for e in ctx.errors)
    @test msgs["$MM_PATH.\"Yield Stress\""] == "missing (required by Correspondence Plastic)"
    @test msgs["$MM_PATH.\"Yeild Stress\""] == "unknown key — did you mean \"Yield Stress\"?"
end

@testset "unknown material model" begin
    input, ctx = MMID.read_input(ut_material_deck(Dict{String,Any}("Mat" => Dict{String,Any}("Material Model" => "Bond based Elastic"))))
    @test input === nothing
    @test ctx.errors[1].path == "$MM_PATH.\"Material Model\""
    @test ctx.errors[1].message ==
          "model \"Bond based Elastic\" not found — did you mean \"Bond-based Elastic\"?"
end

@testset "malformed Material Models sections" begin
    input, ctx = MMID.read_input(ut_material_deck(Dict{String,Any}("Mat" => 5)))
    @test input === nothing
    @test ctx.errors[1].path == MM_PATH
    @test ctx.errors[1].message == "expected a section of `key: value` entries, got 5"
    deck = ut_material_deck(Dict{String,Any}())
    deck["Models"]["Material Models"] = "x"
    input, ctx = MMID.read_input(deck)
    @test input === nothing
    @test ctx.errors[1].path == "Models.\"Material Models\""
end

@testset "non-strict mode downgrades unknown material keys" begin
    input, ctx = MMID.read_input(ut_material_deck(Dict{String,Any}("Mat" => Dict{String,Any}("Material Model" => "PD Solid Elastic",
                                                                                             "my new parameter" => 1))),
                                 strict = false)
    @test input isa MMID.PeriLabInput
    @test only(ctx.errors).severity == :warning
end
