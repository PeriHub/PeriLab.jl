# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

const DMID = PeriLab.InputDeck
const DM_PATH = "Models.\"Damage Models\".Dam"

function ut_damage_deck(damages)
    return Dict{String,Any}("Discretization" => Dict{String,Any}("Type" => "Text File",
                                                                 "Input Mesh File" => "m.txt"),
                            "Blocks" => Dict{String,Any}("block_1" => Dict{String,Any}("Block ID" => 1,
                                                                                       "Density" => 1.0,
                                                                                       "Horizon" => 1.0,
                                                                                       "Damage Model" => "Dam")),
                            "Models" => Dict{String,Any}("Damage Models" => damages),
                            "Solver" => Dict{String,Any}("Initial Time" => 0.0,
                                                         "Final Time" => 1.0,
                                                         "Verlet" => Dict{String,Any}()))
end
ut_damage(entry) = DMID.read_input(ut_damage_deck(Dict{String,Any}("Dam" => Dict{String,Any}(entry))))
ut_messages(ctx) = Dict(e.path => e.message for e in ctx.errors)

@testset "damage models are parsed into typed structs" begin
    input, ctx = ut_damage(Dict("Damage Model" => "Critical Stretch", "Critical Value" => 0.1))
    @test isempty(ctx.errors)
    d = input.damages["Dam"]
    @test d isa PeriLab.ParameterSpec.WithBase
    @test PeriLab.ParameterSpec.value(d.base.critical_value, 1) == 0.1
    @test d.model.only_tension
    @test d.base.interblock_damage === nothing && d.base.local_damping === nothing
    input, _ = ut_damage(Dict("Damage Model" => "Critical Energy", "Critical Value" => 2.0,
                              "Thickness" => 0.5, "Only Tension" => false))
    @test input.damages["Dam"].model.thickness == 0.5
    @test !input.damages["Dam"].model.only_tension
    input, _ = ut_damage(Dict("Damage Model" => "Critical Energy Anisotropic",
                              "Critical Value" => 2.0,
                              "Anisotropic Damage" => Dict("Critical Value X" => 1.0,
                                                           "Critical Value Y" => 2.0)))
    d = input.damages["Dam"]
    @test !d.model.only_tension                      # legacy default of this model
    @test d.model.thickness == 1.0
    @test d.base.anisotropic_damage.critical_value_y == 2.0
    @test d.base.anisotropic_damage.critical_value_z === nothing
end

@testset "damage input errors" begin
    _, ctx = ut_damage(Dict("Damage Model" => "Critical Stretch"))
    @test ut_messages(ctx)["$DM_PATH.\"Critical Value\""] ==
          "missing (required by every damage model)"
    _, ctx = ut_damage(Dict("Damage Model" => "Critical Stretch", "Critical Value" => 0.1,
                            "Critcal Value" => 0.2))
    @test ut_messages(ctx)["$DM_PATH.\"Critcal Value\""] ==
          "unknown key — did you mean \"Critical Value\"?"
    _, ctx = ut_damage(Dict("Damage Model" => "Critical Strech", "Critical Value" => 0.1))
    @test ut_messages(ctx)["$DM_PATH.\"Damage Model\""] ==
          "model \"Critical Strech\" not found — did you mean \"Critical Stretch\"?"
    _, ctx = ut_damage(Dict("Damage Model" => "Critical Stretch + Critical Energy",
                            "Critical Value" => 0.1))
    @test ut_messages(ctx)["$DM_PATH.\"Damage Model\""] ==
          "damage models cannot be combined with +"
end

@testset "interblock damage keys" begin
    input, ctx = ut_damage(Dict("Damage Model" => "Critical Stretch", "Critical Value" => 0.1,
                                "Interblock Damage" => Dict("Interblock Critical Value 1_2" => 0.2)))
    @test isempty(ctx.errors)
    @test input.damages["Dam"].base.interblock_damage["Interblock Critical Value 1_2"] == 0.2
    _, ctx = ut_damage(Dict("Damage Model" => "Critical Stretch", "Critical Value" => 0.1,
                            "Interblock Damage" => Dict("Interblock Value 1_2" => 0.2)))
    @test ut_messages(ctx)["$DM_PATH.\"Interblock Damage\".\"Interblock Value 1_2\""] ==
          "unknown key — expected \"Interblock Critical Value <block>_<block>\""
    file = joinpath(mktempdir(), "crit.txt")
    write(file, "header: Temperature Critical_Value\n0.0 1.0\n100.0 2.0\n")
    _, ctx = ut_damage(Dict("Damage Model" => "Critical Stretch", "Critical Value" => file,
                            "Interblock Damage" => Dict("Interblock Critical Value 1_2" => 0.2)))
    @test ut_messages(ctx)["$DM_PATH.\"Critical Value\""] ==
          "must be a number when Interblock Damage is used"
end

@testset "local damping keys are required" begin
    input, ctx = ut_damage(Dict("Damage Model" => "Critical Stretch", "Critical Value" => 0.1,
                                "Local Damping" => Dict("Representative Young's modulus" => 70.0e9,
                                                        "Damping coefficient" => 0.5)))
    @test isempty(ctx.errors)
    @test input.damages["Dam"].base.local_damping.damping_coefficient == 0.5
    _, ctx = ut_damage(Dict("Damage Model" => "Critical Stretch", "Critical Value" => 0.1,
                            "Local Damping" => Dict("Representative Young's modulus" => 70.0e9)))
    @test startswith(ut_messages(ctx)["$DM_PATH.\"Local Damping\".\"Damping coefficient\""],
                     "missing")
end
