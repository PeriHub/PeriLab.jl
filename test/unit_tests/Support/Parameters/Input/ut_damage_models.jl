# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

const DMID = PeriLab.InputDeck
const DM_PATH = "Models.\"Damage Models\".Dam"

function ut_damage_deck(damages)
    return Dict{String,Any}("Discretization" => Dict{String,Any}("Type" => "Text File",
                                             "Input Mesh File" => "m.txt"),
                            "Blocks" => Dict{String,Any}("block_1" =>
                                                 Dict{String,Any}("Block ID" => 1,
                                                                  "Density" => 1.0,
                                                                  "Horizon" => 1.0,
                                                                  "Damage Model" => "Dam")),
                            "Models" => Dict{String,Any}("Damage Models" => damages),
                            "Solver" => Dict{String,Any}("Initial Time" => 0.0,
                                             "Final Time" => 1.0,
                                             "Verlet" => Dict{String,Any}()))
end
function ut_damage(entry)
    DMID.read_input(ut_damage_deck(Dict{String,Any}("Dam" => Dict{String,Any}(entry))))
end
ut_messages(ctx) = Dict(e.path => e.message for e in ctx.errors)

@testset "damage models are parsed into typed structs" begin
    input, ctx = ut_damage(Dict("Damage Model" => "Critical Stretch",
                                "Critical Value" => 0.1))
    @test isempty(ctx.errors)
    d = input.damages["Dam"]
    @test d isa PeriLab.ParameterSpec.WithBase
    @test PeriLab.ParameterSpec.value(d.base.critical_value, 1) == 0.1
    @test d.model.compression == "none"
    @test d.base.interblock_damage === nothing && d.base.local_damping === nothing
    input, _ = ut_damage(Dict("Damage Model" => "Critical Energy", "Critical Value" => 2.0,
                              "Thickness" => 0.5, "Compression" => "shear"))
    @test input.damages["Dam"].model.thickness == 0.5
    @test input.damages["Dam"].model.compression == "shear"
    input, _ = ut_damage(Dict("Damage Model" => "Critical Energy Anisotropic",
                              "Critical Value" => 2.0,
                              "Anisotropic Damage" => Dict("Critical Value X" => 1.0,
                                   "Critical Value Y" => 2.0)))
    d = input.damages["Dam"]
    @test d.model.compression == "none"
    @test d.model.thickness == 1.0
    @test d.base.anisotropic_damage.critical_value_y == 2.0
    @test d.base.anisotropic_damage.critical_value_z === nothing
    @test d.base.anisotropic_damage.bond_criterion == "Axis"
    @test d.base.anisotropic_damage.force_filter
    @test d.base.anisotropic_damage.exponent == 1.0
end

@testset "anisotropic damage block of the CT50 paper deck" begin
    input, ctx = ut_damage(Dict("Damage Model" => "Critical Energy Anisotropic",
                                "Critical Value" => 20.3797,
                                "Anisotropic Damage" =>
                                    Dict("Critical Value X" => 34.1385,
                                         "Critical Value Y" => 20.3797,
                                         "Bond Criterion" => "Axis Force",
                                         "Force Filter" => false)))
    @test isempty(ctx.errors)
    a = input.damages["Dam"].base.anisotropic_damage
    @test a.bond_criterion == "Axis Force" && !a.force_filter
    @test a.critical_value_x == 34.1385
    input, ctx = ut_damage(Dict("Damage Model" => "Critical Energy Anisotropic",
                                "Critical Value" => 0.01,
                                "Anisotropic Damage" => Dict("Critical Value X" => 0.02,
                                     "Critical Value Y" => 0.01,
                                     "Bond Criterion" => "Directional Stretch",
                                     "Exponent" => 8.0)))
    @test isempty(ctx.errors)
    @test input.damages["Dam"].base.anisotropic_damage.exponent == 8.0
    _, ctx = ut_damage(Dict("Damage Model" => "Critical Energy Anisotropic",
                            "Critical Value" => 1.0,
                            "Anisotropic Damage" => Dict("Critical Value X" => 1.0,
                                 "Critical Value Y" => 1.0,
                                 "Bond Criterion" => "Axes")))
    @test !isempty(ctx.errors)
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
    _, ctx = ut_damage(Dict("Damage Model" => "Critical Energy", "Critical Value" => 0.1,
                            "Compression" => "partial"))
    @test ut_messages(ctx)["$DM_PATH.Compression"] ==
          "\"partial\" is not one of: \"full\", \"none\", \"shear\""
end

@testset "Only Tension was removed" begin
    message = "\"Only Tension\" was removed; use \"Compression\" (Only Tension: true → Compression: none, false → Compression: full)"
    for strict in (true, false)
        deck = ut_damage_deck(Dict{String,Any}("Dam" =>
                                                   Dict{String,Any}("Damage Model" => "Critical Stretch",
                                                                    "Critical Value" => 0.1,
                                                                    "Only Tension" => true)))
        _, ctx = DMID.read_input(deck; strict = strict)
        @test ut_messages(ctx) == Dict("$DM_PATH.\"Only Tension\"" => message)
    end
end

@testset "interblock damage keys" begin
    input, ctx = ut_damage(Dict("Damage Model" => "Critical Stretch",
                                "Critical Value" => 0.1,
                                "Interblock Damage" =>
                                    Dict("Interblock Critical Value 1_2" => 0.2)))
    @test isempty(ctx.errors)
    @test input.damages["Dam"].base.interblock_damage["Interblock Critical Value 1_2"] ==
          0.2
    _, ctx = ut_damage(Dict("Damage Model" => "Critical Stretch", "Critical Value" => 0.1,
                            "Interblock Damage" => Dict("Interblock Value 1_2" => 0.2)))
    @test ut_messages(ctx)["$DM_PATH.\"Interblock Damage\".\"Interblock Value 1_2\""] ==
          "unknown key — expected \"Interblock Critical Value <block>_<block>\""
    file = joinpath(mktempdir(), "crit.txt")
    write(file, "header: Temperature Critical_Value\n0.0 1.0\n100.0 2.0\n")
    _, ctx = ut_damage(Dict("Damage Model" => "Critical Stretch", "Critical Value" => file,
                            "Interblock Damage" =>
                                Dict("Interblock Critical Value 1_2" => 0.2)))
    @test ut_messages(ctx)["$DM_PATH.\"Critical Value\""] ==
          "must be a number when Interblock Damage is used"
end

@testset "local damping keys are required" begin
    input, ctx = ut_damage(Dict("Damage Model" => "Critical Stretch",
                                "Critical Value" => 0.1,
                                "Local Damping" =>
                                    Dict("Representative Young's modulus" => 70.0e9,
                                         "Damping coefficient" => 0.5)))
    @test isempty(ctx.errors)
    @test input.damages["Dam"].base.local_damping.damping_coefficient == 0.5
    _, ctx = ut_damage(Dict("Damage Model" => "Critical Stretch", "Critical Value" => 0.1,
                            "Local Damping" =>
                                Dict("Representative Young's modulus" => 70.0e9)))
    @test startswith(ut_messages(ctx)["$DM_PATH.\"Local Damping\".\"Damping coefficient\""],
                     "missing")
end
