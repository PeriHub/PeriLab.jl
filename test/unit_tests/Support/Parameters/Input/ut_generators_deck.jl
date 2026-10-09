# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
import JSON3
const GID = PeriLab.InputDeck

const UT_SCHEMA = PeriLab.to_json_schema(GID.PeriLabInput)
const UT_DECK = UT_SCHEMA["properties"]["PeriLab"]

@testset "deck schema" begin
    @test UT_SCHEMA["\$schema"] == "https://json-schema.org/draft/2020-12/schema"
    @test UT_SCHEMA["required"] == ["PeriLab"]
    @test issubset(["Blocks", "Discretization", "Models"], UT_DECK["required"])
    @test haskey(UT_DECK["properties"], "Solver") && haskey(UT_DECK["properties"], "Contact")
    @test UT_DECK["properties"]["Globals"] == Dict("type" => "object")
    @test JSON3.read(JSON3.write(UT_SCHEMA)) isa JSON3.Object       # serialisable
end

@testset "model entries in the schema" begin
    models = UT_DECK["properties"]["Models"]["properties"]
    entry = models["Material Models"]["additionalProperties"]
    @test "Material Model" in entry["required"]
    @test haskey(entry["properties"], "Symmetry")                  # base part
    plastic = only(filter(c -> c["if"]["properties"]["Material Model"]["const"] ==
                               "Correspondence Plastic", entry["allOf"]))
    @test haskey(plastic["then"]["properties"], "Yield Stress")
    @test "Yield Stress" in plastic["then"]["required"]
    umat = only(filter(c -> c["if"]["properties"]["Material Model"]["const"] ==
                            "Correspondence UMAT", entry["allOf"]))
    @test haskey(umat["then"]["patternProperties"], "^Property_\\d+\$")
    @test haskey(models["Damage Models"]["additionalProperties"]["properties"], "Critical Value")
    switches = models["Pre Calculation Global"]
    @test switches["properties"]["Shape Tensor"] == Dict("type" => "boolean")
    @test haskey(switches["properties"], "Bond Associated Deformation Gradient")
    @test switches["additionalProperties"] == false
    @test models["Pre Calculation Models"]["additionalProperties"] == switches
    @test UT_DECK["properties"]["Models"]["additionalProperties"] == false
end

# every key a shipped deck uses at top level, under Models and in blocks is a
# declared property of the schema
function ut_schema_keys(deck)
    keys_ok = String[]
    bad = String[]
    allowed(props, key) = haskey(props, key)
    for k in keys(deck)
        allowed(UT_DECK["properties"], string(k)) || push!(bad, "PeriLab.$k")
    end
    models = get(deck, "Models", Dict())
    if models isa AbstractDict
        for k in keys(models)
            allowed(UT_DECK["properties"]["Models"]["properties"], string(k)) ||
                push!(bad, "Models.$k")
        end
    end
    block_props = UT_DECK["properties"]["Blocks"]["additionalProperties"]["properties"]
    for (name, block) in get(deck, "Blocks", Dict())
        block isa AbstractDict || continue
        for k in keys(block)
            allowed(block_props, string(k)) || push!(bad, "Blocks.$name.$k")
        end
    end
    return bad
end

@testset "shipped decks fit the schema keys" begin
    root = normpath(joinpath(@__DIR__, "..", "..", "..", "..", ".."))
    for dir in ("examples", "test"), (path, _, files) in walkdir(joinpath(root, dir)),
        file in files

        endswith(file, ".yaml") || continue
        raw = PeriLab.IO.read_input(joinpath(path, file))
        (raw isa AbstractDict && haskey(raw, "PeriLab")) || continue
        relpath(joinpath(path, file), root) ==
        "test/fullscale_tests/test_PD_Solid_Elastic/strain_xx_external.yaml" && continue
        bad = ut_schema_keys(raw["PeriLab"])
        isempty(bad) || @error "$(relpath(joinpath(path, file), root)): $bad"
        @test isempty(bad)
    end
end
