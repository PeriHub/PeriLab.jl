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

@testset "contact models in the schema" begin
    contact = UT_DECK["properties"]["Contact"]
    @test haskey(contact["properties"], "Globals")
    entry = contact["additionalProperties"]
    @test "Type" in entry["required"]
    @test haskey(entry["properties"], "Contact Groups")             # shared keys
    penalty = only(filter(c -> c["if"]["properties"]["Type"]["const"] == "Penalty Contact",
                          entry["allOf"]))
    @test haskey(penalty["then"]["properties"], "Contact Stiffness")
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

ut_describe(name; kw...) = sprint(io -> PeriLab.describe(io, name; kw...))

@testset "describe a model" begin
    text = ut_describe("Correspondence Plastic")
    @test startswith(text, "Correspondence Plastic (material model)")
    @test occursin("Yield Stress", text) && occursin("number or data file", text)
    @test occursin("Shared material keys", text) && occursin("Symmetry", text)
end

@testset "describe a section" begin
    text = ut_describe("Solver")
    @test startswith(text, "Solver (section)")
    @test occursin("Final Time", text)
end

@testset "describe template" begin
    yaml = ut_describe("Correspondence Plastic"; template = true)
    @test occursin("Material Model: \"Correspondence Plastic\"", yaml)
    @test occursin(r"\n  Yield Stress: .*# required", yaml)
    @test occursin(r"\n  # Symmetry:", yaml)                       # optional keys commented
    file = joinpath(mktempdir(), "t.yaml")
    write(file, replace(yaml, r"<[^>]*>" => "1"))
    @test PeriLab.IO.read_input(file) isa AbstractDict
end

@testset "describe unknown name" begin
    @test_logs (:error, r"did you mean \"Correspondence Plastic\"") match_mode=:any @test_throws PeriLab.PeriLabError ut_describe("Correspondence Plastik")
end

@testset "parameter docs" begin
    dir = mktempdir()
    files = PeriLab.generate_parameter_docs(joinpath(dir, "generated"))
    @test basename.(files) == ["input_sections.md", "input_models.md"]
    sections = read(files[1], String)
    @test startswith(sections, "<!-- generated by PeriLab.generate_parameter_docs")
    @test occursin("## Solver", sections) && occursin("| Final Time |", sections)
    @test occursin("## Contact", sections)
    models = read(files[2], String)
    @test occursin("## Material Models", models)
    @test occursin("### Correspondence Plastic", models)
    @test occursin("| Yield Stress | number or data file | required |", models)
    @test occursin("### Shared keys", models)
    @test occursin("## Contact\n", models) && occursin("### Penalty Contact", models)
    @test occursin("| Contact Stiffness |", models)
    @test occursin("## Pre Calculation", models) && occursin("Shape Tensor", models)
    @test occursin("in your consistent unit system", models)
end

# a small JSON Schema validator for the keywords the generator emits
ut_is_type(t, x) = t == "object" ? x isa AbstractDict :
                   t == "array" ? x isa AbstractVector :
                   t == "string" ? x isa AbstractString :
                   t == "boolean" ? x isa Bool :
                   t == "integer" ?
                   ((x isa Integer && !(x isa Bool)) || (x isa AbstractFloat && isinteger(x))) :
                   t == "number" ? (x isa Real && !(x isa Bool)) : false

function ut_valid(schema::AbstractDict, x)::Bool
    if haskey(schema, "oneOf")
        count(s -> ut_valid(s, x), schema["oneOf"]) == 1 || return false
    end
    if haskey(schema, "allOf")
        all(s -> ut_valid(s, x), schema["allOf"]) || return false
    end
    if haskey(schema, "if") && ut_valid(schema["if"], x) && haskey(schema, "then")
        ut_valid(schema["then"], x) || return false
    end
    haskey(schema, "const") && x != schema["const"] && return false
    haskey(schema, "enum") && !(x in schema["enum"]) && return false
    if haskey(schema, "type")
        types = schema["type"] isa AbstractVector ? schema["type"] : [schema["type"]]
        any(t -> ut_is_type(t, x), types) || return false
    end
    if x isa Real && !(x isa Bool)
        haskey(schema, "minimum") && x < schema["minimum"] && return false
        haskey(schema, "maximum") && x > schema["maximum"] && return false
    end
    if x isa AbstractVector && haskey(schema, "items")
        all(v -> ut_valid(schema["items"], v), x) || return false
    end
    if x isa AbstractDict
        props = get(schema, "properties", Dict())
        patterns = get(schema, "patternProperties", Dict())
        all(k -> haskey(x, k), get(schema, "required", String[])) || return false
        for (k, v) in x
            key = string(k)
            matched = false
            if haskey(props, key)
                matched = true
                ut_valid(props[key], v) || return false
            end
            for (p, s) in patterns
                occursin(Regex(p), key) || continue
                matched = true
                ut_valid(s, v) || return false
            end
            matched && continue
            extra = get(schema, "additionalProperties", true)
            extra === false && return false
            extra isa AbstractDict && !ut_valid(extra, v) && return false
        end
    end
    return true
end

@testset "single models reject undeclared keys" begin
    entry = UT_DECK["properties"]["Models"]["properties"]["Material Models"]["additionalProperties"]
    plastic = Dict{String,Any}("Material Model" => "Correspondence Plastic",
                               "Symmetry" => "isotropic", "Bulk Modulus" => 1.0,
                               "Shear Modulus" => 1.0, "Yield Stress" => 1.0)
    @test ut_valid(entry, plastic)
    @test !ut_valid(entry, merge(plastic, Dict{String,Any}("Bogus" => 1)))
    @test !ut_valid(entry, merge(plastic, Dict{String,Any}("Yeild Stress" => 1)))
    umat = Dict{String,Any}("Material Model" => "Correspondence UMAT", "File" => "u.so",
                            "Number of Properties" => 1, "Property_1" => 2.0)
    @test ut_valid(entry, umat)
    composite = merge(plastic,
                      Dict{String,Any}("Material Model" => "Correspondence Elastic + Correspondence Plastic"))
    @test ut_valid(entry, composite)                          # composites: shared keys only
    contact = UT_DECK["properties"]["Contact"]["additionalProperties"]
    penalty = Dict{String,Any}("Type" => "Penalty Contact", "Contact Radius" => 0.1,
                               "Contact Stiffness" => 1.0, "Contact Groups" => Dict{String,Any}())
    @test ut_valid(contact, penalty)
    @test !ut_valid(contact, merge(penalty, Dict{String,Any}("Bogus" => 1)))
end

@testset "shipped decks validate against the schema" begin
    root = normpath(joinpath(@__DIR__, "..", "..", "..", "..", ".."))
    for dir in ("examples", "test"), (path, _, files) in walkdir(joinpath(root, dir)),
        file in files

        endswith(file, ".yaml") || continue
        rel = relpath(joinpath(path, file), root)
        rel == "test/fullscale_tests/test_PD_Solid_Elastic/strain_xx_external.yaml" && continue
        raw = PeriLab.IO.read_input(joinpath(path, file))
        (raw isa AbstractDict && haskey(raw, "PeriLab")) || continue
        ok = ut_valid(UT_SCHEMA, raw)
        ok || @error "$rel does not validate against the schema"
        @test ok
    end
end

ut_template_parses(yaml) = begin
    file = joinpath(mktempdir(), "t.yaml")
    write(file, replace(yaml, r"<[^>]*>" => "1"))
    PeriLab.IO.read_input(file)
end

@testset "templates of every section and model" begin
    names = [[fs.alias for fs in PeriLab.ParameterSpec.parameter_spec(GID.PeriLabSections)
              if PeriLab.ParameterSpec._params_type(fs.type) !== nothing];
             "Contact";
             [first(m) for (_, c, _) in GID.MODEL_SECTIONS
              for m in PeriLab.ParameterSpec.registered_models(c)];
             first.(PeriLab.ParameterSpec.registered_models(:contact))]
    for name in names
        yaml = ut_describe(name; template = true)
        @test ut_template_parses(yaml) isa AbstractDict
    end
    blocks = ut_template_parses(ut_describe("Blocks"; template = true))["Blocks"]
    @test length(blocks) == 1 && haskey(only(values(blocks)), "Block ID")
    contact = ut_template_parses(ut_describe("Contact"; template = true))["Contact"]
    @test haskey(only(values(contact)), "Contact Radius")
    @test only(values(contact))["Type"] == "Penalty Contact"
    penalty = ut_describe("Penalty Contact")
    @test startswith(penalty, "Penalty Contact (contact model)")
    @test occursin("Contact Stiffness", penalty) && occursin("Shared contact keys:", penalty)
    @test occursin("# Maximum Damage: .inf", ut_describe("Solver"; template = true))
end
