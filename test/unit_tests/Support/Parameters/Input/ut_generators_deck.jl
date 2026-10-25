# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
import JSON3
const GID = PeriLab.InputDeck

const UT_SCHEMA = PeriLab.to_json_schema(GID.PeriLabInput)
const UT_DECK = UT_SCHEMA["properties"]["PeriLab"]

# the `allOf` case for `key` = `name` (model cases and conditional keys share allOf)
function ut_case(cases, key, name)
    is_case(c) = get(get(c["if"]["properties"], key, Dict()), "const", nothing) == name
    return only(filter(is_case, cases))
end

@testset "deck schema" begin
    @test UT_SCHEMA["\$schema"] == "https://json-schema.org/draft/2020-12/schema"
    @test UT_SCHEMA["required"] == ["PeriLab"]
    @test issubset(["Blocks", "Discretization", "Models"], UT_DECK["required"])
    @test haskey(UT_DECK["properties"], "Solver") &&
          haskey(UT_DECK["properties"], "Contact")
    @test UT_DECK["properties"]["Globals"] == Dict("type" => "object")
    @test JSON3.read(JSON3.write(UT_SCHEMA)) isa JSON3.Object       # serialisable
end

@testset "model entries in the schema" begin
    models = UT_DECK["properties"]["Models"]["properties"]
    entry = models["Material Models"]["additionalProperties"]
    @test "Material Model" in entry["required"]
    @test haskey(entry["properties"], "Symmetry")                  # base part
    plastic = ut_case(entry["allOf"], "Material Model", "Correspondence Plastic")
    @test haskey(plastic["then"]["properties"], "Yield Stress")
    @test "Yield Stress" in plastic["then"]["required"]
    umat = ut_case(entry["allOf"], "Material Model", "Correspondence UMAT")
    @test haskey(umat["then"]["patternProperties"], "^Property_\\d+\$")
    @test haskey(models["Damage Models"]["additionalProperties"]["properties"],
                 "Critical Value")
    switches = models["Pre Calculation Global"]
    @test switches["properties"]["Shape Tensor"]["type"] == "boolean"
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
    penalty = ut_case(entry["allOf"], "Type", "Penalty Contact")
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
    for file in ut_all_decks()
        rel = relpath(file, UT_REPO_ROOT)
        haskey(UT_DECK_ALLOWLIST, rel) && continue
        raw = PeriLab.IO.read_input(file)
        (raw isa AbstractDict && haskey(raw, "PeriLab")) || continue
        bad = ut_schema_keys(raw["PeriLab"])
        isempty(bad) || @error "$rel: $bad"
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
function ut_is_type(t, x)
    t == "object" ? x isa AbstractDict :
    t == "array" ? x isa AbstractVector :
    t == "string" ? x isa AbstractString :
    t == "boolean" ? x isa Bool :
    t == "null" ? x === nothing :
    t == "integer" ?
    ((x isa Integer && !(x isa Bool)) || (x isa AbstractFloat && isinteger(x))) :
    t == "number" ? (x isa Real && !(x isa Bool)) : false
end

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
                               "Contact Stiffness" => 1.0,
                               "Contact Groups" => Dict{String,Any}())
    @test ut_valid(contact, penalty)
    @test !ut_valid(contact, merge(penalty, Dict{String,Any}("Bogus" => 1)))
end

@testset "shipped decks validate against the schema" begin
    for file in ut_all_decks()
        rel = relpath(file, UT_REPO_ROOT)
        haskey(UT_DECK_ALLOWLIST, rel) && continue
        raw = PeriLab.IO.read_input(file)
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
    names = [[fs.alias
              for fs in PeriLab.ParameterSpec.parameter_spec(GID.PeriLabSections)
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
    @test occursin("Contact Stiffness", penalty) &&
          occursin("Shared contact keys:", penalty)
    @test occursin("# Maximum Damage: .inf", ut_describe("Solver"; template = true))
end

@testset "solver choice in schema, describe and template" begin
    solver = UT_DECK["properties"]["Solver"]
    @test haskey(solver["properties"], "Verlet") && haskey(solver["properties"], "Newmark")
    @test Dict{String,Any}("required" => ["Static"]) in solver["oneOf"]
    text = ut_describe("Solver")
    @test occursin("<solver>", text) && occursin("one of: Linear Static Matrix Based", text)
    @test occursin("Verlet:", text) && occursin("Safety Factor", text)
    @test occursin("# exactly one of: Linear Static Matrix Based, Newmark, Static, Verlet, Verlet Matrix Based",
                   ut_describe("Solver"; template = true))
    sections = read(first(PeriLab.generate_parameter_docs(mktempdir())), String)
    @test occursin("## Solver → Verlet", sections) &&
          occursin("## Solver → Newmark", sections)
    @test occursin("| <solver> | one of: Linear Static Matrix Based", sections)
end

@testset "evaluated keys are marked" begin
    p = UT_DECK["properties"]
    @test p["Boundary Conditions"]["additionalProperties"]["properties"]["Value"]["x-perilab-eval"] == true
    @test p["Compute Class Parameters"]["additionalProperties"]["properties"]["Equation"]["x-perilab-eval"] ==
          true
    d = p["Discretization"]["properties"]
    @test d["Influence Function"]["x-perilab-eval"] == true
    @test d["Node Sets"]["x-perilab-eval"] == true
    @test d["Gcode"]["properties"]["Blocks"]["x-perilab-eval"] == true
    @test !haskey(d["Input Mesh File"], "x-perilab-eval")
    thermal = p["Models"]["properties"]["Thermal Models"]["additionalProperties"]["allOf"]
    ht = ut_case(thermal, "Thermal Model", "Heat Transfer")
    @test ht["then"]["properties"]["Environmental Temperature"]["x-perilab-eval"] == true
end

@testset "print_json_schema" begin
    io = IOBuffer()
    PeriLab.print_json_schema(io)
    schema = JSON3.read(String(take!(io)))
    @test haskey(schema["properties"], "PeriLab")
    solver = schema["properties"]["PeriLab"]["properties"]["Solver"]
    @test !isempty(solver["description"])
    @test !isempty(solver["properties"]["Verlet"]["description"])
end

# the description without the generated "quantity: ..." and "one of: ..." parts
function ut_own_description(s)
    parts = split(get(s, "description", ""), "; ")
    return strip(join(filter(p -> !startswith(p, "quantity: ") && !startswith(p, "one of: "),
                             parts), "; "))
end

# paths of sections without a description and of keys without an own description
function ut_undescribed!(bad, s, path)
    s isa AbstractDict || return bad
    if haskey(s, "properties")
        isempty(strip(get(s, "description", ""))) && push!(bad, path)
        for (key, p) in s["properties"]
            key == "Globals" && continue
            haskey(p, "properties") || isempty(ut_own_description(p)) &&
                                       push!(bad, path * " → " * key)
            ut_undescribed!(bad, p, path * " → " * key)
        end
    end
    ut_undescribed!(bad, get(s, "additionalProperties", nothing), path * " → <entry>")
    ut_undescribed!(bad, get(s, "items", nothing), path * " → <item>")
    for case in get(s, "allOf", Any[])
        haskey(case["then"], "properties") || continue        # conditional required keys
        name = only(values(case["if"]["properties"]))["const"]
        ut_undescribed!(bad, case["then"], path * " → <" * name * ">")
    end
    return bad
end

@testset "every section and key has a description" begin
    bad = unique!(ut_undescribed!(String[], UT_SCHEMA, "PeriLab input deck"))
    isempty(bad) || @info "undescribed in the JSON schema:\n" * join(bad, "\n")
    @test isempty(bad)
    p = UT_DECK["properties"]
    @test !isempty(p["Blocks"]["additionalProperties"]["description"])
    @test !isempty(p["Solver"]["properties"]["Static"]["description"])
    plastic = ut_case(p["Models"]["properties"]["Material Models"]["additionalProperties"]["allOf"],
                      "Material Model", "Correspondence Plastic")
    @test !isempty(plastic["then"]["description"])
    @test occursin(PeriLab.ParameterSpec.params_description(GID.SolverParams),
                   ut_describe("Solver"))
end

# every object schema in `s` (sections, entries, items, allOf branches) with its path
function ut_objects!(found, s, path)
    s isa AbstractDict || return found
    haskey(s, "properties") && push!(found, path => s)
    for (key, p) in get(s, "properties", Dict())
        ut_objects!(found, p, path * " → " * key)
    end
    ut_objects!(found, get(s, "additionalProperties", nothing), path * " → <entry>")
    ut_objects!(found, get(s, "items", nothing), path * " → <item>")
    for (i, case) in enumerate(get(s, "allOf", Any[]))
        ut_objects!(found, case["then"], path * " → <case $i>")
    end
    return found
end

const UT_OBJECTS = ut_objects!(Pair{String,Any}[], UT_SCHEMA, "PeriLab input deck")

@testset "no key is declared only to be ignored" begin
    bad = [path * " → " * key for (path, s) in UT_OBJECTS
           for (key, p) in s["properties"]
           if occursin(r"not used"i, get(p, "description", ""))]
    isempty(bad) || @info "keys described as not used:\n" * join(bad, "\n")
    @test isempty(bad)
end

@testset "x-perilab-when refers to sibling choices" begin
    solvers = PeriLab.ParameterSpec.registered_names(:solver)
    # the values a condition key can take: its enum, or the spellings of a Julia enum
    choices(p) = haskey(p, "enum") ? p["enum"] :
                 haskey(p, "examples") && get(p, "type", "") == "string" ? p["examples"] :
                 nothing
    n = 0
    for (path, s) in UT_OBJECTS, (key, p) in s["properties"]
        haskey(p, "x-perilab-when") || continue
        for (sibling, values) in p["x-perilab-when"]
            n += 1
            allowed = sibling == "\$solver" ? solvers :
                      choices(get(s["properties"], sibling, Dict()))
            ok = allowed !== nothing &&
                 all(v -> any(a -> lowercase(a) == lowercase(v), allowed), values)
            ok ||
                @error "$path → $key: x-perilab-when $sibling => $values, choices $allowed"
            @test ok
        end
    end
    @test n > 40
    solver = UT_DECK["properties"]["Solver"]["properties"]
    @test solver["Number of Steps"]["x-perilab-when"]["\$solver"] ==
          ["Static", "Newmark", "Linear Static Matrix Based", "Verlet Matrix Based"]
    bond_filter = UT_DECK["properties"]["Discretization"]["properties"]["Bond Filters"]["additionalProperties"]
    @test bond_filter["properties"]["Type"]["enum"] == ["Disk", "Rectangular_Plane"]
    @test bond_filter["properties"]["Radius"]["x-perilab-when"] == Dict("Type" => ["Disk"])
    @test bond_filter["properties"]["Side Length"]["x-perilab-when"] ==
          Dict("Type" => ["Rectangular_Plane"])
end

@testset "properties follow the declaration order" begin
    PS = PeriLab.ParameterSpec
    structs = [m.sig.parameters[2].var.ub
               for m in methods(PS.is_params)
               if m.sig.parameters[2] isa UnionAll && m.sig.parameters[2].var.ub !== Any]
    @test length(structs) > 50
    for T in structs
        expected = String[]
        for fs in PS.parameter_spec(T)
            fs.category === nothing ? push!(expected, fs.alias) :
            append!(expected, PS.registered_names(fs.category))
        end
        keys_ = filter(!=("Globals"), collect(keys(PS.json_schema(T)["properties"])))
        keys_ == expected || @error "$T: $keys_ != $expected"
        @test keys_ == expected
    end
    @test collect(keys(UT_DECK["properties"]))[1:4] ==
          ["Discretization", "Blocks", "Models", "Contact"]
    material = UT_DECK["properties"]["Models"]["properties"]["Material Models"]["additionalProperties"]
    @test first(keys(material["properties"])) == "Material Model"
end

@testset "keys required in some location are required there" begin
    # keys an object (or one of its allOf branches) requires
    required(s) = union(get(s, "required", String[]),
                        [k for c in get(s, "allOf", Any[])
                         for k in get(c["then"], "required", String[])])
    bad = String[]
    for (path, s) in UT_OBJECTS, (key, p) in s["properties"]
        occursin(r"required (in|there)", get(p, "description", "")) || continue
        any(o -> haskey(last(o)["properties"], key) && key in required(last(o)),
            UT_OBJECTS) ||
            push!(bad, path * " → " * key)
    end
    isempty(bad) || @info "required in the text only:\n" * join(bad, "\n")
    @test isempty(bad)
    p = UT_DECK["properties"]
    @test p["Solver"]["required"] == ["Final Time", "Initial Time"]
    @test p["Multistep Solver"]["additionalProperties"]["required"] == ["Step ID"]
    bond_filter = p["Discretization"]["properties"]["Bond Filters"]["additionalProperties"]
    rectangle = ut_case(bond_filter["allOf"], "Type", "Rectangular_Plane")
    @test "Side Length" in rectangle["then"]["required"]
    @test ut_valid(bond_filter,
                   Dict{String,Any}("Type" => "Disk", "Normal X" => 0.0, "Normal Y" => 1.0,
                                    "Normal Z" => 0.0, "Center X" => 0.0, "Center Y" => 0.0,
                                    "Center Z" => 0.0, "Radius" => 1.0))
    @test !ut_valid(bond_filter,
                    Dict{String,Any}("Type" => "Rectangular_Plane", "Normal X" => 0.0,
                                     "Normal Y" => 1.0))
    compute = p["Compute Class Parameters"]["additionalProperties"]
    @test !ut_valid(compute,
                    Dict{String,Any}("Compute Class" => "Block_Data",
                                     "Variable" => "Forces"))
    material = p["Models"]["properties"]["Material Models"]["additionalProperties"]
    elastic = Dict{String,Any}("Material Model" => "Correspondence Elastic",
                               "Symmetry" => "orthotropic", "Young's Modulus X" => 1.0)
    @test !ut_valid(material, elastic)
end
