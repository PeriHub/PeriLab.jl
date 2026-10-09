# Typed Input Parameters — Phase 5: Generators

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** The registered `@params` declarations generate three things:
- a JSON Schema of the whole input deck (`PeriLab.to_json_schema(PeriLab.InputDeck.PeriLabInput)`);
- Documenter reference pages (`PeriLab.generate_parameter_docs(dir)`, built with the docs);
- a `PeriLab.describe(name; template)` that prints a section's or model's parameters, or a ready-to-fill YAML block.

The stale prose about the deleted `parameter_handling.jl` is replaced.

**Architecture:**
- **ParameterSpec** gains pure per-struct generators, all from `parameter_spec` metadata (type, required, default, min/max, allowed, quantity, description):
  - `json_schema(T)` and `type_schema(T)`;
  - `parameter_rows(T)`, a table of strings;
  - `markdown_table(rows)`.
- **InputDeck** assembles the deck:
  - the fixed sections from `PeriLabSections`;
  - `Contact`;
  - the model sections from the registry (base part plus one `if/then` per registered model);
  - the pre-calculation switches;
  - `Globals`.

  `describe` and `generate_parameter_docs` live in InputDeck as well.
- **What the generators show.** They read the registry at call time, so they show exactly the installed (and licensed) modules (spec §4). PeriLab re-exposes the three entry points (`PeriLab.describe` etc.) without exporting them (`describe` would clash with `DataFrames.describe`).

**Tech Stack:** Julia 1.12, `PeriLab.ParameterSpec`, `PeriLab.InputDeck`, JSON3 (already a dependency), Documenter.

**Spec:** `docs/superpowers/specs/2026-10-03-typed-input-parameters-design.md` §4 Generators:
- "`to_json_schema(PeriLabInput)` — JSON Schema document for PeriHub and editors: types, required, defaults, `minimum` / `maximum`, `enum`, `description`; `Dependent` becomes `number | string (data file path)`; composites and user-named entries map to `additionalProperties`";
- "`generate_parameter_docs(dir)` — Documenter pages under `docs/src/` listing, per section and per model: YAML key, type, required, default, range, quantity ("in your consistent unit system"), description";
- "`describe(name)` — prints a model's or section's parameters; `describe(name; template = true)` prints a commented YAML block".

## Global Constraints

- **No behaviour change of the solver or the input reader.** The full suite stays green.
- **Pure functions.** Generators do not change state and do not read files; only `generate_parameter_docs` writes, and only into the directory it is given.
- **Matching of enum values.** Enum YAML values match case-, space- and punctuation-insensitively (`ParameterSpec._normalize`). The schema therefore lists enum spellings in `description` and does NOT emit a JSON `enum` for `@enum` fields, so it never rejects a spelling the reader accepts. A field's `allowed` list is matched exactly by the reader and becomes a JSON `enum`.
- **`Globals`** is accepted in every object (`check_unknown!` never reports it), so every object schema lists `"Globals": {"type": "object"}`.
- **Generated pages are build artifacts.** They are written to `docs/src/generated/` by `docs/make.jl` and git-ignored. Nothing generated is committed.
- **Memory.** The machine has ~8 GB. Before every full-suite run, kill orphaned MPI ranks (`PeriLab.run`, `hydra_pmi_proxy`, `mpiexec` in `/proc/*/cmdline`) with `kill -9`, and delete stale `/dev/shm/mpich_shm_*` files.
- **Git.**
  - Branch `feature/typed-input-parameters`; one commit per task; never merge.
  - Add files by path. Never `git add test` or `git add .`: the user keeps untracked files in `test/`.
  - Never run `git checkout`/`git stash` in the main tree for comparisons; use `git show <rev>:<path>` or a separate worktree.

## Review Focus

1. **A model-specific key** (e.g. `Yield Stress` of Correspondence Plastic) is allowed by the schema for that model. A key no model or base declares is not allowed for a single model. Pinned in Task 2 (`model entries in the schema`).
2. **A `Dependent` field** accepts both a number and a data-file string in the schema. A number-or-list field accepts both forms. Pinned in Task 1 (`type schemas`).
3. **An enum field written in another spelling** ("bond based" for `BondBased`) is not rejected by the schema. Pinned in Task 1 (`enum fields are not strict in the schema`).
4. **Every shipped deck's keys** (top level, `Models` sections, block keys) are allowed by the schema's `properties`. Pinned in Task 2 (`shipped decks fit the schema keys`).
5. **`describe` of an unknown name** aborts with a suggestion instead of printing nothing. Pinned in Task 3 (`describe unknown name`).

## Shared test commands

Workspace runner `<workspace>/unit.jl` (create it once):

```julia
# usage: JULIA_PROJECT=/home/PeriLab.jl julia unit.jl <path under test/, without .jl> ...
using Test, Logging, MPI
import PeriLab
Logging.disable_logging(Logging.Warn)
const TESTDIR = "/home/PeriLab.jl/test"
cd(TESTDIR)
include(joinpath(TESTDIR, "helper.jl"))
MPI.Init()
@testset "t" begin
    for p in ARGS
        include(joinpath(TESTDIR, p * ".jl"))
    end
end
```

- Single files: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/<path> [...]`
- Full suite (~35 min, background, read the tail): `julia --project=. -e 'using Pkg; Pkg.test()' > <workspace>/suite.log 2>&1`

---

### Task 1: Per-struct generators in ParameterSpec

**Files:**
- Create: `src/Support/Parameters/Spec/generators.jl`, included last in `src/Support/Parameters/Spec/ParameterSpec.jl`
- Modify: `src/Support/Parameters/Spec/registry.jl` (`registered_models`)
- Create: `test/unit_tests/Support/Parameters/Spec/ut_generators.jl`, included where the other Spec tests are included (`grep -n 'ut_convert' test/runtests.jl test/unit_tests/Support/Parameters/Spec/*.jl`)

**Interfaces:**
- Produces, in `ParameterSpec`:
  - `type_schema(T)::Dict{String,Any}`;
  - `json_schema(T)::Dict{String,Any}` for an `@params` struct (also a `UnionAll` one);
  - `type_label(T)::String`;
  - `parameter_rows(T)::Vector{NamedTuple{(:key, :type, :required, :default, :range, :quantity, :description)}}` with `String` fields;
  - `nested_params(T)::Vector{Pair{String,Any}}`: alias ⇒ nested `@params` type, in field order;
  - `markdown_table(rows)::String`;
  - `registered_models(category)::Vector{Pair{String,Any}}`: registered names ⇒ types, sorted by name, unavailable stubs left out.

- [ ] **Step 1: Write the failing tests**

Create `test/unit_tests/Support/Parameters/Spec/ut_generators.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const GPS = PeriLab.ParameterSpec

@enum UTGenMode UTGenFast UTGenSlow

GPS.@params struct UTGenInner
    depth::Int64 = req("Depth"; min = 1)
end

GPS.@params struct UTGenParams
    stiffness::Float64 = req("Stiffness"; min = 0, quantity = :stress,
                             description = "elastic stiffness")
    table::GPS.Dependent = opt("Table Value"; default = 1.0)
    alpha::Union{Float64,Vector{Float64}} = opt("Alpha"; default = 0.5, min = 0)
    kind::String = opt("Kind"; default = "a", allowed = ["a", "b"])
    mode::UTGenMode = opt("Mode"; default = UTGenFast)
    label::Union{Nothing,String} = opt("Label"; default = nothing)
    id_or_name::Union{Int64,String} = opt("Id"; default = 1)
    inner::Union{Nothing,UTGenInner} = opt("Inner"; default = nothing)
    named::Dict{String,UTGenInner} = opt("Named"; default = Dict{String,UTGenInner}())
end

@testset "type schemas" begin
    @test GPS.type_schema(Float64) == Dict("type" => "number")
    @test GPS.type_schema(Union{Nothing,Int64}) == Dict("type" => "integer")
    @test GPS.type_schema(Union{Int64,String})["type"] == ["integer", "string"]
    dep = GPS.type_schema(GPS.Dependent)["oneOf"]
    @test Dict("type" => "number") in dep
    @test any(s -> get(s, "type", "") == "string", dep)
    nl = GPS.type_schema(Union{Float64,Vector{Float64}})["oneOf"]
    @test Dict("type" => "number") in nl
    @test Dict("type" => "array", "items" => Dict("type" => "number")) in nl
    @test GPS.type_schema(Vector{Int64}) ==
          Dict("type" => "array", "items" => Dict("type" => "integer"))
end

@testset "struct schema" begin
    s = GPS.json_schema(UTGenParams)
    p = s["properties"]
    @test s["type"] == "object" && s["additionalProperties"] == false
    @test s["required"] == ["Stiffness"]
    @test p["Stiffness"]["minimum"] == 0.0
    @test occursin("elastic stiffness", p["Stiffness"]["description"])
    @test occursin("stress", p["Stiffness"]["description"])
    @test p["Kind"]["enum"] == ["a", "b"] && p["Kind"]["default"] == "a"
    @test all(b -> get(b, "minimum", 0.0) == 0.0 &&
                   get(get(b, "items", Dict()), "minimum", 0.0) == 0.0, p["Alpha"]["oneOf"])
    @test p["Inner"]["properties"]["Depth"]["minimum"] == 1.0
    @test p["Named"]["additionalProperties"]["required"] == ["Depth"]
    @test p["Globals"] == Dict("type" => "object")
    @test !haskey(p["Label"], "default")                 # `nothing` has no JSON default
end

@testset "enum fields are not strict in the schema" begin
    m = GPS.json_schema(UTGenParams)["properties"]["Mode"]
    @test m["type"] == "string"
    @test !haskey(m, "enum")
    @test occursin("UTGenFast", m["description"]) && occursin("UTGenSlow", m["description"])
    @test m["default"] == "UTGenFast"
end

@testset "parameter rows and markdown" begin
    rows = GPS.parameter_rows(UTGenParams)
    stiff = only(filter(r -> r.key == "Stiffness", rows))
    @test stiff.type == "number" && stiff.required == "required" && stiff.range == "≥ 0"
    @test stiff.quantity == "stress" && stiff.description == "elastic stiffness"
    @test only(filter(r -> r.key == "Table Value", rows)).type == "number or data file"
    @test only(filter(r -> r.key == "Kind", rows)).range == "one of: a, b"
    @test only(filter(r -> r.key == "Label", rows)).default == "—"
    @test first.(GPS.nested_params(UTGenParams)) == ["Inner", "Named"]
    md = GPS.markdown_table(rows)
    @test startswith(md, "| YAML key | Type | Required | Default | Range | Quantity | Description |")
    @test occursin("| Stiffness | number | required |", md)
end

@testset "registered models" begin
    names = first.(GPS.registered_models(:material))
    @test "Correspondence Plastic" in names
    @test issorted(names)
    @test all(T -> GPS.is_params(T), last.(GPS.registered_models(:material)))
end
```

- [ ] **Step 2: Run to verify they fail**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Support/Parameters/Spec/ut_generators`
Expected: `UndefVarError: type_schema` (and the other new names).

- [ ] **Step 3: Implement**

`registry.jl`, below `registered_names`:

```julia
"Registered models of `category` as name => parameter type, sorted by name."
function registered_models(category::Symbol)
    models = get(REGISTRY, category, Dict{String,Any}())
    return sort!([name => T for (name, T) in models if T isa Type]; by = first)
end
```

`generators.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

export json_schema, type_schema, type_label, parameter_rows, nested_params, markdown_table

const _GLOBALS_SCHEMA = Dict{String,Any}("type" => "object")

"""
    type_schema(T)

JSON Schema of a declared field type. `Dependent` is a number or a data file
path; enums are strings (their spellings match loosely, so no JSON `enum`).
"""
function type_schema(T)
    T === Float64 && return Dict{String,Any}("type" => "number")
    T === Int64 && return Dict{String,Any}("type" => "integer")
    T === Bool && return Dict{String,Any}("type" => "boolean")
    T === String && return Dict{String,Any}("type" => "string")
    if T === Dependent
        return Dict{String,Any}("oneOf" => Any[Dict{String,Any}("type" => "number"),
                                               Dict{String,Any}("type" => "string",
                                                                "description" => "data file path")])
    end
    if T === Union{Float64,Vector{Float64}}
        return Dict{String,Any}("oneOf" => Any[Dict{String,Any}("type" => "number"),
                                               Dict{String,Any}("type" => "array",
                                                                "items" => Dict{String,Any}("type" => "number"))])
    end
    if _is_scalar_union(T)
        types = [type_schema(m)["type"] for m in Base.uniontypes(T) if m !== Nothing]
        return Dict{String,Any}("type" => length(types) == 1 ? only(types) : types)
    end
    T isa Union && return type_schema(_nonnothing(T))
    if T isa DataType && T <: Enum
        spellings = vcat(string.(instances(T)), collect(keys(enum_aliases(T))))
        return Dict{String,Any}("type" => "string",
                                "description" => "one of: " * join(spellings, ", ") *
                                                 " (case, spaces and punctuation are ignored)")
    end
    if T isa DataType && T <: Vector
        return Dict{String,Any}("type" => "array", "items" => type_schema(eltype(T)))
    end
    if T isa DataType && T <: Dict
        return Dict{String,Any}("type" => "object", "additionalProperties" => type_schema(valtype(T)))
    end
    is_params(T) && return json_schema(T)
    throw(ArgumentError("type_schema: unsupported type $T"))
end

# min / max on the numbers of a schema (also inside arrays and oneOf branches)
function _bounds!(s::Dict{String,Any}, lo, hi)
    if haskey(s, "oneOf")
        foreach(b -> _bounds!(b, lo, hi), s["oneOf"])
    elseif get(s, "type", "") == "array"
        _bounds!(s["items"], lo, hi)
    elseif get(s, "type", "") in ("number", "integer")
        lo === nothing || (s["minimum"] = lo)
        hi === nothing || (s["maximum"] = hi)
    end
    return s
end

_json_value(x::Union{Real,String}) = x
_json_value(x::Enum) = string(x)
_json_value(x::AbstractVector) = all(v -> _json_value(v) !== nothing, x) ?
                                 Any[_json_value(v) for v in x] : nothing
_json_value(x) = nothing

function _quantity_text(q)
    q === nothing && return ""
    return "quantity: $q (in your consistent unit system)"
end

function field_schema(fs::FieldSpec)
    s = type_schema(fs.type)
    _bounds!(s, fs.min, fs.max)
    fs.allowed === nothing || (s["enum"] = Any[_json_value(v) for v in fs.allowed])
    if !fs.required
        v = _json_value(fs.default)
        v === nothing || (s["default"] = v)
    end
    text = filter(!isempty, [fs.description, get(s, "description", ""), _quantity_text(fs.quantity)])
    isempty(text) || (s["description"] = join(text, "; "))
    return s
end

"""
    json_schema(T)

JSON Schema object of the `@params` struct `T`: its YAML keys with types,
defaults, bounds, allowed values and descriptions; `required` lists the
required keys; indexed keys (`key_patterns`) become `patternProperties`.
"""
function json_schema(T)
    props = Dict{String,Any}("Globals" => copy(_GLOBALS_SCHEMA))
    required = String[]
    for fs in parameter_spec(T)
        props[fs.alias] = field_schema(fs)
        fs.required && push!(required, fs.alias)
    end
    schema = Dict{String,Any}("type" => "object", "properties" => props,
                              "additionalProperties" => false)
    isempty(required) || (schema["required"] = sort!(required))
    patterns = key_patterns(T)
    if !isempty(patterns)
        schema["patternProperties"] = Dict{String,Any}(first(p).pattern => type_schema(last(p))
                                                       for p in patterns)
    end
    return schema
end

"Readable name of a declared field type (docs and `describe`)."
function type_label(T)
    T === Float64 && return "number"
    T === Int64 && return "integer"
    T === Bool && return "true/false"
    T === String && return "text"
    T === Dependent && return "number or data file"
    T === Union{Float64,Vector{Float64}} && return "number or list of numbers"
    if _is_scalar_union(T)
        return join([type_label(m) for m in Base.uniontypes(T) if m !== Nothing], " or ")
    end
    T isa Union && return type_label(_nonnothing(T))
    T isa DataType && T <: Enum && return "one of: " * join(string.(instances(T)), ", ")
    T isa DataType && T <: Vector && return "list of " * type_label(eltype(T))
    T isa DataType && T <: Dict && return "named entries of " * type_label(valtype(T))
    is_params(T) && return "section"
    return string(T)
end

_fmt_value(x::AbstractFloat) = isinteger(x) && abs(x) < 1e15 ? string(Int64(x)) : string(x)
_fmt_value(x) = string(x)

function _range_text(fs::FieldSpec)
    fs.allowed === nothing || return "one of: " * join(_fmt_value.(fs.allowed), ", ")
    fs.min !== nothing && fs.max !== nothing &&
        return _fmt_value(fs.min) * " … " * _fmt_value(fs.max)
    fs.min !== nothing && return "≥ " * _fmt_value(fs.min)
    fs.max !== nothing && return "≤ " * _fmt_value(fs.max)
    return ""
end

function _default_text(fs::FieldSpec)
    fs.required && return ""
    fs.default === nothing && return "—"
    fs.default isa AbstractDict && isempty(fs.default) && return "—"
    return _fmt_value(fs.default)
end

"""
    parameter_rows(T)

One row of strings per YAML key of `T`: key, type, required, default, range,
quantity, description.
"""
function parameter_rows(T)
    return [(key = fs.alias, type = type_label(fs.type),
             required = fs.required ? "required" : "optional",
             default = _default_text(fs), range = _range_text(fs),
             quantity = fs.quantity === nothing ? "" : string(fs.quantity),
             description = fs.description) for fs in parameter_spec(T)]
end

function _params_type(T)
    T isa Union && return _params_type(_nonnothing(T))
    T isa DataType && T <: Dict && return _params_type(valtype(T))
    T isa DataType && T <: Vector && return _params_type(eltype(T))
    return is_params(T) ? T : nothing
end

"The nested `@params` sections of `T` as alias => type, in field order."
function nested_params(T)
    nested = Pair{String,Any}[]
    for fs in parameter_spec(T)
        P = _params_type(fs.type)
        P === nothing || push!(nested, fs.alias => P)
    end
    return nested
end

_md_cell(s::AbstractString) = replace(s, "|" => "\\|", "\n" => " ")

"Markdown table of `parameter_rows`."
function markdown_table(rows)
    lines = ["| YAML key | Type | Required | Default | Range | Quantity | Description |",
             "|---|---|---|---|---|---|---|"]
    for r in rows
        push!(lines,
              "| " * join(_md_cell.([r.key, r.type, r.required, r.default, r.range,
                                     r.quantity, r.description]), " | ") * " |")
    end
    return join(lines, "\n") * "\n"
end
```

Check that `_is_scalar_union`, `_nonnothing`, `enum_aliases`, `key_patterns`, `is_params`, `parameter_spec` and `FieldSpec` are visible in `generators.jl`; they are defined in earlier `ParameterSpec` files. In `ParameterSpec.jl`, add `include("generators.jl")` after the last existing `include`.

If `parameter_spec(T)` needs a concrete type for a parametric `@params` struct, look at how `build` calls it and do the same. Registered types may be `UnionAll`.

- [ ] **Step 4: Run to verify they pass**

Run the Step 2 command.
Expected: all pass.

- [ ] **Step 5: Full suite, then commit**

```bash
git add src/Support/Parameters/Spec/generators.jl src/Support/Parameters/Spec/ParameterSpec.jl src/Support/Parameters/Spec/registry.jl test/unit_tests/Support/Parameters/Spec/ut_generators.jl test/runtests.jl
git commit -m "ParameterSpec generators: JSON Schema, readable rows and Markdown tables per @params struct"
```

Add the file that includes the Spec tests instead of `test/runtests.jl` if that is where the include went.

---

### Task 2: `to_json_schema(PeriLabInput)` — the whole deck

**Files:**
- Create: `src/Support/Parameters/Input/generators.jl`, included last in `InputDeck.jl` (after `input.jl`)
- Modify: `src/Support/Parameters/Input/input.jl` (`MODEL_SECTIONS` constant, also used for `_MODEL_SECTIONS`)
- Modify: `src/PeriLab.jl` (`using .InputDeck: to_json_schema, describe, generate_parameter_docs`, without exporting them; `describe` and `generate_parameter_docs` arrive in Tasks 3–4, so import only `to_json_schema` now)
- Create: `test/unit_tests/Support/Parameters/Input/ut_generators_deck.jl`, added to `input_tests.jl`

**Interfaces:**
- Consumes (Task 1): `json_schema`, `type_schema`, `registered_models`, `registered_names`, `base_model`.
- Produces:
  - `InputDeck.MODEL_SECTIONS::Tuple`, holding `(section, category, name_key)` per model category, e.g. `("Material Models", :material, "Material Model")`;
  - `InputDeck.to_json_schema(::Type{PeriLabInput})::Dict{String,Any}`;
  - `PeriLab.to_json_schema`.

- [ ] **Step 1: Write the failing tests**

Create `ut_generators_deck.jl`:

```julia
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
```

Add `"ut_generators_deck.jl"` to the file list in `input_tests.jl`.

- [ ] **Step 2: Run to verify they fail**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Support/Parameters/Input/ut_generators_deck`
Expected: `UndefVarError: to_json_schema not defined in PeriLab`.

- [ ] **Step 3: Implement**

`input.jl`: replace the `_MODEL_SECTIONS` definition with:

```julia
# model categories: section under `Models`, registry category, key naming the model
const MODEL_SECTIONS = (("Material Models", :material, "Material Model"),
                        ("Damage Models", :damage, "Damage Model"),
                        ("Thermal Models", :thermal, "Thermal Model"),
                        ("Additive Models", :additive, "Additive Model"),
                        ("Degradation Models", :degradation, "Degradation Model"))

# sections of `Models`; any other key is reported like an unknown key
const _MODEL_SECTIONS = Set([first.(MODEL_SECTIONS)...; "Pre Calculation Global";
                             "Pre Calculation Models"])
```

`generators.jl` (InputDeck):

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

export to_json_schema

# a named model entry: the category's base keys, the key naming the model, and
# per registered model (`if` the name is exactly that model) its own keys
function _model_entry_schema(category::Symbol, name_key::String)
    base = ParameterSpec.base_model(category)
    entry = base === nothing ?
            Dict{String,Any}("type" => "object",
                             "properties" => Dict{String,Any}("Globals" => Dict{String,Any}("type" => "object"))) :
            ParameterSpec.json_schema(base)
    delete!(entry, "additionalProperties")      # model keys come from the matching model
    names = ParameterSpec.registered_names(category)
    entry["properties"][name_key] = Dict{String,Any}("type" => "string",
                                                     "description" => "one of: " *
                                                                      join(names, ", ") *
                                                                      "; models can be combined with +")
    entry["required"] = sort!(unique([get(entry, "required", String[]); name_key]))
    cases = Any[]
    for (name, T) in ParameterSpec.registered_models(category)
        model = ParameterSpec.json_schema(T)
        then = Dict{String,Any}("properties" => model["properties"])
        haskey(model, "required") && (then["required"] = model["required"])
        haskey(model, "patternProperties") &&
            (then["patternProperties"] = model["patternProperties"])
        push!(cases,
              Dict{String,Any}("if" => Dict{String,Any}("properties" => Dict{String,Any}(name_key => Dict{String,Any}("const" => name))),
                               "then" => then))
    end
    isempty(cases) || (entry["allOf"] = cases)
    return entry
end

function _pre_calculation_switches_schema()
    props = Dict{String,Any}(name => Dict{String,Any}("type" => "boolean")
                             for name in ParameterSpec.registered_names(:pre_calculation))
    for (old, new) in DEPRECATED_PRE_CALCULATIONS
        props[old] = Dict{String,Any}("type" => "boolean", "const" => false,
                                      "description" => "deprecated, use \"$new\"")
    end
    return Dict{String,Any}("type" => "object", "properties" => props,
                            "additionalProperties" => false)
end

"""
    to_json_schema(PeriLabInput)

JSON Schema (draft 2020-12) of an input deck, generated from the `@params`
declarations and the registered models (installed and licensed modules). Model
entries list the category's shared keys; the keys of a single named model are
added for exactly that name (models combined with `+` are checked only for the
shared keys).
"""
function to_json_schema(::Type{PeriLabInput})
    deck = ParameterSpec.json_schema(PeriLabSections)
    props = deck["properties"]
    models = Dict{String,Any}()
    for (section, category, name_key) in MODEL_SECTIONS
        models[section] = Dict{String,Any}("type" => "object",
                                           "additionalProperties" => _model_entry_schema(category,
                                                                                         name_key))
    end
    switches = _pre_calculation_switches_schema()
    models["Pre Calculation Global"] = switches
    models["Pre Calculation Models"] = Dict{String,Any}("type" => "object",
                                                        "additionalProperties" => switches)
    props["Models"] = Dict{String,Any}("type" => "object", "properties" => models,
                                       "additionalProperties" => false)
    contact = ParameterSpec.json_schema(ContactModelParams)
    props["Contact"] = Dict{String,Any}("type" => "object",
                                        "properties" => Dict{String,Any}("Globals" => ParameterSpec.json_schema(ContactGlobalsParams)),
                                        "additionalProperties" => contact)
    props["Globals"] = Dict{String,Any}("type" => "object")
    deck["required"] = sort!(unique([get(deck, "required", String[]); "Models"]))
    return Dict{String,Any}("\$schema" => "https://json-schema.org/draft/2020-12/schema",
                            "title" => "PeriLab input deck", "type" => "object",
                            "properties" => Dict{String,Any}("PeriLab" => deck),
                            "required" => ["PeriLab"])
end
```

Check that `DEPRECATED_PRE_CALCULATIONS` is defined in `input.jl` before `generators.jl` is included. In `InputDeck.jl`, add `include("generators.jl")` after `include("input.jl")`. In `PeriLab.jl`, add `using .InputDeck: to_json_schema` next to the other `using` lines, without exporting it.

If "shipped decks fit the schema keys" reports a key, it is one of two things:
- the schema misses a declared key (a bug in the generator: fix it);
- the deck uses a key the reader rejects in strict mode. That cannot happen, because the golden decks parse. Investigate before changing anything.

- [ ] **Step 4: Run to verify they pass**

Run the Step 2 command, plus `unit_tests/Support/Parameters/Input/input_tests`.
Expected: all pass.

- [ ] **Step 5: Full suite, then commit**

```bash
git add src/Support/Parameters/Input/generators.jl src/Support/Parameters/Input/InputDeck.jl src/Support/Parameters/Input/input.jl src/PeriLab.jl test/unit_tests/Support/Parameters/Input/ut_generators_deck.jl test/unit_tests/Support/Parameters/Input/input_tests.jl
git commit -m "to_json_schema: JSON Schema of the whole input deck from the typed declarations"
```

---

### Task 3: `describe(name; template)`

**Files:**
- Modify: `src/Support/Parameters/Input/generators.jl` (`describe`)
- Modify: `src/PeriLab.jl` (import `describe`)
- Test: append to `ut_generators_deck.jl`

**Interfaces:**
- Consumes (Tasks 1–2): `parameter_rows`, `nested_params`, `registered_models`, `MODEL_SECTIONS`, `base_model`.
- Produces:
  - `InputDeck.describe(io::IO, name::AbstractString; template::Bool = false)`;
  - `InputDeck.describe(name; template = false)` (prints to `stdout`);
  - `PeriLab.describe`.
- Names `describe` knows:
  - every top-level section alias of `PeriLabSections` ("Discretization", "Blocks", "Solver", "Outputs", …), plus "Contact";
  - every registered model name of the five model categories.

- [ ] **Step 1: Write the failing tests**

Append to `ut_generators_deck.jl`:

```julia
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
    @test PeriLab.IO.YAML.load(replace(yaml, r"<[^>]*>" => "1")) isa AbstractDict
end

@testset "describe unknown name" begin
    @test_logs (:error, r"did you mean \"Correspondence Plastic\"") match_mode=:any @test_throws PeriLab.PeriLabError ut_describe("Correspondence Plastik")
end
```

If `YAML` is not reachable as `PeriLab.IO.YAML` (the IO file imports only `load_file, ParserError`), use `import YAML` in the test only if the test environment has it. Otherwise check with `PeriLab.IO.read_input` on a temporary file holding the template:

```julia
file = joinpath(mktempdir(), "t.yaml")
write(file, replace(yaml, r"<[^>]*>" => "1"))
@test PeriLab.IO.read_input(file) isa AbstractDict
```

- [ ] **Step 2: Run to verify they fail**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Support/Parameters/Input/ut_generators_deck`
Expected: `UndefVarError: describe not defined in PeriLab`.

- [ ] **Step 3: Implement**

Append to the InputDeck `generators.jl`:

```julia
export describe

# (title, struct, extra rows' struct with its title) for a section or model name
function _described(name::AbstractString)
    for fs in ParameterSpec.parameter_spec(PeriLabSections)
        fs.alias == name || continue
        T = ParameterSpec._params_type(fs.type)
        T === nothing || return ("$name (section)", T, nothing, nothing)
    end
    name == "Contact" && return ("Contact (section, one entry per contact model)",
                                 ContactModelParams, nothing, nothing)
    for (section, category, name_key) in MODEL_SECTIONS
        for (model, T) in ParameterSpec.registered_models(category)
            model == name || continue
            kind = lowercase(replace(section, " Models" => ""))
            return ("$name ($kind model)", T, ParameterSpec.base_model(category),
                    (name_key, "Shared $kind keys"))
        end
    end
    candidates = [[fs.alias for fs in ParameterSpec.parameter_spec(PeriLabSections)];
                  "Contact";
                  [first(m) for (_, c, _) in MODEL_SECTIONS
                   for m in ParameterSpec.registered_models(c)]]
    suggestion = ParameterSpec.suggest(name, candidates)
    @abort "unknown section or model \"$name\"" *
           (suggestion === nothing ? "" : " — did you mean \"$suggestion\"?")
end

function _describe_table(io::IO, T, indent::String)
    rows = ParameterSpec.parameter_rows(T)
    isempty(rows) && return println(io, indent, "(no keys)")
    width = maximum(length(r.key) for r in rows)
    for r in rows
        info = filter(!isempty, [r.type, r.required == "required" ? "required" : "default " * r.default,
                                 r.range, r.quantity, r.description])
        println(io, indent, rpad(r.key, width), "  ", join(info, "; "))
    end
    for (alias, P) in ParameterSpec.nested_params(T)
        println(io, indent, alias, ":")
        _describe_table(io, P, indent * "  ")
    end
end

function _template(io::IO, T, indent::String)
    for fs in ParameterSpec.parameter_spec(T)
        P = ParameterSpec._params_type(fs.type)
        prefix = fs.required ? "" : "# "
        if P !== nothing && P === ParameterSpec._nonnothing(fs.type)
            println(io, indent, prefix, fs.alias, ":")
            _template(io, P, indent * "  " * prefix)
            continue
        end
        value = fs.required ? "<" * ParameterSpec.type_label(fs.type) * ">" :
                ParameterSpec._default_text(ParameterSpec.FieldSpec(fs.name, fs.type, fs.alias,
                                                                    fs.required, fs.default,
                                                                    fs.min, fs.max, fs.allowed,
                                                                    fs.quantity, fs.description))
        value == "—" && (value = "<" * ParameterSpec.type_label(fs.type) * ">")
        note = filter(!isempty, [fs.required ? "required" : "optional",
                                 ParameterSpec._range_text(fs),
                                 fs.quantity === nothing ? "" : string(fs.quantity),
                                 fs.description])
        println(io, indent, prefix, fs.alias, ": ", value, "  # ", join(note, "; "))
    end
end

"""
    describe(name; template = false)
    describe(io, name; template = false)

Prints the parameters of a section (e.g. "Solver") or a registered model (e.g.
"Correspondence Plastic"). With `template = true`, prints a YAML block to copy
into an input deck: required keys with a placeholder, optional keys commented
out with their default.
"""
function describe(io::IO, name::AbstractString; template::Bool = false)
    title, T, base, model_key = _described(name)
    if template
        println(io, model_key === nothing ? "$name:" : "My $(lowercase(first(model_key))):")
        model_key === nothing || println(io, "  ", first(model_key), ": \"", name, "\"")
        _template(io, T, "  ")
        base === nothing || _template(io, base, "  ")
        return nothing
    end
    println(io, title)
    _describe_table(io, T, "  ")
    if base !== nothing
        println(io, last(model_key), ":")
        _describe_table(io, base, "  ")
    end
    return nothing
end

describe(name::AbstractString; template::Bool = false) = describe(stdout, name;
                                                                  template = template)
```

Notes:
- `_default_text` and `_range_text` take a `FieldSpec`; `fs` already is one, so call them with `fs` directly and drop the reconstruction shown above if it is unnecessary.
- `@abort` must be available in `InputDeck` (it is: `using ..PeriLabExceptions: @abort`).
- The template line for a model with a base part lists the base's optional keys commented out; that is long for materials (C11…C66) but complete. Keep it.
- In `PeriLab.jl`, extend the `using .InputDeck:` line with `describe`.

- [ ] **Step 4: Run to verify they pass**

Run the Step 2 command.
Expected: all pass. Also run it by hand: `JULIA_PROJECT=/home/PeriLab.jl julia -e 'import PeriLab; PeriLab.describe("Critical Stretch"); PeriLab.describe("Critical Stretch"; template = true)'`. Read the output: it should be readable, and the template should be valid YAML.

- [ ] **Step 5: Full suite, then commit**

```bash
git add src/Support/Parameters/Input/generators.jl src/PeriLab.jl test/unit_tests/Support/Parameters/Input/ut_generators_deck.jl
git commit -m "describe: a section's or model's parameters, or a YAML template"
```

---

### Task 4: `generate_parameter_docs(dir)` and the docs build

**Files:**
- Modify: `src/Support/Parameters/Input/generators.jl` (`generate_parameter_docs`)
- Modify: `src/PeriLab.jl` (import `generate_parameter_docs`)
- Modify: `docs/make.jl`:
  - call `PeriLab.generate_parameter_docs(joinpath(@__DIR__, "src", "generated"))` before `makedocs`;
  - add the pages `"Input Reference" => Any["Sections" => "generated/input_sections.md", "Models" => "generated/input_models.md"]` under "User Guide", after "Input File".
- Modify: `.gitignore` (add `docs/src/generated/`)
- Test: append to `ut_generators_deck.jl`

**Interfaces:**
- Produces:
  - `InputDeck.generate_parameter_docs(dir::AbstractString)::Vector{String}`, which writes `input_sections.md` and `input_models.md` into `dir` (created if missing) and returns their paths;
  - `PeriLab.generate_parameter_docs`.

- [ ] **Step 1: Write the failing test**

```julia
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
    @test occursin("## Pre Calculation", models) && occursin("Shape Tensor", models)
    @test occursin("in your consistent unit system", models)
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Support/Parameters/Input/ut_generators_deck`
Expected: `UndefVarError: generate_parameter_docs not defined in PeriLab`.

- [ ] **Step 3: Implement**

Append to the InputDeck `generators.jl`:

```julia
export generate_parameter_docs

const _DOCS_HEADER = "<!-- generated by PeriLab.generate_parameter_docs; do not edit -->\n\n"
const _QUANTITY_NOTE = "Quantities name what a value measures; PeriLab has no fixed units, use them in your consistent unit system.\n\n"

function _docs_section(io::IO, title::String, T, level::String)
    println(io, level, " ", title, "\n")
    print(io, ParameterSpec.markdown_table(ParameterSpec.parameter_rows(T)), "\n")
    for (alias, P) in ParameterSpec.nested_params(T)
        _docs_section(io, title * " → " * alias, P, level * "#")
    end
end

"""
    generate_parameter_docs(dir)

Writes the input reference pages `input_sections.md` (the deck's sections) and
`input_models.md` (shared keys and every registered model per category, plus the
pre-calculation switches) into `dir` and returns their paths.
"""
function generate_parameter_docs(dir::AbstractString)
    mkpath(dir)
    sections_file = joinpath(dir, "input_sections.md")
    open(sections_file, "w") do io
        print(io, _DOCS_HEADER, "# Input Sections\n\n", _QUANTITY_NOTE)
        for fs in ParameterSpec.parameter_spec(PeriLabSections)
            T = ParameterSpec._params_type(fs.type)
            T === nothing || _docs_section(io, fs.alias, T, "##")
        end
        _docs_section(io, "Contact", ContactModelParams, "##")
        _docs_section(io, "Contact → Globals", ContactGlobalsParams, "##")
    end
    models_file = joinpath(dir, "input_models.md")
    open(models_file, "w") do io
        print(io, _DOCS_HEADER, "# Input Models\n\n", _QUANTITY_NOTE)
        for (section, category, name_key) in MODEL_SECTIONS
            println(io, "## ", section, "\n")
            println(io, "Each entry names its model in `", name_key, "`.\n")
            base = ParameterSpec.base_model(category)
            base === nothing || _docs_section(io, "Shared keys", base, "###")
            for (name, T) in ParameterSpec.registered_models(category)
                _docs_section(io, name, T, "###")
            end
        end
        println(io, "## Pre Calculation\n")
        println(io, "`Pre Calculation Global` and each entry of `Pre Calculation Models` switch these on (`true`) or off (`false`):\n")
        for name in ParameterSpec.registered_names(:pre_calculation)
            println(io, "- ", name)
        end
    end
    return [sections_file, models_file]
end
```

`_docs_section` puts its own table under its own heading: a nested section gets a heading one level deeper.

`docs/make.jl`: before `makedocs(`, add

```julia
# input reference pages from the typed parameter declarations
PeriLab.generate_parameter_docs(joinpath(@__DIR__, "src", "generated"))
```

Add the pages as listed above.

`.gitignore`: add a line `docs/src/generated/`.

If Documenter is installed in this environment (`ls docs/Manifest.toml` or `julia --project=docs -e 'using Documenter'`), build the docs once: `julia --project=docs docs/make.jl > <workspace>/docs.log 2>&1`, and check the tail for errors. If it is not installed (no network), ledger that the docs build was not run and rely on the generator test.

- [ ] **Step 4: Run to verify it passes**

Run the Step 2 command.
Expected: all pass.

- [ ] **Step 5: Full suite, then commit**

```bash
git add src/Support/Parameters/Input/generators.jl src/PeriLab.jl docs/make.jl .gitignore test/unit_tests/Support/Parameters/Input/ut_generators_deck.jl
git commit -m "generate_parameter_docs: input reference pages built with the docs"
```

---

### Task 5: Stale prose and `read_input_deck` returns only the typed input

**Files:**
- Modify: `src/IO/read_inputdeck.jl`:
  - `read_input_deck` returns only the `PeriLabInput`;
  - the docstring says so;
  - `validate_input` keeps returning `(deck, input)`, because its tests use the deck.
- Modify: `src/IO/IO.jl` (`input = read_input_deck(...)` in `initialize_data`)
- Modify tests: `test/unit_tests/IO/ut_read_inputdeck.jl` and `test/unit_tests/Support/Parameters/Input/ut_validate_yaml.jl` (callers of `read_input_deck`)
- Modify:
  - `README.md:154`: link the generated input reference (the "Input Reference" docs pages) instead of `parameter_handling.jl`;
  - `CONTRIBUTING.md:58`: "Declare new input values as `@params` fields of your module (see the templates and `docs/src/man/dev/parameters.md`)";
  - `docs/src/man/dev/parameters.md`: rewrite (content below);
  - `docs/src/man/dev/module_overview.md`: replace `Parameter_Handling` with `ParameterSpec` and `InputDeck`.

**Interfaces:**
- Produces: `IO.read_input_deck(filename; directory, no_strict)::PeriLabInput`.

- [ ] **Step 1: Write the failing test**

In `ut_read_inputdeck.jl`, the success case becomes:

```julia
    input = PeriLab.IO.read_input_deck(filename; no_strict = true)  # Models keys d, a are placeholders
    @test input isa PeriLab.InputDeck.PeriLabInput
    @test input.sections.discretization.input_mesh_file == "test"
    @test input.sections.solver.final_time == 1.0
```

Keep the other assertions only if they read the typed input. The removed raw-dict assertions (`dict["Models"]["d"] == 3`, …) tested YAML loading. That is still tested by `PeriLab.IO.read_input`.

In `ut_validate_yaml.jl`, a `deck, input = PeriLab.IO.read_input_deck(file)` becomes `input = PeriLab.IO.read_input_deck(file)`. Drop or rewrite any assertion on `deck`.

- [ ] **Step 2: Run to verify it fails**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/IO/ut_read_inputdeck unit_tests/Support/Parameters/Input/ut_validate_yaml`
Expected: failures: `read_input_deck` still returns a tuple.

- [ ] **Step 3: Implement**

`read_input_deck`: `return last(validate_input(read_input(filename); directory = directory, no_strict = no_strict))`. Docstring: "Reads and validates the input deck; returns the typed `PeriLabInput`." In `IO.initialize_data`, `_, input = read_input_deck(...)` becomes `input = read_input_deck(...)`.

`docs/src/man/dev/parameters.md`, new content. Keep the file's existing license header lines, if any, and keep its CRLF/LF line endings: check with `file docs/src/man/dev/parameters.md`.

```markdown
# Parameters

The input deck is read into typed structs. Every section and every model declares
its YAML keys as an `@params` struct (`src/Support/Parameters/Spec`): type,
required or default, `min` / `max`, allowed values, a quantity (documentation only;
PeriLab has no fixed units) and a description.

`PeriLab.InputDeck.read_input` validates a deck against these declarations and
reports every problem at once (unknown keys with a suggestion, wrong types, values
out of range). With `Strict Validation: false` (or `--no_strict`) unknown keys are
warnings instead of errors.

A model module declares its own keys and registers them under its model name, e.g.

    @params struct MyMaterialParams
        yield_stress::Dependent = req("Yield Stress"; min = 0, quantity = :stress)
    end
    __init__() = register_material("My Material", MyMaterialParams)

See the templates under `src/Models/*/…_template` for every category.

The [input reference](@ref "Input Sections") lists every section and model with its
keys; it is generated from the declarations (`PeriLab.generate_parameter_docs`). In a
Julia session, `PeriLab.describe("Correspondence Plastic")` prints the keys of one
model and `PeriLab.describe("Correspondence Plastic"; template = true)` a YAML block
to start from. `PeriLab.to_json_schema(PeriLab.InputDeck.PeriLabInput)` returns a
JSON Schema of the whole deck for editors.

!!! note "Good start"
    Please check some of the full scale tests. There are several yaml files with parameter definitions.
```

The `@ref "Input Sections"` must match the generated page's first heading (`# Input Sections`). If the docs are not built in this environment, use a plain relative link `[input reference](../../generated/input_sections.md)` instead and ledger it.

Run `grep -rn 'parameter_handling\|Parameter_Handling\|expected_structure\|validate_yaml' README.md CONTRIBUTING.md docs/src src --include=*.md --include=*.jl`.
Expected: no matches, except historical plan/spec files under `docs/superpowers/`; leave those as they are.

- [ ] **Step 4: Run to verify it passes**

Run the Step 2 command, plus `unit_tests/ut_docs_references`.
Expected: all pass.

- [ ] **Step 5: Full suite, then commit**

```bash
git add src/IO/read_inputdeck.jl src/IO/IO.jl test/unit_tests/IO/ut_read_inputdeck.jl test/unit_tests/Support/Parameters/Input/ut_validate_yaml.jl README.md CONTRIBUTING.md docs/src/man/dev/parameters.md docs/src/man/dev/module_overview.md
git commit -m "Docs describe the typed input; read_input_deck returns the typed input"
```
