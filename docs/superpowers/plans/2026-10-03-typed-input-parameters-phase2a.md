# Typed Input Parameters — Phase 2a (Fixed Sections, Validation in the Live Path) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Declare every fixed input section (Discretization, Blocks, FEM, Surface Correction, Solver / Multistep Solver, Outputs, Boundary Conditions, Compute Class Parameters, Contact) as `@params` structs, read the whole deck into a typed `PeriLabInput`, and make it the validator of the live input path — with strict mode, a `--no_strict` command line flag, and a golden test proving every shipped deck still parses.

**Architecture:** Phase 1's `ParameterSpec` gets four small, generic extensions (scalar unions, dicts of scalars, a `check!` cross-field hook, `Globals` skipped in named sections). A new module `PeriLab.InputDeck` (`src/Support/Parameters/Input/`) holds the PeriLab-specific section structs and `read_input`. `Parameter_Handling.validate_yaml` calls `read_input` and reports all errors at once; model parameters (`Models`) keep the legacy validator until phase 3. Consumers still receive the validated `Dict` unchanged (the spec's bridge) — switching them to struct fields is phase 2b.

**Tech Stack:** Julia 1.12, existing dependencies only, Test stdlib.

**Spec:** `docs/superpowers/specs/2026-10-03-typed-input-parameters-design.md` — this plan implements spec §5 phase 2 *except* moving consumers to struct fields (phase 2b, separate plan), plus the strict-mode flag of §2.3 and the golden compatibility test of §6. Phase 1 code: `src/Support/Parameters/Spec/` (module `PeriLab.ParameterSpec`).

## Global Constraints

- Julia `1.12`; no new packages in `Project.toml`.
- Work on branch `feature/typed-input-parameters`; commit after each task (message ends with the `Co-Authored-By` line the session prescribes). Never merge.
- Full backward compatibility of YAML keys: aliases are copied exactly from existing decks/code (spaces, capitalisation, `Step_X`, `m`, `Maximum number of iterations`, …). No YAML format change.
- Consumers keep receiving `params["PeriLab"]` as a `Dict` exactly as today; phase 2a must not change simulation results.
- Strict mode is the default. Strict off via deck key `Strict Validation: false` or command line flag `--no_strict` (underscore, matching the existing `--dry_run` / `--output_dir` flags). A key that matches a declared key apart from case/spaces/punctuation is an error even in non-strict mode (phase 1 behaviour).
- `Globals` at any level is an escape hatch and never validated — except `Contact.Globals`, which is real contact configuration (see Task 4).
- Models (`PeriLab.Models`) are validated by the legacy validator in this phase (phase 3 replaces it).
- Every new source file starts with the SPDX header:
  ```
  # SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
  #
  # SPDX-License-Identifier: BSD-3-Clause
  ```
- Optional keys: if the code reading the key has a default (`get(params, key, default)`), declare that default; otherwise declare `Union{Nothing,T}` with `default = nothing`. Required flags follow the legacy `expected_structure` (it already enforced them).
- Test commands (from repository root):
  - Phase 1 spec tests: `julia --project=. test/unit_tests/Support/Parameters/Spec/run_spec_tests.jl`
  - Input tests (new): `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl`
  - Full suite: `julia --project=. -e 'using Pkg; Pkg.test()'` (≈30 min; run in the background, write output to a file).
  Test files use `const PS = PeriLab.ParameterSpec` / `const ID = PeriLab.InputDeck` and never `using` their exports. `test/runtests.jl` disables warn-level logging, so tests must not rely on `@test_logs` for warnings.

## Review Focus

1. A deck giving both `Output Frequency` and `Number of Output Steps` must stay valid (today: warning, first wins), not become an error — test in Task 4 ("output frequency rules").
2. A misspelled top-level section (`Boundary Condition` for `Boundary Conditions`) must be an error with a suggestion, not silently ignored — test in Task 5 ("unknown top-level section").
3. `--no_strict` must let a deck with a genuinely unknown key run (warning), while a case-only typo still stops it — test in Task 7 ("no_strict flag").
4. A boundary condition `Value` written as a quoted number (`"0.0"`) or an expression (`"0*t"`) must be accepted as text, not rejected as "expected a number" — test in Task 4 ("boundary condition values").
5. Node set entries given as an integer (`5`), a list string (`"2 3 4"`) or a file name (`"ns10.txt"`) must all be accepted — test in Task 2 ("node set entries").

---

## File Structure

| File | Responsibility |
|---|---|
| `src/Support/Parameters/Spec/params_macro.jl` (modify) | `supported_type` for scalar unions and dicts of scalars; `check!` hook |
| `src/Support/Parameters/Spec/convert.jl` (modify) | conversion of scalar unions; named sections of any supported value type; skip `Globals` |
| `src/Support/Parameters/Spec/build.jl`, `model.jl` (modify) | call `check!` after building |
| `src/Support/Parameters/Input/InputDeck.jl` | module `InputDeck`, includes, exports |
| `src/Support/Parameters/Input/discretization.jl` | Discretization, bond filters, Gcode, surface extrusion, external topology |
| `src/Support/Parameters/Input/blocks.jl` | Blocks, FEM, Surface Correction |
| `src/Support/Parameters/Input/solver.jl` | Solver / Multistep step and solver option sections |
| `src/Support/Parameters/Input/outputs.jl` | Outputs |
| `src/Support/Parameters/Input/conditions.jl` | Boundary Conditions, Compute Class Parameters |
| `src/Support/Parameters/Input/contact.jl` | Contact (custom reader, `Globals` is real config) |
| `src/Support/Parameters/Input/input.jl` | `PeriLabSections`, `PeriLabInput`, `read_input` |
| `src/Support/Parameters/parameter_handling.jl` (modify) | `validate_yaml` uses `read_input`; legacy `validate_models` |
| `src/IO/read_inputdeck.jl`, `src/IO/IO.jl`, `src/PeriLab.jl` (modify) | `--no_strict` plumbing, include `InputDeck` |
| `test/unit_tests/Support/Parameters/Spec/ut_extensions.jl` | Task 1 tests |
| `test/unit_tests/Support/Parameters/Input/*` | runner, list, one test file per task |
| `test/unit_tests/Support/Parameters/ut_parameter_handling.jl` (modify) | old `validate_yaml` tests updated to the new contract |
| 6 YAML decks (modify) | remove dead keys (Task 6) |

---

### Task 1: ParameterSpec extensions — scalar unions, dicts of scalars, `check!`, `Globals` in named sections

**Files:**
- Modify: `src/Support/Parameters/Spec/params_macro.jl`
- Modify: `src/Support/Parameters/Spec/convert.jl`
- Modify: `src/Support/Parameters/Spec/build.jl`
- Modify: `src/Support/Parameters/Spec/model.jl`
- Create: `test/unit_tests/Support/Parameters/Spec/ut_extensions.jl`
- Modify: `test/unit_tests/Support/Parameters/Spec/spec_tests.jl`

**Interfaces:**
- Consumes (phase 1): `supported_type`, `convert_value`, `_convert_named_sections`, `_parse_section`, `_parse_model`, `ParseContext`, `add_error!`, `_describe`, `FAILED`.
- Produces:
  - `_is_scalar_union(T)::Bool` — `T` is a `Union` whose members are all in `(Nothing, Int64, Float64, String, Bool)`, not both `Int64` and `Float64`, at least one non-`Nothing` member.
  - `supported_type` accepts scalar unions (e.g. `Union{Int64,String}`, `Union{Nothing,Int64,String}`, `Union{Float64,String}`) and `Dict{String,V}` for any `V` that is a scalar type, a scalar union, or a non-parametric `@params` struct.
  - `convert_value` converts scalar unions by the raw value's kind (Bool→Bool; integer→Int64, else Float64; float→Float64, else Int64 if integral; string→String).
  - `check!(p, path::String, ctx::ParseContext)` — exported hook, default `nothing`; called by `parse_section` and `parse_model` on every successfully built struct, before `derive`.
  - Named sections (`Dict{String,V}` fields) skip a `Globals` entry.

- [ ] **Step 1: Write the failing test**

`test/unit_tests/Support/Parameters/Spec/ut_extensions.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const PS = PeriLab.ParameterSpec

PS.@params struct UTUnions
    degree::Union{Int64,String} = req("Degree")
    value::Union{Float64,String} = req("Value")
    step::Union{Nothing,Int64,String} = opt("Step ID"; default = nothing)
    sets::Dict{String,Union{Int64,String}} = opt("Node Sets";
                                                 default = Dict{String,Union{Int64,String}}())
    flags::Dict{String,Bool} = opt("Flags"; default = Dict{String,Bool}())
end

PS.@params struct UTChecked
    low::Float64 = req("Low")
    high::Float64 = req("High")
end
function PS.check!(p::UTChecked, path::String, ctx::PS.ParseContext)
    p.low <= p.high || PS.add_error!(ctx, path, "\"Low\" must not exceed \"High\"")
    return nothing
end

PS.@params struct UTCheckedHolder
    items::Dict{String,UTChecked} = req("Items")
end

function ut_ext_conv(T, raw)
    ctx = PS.ParseContext()
    return PS.convert_value(T, raw, "k", ctx), ctx
end

@testset "scalar unions are supported" begin
    @test PS.supported_type(Union{Int64,String})
    @test PS.supported_type(Union{Nothing,Int64,String})
    @test PS.supported_type(Union{Float64,String})
    @test !PS.supported_type(Union{Int64,Float64})
    @test PS.supported_type(Dict{String,Union{Int64,String}})
    @test PS.supported_type(Dict{String,Bool})
end

@testset "scalar union conversion" begin
    @test ut_ext_conv(Union{Int64,String}, 5)[1] === 5
    @test ut_ext_conv(Union{Int64,String}, "1 1")[1] == "1 1"
    @test ut_ext_conv(Union{Int64,String}, 2.0)[1] === 2
    @test ut_ext_conv(Union{Float64,String}, -25)[1] === -25.0
    @test ut_ext_conv(Union{Float64,String}, "0*t")[1] == "0*t"
    @test ut_ext_conv(Union{Nothing,Int64,String}, nothing)[1] === nothing
    v, ctx = ut_ext_conv(Union{Int64,String}, 2.5)
    @test v === PS.FAILED
    @test ctx.errors[1].message == "expected an integer or text, got 2.5"
    v, ctx = ut_ext_conv(Union{Float64,String}, true)
    @test ctx.errors[1].message == "expected a number or text, got true"
end

@testset "dicts of scalars, Globals skipped in named sections" begin
    ctx = PS.ParseContext()
    p = PS.parse_section(UTUnions,
                         Dict{String,Any}("Degree" => "1 1", "Value" => 0,
                                          "Node Sets" => Dict{String,Any}("a" => 5,
                                                                          "b" => "2 3 4",
                                                                          "Globals" => Dict{String,Any}()),
                                          "Flags" => Dict{String,Any}("x" => true)),
                         "S", ctx)
    @test isempty(ctx.errors)
    @test p.degree == "1 1" && p.value === 0.0 && p.step === nothing
    @test p.sets == Dict{String,Union{Int64,String}}("a" => 5, "b" => "2 3 4")
    @test p.flags == Dict("x" => true)
    ctx = PS.ParseContext()
    PS.parse_section(UTUnions,
                     Dict{String,Any}("Degree" => 1, "Value" => 1,
                                      "Flags" => Dict{String,Any}("x" => "yes")), "S", ctx)
    @test ctx.errors[1].path == "S.Flags.x"
    @test ctx.errors[1].message == "expected true or false, got \"yes\""
end

@testset "check! runs after building, also for nested sections" begin
    ctx = PS.ParseContext()
    @test PS.parse_section(UTChecked, Dict{String,Any}("Low" => 2.0, "High" => 1.0), "C",
                           ctx) isa UTChecked
    @test ctx.errors[1].path == "C"
    @test ctx.errors[1].message == "\"Low\" must not exceed \"High\""
    ctx = PS.ParseContext()
    PS.parse_section(UTCheckedHolder,
                     Dict{String,Any}("Items" => Dict{String,Any}("a" => Dict{String,Any}("Low" => 3.0,
                                                                                         "High" => 1.0))),
                     "H", ctx)
    @test ctx.errors[1].path == "H.Items.a"
    ctx = PS.ParseContext()
    PS.parse_section(UTChecked, Dict{String,Any}("Low" => 1.0), "C", ctx)
    @test length(ctx.errors) == 1           # check! is not run on a struct that failed to build
end
```

Add `"ut_extensions.jl"` at the end of the list in `test/unit_tests/Support/Parameters/Spec/spec_tests.jl`:

```julia
for file in ["ut_errors.jl", "ut_dependent.jl", "ut_convert.jl", "ut_params.jl", "ut_model.jl",
             "ut_end_to_end.jl", "ut_extensions.jl"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. test/unit_tests/Support/Parameters/Spec/run_spec_tests.jl`
Expected: FAIL in `ut_extensions.jl` with `ParamsDefinitionError: UTUnions.degree: unsupported field type Union{Int64, String}`.

- [ ] **Step 3: Write the implementation**

In `params_macro.jl`, export the hook — replace

```julia
export @params, derive
```

with

```julia
export @params, derive, check!
```

and after the `derive(p) = p` definition add:

```julia
"""
    check!(p, path, ctx)

Hook for cross-field rules of an `@params` struct (e.g. "exactly one solver",
"Low must not exceed High"). Runs once after `p` was built successfully and
before `derive`; add problems with `add_error!(ctx, path, msg)`. The default
does nothing.
"""
check!(p, path::String, ctx::ParseContext) = nothing

const _SCALAR_UNION_MEMBERS = (Nothing, Int64, Float64, String, Bool)

function _is_scalar_union(T)
    T isa Union || return false
    members = Base.uniontypes(T)
    return all(m -> m in _SCALAR_UNION_MEMBERS, members) &&
           !(Int64 in members && Float64 in members) &&
           any(m -> m !== Nothing, members)
end
```

Still in `params_macro.jl`, replace the body of `supported_type`:

```julia
function supported_type(T)
    if T isa Union
        S = _nonnothing(T)
        return Nothing <: T && !(S isa Union) && S !== Dependent && supported_type(S)
    end
    (T in _SCALAR_TYPES || T in _VECTOR_TYPES || T === Dependent) && return true
    T isa DataType || return false
    T <: Enum && return true
    if T <: Dict && T.parameters[1] === String
        V = T.parameters[2]
        return V isa DataType && is_params(V)
    end
    return is_params(T)
end
```

with

```julia
function supported_type(T)
    _is_scalar_union(T) && return true
    if T isa Union
        S = _nonnothing(T)
        return Nothing <: T && !(S isa Union) && S !== Dependent && supported_type(S)
    end
    (T in _SCALAR_TYPES || T in _VECTOR_TYPES || T === Dependent) && return true
    T isa DataType || return false
    T <: Enum && return true
    if T <: Dict && T.parameters[1] === String
        V = T.parameters[2]
        return V in _SCALAR_TYPES || _is_scalar_union(V) || (V isa DataType && is_params(V))
    end
    return is_params(T)
end
```

In `convert.jl`, replace

```julia
    raw isa T && return _fresh(raw)
    if T isa Union
        return convert_value(_nonnothing(T), raw, path, ctx; alias = alias)
```

with

```julia
    raw isa T && return _fresh(raw)
    if _is_scalar_union(T)
        return _convert_scalar_union(T, raw, path, ctx)
    elseif T isa Union
        return convert_value(_nonnothing(T), raw, path, ctx; alias = alias)
```

replace

```julia
    elseif T isa DataType && T <: Dict && T.parameters[1] === String &&
           is_params(T.parameters[2])
        return _convert_named_sections(T, raw, path, ctx)
```

with

```julia
    elseif T isa DataType && T <: Dict && T.parameters[1] === String &&
           supported_type(T)
        return _convert_named_sections(T, raw, path, ctx)
```

in `_convert_named_sections`, replace

```julia
    for (k, item) in raw
        v = convert_value(V, item, join_path(path, string(k)), ctx)
```

with

```julia
    for (k, item) in raw
        string(k) == "Globals" && continue
        v = convert_value(V, item, join_path(path, string(k)), ctx)
```

and append at the end of `convert.jl`:

```julia
const _SCALAR_KIND_NAMES = Dict{Any,String}(Int64 => "an integer", Float64 => "a number",
                                            String => "text", Bool => "true or false")

function _convert_scalar_union(T, raw, path::String, ctx::ParseContext)
    members = Base.uniontypes(T)
    if raw isa Bool
        Bool in members && return raw
    elseif raw isa Integer
        Int64 in members && return Int64(raw)
        Float64 in members && return Float64(raw)
    elseif raw isa AbstractFloat
        Float64 in members && return Float64(raw)
        Int64 in members && isinteger(raw) && return Int64(raw)
    elseif raw isa AbstractString
        String in members && return String(raw)
    end
    # fixed order, independent of how Julia orders union members
    expected = join([_SCALAR_KIND_NAMES[m] for m in (Int64, Float64, String, Bool)
                     if m in members], " or ")
    return _fail(ctx, path, "expected $expected, got $(_describe(raw))")
end
```

In `build.jl`, inside `_parse_section`, replace

```julia
    p = build(T, dict, path, ctx)
    check_unknown!(dict, aliases(T), path, ctx)
    return p === nothing ? nothing : derive(p)
```

with

```julia
    p = build(T, dict, path, ctx)
    check_unknown!(dict, aliases(T), path, ctx)
    p === nothing && return nothing
    check!(p, path, ctx)
    return derive(p)
```

In `model.jl`, inside `_parse_model`, replace

```julia
        part = build(T, dict, path, ctx; owner = name)
        push!(parts, part === nothing ? nothing : derive(part))
```

with

```julia
        part = build(T, dict, path, ctx; owner = name)
        part === nothing || check!(part, path, ctx)
        push!(parts, part === nothing ? nothing : derive(part))
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. test/unit_tests/Support/Parameters/Spec/run_spec_tests.jl`
Expected: PASS (all phase 1 tests plus `ut_extensions.jl`). The phase 1 test `ut_convert.jl` "Union{Nothing,T}" expects `"expected a number, got \"x\""` — the scalar-union path produces exactly that.

- [ ] **Step 5: Commit**

```bash
git add src/Support/Parameters/Spec test/unit_tests/Support/Parameters/Spec
git commit -m "ParameterSpec: scalar unions, dicts of scalars, check! hook"
```

---

### Task 2: `InputDeck` module — Discretization, Blocks, FEM, Surface Correction

**Files:**
- Create: `src/Support/Parameters/Input/InputDeck.jl`
- Create: `src/Support/Parameters/Input/discretization.jl`
- Create: `src/Support/Parameters/Input/blocks.jl`
- Modify: `src/PeriLab.jl` (include after the `ParameterSpec` include)
- Create: `test/unit_tests/Support/Parameters/Input/run_input_tests.jl`
- Create: `test/unit_tests/Support/Parameters/Input/input_tests.jl`
- Create: `test/unit_tests/Support/Parameters/Input/ut_mesh_blocks.jl`
- Modify: `test/runtests.jl` (Support → Parameters testset)

**Interfaces:**
- Consumes: Task 1 and phase 1 (`@params`, `req`, `opt`, `parse_section`, `ParseContext`).
- Produces (all `@params`, field names are the snake_case of the YAML key):
  - `ExternalTopologyParams(file::String, add_neighbor_search::Union{Nothing,Bool})`
  - `SurfaceExtrusionParams(direction::String, step_x, step_y, step_z::Float64, number::Int64)`
  - `GcodeParams(overwrite_mesh::Bool, sampling, width, height::Float64, scale::Float64, start_command, stop_command, end_command::Union{Nothing,String}, blocks::Union{Nothing,Dict{String,String}})`
  - `BondFilterParams` (type, normals, corners, unit vectors, centers, radius, lengths, `allow_contact::Bool`)
  - `DiscretizationParams` (fields listed in Step 3)
  - `BlockParams`, `FEMCouplingParams`, `FEMParams`, `SurfaceCorrectionParams`

- [ ] **Step 1: Create runner, list, hook into runtests, write the failing test**

`test/unit_tests/Support/Parameters/Input/run_input_tests.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# Standalone runner for the InputDeck unit tests:
#   julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl
using Test
import PeriLab

@testset "InputDeck" begin
    include(joinpath(@__DIR__, "input_tests.jl"))
end
```

`test/unit_tests/Support/Parameters/Input/input_tests.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

for file in ["ut_mesh_blocks.jl"]
    @testset "$file" begin
        include(joinpath(@__DIR__, file))
    end
end
```

In `test/runtests.jl`, replace

```julia
                @testset "ParameterSpec" begin
                    include("unit_tests/Support/Parameters/Spec/spec_tests.jl")
                end
```

with

```julia
                @testset "ParameterSpec" begin
                    include("unit_tests/Support/Parameters/Spec/spec_tests.jl")
                end
                @testset "InputDeck" begin
                    include("unit_tests/Support/Parameters/Input/input_tests.jl")
                end
```

`test/unit_tests/Support/Parameters/Input/ut_mesh_blocks.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const PS = PeriLab.ParameterSpec
const ID = PeriLab.InputDeck

function ut_section(T, dict)
    ctx = PS.ParseContext()
    return PS.parse_section(T, dict, "X", ctx), ctx
end

@testset "Discretization: minimal and full" begin
    d, ctx = ut_section(ID.DiscretizationParams,
                        Dict{String,Any}("Type" => "Text File", "Input Mesh File" => "mesh.txt"))
    @test isempty(ctx.errors)
    @test d.type == "Text File" && d.input_mesh_file == "mesh.txt"
    @test isempty(d.node_sets) && isempty(d.bond_filters)
    @test d.gcode === nothing && d.surface_extrusion === nothing
    @test d.influence_function === nothing
    full = Dict{String,Any}("Type" => "Exodus", "Input Mesh File" => "m.g",
                            "Distribution Type" => "Neighbor based",
                            "Influence Function" => "1/xi^2",
                            "Horizon Mesh Scaling X" => 1.5,
                            "Input External Topology" => Dict{String,Any}("File" => "t.txt",
                                                                          "Add Neighbor Search" => true),
                            "Surface Extrusion" => Dict{String,Any}("Direction" => "X",
                                                                    "Step_X" => 0.1,
                                                                    "Step_Y" => 0.1,
                                                                    "Step_Z" => 0,
                                                                    "Number" => 3),
                            "Gcode" => Dict{String,Any}("Overwrite Mesh" => true,
                                                        "Sampling" => 0.5, "Width" => 1,
                                                        "Height" => 0.2,
                                                        "Blocks" => Dict{String,Any}("1" => "block_1")),
                            "Bond Filters" => Dict{String,Any}("bf_1" => Dict{String,Any}("Type" => "Rectangular_Plane",
                                                                                          "Normal X" => 0.0,
                                                                                          "Normal Y" => 1.0,
                                                                                          "Allow Contact" => true)))
    d, ctx = ut_section(ID.DiscretizationParams, full)
    @test isempty(ctx.errors)
    @test d.influence_function == "1/xi^2"
    @test d.horizon_mesh_scaling_x === 1.5 && d.horizon_mesh_scaling_y === nothing
    @test d.input_external_topology.add_neighbor_search === true
    @test d.surface_extrusion.step_z === 0.0 && d.surface_extrusion.number === 3
    @test d.gcode.scale === 1.0 && d.gcode.blocks == Dict("1" => "block_1")
    @test d.bond_filters["bf_1"].allow_contact && d.bond_filters["bf_1"].normal_z === nothing
end

@testset "node set entries" begin
    d, ctx = ut_section(ID.DiscretizationParams,
                        Dict{String,Any}("Type" => "Text File", "Input Mesh File" => "m.txt",
                                         "Node Sets" => Dict{String,Any}("a" => 5,
                                                                         "b" => "2 3 4",
                                                                         "c" => "ns10.txt")))
    @test isempty(ctx.errors)
    @test d.node_sets["a"] === 5
    @test d.node_sets["b"] == "2 3 4"
    @test d.node_sets["c"] == "ns10.txt"
end

@testset "Discretization errors" begin
    d, ctx = ut_section(ID.DiscretizationParams,
                        Dict{String,Any}("Type" => "Text File",
                                         "Surface Extrusion" => Dict{String,Any}("Direction" => "W",
                                                                                 "Step_X" => 1,
                                                                                 "Step_Y" => 1,
                                                                                 "Step_Z" => 1,
                                                                                 "Number" => 1)))
    @test d === nothing
    msgs = Dict(e.path => e.message for e in ctx.errors)
    @test msgs["X.\"Input Mesh File\""] == "missing (required by DiscretizationParams)"
    @test msgs["X.\"Surface Extrusion\".Direction"] ==
          "\"W\" is not one of: \"X\", \"Y\", \"Z\""
end

@testset "Blocks" begin
    b, ctx = ut_section(ID.BlockParams,
                        Dict{String,Any}("Block ID" => 1, "Density" => 2700,
                                         "Horizon" => 2, "Material Model" => "Steel",
                                         "Degradation Model" => "deg", "FEM" => true,
                                         "Step ID" => "1,2"))
    @test isempty(ctx.errors)
    @test b.block_id === 1 && b.density === 2700.0 && b.horizon === 2.0
    @test b.material_model == "Steel" && b.damage_model === nothing
    @test b.degradation_model == "deg" && b.fem === true && b.step_id == "1,2"
    b, ctx = ut_section(ID.BlockParams,
                        Dict{String,Any}("Block ID" => 1, "Density" => 1.0, "Horizon" => -1.0))
    @test ctx.errors[1].message == "-1.0 is below minimum 0"
end

@testset "FEM and Surface Correction" begin
    f, ctx = ut_section(ID.FEMParams,
                        Dict{String,Any}("Element Type" => "Lagrange", "Degree" => "1 1",
                                         "Material Model" => "Elastic",
                                         "Coupling" => Dict{String,Any}("Coupling Type" => "Arlequin",
                                                                        "PD Weight" => 0.5,
                                                                        "Coupling Block" => 2)))
    @test isempty(ctx.errors)
    @test f.degree == "1 1" && f.coupling.pd_weight === 0.5 && f.coupling.coupling_block === 2
    f, ctx = ut_section(ID.FEMParams,
                        Dict{String,Any}("Element Type" => "Lagrange", "Degree" => 1,
                                         "Material Model" => "Elastic"))
    @test f.degree === 1 && f.coupling === nothing
    s, ctx = ut_section(ID.SurfaceCorrectionParams,
                        Dict{String,Any}("Type" => "Volume Correction"))
    @test isempty(ctx.errors) && s.update === false
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl`
Expected: FAIL with `UndefVarError: InputDeck not defined in PeriLab`.

- [ ] **Step 3: Write the implementation**

`src/Support/Parameters/Input/InputDeck.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

"""
    InputDeck

Typed declarations of PeriLab's fixed input sections and `read_input`, which
turns the `PeriLab:` part of a YAML deck into a validated `PeriLabInput`.
Model parameters (`Models`) are declared by the model modules (phase 3).
"""
module InputDeck

using ..ParameterSpec: @params, req, opt, ParseContext, add_error!, join_path,
                       parse_section
import ..ParameterSpec: check!

include("discretization.jl")
include("blocks.jl")

end
```

`src/Support/Parameters/Input/discretization.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

@params struct ExternalTopologyParams
    file::String = req("File"; description = "File with the external element topology")
    add_neighbor_search::Union{Nothing,Bool} = opt("Add Neighbor Search"; default = nothing)
end

@params struct SurfaceExtrusionParams
    direction::String = req("Direction"; allowed = ["X", "Y", "Z"])
    step_x::Float64 = req("Step_X"; quantity = :length)
    step_y::Float64 = req("Step_Y"; quantity = :length)
    step_z::Float64 = req("Step_Z"; quantity = :length)
    number::Int64 = req("Number"; min = 0)
end

@params struct GcodeParams
    overwrite_mesh::Bool = req("Overwrite Mesh")
    sampling::Float64 = req("Sampling"; min = 0, quantity = :length)
    width::Float64 = req("Width"; min = 0, quantity = :length)
    height::Float64 = req("Height"; min = 0, quantity = :length)
    scale::Float64 = opt("Scale"; default = 1.0)
    start_command::Union{Nothing,String} = opt("Start Command"; default = nothing)
    stop_command::Union{Nothing,String} = opt("Stop Command"; default = nothing)
    end_command::Union{Nothing,String} = opt("End Command"; default = nothing)
    blocks::Union{Nothing,Dict{String,String}} = opt("Blocks"; default = nothing)
end

@params struct BondFilterParams
    type::String = req("Type")
    normal_x::Float64 = req("Normal X")
    normal_y::Float64 = req("Normal Y")
    normal_z::Union{Nothing,Float64} = opt("Normal Z"; default = nothing)
    lower_left_corner_x::Union{Nothing,Float64} = opt("Lower Left Corner X"; default = nothing)
    lower_left_corner_y::Union{Nothing,Float64} = opt("Lower Left Corner Y"; default = nothing)
    lower_left_corner_z::Union{Nothing,Float64} = opt("Lower Left Corner Z"; default = nothing)
    bottom_unit_vector_x::Union{Nothing,Float64} = opt("Bottom Unit Vector X"; default = nothing)
    bottom_unit_vector_y::Union{Nothing,Float64} = opt("Bottom Unit Vector Y"; default = nothing)
    bottom_unit_vector_z::Union{Nothing,Float64} = opt("Bottom Unit Vector Z"; default = nothing)
    center_x::Union{Nothing,Float64} = opt("Center X"; default = nothing)
    center_y::Union{Nothing,Float64} = opt("Center Y"; default = nothing)
    center_z::Union{Nothing,Float64} = opt("Center Z"; default = nothing)
    radius::Union{Nothing,Float64} = opt("Radius"; default = nothing, min = 0)
    bottom_length::Union{Nothing,Float64} = opt("Bottom Length"; default = nothing, min = 0)
    side_length::Union{Nothing,Float64} = opt("Side Length"; default = nothing, min = 0)
    allow_contact::Bool = opt("Allow Contact"; default = false)
end

@params struct DiscretizationParams
    type::String = req("Type"; description = "Mesh format, e.g. \"Text File\" or \"Exodus\"")
    input_mesh_file::String = req("Input Mesh File")
    input_external_topology::Union{Nothing,ExternalTopologyParams} = opt("Input External Topology";
                                                                         default = nothing)
    node_sets::Dict{String,Union{Int64,String}} = opt("Node Sets";
                                                      default = Dict{String,
                                                                     Union{Int64,String}}(),
                                                      description = "Node id, list of ids, or file")
    distribution_type::Union{Nothing,String} = opt("Distribution Type"; default = nothing)
    influence_function::Union{Nothing,String} = opt("Influence Function"; default = nothing)
    surface_extrusion::Union{Nothing,SurfaceExtrusionParams} = opt("Surface Extrusion";
                                                                   default = nothing)
    bond_filters::Dict{String,BondFilterParams} = opt("Bond Filters";
                                                      default = Dict{String,BondFilterParams}())
    horizon_mesh_scaling_x::Union{Nothing,Float64} = opt("Horizon Mesh Scaling X";
                                                         default = nothing)
    horizon_mesh_scaling_y::Union{Nothing,Float64} = opt("Horizon Mesh Scaling Y";
                                                         default = nothing)
    horizon_mesh_scaling_z::Union{Nothing,Float64} = opt("Horizon Mesh Scaling Z";
                                                         default = nothing)
    gcode::Union{Nothing,GcodeParams} = opt("Gcode"; default = nothing)
end
```

`src/Support/Parameters/Input/blocks.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

@params struct BlockParams
    block_id::Int64 = req("Block ID"; min = 1)
    density::Float64 = req("Density"; min = 0, quantity = :density)
    horizon::Float64 = req("Horizon"; min = 0, quantity = :length)
    specific_heat_capacity::Union{Nothing,Float64} = opt("Specific Heat Capacity";
                                                         default = nothing, min = 0)
    material_model::Union{Nothing,String} = opt("Material Model"; default = nothing)
    damage_model::Union{Nothing,String} = opt("Damage Model"; default = nothing)
    thermal_model::Union{Nothing,String} = opt("Thermal Model"; default = nothing)
    additive_model::Union{Nothing,String} = opt("Additive Model"; default = nothing)
    pre_calculation_model::Union{Nothing,String} = opt("Pre Calculation Model";
                                                       default = nothing)
    degradation_model::Union{Nothing,String} = opt("Degradation Model"; default = nothing)
    angle_x::Union{Nothing,Float64} = opt("Angle X"; default = nothing)
    angle_y::Union{Nothing,Float64} = opt("Angle Y"; default = nothing)
    angle_z::Union{Nothing,Float64} = opt("Angle Z"; default = nothing)
    step_id::Union{Nothing,Int64,String} = opt("Step ID"; default = nothing)
    fem::Union{Nothing,Bool} = opt("FEM"; default = nothing)
end

@params struct FEMCouplingParams
    coupling_type::String = req("Coupling Type")
    pd_weight::Union{Nothing,Float64} = opt("PD Weight"; default = nothing)
    kappa::Union{Nothing,Float64} = opt("Kappa"; default = nothing)
    coupling_block::Union{Nothing,Int64} = opt("Coupling Block"; default = nothing)
end

@params struct FEMParams
    element_type::String = req("Element Type")
    degree::Union{Int64,String} = req("Degree")
    material_model::String = req("Material Model")
    coupling::Union{Nothing,FEMCouplingParams} = opt("Coupling"; default = nothing)
end

@params struct SurfaceCorrectionParams
    type::String = req("Type")
    update::Bool = opt("Update"; default = false)
end
```

In `src/PeriLab.jl`, replace

```julia
include("./Support/Parameters/Spec/ParameterSpec.jl")
```

with

```julia
include("./Support/Parameters/Spec/ParameterSpec.jl")
include("./Support/Parameters/Input/InputDeck.jl")
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/Support/Parameters/Input src/PeriLab.jl test/runtests.jl test/unit_tests/Support/Parameters/Input
git commit -m "InputDeck: Discretization, Blocks, FEM, Surface Correction sections"
```

---

### Task 3: Solver and Multistep step sections

**Files:**
- Create: `src/Support/Parameters/Input/solver.jl`
- Modify: `src/Support/Parameters/Input/InputDeck.jl`
- Create: `test/unit_tests/Support/Parameters/Input/ut_solver.jl`
- Modify: `test/unit_tests/Support/Parameters/Input/input_tests.jl`

**Interfaces:**
- Consumes: Task 1 (`check!`), Task 2 (module).
- Produces:
  - Option sections, each with `safety_factor::Float64 = 1.0`, `fixed_dt::Float64 = -1.0`, `numerical_damping::Float64 = 0.0` (the solver getters read these three from whichever solver section is active): `VerletParams`, `StaticParams` (+ `solution_tolerance`, `residual_tolerance`, `maximum_number_of_iterations`, `show_solver_iteration`, `residual_scaling`, `m`, `linear_start_value`, `nlsolve`, `solver_type`), `LinearStaticMatrixParams` (+ `matrix_update`), `ModelReductionParams`, `VerletMatrixParams` (+ `model_reduction`), `NewmarkParams` (+ `matrix_update`, `newmark_delta`, `newmark_alpha`).
  - `SolverParams` — used for `Solver` and for every `Multistep Solver` step; fields `initial_time`, `final_time`, `additional_time::Union{Nothing,Float64}`, `number_of_steps::Int64`, `maximum_damage::Float64`, `step_id::Union{Nothing,Int64}`, the six model flags, the three calculation flags, and `verlet`, `static`, `linear_static_matrix_based`, `verlet_matrix_based`, `newmark` (each `Union{Nothing,…}`).
  - `check!(::SolverParams)`: exactly one solver section; `Final Time` or `Additional Time` present.
  - `const SOLVER_NAMES` (the five YAML names, in the order above).

- [ ] **Step 1: Write the failing test**

`test/unit_tests/Support/Parameters/Input/ut_solver.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const PS = PeriLab.ParameterSpec
const ID = PeriLab.InputDeck

function ut_solver(dict)
    ctx = PS.ParseContext()
    return PS.parse_section(ID.SolverParams, dict, "Solver", ctx), ctx
end

@testset "Verlet solver with defaults" begin
    s, ctx = ut_solver(Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1,
                                        "Verlet" => Dict{String,Any}("Safety Factor" => 0.9)))
    @test isempty(ctx.errors)
    @test s.final_time === 1.0 && s.number_of_steps === 1 && s.maximum_damage === Inf
    @test s.material_models && s.pre_calculation_models && !s.damage_models
    @test s.verlet.safety_factor === 0.9 && s.verlet.fixed_dt === -1.0
    @test s.verlet.numerical_damping === 0.0
    @test s.static === nothing && s.newmark === nothing
end

@testset "Static, matrix based, Newmark and model reduction options" begin
    s, ctx = ut_solver(Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                        "Static" => Dict{String,Any}("NLSolve" => true,
                                                                     "Fixed dt" => 0.1,
                                                                     "Residual scaling" => 70000)))
    @test isempty(ctx.errors)
    @test s.static.residual_scaling === 70000.0 && s.static.m === 15
    @test s.static.maximum_number_of_iterations === 100 && s.static.nlsolve === true
    @test s.static.fixed_dt === 0.1
    s, ctx = ut_solver(Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                        "Verlet Matrix Based" => Dict{String,Any}("Model Reduction" => Dict{String,Any}("Type" => "Craig Bampton",
                                                                                                                     "Number of Modes" => 5))))
    @test isempty(ctx.errors)
    @test s.verlet_matrix_based.model_reduction.number_of_modes === 5
    @test s.verlet_matrix_based.model_reduction.material_point_region === true
    s, ctx = ut_solver(Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                        "Newmark" => Dict{String,Any}("Matrix Update" => true)))
    @test isempty(ctx.errors)
    @test s.newmark.matrix_update && s.newmark.newmark_delta === 0.5
    @test s.newmark.newmark_alpha === nothing
end

@testset "solver rules" begin
    s, ctx = ut_solver(Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0))
    @test ctx.errors[1].path == "Solver"
    @test ctx.errors[1].message ==
          "one solver is required: Verlet, Static, Linear Static Matrix Based, Verlet Matrix Based or Newmark"
    s, ctx = ut_solver(Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                        "Verlet" => Dict{String,Any}(),
                                        "Static" => Dict{String,Any}()))
    @test ctx.errors[1].message == "only one solver may be given, found: Verlet, Static"
    s, ctx = ut_solver(Dict{String,Any}("Initial Time" => 0.0, "Verlet" => Dict{String,Any}()))
    @test ctx.errors[1].message == "\"Final Time\" or \"Additional Time\" is required"
    s, ctx = ut_solver(Dict{String,Any}("Additional Time" => 1.0, "Step ID" => 2,
                                        "Verlet" => Dict{String,Any}()))
    @test isempty(ctx.errors) && s.step_id === 2
    s, ctx = ut_solver(Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                        "External" => Dict{String,Any}(),
                                        "Verlet" => Dict{String,Any}()))
    @test ctx.errors[1].path == "Solver.External"
    @test startswith(ctx.errors[1].message, "unknown key")
end
```

Update `input_tests.jl`:

```julia
for file in ["ut_mesh_blocks.jl", "ut_solver.jl"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl`
Expected: FAIL in `ut_solver.jl` with `UndefVarError: SolverParams not defined in PeriLab.InputDeck`.

- [ ] **Step 3: Write the implementation**

In `InputDeck.jl`, replace

```julia
include("blocks.jl")
```

with

```julia
include("blocks.jl")
include("solver.jl")
```

`src/Support/Parameters/Input/solver.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# Every solver option section carries "Safety Factor", "Fixed dt" and
# "Numerical Damping": the solver getters read them from the active section.

@params struct VerletParams
    safety_factor::Float64 = opt("Safety Factor"; default = 1.0, min = 0)
    fixed_dt::Float64 = opt("Fixed dt"; default = -1.0, quantity = :time)
    numerical_damping::Float64 = opt("Numerical Damping"; default = 0.0, min = 0)
end

@params struct StaticParams
    safety_factor::Float64 = opt("Safety Factor"; default = 1.0, min = 0)
    fixed_dt::Float64 = opt("Fixed dt"; default = -1.0, quantity = :time)
    numerical_damping::Float64 = opt("Numerical Damping"; default = 0.0, min = 0)
    solution_tolerance::Float64 = opt("Solution tolerance"; default = 1e-7, min = 0)
    residual_tolerance::Float64 = opt("Residual tolerance"; default = 1e-7, min = 0)
    maximum_number_of_iterations::Int64 = opt("Maximum number of iterations"; default = 100,
                                              min = 1)
    show_solver_iteration::Bool = opt("Show solver iteration"; default = false)
    residual_scaling::Float64 = opt("Residual scaling"; default = 1e6)
    m::Int64 = opt("m"; default = 15, min = 1)
    linear_start_value::Union{Nothing,String} = opt("Linear Start Value"; default = nothing)
    nlsolve::Union{Nothing,Bool} = opt("NLSolve"; default = nothing)
    solver_type::Union{Nothing,String} = opt("Solver Type"; default = nothing)
end

@params struct LinearStaticMatrixParams
    safety_factor::Float64 = opt("Safety Factor"; default = 1.0, min = 0)
    fixed_dt::Float64 = opt("Fixed dt"; default = -1.0, quantity = :time)
    numerical_damping::Float64 = opt("Numerical Damping"; default = 0.0, min = 0)
    matrix_update::Bool = opt("Matrix Update"; default = false)
end

@params struct ModelReductionParams
    type::String = req("Type")
    number_of_modes::Int64 = opt("Number of Modes"; default = 1, min = 1)
    material_point_region::Bool = opt("Material Point Region"; default = true)
    reduction_blocks::Union{Nothing,String} = opt("Reduction Blocks"; default = nothing)
end

@params struct VerletMatrixParams
    safety_factor::Float64 = opt("Safety Factor"; default = 1.0, min = 0)
    fixed_dt::Float64 = opt("Fixed dt"; default = -1.0, quantity = :time)
    numerical_damping::Float64 = opt("Numerical Damping"; default = 0.0, min = 0)
    model_reduction::Union{Nothing,ModelReductionParams} = opt("Model Reduction";
                                                               default = nothing)
end

@params struct NewmarkParams
    safety_factor::Float64 = opt("Safety Factor"; default = 1.0, min = 0)
    fixed_dt::Float64 = opt("Fixed dt"; default = -1.0, quantity = :time)
    numerical_damping::Float64 = opt("Numerical Damping"; default = 0.0, min = 0)
    matrix_update::Bool = opt("Matrix Update"; default = false)
    newmark_delta::Float64 = opt("Newmark Delta"; default = 0.5)
    newmark_alpha::Union{Nothing,Float64} = opt("Newmark Alpha"; default = nothing)
end

"Used for `Solver` and for every step of `Multistep Solver`."
@params struct SolverParams
    initial_time::Union{Nothing,Float64} = opt("Initial Time"; default = nothing,
                                               quantity = :time)
    final_time::Union{Nothing,Float64} = opt("Final Time"; default = nothing, quantity = :time)
    additional_time::Union{Nothing,Float64} = opt("Additional Time"; default = nothing,
                                                  quantity = :time)
    number_of_steps::Int64 = opt("Number of Steps"; default = 1, min = 1)
    maximum_damage::Float64 = opt("Maximum Damage"; default = Inf)
    step_id::Union{Nothing,Int64} = opt("Step ID"; default = nothing)
    additive_models::Bool = opt("Additive Models"; default = false)
    degradation_models::Bool = opt("Degradation Models"; default = false)
    damage_models::Bool = opt("Damage Models"; default = false)
    material_models::Bool = opt("Material Models"; default = true)
    thermal_models::Bool = opt("Thermal Models"; default = false)
    pre_calculation_models::Bool = opt("Pre Calculation Models"; default = true)
    calculate_cauchy::Bool = opt("Calculate Cauchy"; default = false)
    calculate_von_mises_stress::Bool = opt("Calculate von Mises stress"; default = false)
    calculate_strain::Bool = opt("Calculate Strain"; default = false)
    verlet::Union{Nothing,VerletParams} = opt("Verlet"; default = nothing)
    static::Union{Nothing,StaticParams} = opt("Static"; default = nothing)
    linear_static_matrix_based::Union{Nothing,LinearStaticMatrixParams} = opt("Linear Static Matrix Based";
                                                                              default = nothing)
    verlet_matrix_based::Union{Nothing,VerletMatrixParams} = opt("Verlet Matrix Based";
                                                                 default = nothing)
    newmark::Union{Nothing,NewmarkParams} = opt("Newmark"; default = nothing)
end

const SOLVER_NAMES = ("Verlet", "Static", "Linear Static Matrix Based", "Verlet Matrix Based",
                      "Newmark")

function check!(p::SolverParams, path::String, ctx::ParseContext)
    sections = (p.verlet, p.static, p.linear_static_matrix_based, p.verlet_matrix_based,
                p.newmark)
    given = [name for (name, section) in zip(SOLVER_NAMES, sections) if section !== nothing]
    if isempty(given)
        add_error!(ctx, path,
                   "one solver is required: Verlet, Static, Linear Static Matrix Based, Verlet Matrix Based or Newmark")
    elseif length(given) > 1
        add_error!(ctx, path, "only one solver may be given, found: $(join(given, ", "))")
    end
    if p.final_time === nothing && p.additional_time === nothing
        add_error!(ctx, path, "\"Final Time\" or \"Additional Time\" is required")
    end
    return nothing
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/Support/Parameters/Input test/unit_tests/Support/Parameters/Input
git commit -m "InputDeck: Solver and Multistep step sections"
```

---

### Task 4: Outputs, Boundary Conditions, Compute Classes, Contact

**Files:**
- Create: `src/Support/Parameters/Input/outputs.jl`
- Create: `src/Support/Parameters/Input/conditions.jl`
- Create: `src/Support/Parameters/Input/contact.jl`
- Modify: `src/Support/Parameters/Input/InputDeck.jl`
- Create: `test/unit_tests/Support/Parameters/Input/ut_outputs_conditions_contact.jl`
- Modify: `test/unit_tests/Support/Parameters/Input/input_tests.jl`

**Interfaces:**
- Consumes: Tasks 1–3; phase 1 `parse_section`, `join_path`, `add_error!`, `ParameterSpec._describe`.
- Produces:
  - `OutputParams` + `check!`: `Output Frequency` or `Number of Output Steps` required (both allowed).
  - `BoundaryConditionParams`, `ComputeClassParams`.
  - `ContactGlobalsParams`, `ContactGroupParams`, `ContactModelParams` (`@params`), and plain `struct ContactInput; globals::ContactGlobalsParams; models::Dict{String,ContactModelParams}; end`.
  - `parse_contact(raw, path::String, ctx::ParseContext)::Union{Nothing,ContactInput}` — `Contact.Globals` is parsed as `ContactGlobalsParams` (defaults if absent); every other entry is a contact model.

- [ ] **Step 1: Write the failing test**

`test/unit_tests/Support/Parameters/Input/ut_outputs_conditions_contact.jl`:

```julia
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
```

Update `input_tests.jl`:

```julia
for file in ["ut_mesh_blocks.jl", "ut_solver.jl", "ut_outputs_conditions_contact.jl"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl`
Expected: FAIL in `ut_outputs_conditions_contact.jl` with `UndefVarError: OutputParams not defined in PeriLab.InputDeck`.

- [ ] **Step 3: Write the implementation**

In `InputDeck.jl`, replace

```julia
include("solver.jl")
```

with

```julia
include("solver.jl")
include("outputs.jl")
include("conditions.jl")
include("contact.jl")
```

`src/Support/Parameters/Input/outputs.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

@params struct OutputParams
    output_filename::String = req("Output Filename")
    output_file_type::String = opt("Output File Type"; default = "Exodus",
                                   allowed = ["Exodus", "CSV"])
    output_frequency::Union{Nothing,Int64,String} = opt("Output Frequency"; default = nothing,
                                                        description = "Every n-th step; per step as \"n1 n2\"")
    number_of_output_steps::Union{Nothing,Int64,String} = opt("Number of Output Steps";
                                                              default = nothing)
    output_variables::Dict{String,Bool} = req("Output Variables")
    flush_file::Bool = opt("Flush File"; default = true)
    write_after_damage::Bool = opt("Write After Damage"; default = false)
    start_time::Float64 = opt("Start Time"; default = 0.0, quantity = :time)
    end_time::Float64 = opt("End Time"; default = Inf, quantity = :time)
    bond_export::Bool = opt("Bond Export"; default = false)
    bond_blocks::Union{Nothing,Int64} = opt("Bond Blocks"; default = nothing)
end

function check!(p::OutputParams, path::String, ctx::ParseContext)
    if p.output_frequency === nothing && p.number_of_output_steps === nothing
        add_error!(ctx, path, "\"Output Frequency\" or \"Number of Output Steps\" is required")
    end
    return nothing
end
```

`src/Support/Parameters/Input/conditions.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

@params struct BoundaryConditionParams
    variable::String = req("Variable")
    node_set::String = req("Node Set"; description = "Node set name; several joined with +")
    value::Union{Float64,String} = req("Value"; description = "Number or expression, e.g. \"10*t\"")
    type::Union{Nothing,String} = opt("Type"; default = nothing,
                                      allowed = ["Initial", "Dirichlet", "Neumann"])
    coordinate::Union{Nothing,String} = opt("Coordinate"; default = nothing)
    step_id::Union{Nothing,Int64,String} = opt("Step ID"; default = nothing)
end

@params struct ComputeClassParams
    compute_class::String = req("Compute Class")
    variable::String = req("Variable")
    calculation_type::Union{Nothing,String} = opt("Calculation Type"; default = nothing)
    block::Union{Nothing,String} = opt("Block"; default = nothing)
    node_set::Union{Nothing,String} = opt("Node Set"; default = nothing)
    equation::Union{Nothing,String} = opt("Equation"; default = nothing)
    x::Union{Nothing,Float64} = opt("X"; default = nothing)
    y::Union{Nothing,Float64} = opt("Y"; default = nothing)
    z::Union{Nothing,Float64} = opt("Z"; default = nothing)
end
```

`src/Support/Parameters/Input/contact.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# In the Contact section, "Globals" is real configuration (not the generic
# escape hatch); every other entry is a contact model.

@params struct ContactGlobalsParams
    global_search_frequency::Int64 = opt("Global Search Frequency"; default = 1, min = 1)
    only_surface_contact_nodes::Bool = opt("Only Surface Contact Nodes"; default = true)
end

@params struct ContactGroupParams
    master_block_id::Int64 = req("Master Block ID"; min = 1)
    slave_block_id::Int64 = req("Slave Block ID"; min = 1)
    search_radius::Float64 = req("Search Radius"; min = 0, quantity = :length)
    global_search_frequency::Union{Nothing,Int64} = opt("Global Search Frequency";
                                                        default = nothing, min = 1)
end

@params struct ContactModelParams
    type::String = req("Type")
    contact_radius::Float64 = req("Contact Radius"; min = 0, quantity = :length)
    contact_stiffness::Float64 = req("Contact Stiffness"; min = 0)
    friction_coefficient::Union{Nothing,Float64} = opt("Friction Coefficient";
                                                       default = nothing, min = 0)
    symmetry::Union{Nothing,String} = opt("Symmetry"; default = nothing)
    contact_groups::Dict{String,ContactGroupParams} = req("Contact Groups")
end

struct ContactInput
    globals::ContactGlobalsParams
    models::Dict{String,ContactModelParams}
end

function _contact_section(raw, path::String, ctx::ParseContext)
    raw isa AbstractDict && return true
    add_error!(ctx, path,
               "expected a section of `key: value` entries, got $(ParameterSpec._describe(raw))")
    return false
end

"""
    parse_contact(raw, path, ctx)

Reads the `Contact` section: `Globals` (search settings for all contact
models, defaults if absent) and the named contact models.
"""
function parse_contact(raw, path::String, ctx::ParseContext)
    _contact_section(raw, path, ctx) || return nothing
    globals_path = join_path(path, "Globals")
    raw_globals = get(raw, "Globals", Dict{String,Any}())
    globals = _contact_section(raw_globals, globals_path, ctx) ?
              parse_section(ContactGlobalsParams, raw_globals, globals_path, ctx) : nothing
    models = Dict{String,ContactModelParams}()
    ok = globals !== nothing
    for (name, entry) in raw
        key = string(name)
        key == "Globals" && continue
        model_path = join_path(path, key)
        if !_contact_section(entry, model_path, ctx)
            ok = false
            continue
        end
        model = parse_section(ContactModelParams, entry, model_path, ctx)
        if model === nothing
            ok = false
        else
            models[key] = model
        end
    end
    return ok ? ContactInput(globals, models) : nothing
end
```

In `InputDeck.jl`, `parse_contact` needs `ParameterSpec` itself: replace

```julia
using ..ParameterSpec: @params, req, opt, ParseContext, add_error!, join_path,
                       parse_section
```

with

```julia
using ..ParameterSpec
using ..ParameterSpec: @params, req, opt, ParseContext, add_error!, join_path,
                       parse_section
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/Support/Parameters/Input test/unit_tests/Support/Parameters/Input
git commit -m "InputDeck: Outputs, Boundary Conditions, Compute Classes, Contact"
```

---

### Task 5: `PeriLabInput` and `read_input`

**Files:**
- Create: `src/Support/Parameters/Input/input.jl`
- Modify: `src/Support/Parameters/Input/InputDeck.jl`
- Create: `test/unit_tests/Support/Parameters/Input/ut_read_input.jl`
- Modify: `test/unit_tests/Support/Parameters/Input/input_tests.jl`

**Interfaces:**
- Consumes: Tasks 2–4.
- Produces:
  - `PeriLabSections` (`@params`): `discretization`, `blocks::Dict{String,BlockParams}`, `fem`, `surface_correction`, `solver::Union{Nothing,SolverParams}`, `multistep_solver::Dict{String,SolverParams}`, `outputs::Dict{String,OutputParams}`, `boundary_conditions::Dict{String,BoundaryConditionParams}`, `compute_class_parameters::Dict{String,ComputeClassParams}`, `strict_validation::Bool`; `check!`: at least one block; `Solver` or `Multistep Solver`; `Solver` needs `Initial Time` and `Final Time`; every multistep step needs `Step ID`.
  - `struct PeriLabInput; sections::PeriLabSections; contact::Union{Nothing,ContactInput}; models::Dict{String,Any}; globals::Dict{String,Any}; end`
  - `read_input(deck::AbstractDict, directory::AbstractString = ""; strict::Bool = true)::Tuple{Union{Nothing,PeriLabInput},ParseContext}` — `deck` is the content of the `PeriLab:` key. Does not abort; the caller reports.
  - exported: `read_input`, `PeriLabInput`.

- [ ] **Step 1: Write the failing test**

`test/unit_tests/Support/Parameters/Input/ut_read_input.jl`:

```julia
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
                            "Models" => Dict{String,Any}("Material Models" => Dict{String,Any}("Mat" => Dict{String,Any}("Material Model" => "PD Solid Elastic"))),
                            "Solver" => Dict{String,Any}("Initial Time" => 0.0,
                                                         "Final Time" => 1.0,
                                                         "Verlet" => Dict{String,Any}()))
end

@testset "minimal deck" begin
    input, ctx = ID.read_input(ut_deck())
    @test isempty(ctx.errors)
    @test input isa ID.PeriLabInput
    @test input.sections.blocks["block_1"].block_id === 1
    @test input.sections.solver.verlet.safety_factor === 1.0
    @test input.contact === nothing
    @test haskey(input.models, "Material Models")     # models stay raw until phase 3
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
    @test input.contact.models["C"].type == "Penalty Contact"
end
```

Update `input_tests.jl`:

```julia
for file in ["ut_mesh_blocks.jl", "ut_solver.jl", "ut_outputs_conditions_contact.jl",
             "ut_read_input.jl"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl`
Expected: FAIL in `ut_read_input.jl` with `UndefVarError: read_input not defined in PeriLab.InputDeck`.

- [ ] **Step 3: Write the implementation**

In `InputDeck.jl`, replace

```julia
include("contact.jl")
```

with

```julia
include("contact.jl")
include("input.jl")

export read_input, PeriLabInput
```

`src/Support/Parameters/Input/input.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

@params struct PeriLabSections
    discretization::DiscretizationParams = req("Discretization")
    blocks::Dict{String,BlockParams} = req("Blocks")
    fem::Union{Nothing,FEMParams} = opt("FEM"; default = nothing)
    surface_correction::Union{Nothing,SurfaceCorrectionParams} = opt("Surface Correction";
                                                                     default = nothing)
    solver::Union{Nothing,SolverParams} = opt("Solver"; default = nothing)
    multistep_solver::Dict{String,SolverParams} = opt("Multistep Solver";
                                                      default = Dict{String,SolverParams}())
    outputs::Dict{String,OutputParams} = opt("Outputs"; default = Dict{String,OutputParams}())
    boundary_conditions::Dict{String,BoundaryConditionParams} = opt("Boundary Conditions";
                                                                    default = Dict{String,
                                                                                   BoundaryConditionParams}())
    compute_class_parameters::Dict{String,ComputeClassParams} = opt("Compute Class Parameters";
                                                                    default = Dict{String,
                                                                                   ComputeClassParams}())
    strict_validation::Bool = opt("Strict Validation"; default = true,
                                  description = "false reports unknown keys as warnings")
end

function check!(p::PeriLabSections, path::String, ctx::ParseContext)
    isempty(p.blocks) && add_error!(ctx, join_path(path, "Blocks"),
                                    "at least one block is required")
    if p.solver === nothing && isempty(p.multistep_solver)
        add_error!(ctx, path, "\"Solver\" or \"Multistep Solver\" is required")
    end
    if p.solver !== nothing
        solver_path = join_path(path, "Solver")
        for (key, value) in (("Initial Time", p.solver.initial_time),
                             ("Final Time", p.solver.final_time))
            value === nothing &&
                add_error!(ctx, join_path(solver_path, key), "missing (required for \"Solver\")")
        end
    end
    for (name, step) in p.multistep_solver
        if step.step_id === nothing
            add_error!(ctx,
                       join_path(join_path(join_path(path, "Multistep Solver"), name), "Step ID"),
                       "missing (required in a multistep solver step)")
        end
    end
    return nothing
end

"""
    PeriLabInput

A validated input deck. `models` stays the raw `Models` dict until the model
modules declare their parameters (phase 3); `globals` is the unvalidated
`Globals` escape hatch.
"""
struct PeriLabInput
    sections::PeriLabSections
    contact::Union{Nothing,ContactInput}
    models::Dict{String,Any}
    globals::Dict{String,Any}
end

const _SPECIAL_KEYS = ("Models", "Contact", "Globals")

"""
    read_input(deck, directory = ""; strict = true) -> (input, ctx)

Reads the content of the `PeriLab:` key of an input deck. `directory` is the
deck's directory (for relative data files). All problems are collected in
`ctx`; `input` is `nothing` if any of them is an error. Does not abort —
call `ParameterSpec.report!(ctx)` for that.
"""
function read_input(deck::AbstractDict, directory::AbstractString = ""; strict::Bool = true)
    ctx = ParseContext(directory = directory, strict = strict)
    plain = Dict{String,Any}(string(k) => v for (k, v) in deck
                             if !(string(k) in _SPECIAL_KEYS))
    sections = parse_section(PeriLabSections, plain, "", ctx)
    contact = haskey(deck, "Contact") ? parse_contact(deck["Contact"], "Contact", ctx) :
              nothing
    models = get(deck, "Models", nothing)
    if models === nothing
        add_error!(ctx, "Models", "missing (required)")
    elseif !(models isa AbstractDict)
        add_error!(ctx, "Models",
                   "expected a section of `key: value` entries, got $(ParameterSpec._describe(models))")
    end
    globals = get(deck, "Globals", Dict{String,Any}())
    if ParameterSpec.has_errors(ctx) || sections === nothing
        return nothing, ctx
    end
    input = PeriLabInput(sections, contact, Dict{String,Any}(string(k) => v for (k, v) in models),
                         globals isa AbstractDict ?
                         Dict{String,Any}(string(k) => v for (k, v) in globals) :
                         Dict{String,Any}())
    return input, ctx
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/Support/Parameters/Input test/unit_tests/Support/Parameters/Input
git commit -m "InputDeck: PeriLabInput and read_input"
```

---

### Task 6: Golden compatibility test over all shipped decks; remove dead keys

**Files:**
- Create: `test/unit_tests/Support/Parameters/Input/ut_golden_decks.jl`
- Modify: `test/unit_tests/Support/Parameters/Input/input_tests.jl`
- Modify (remove dead keys, nothing else):
  - `examples/Dogbone/Dogbone_plastic_corr_static.yaml` — `NLsolver: true` (line 108)
  - `test/fullscale_tests/test_Craig_Bampton/cb_mat_point.yaml` — all `Compute Property:` lines
  - `test/fullscale_tests/test_Craig_Bampton/cb_without_mat_point.yaml` — all `Compute Property:` lines
  - `test/fullscale_tests/test_DCB/DCBmodel_PD_solid_switch.yaml` — `Solve For Displacement: True` (line 101) and the step-level `Numerical Damping: 5.0e-06` (line 117, 6-space indent)
  - `test/fullscale_tests/test_thermal_expansion/thermal_expansion_lin_mat_const_T.yaml` and `thermal_expansion_lin_mat_var_T.yaml` — `Update Matrix:` lines

**Interfaces:**
- Consumes: Task 5 `read_input`; `PeriLab.IO.read_input(filename)` (existing, returns the raw YAML dict).
- Produces: a test that fails whenever a shipped deck stops parsing in strict mode.

Background (verified while writing this plan): no code reads `NLsolver`, `Compute Property`, `Solve For Displacement`, `Update Matrix` (code reads `Matrix Update`), or step-level `Numerical Damping` (code reads it from the solver section). They are ignored today with a "Key not known" warning, so removing them does not change results. `test/fullscale_tests/test_PD_Solid_Elastic/strain_xx_external.yaml` uses an `External` solver that does not exist and is not run by any test; it is allowlisted, not edited.

- [ ] **Step 1: Write the failing test**

`test/unit_tests/Support/Parameters/Input/ut_golden_decks.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const ID = PeriLab.InputDeck

const UT_REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", "..", ".."))

# Decks that are known not to parse, with the reason. Fix or document them;
# never add a deck here to hide a validator bug.
const UT_DECK_ALLOWLIST = Dict("test/fullscale_tests/test_PD_Solid_Elastic/strain_xx_external.yaml" => "uses a non-existent \"External\" solver; not run by any test")

function ut_all_decks()
    decks = String[]
    for root in ("examples", "test"), (dir, _, files) in walkdir(joinpath(UT_REPO_ROOT, root))
        for file in files
            endswith(file, ".yaml") && push!(decks, joinpath(dir, file))
        end
    end
    return sort!(decks)
end

@testset "every shipped deck parses in strict mode" begin
    decks = ut_all_decks()
    @test length(decks) > 100
    for file in decks
        relative = relpath(file, UT_REPO_ROOT)
        haskey(UT_DECK_ALLOWLIST, relative) && continue
        raw = PeriLab.IO.read_input(file)
        if !(raw isa AbstractDict && haskey(raw, "PeriLab"))
            continue                        # not an input deck (e.g. a config file)
        end
        _, ctx = ID.read_input(raw["PeriLab"], dirname(file); strict = true)
        errors = filter(e -> e.severity == :error, ctx.errors)
        if !isempty(errors)
            @error "$relative\n" * PeriLab.ParameterSpec.format_errors(errors)
        end
        @test isempty(errors)
    end
end
```

Update `input_tests.jl`:

```julia
for file in ["ut_mesh_blocks.jl", "ut_solver.jl", "ut_outputs_conditions_contact.jl",
             "ut_read_input.jl", "ut_golden_decks.jl"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl > /tmp/golden.log 2>&1; grep -A4 "Error: .*yaml" /tmp/golden.log`
Expected: FAIL, with errors for exactly the six decks listed under **Files** (unknown keys `NLsolver`, `Compute Property`, `Solve For Displacement`, `Numerical Damping` under a multistep step, `Update Matrix`). If any other deck fails, classify each error before continuing: (a) the key is read by code → add it to the matching struct (find the reader with `grep -rn '"<key>"' src`), with the reader's default; (b) no code reads it → delete the line from the deck; (c) the deck cannot run at all → allowlist it with the reason. Record each such decision as a ruling.

- [ ] **Step 3: Remove the dead keys**

```bash
sed -i '108{/NLsolver: true/d}' examples/Dogbone/Dogbone_plastic_corr_static.yaml
sed -i '/^\s*Compute Property:/d' test/fullscale_tests/test_Craig_Bampton/cb_mat_point.yaml test/fullscale_tests/test_Craig_Bampton/cb_without_mat_point.yaml
sed -i -e '101{/Solve For Displacement/d}' -e '117{/^      Numerical Damping:/d}' test/fullscale_tests/test_DCB/DCBmodel_PD_solid_switch.yaml
sed -i '/^\s*Update Matrix:/d' test/fullscale_tests/test_thermal_expansion/thermal_expansion_lin_mat_const_T.yaml test/fullscale_tests/test_thermal_expansion/thermal_expansion_lin_mat_var_T.yaml
git diff --stat -- examples test/fullscale_tests
```

Expected `git diff --stat`: 6 files changed, only deletions (`Compute Property` lines: 12 per Craig-Bampton deck).

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add examples test/fullscale_tests test/unit_tests/Support/Parameters/Input
git commit -m "Golden test: all shipped decks parse in strict mode; remove dead keys"
```

---

### Task 7: Validate the live input path with `read_input`; `--no_strict` flag

**Files:**
- Modify: `src/Support/Parameters/parameter_handling.jl` (`validate_yaml`, new `validate_models`)
- Modify: `src/IO/read_inputdeck.jl` (`read_input_file`)
- Modify: `src/IO/IO.jl` (`initialize_data`)
- Modify: `src/PeriLab.jl` (`parse_commandline`, `main`, `run`)
- Modify: `test/unit_tests/Support/Parameters/ut_parameter_handling.jl` (`ut_validate_yaml`)
- Modify: `test/unit_tests/ut_perilab.jl` (`ut_parse_commandline`)
- Create: `test/unit_tests/Support/Parameters/Input/ut_validate_yaml.jl`
- Modify: `test/unit_tests/Support/Parameters/Input/input_tests.jl`

**Interfaces:**
- Consumes: Task 5 `read_input`; phase 1 `ParameterSpec.report!`, `ParameterSpec.strict_mode`, `ParameterSpec.add_error!`; legacy `validate_structure_recursive`, `get_all_keys`, `expected_structure`.
- Produces:
  - `validate_yaml(params::Dict; directory::AbstractString = "", no_strict::Bool = false)` — returns `params["PeriLab"]` (unchanged dict) or aborts with all errors.
  - `validate_models(deck::AbstractDict)::Bool` — legacy validation of `Models` only (warnings for unknown keys as before).
  - `read_input_file(filename::String; directory::AbstractString = dirname(filename), no_strict::Bool = false)`
  - `IO.initialize_data(filename, filedirectory, comm; no_strict::Bool = false)`
  - `PeriLab.run(filename; ..., no_strict::Bool = false)`; command line flag `--no_strict`.

- [ ] **Step 1: Write the failing tests**

`test/unit_tests/Support/Parameters/Input/ut_validate_yaml.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test

function ut_valid_params()
    return Dict{Any,Any}("PeriLab" => Dict{Any,Any}("Discretization" => Dict{Any,Any}("Type" => "Text File",
                                                                                     "Input Mesh File" => "m.txt"),
                                                    "Blocks" => Dict{Any,Any}("block_1" => Dict{Any,Any}("Block ID" => 1,
                                                                                                         "Density" => 1.0,
                                                                                                         "Horizon" => 1.0)),
                                                    "Models" => Dict{Any,Any}("Material Models" => Dict{Any,Any}("mat_1" => Dict{Any,Any}("Material Model" => "a"))),
                                                    "Solver" => Dict{Any,Any}("Initial Time" => 0.0,
                                                                              "Final Time" => 1.0,
                                                                              "Verlet" => Dict{Any,Any}())))
end

@testset "validate_yaml returns the deck dict unchanged" begin
    params = ut_valid_params()
    @test PeriLab.Parameter_Handling.validate_yaml(params) === params["PeriLab"]
end

@testset "validate_yaml aborts with all input errors" begin
    params = ut_valid_params()
    params["PeriLab"]["Blocks"]["block_1"]["Horizon"] = "1.0"
    params["PeriLab"]["Bocks"] = Dict{Any,Any}()
    @test_throws PeriLab.PeriLabError PeriLab.Parameter_Handling.validate_yaml(params)
end

@testset "no_strict flag" begin
    params = ut_valid_params()
    params["PeriLab"]["Solver"]["Unused Option"] = 1
    @test_throws PeriLab.PeriLabError PeriLab.Parameter_Handling.validate_yaml(params)
    @test PeriLab.Parameter_Handling.validate_yaml(params; no_strict = true) ===
          params["PeriLab"]
    params["PeriLab"]["Solver"]["final time"] = 2.0       # case-only typo: always an error
    @test_throws PeriLab.PeriLabError PeriLab.Parameter_Handling.validate_yaml(params;
                                                                              no_strict = true)
    @test PeriLab.parse_commandline(["--no_strict", "a.yaml"])["no_strict"] == true
    @test PeriLab.parse_commandline(["a.yaml"])["no_strict"] == false
end

@testset "legacy model validation still applies" begin
    params = ut_valid_params()
    params["PeriLab"]["Models"]["Material Models"]["mat_1"]["Material Model"] = 5
    @test_throws PeriLab.PeriLabError PeriLab.Parameter_Handling.validate_yaml(params)
end
```

Update `input_tests.jl`:

```julia
for file in ["ut_mesh_blocks.jl", "ut_solver.jl", "ut_outputs_conditions_contact.jl",
             "ut_read_input.jl", "ut_golden_decks.jl", "ut_validate_yaml.jl"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl`
Expected: FAIL in `ut_validate_yaml.jl` — "aborts with all input errors" does not throw for `"Bocks"` (the old validator only warns), and `validate_yaml(params; no_strict = true)` errors with `MethodError: no method matching validate_yaml(...; no_strict::Bool)`.

- [ ] **Step 3: Write the implementation**

In `src/Support/Parameters/parameter_handling.jl`, replace

```julia
using ...PeriLabExceptions: @abort
```

with

```julia
using ...PeriLabExceptions: @abort
using ..ParameterSpec: report!, strict_mode, add_error!
using ..InputDeck: read_input
```

and replace the whole `function validate_yaml(params::Dict) … end` (the last function before the module's closing `end`) with:

```julia
"""
    validate_models(deck)

Legacy validation of the `Models` section (types and required keys from
`expected_structure`; unknown keys are warnings). Replaced in phase 3, when
model modules declare their parameters.
"""
function validate_models(deck::AbstractDict)
    models = get(deck, "Models", nothing)
    models isa Dict || return true            # missing / malformed: reported by read_input
    checked_keys = []
    valid = try
        first(validate_structure_recursive(expected_structure["PeriLab"][1]["Models"][1],
                                           models, true, checked_keys, "Models"))
    catch
        false
    end
    for key in get_all_keys(models)
        key isa Int64 && continue
        if !(key in checked_keys) && !contains(key, "Property_")
            @warn "Key not known - $key, going to ignore it"
        end
    end
    return valid
end

"""
    validate_yaml(params; directory = "", no_strict = false)

Validates a loaded input deck against the typed input declarations
(`InputDeck.read_input`) and, for `Models`, the legacy structure. Reports
every problem at once and aborts if there is an error; otherwise returns
`params["PeriLab"]` unchanged.
"""
function validate_yaml(params::Dict; directory::AbstractString = "", no_strict::Bool = false)
    if !haskey(params, "PeriLab") || !(params["PeriLab"] isa AbstractDict) ||
       length(params["PeriLab"]) < 2
        @abort "Yaml file is not valid."
        return
    end
    deck = params["PeriLab"]
    _, ctx = read_input(deck, directory; strict = strict_mode(deck; no_strict_flag = no_strict))
    validate_models(deck) ||
        add_error!(ctx, "Models", "invalid model parameters (see the warnings above)")
    report!(ctx)
    return deck
end
```

In `src/IO/read_inputdeck.jl`, replace

```julia
function read_input_file(filename::String)
```

with

```julia
function read_input_file(filename::String; directory::AbstractString = dirname(filename),
                         no_strict::Bool = false)
```

and, in the same function, replace

```julia
    return validate_yaml(read_input(filename))
```

with

```julia
    return validate_yaml(read_input(filename); directory = directory, no_strict = no_strict)
```

In `src/IO/IO.jl`, replace

```julia
function initialize_data(filename::String,
                         filedirectory::String,
                         comm::MPI.Comm)
```

with

```julia
function initialize_data(filename::String,
                         filedirectory::String,
                         comm::MPI.Comm;
                         no_strict::Bool = false)
```

and, in the same function, replace

```julia
    @timeit "init_data" params=init_data(read_input_file(filename),
                                         filedirectory, comm)
```

with

```julia
    @timeit "init_data" params=init_data(read_input_file(filename;
                                                         directory = filedirectory,
                                                         no_strict = no_strict),
                                         filedirectory, comm)
```

In `src/PeriLab.jl`, in `parse_commandline`, replace

```julia
        "--reload", "-r"
        help = "reload"
        action = :store_true
```

with

```julia
        "--reload", "-r"
        help = "reload"
        action = :store_true
        "--no_strict"
        help = "report unknown input keys as warnings instead of errors"
        action = :store_true
```

in `main`, replace

```julia
                         reload = parsed_args["reload"])
```

with

```julia
                         reload = parsed_args["reload"],
                         no_strict = parsed_args["no_strict"])
```

in `run`, replace

```julia
             silent::Bool = false,
             reload::Bool = false,)
```

with

```julia
             silent::Bool = false,
             reload::Bool = false,
             no_strict::Bool = false)
```

and in the same function replace

```julia
            @timeit "IO.initialize_data" params,
                                         steps=IO.initialize_data(filename,
                                                                  filedirectory,
                                                                  comm)
```

with

```julia
            @timeit "IO.initialize_data" params,
                                         steps=IO.initialize_data(filename,
                                                                  filedirectory,
                                                                  comm;
                                                                  no_strict = no_strict)
```

Update the docstring of `run` (the `- reload::Bool=false` bullet list) by adding after the `reload` line:

```julia
- `no_strict::Bool=false`: Report unknown input keys as warnings instead of errors.
```

Update the old tests in `test/unit_tests/Support/Parameters/ut_parameter_handling.jl` (`@testset "ut_validate_yaml"`) to the new contract:
- The three cases whose `PeriLab` dict has fewer than two keys (empty `params`; only `Models`; only `Blocks`) keep `(:error, "Yaml file is not valid.")`.
- Every other `@test_logs (:error, "Yaml file is not valid.") @test_throws …` becomes `@test_logs (:error, r"^Input errors") match_mode=:any @test_throws …`.
- In the final, expected-valid case, delete `"Block Names" => "Block_1",` from `Block_1` and add `"Verlet" => Dict{Any,Any}()` to `"Solver"` (a deck without a solver type and with an unknown key is no longer valid).

Update `test/unit_tests/ut_perilab.jl` (`@testset "ut_parse_commandline"`), which compares the complete parsed-argument dict: add `"no_strict" => false` to both expected `Dict(...)` literals.

- [ ] **Step 4: Run tests to verify they pass**

Run: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl`
Expected: PASS.
Run: `julia --project=. test/unit_tests/Support/Parameters/Spec/run_spec_tests.jl`
Expected: PASS.
Run (background, ≈30 min): `julia --project=. -e 'using Pkg; Pkg.test()' > /tmp/full.log 2>&1; tail -5 /tmp/full.log`
Expected: `Testing PeriLab tests passed`; every fullscale deck now goes through strict validation.

- [ ] **Step 5: Commit**

```bash
git add src test examples
git commit -m "Validate input decks with typed sections; add --no_strict flag"
```
