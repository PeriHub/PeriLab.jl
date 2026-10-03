# Typed Input Parameters — Phase 1 (Core Machinery) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the self-contained `ParameterSpec` module: the `@params` declaration macro, typed value conversion with constraints, dependent values (`Constant` / `Table1D`), the model registry, `"A + B"` composition, field binding, and error collection with suggestions — without changing PeriLab's existing input path.

**Architecture:** One new Julia submodule `PeriLab.ParameterSpec` under `src/Support/Parameters/Spec/`, split into focused files (errors, dependent values, field declarations and conversion, macro, building, registry and models, binding). Module authors declare parameters with `@params struct ... end`; the macro emits a plain (possibly parametric) struct plus a `parameter_spec` method. Everything else (validation, building from YAML dicts, binding) is ordinary functions reading that spec. Nothing outside the new module and the test runner is touched except one `include` line in `src/PeriLab.jl`.

**Tech Stack:** Julia 1.12, Dierckx (existing dependency, splines), Test stdlib. No new dependencies.

**Spec:** `docs/superpowers/specs/2026-10-03-typed-input-parameters-design.md` (this plan implements spec §5 phase 1: §2.1, §2.2, §2.3 without the top-level `PeriLabInput`, §2.4, §2.6, and the module-facing parts of §2.7). Phases 2–5 get their own plans once this code exists.

## Global Constraints

- Julia `1.12` (as in `Project.toml` `[compat] julia = "1.12"`). No new packages in `Project.toml`.
- **Do not commit or stage anything.** The user asked for no commits; leave all changes in the working tree.
- Phase 1 must not change behaviour of the existing input path: `parameter_handling*.jl`, `read_inputdeck.jl`, `Data_manager.jl` and all model code stay untouched.
- YAML keys are matched exactly as written in existing decks (spaces, apostrophes, capitalisation); no YAML format change.
- No unit handling: `quantity` is documentation only, never converted or checked.
- Model registration must happen at runtime, never during precompilation: package modules call `register_*` from their `__init__()` function (a top-level call in a precompiled module would be lost). Runtime-loaded (licensed) modules may do the same.
- Every new source file starts with the repository's SPDX header:
  ```
  # SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
  #
  # SPDX-License-Identifier: BSD-3-Clause
  ```
- Error path format: segments joined with `.`; a segment that is not purely `[A-Za-z0-9_]` is wrapped in double quotes, e.g. `Models."Material Models".Steel."Young's Modulus"`.
- Test command (from repository root), used by every task:
  `julia --project=. test/unit_tests/Support/Parameters/Spec/run_spec_tests.jl`
  Test files access the module as `const PS = PeriLab.ParameterSpec` and never `using` its exports (avoids name clashes in `Main` with other test files). Note that `test/runtests.jl` disables warn-level logging globally, so tests must not rely on `@test_logs`.

## Review Focus

1. A YAML key present with no value (`Horizon:` → `nothing`) must give "expected a number, got an empty value", not a crash or a "missing" message — test in Task 4 ("missing and empty values").
2. Integral numbers across numeric types (`100.0` for an `Int64` field, `1` for a `Float64` field) must be accepted — tests in Task 3 (Int64 / Float64) and Task 4 ("Int field accepts integral float").
3. Data files written on Windows or by hand (CRLF line endings, tabs, blank lines) must parse — test in Task 2 ("read_table tolerates CRLF, tabs and blank lines").
4. A model name with an empty `+` part (`"Elastic + "`) must give a clear error, not look up `""` — test in Task 5 ("model name errors").
5. A key differing only in case or punctuation (`"radius"` vs `"Radius"`) must be rejected with a suggestion, never silently accepted — test in Task 4 ("case variants are not accepted silently").

---

## File Structure

| File | Responsibility |
|---|---|
| `src/Support/Parameters/Spec/ParameterSpec.jl` | Module definition, imports, `include`s, exports |
| `src/Support/Parameters/Spec/errors.jl` | `ParamsDefinitionError`, `InputError`, `ParseContext`, path joining, suggestions, error report, `strict_mode` |
| `src/Support/Parameters/Spec/dependent.jl` | `Dependent`, `Constant`, `Table1D`, `value`, `bind_table!`, `read_table`, `combine` |
| `src/Support/Parameters/Spec/field_spec.jl` | `req` / `opt` → `FieldDecl`, `FieldSpec` |
| `src/Support/Parameters/Spec/convert.jl` | `convert_value` (YAML value → declared type), `enum_aliases`, `check_constraints!` |
| `src/Support/Parameters/Spec/params_macro.jl` | `@params`, `build_spec`, `supported_type`, `is_params`, `parameter_spec`, `derive` |
| `src/Support/Parameters/Spec/build.jl` | `build`, `parse_section`, `check_unknown!`, `aliases` |
| `src/Support/Parameters/Spec/registry.jl` | Model registry, `register_*`, `lookup_model`, `registered_names` |
| `src/Support/Parameters/Spec/model.jl` | `NoModel`, `Composite`, `parse_model` |
| `src/Support/Parameters/Spec/bind.jl` | `bind_dependents!` |
| `src/PeriLab.jl` (modify) | one `include` line |
| `test/unit_tests/Support/Parameters/Spec/run_spec_tests.jl` | standalone runner |
| `test/unit_tests/Support/Parameters/Spec/spec_tests.jl` | list of spec test files (shared by runner and `test/runtests.jl`) |
| `test/unit_tests/Support/Parameters/Spec/ut_*.jl` | one test file per task |
| `test/runtests.jl` (modify) | include `spec_tests.jl` |

---

### Task 1: Module skeleton, error collection, suggestions, strict mode

**Files:**
- Create: `src/Support/Parameters/Spec/ParameterSpec.jl`
- Create: `src/Support/Parameters/Spec/errors.jl`
- Modify: `src/PeriLab.jl:39` (add include after `include("./IO/exceptions.jl")`)
- Create: `test/unit_tests/Support/Parameters/Spec/run_spec_tests.jl`
- Create: `test/unit_tests/Support/Parameters/Spec/spec_tests.jl`
- Create: `test/unit_tests/Support/Parameters/Spec/ut_errors.jl`
- Modify: `test/runtests.jl:49-53` (Support → Parameters testset)

**Interfaces:**
- Consumes: `PeriLab.PeriLabExceptions.@abort`, `PeriLabError` (existing, `src/IO/exceptions.jl`).
- Produces:
  - `struct ParamsDefinitionError <: Exception; msg::String; end`
  - `struct InputError; path::String; message::String; severity::Symbol; end` (`:error` / `:warning`)
  - `mutable struct ParseContext; directory::String; strict::Bool; errors::Vector{InputError}; end`, constructor `ParseContext(; directory = "", strict = true)`
  - `add_error!(ctx, path, msg)`, `add_warning!(ctx, path, msg)`, `has_errors(ctx)::Bool`
  - `join_path(path::AbstractString, key::AbstractString)::String`
  - `_normalize(s)::String`, `levenshtein(a, b)::Int`, `suggest(key, candidates)::Union{Nothing,String}`, `unknown_key_message(key, candidates)::String`
  - `format_errors(errors::Vector{InputError})::String`, `report!(ctx)::Nothing` (aborts via `@abort` if any error)
  - `strict_mode(input::AbstractDict; no_strict_flag::Bool = false)::Bool`

- [ ] **Step 1: Create the test runner and test list**

`test/unit_tests/Support/Parameters/Spec/run_spec_tests.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# Standalone runner for the ParameterSpec unit tests:
#   julia --project=. test/unit_tests/Support/Parameters/Spec/run_spec_tests.jl
using Test
import PeriLab

@testset "ParameterSpec" begin
    include(joinpath(@__DIR__, "spec_tests.jl"))
end
```

`test/unit_tests/Support/Parameters/Spec/spec_tests.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

for file in ["ut_errors.jl"]
    @testset "$file" begin
        include(joinpath(@__DIR__, file))
    end
end
```

In `test/runtests.jl`, replace

```julia
            @testset "Parameters" begin
                @testset "ut_parameter_handling" begin
                    include("unit_tests/Support/Parameters/ut_parameter_handling.jl")
                end
            end
```

with

```julia
            @testset "Parameters" begin
                @testset "ut_parameter_handling" begin
                    include("unit_tests/Support/Parameters/ut_parameter_handling.jl")
                end
                @testset "ParameterSpec" begin
                    include("unit_tests/Support/Parameters/Spec/spec_tests.jl")
                end
            end
```

- [ ] **Step 2: Write the failing test**

`test/unit_tests/Support/Parameters/Spec/ut_errors.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const PS = PeriLab.ParameterSpec

@testset "join_path" begin
    @test PS.join_path("", "Blocks") == "Blocks"
    @test PS.join_path("Blocks", "block_1") == "Blocks.block_1"
    @test PS.join_path("Models", "Material Models") == "Models.\"Material Models\""
    @test PS.join_path("M", "Young's Modulus") == "M.\"Young's Modulus\""
end

@testset "levenshtein" begin
    @test PS.levenshtein("kitten", "sitting") == 3
    @test PS.levenshtein("", "abc") == 3
    @test PS.levenshtein("same", "same") == 0
end

@testset "suggest" begin
    candidates = ["Young's Modulus", "Poisson's Ratio", "Density"]
    @test PS.suggest("Poissons Ratio", candidates) == "Poisson's Ratio"
    @test PS.suggest("young's modulus", candidates) == "Young's Modulus"
    @test PS.suggest("Youngs Modulos", candidates) == "Young's Modulus"
    @test PS.suggest("Horizon", candidates) === nothing
    @test PS.suggest("x", String[]) === nothing
    @test PS.unknown_key_message("Densty", candidates) ==
          "unknown key — did you mean \"Density\"?"
    @test PS.unknown_key_message("Horizon", candidates) == "unknown key"
end

@testset "ParseContext and report!" begin
    ctx = PS.ParseContext(directory = "/tmp", strict = true)
    @test ctx.directory == "/tmp"
    @test ctx.strict
    @test !PS.has_errors(ctx)
    PS.add_warning!(ctx, "a", "only a warning")
    @test !PS.has_errors(ctx)
    @test PS.report!(ctx) === nothing
    PS.add_error!(ctx, "Blocks.block_1.Horizon", "-0.1 is below minimum 0")
    PS.add_error!(ctx, "x", "missing")
    @test PS.has_errors(ctx)
    errors = filter(e -> e.severity == :error, ctx.errors)
    @test PS.format_errors(errors) ==
          "Input errors (2):\n  Blocks.block_1.Horizon: -0.1 is below minimum 0\n  x: missing"
    @test_throws PeriLab.PeriLabExceptions.PeriLabError PS.report!(ctx)
end

@testset "strict_mode" begin
    @test PS.strict_mode(Dict{String,Any}()) == true
    @test PS.strict_mode(Dict{String,Any}("Strict Validation" => false)) == false
    @test PS.strict_mode(Dict{String,Any}("Strict Validation" => true);
                         no_strict_flag = true) == false
    @test_throws PeriLab.PeriLabExceptions.PeriLabError PS.strict_mode(Dict{String,Any}("Strict Validation" => "no"))
end

@testset "ParamsDefinitionError" begin
    e = PS.ParamsDefinitionError("X.y: bad")
    @test sprint(showerror, e) == "ParamsDefinitionError: X.y: bad"
end
```

- [ ] **Step 3: Run test to verify it fails**

Run: `julia --project=. test/unit_tests/Support/Parameters/Spec/run_spec_tests.jl`
Expected: FAIL with `UndefVarError: ParameterSpec not defined in PeriLab`.

- [ ] **Step 4: Write the implementation**

`src/Support/Parameters/Spec/ParameterSpec.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

"""
    ParameterSpec

Typed, validated input parameters. Modules declare their parameters with
`@params`; YAML input is converted into those structs and validated against
the declarations. See
`docs/superpowers/specs/2026-10-03-typed-input-parameters-design.md`.
"""
module ParameterSpec

using ..PeriLabExceptions: @abort

include("errors.jl")

end
```

`src/Support/Parameters/Spec/errors.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

"""
    ParamsDefinitionError(msg)

Raised when an `@params` declaration or a registration is invalid. This is a
programming error in a module, not a problem in the user's input deck.
"""
struct ParamsDefinitionError <: Exception
    msg::String
end
Base.showerror(io::IO, e::ParamsDefinitionError) = print(io, "ParamsDefinitionError: ", e.msg)

"""
    InputError(path, message, severity)

One problem found in the input deck. `severity` is `:error` or `:warning`.
"""
struct InputError
    path::String
    message::String
    severity::Symbol
end

"""
    ParseContext(; directory = "", strict = true)

State shared while reading an input deck: the deck's directory (for relative
data file paths), strict mode, and all problems found so far.
"""
mutable struct ParseContext
    directory::String
    strict::Bool
    errors::Vector{InputError}
end
ParseContext(; directory::AbstractString = "", strict::Bool = true) = ParseContext(String(directory),
                                                                                    strict,
                                                                                    InputError[])

add_error!(ctx::ParseContext, path::AbstractString, msg::AbstractString) = push!(ctx.errors,
                                                                                  InputError(path,
                                                                                             msg,
                                                                                             :error))
add_warning!(ctx::ParseContext, path::AbstractString, msg::AbstractString) = push!(ctx.errors,
                                                                                    InputError(path,
                                                                                               msg,
                                                                                               :warning))
has_errors(ctx::ParseContext) = any(e -> e.severity == :error, ctx.errors)

"""
    join_path(path, key)

Appends `key` to an error path. Keys that are not plain identifiers are quoted.
"""
function join_path(path::AbstractString, key::AbstractString)
    segment = occursin(r"^[A-Za-z0-9_]+$", key) ? String(key) : "\"$key\""
    return isempty(path) ? segment : "$path.$segment"
end

_normalize(s::AbstractString) = lowercase(filter(c -> isletter(c) || isdigit(c), s))

function levenshtein(a::AbstractString, b::AbstractString)
    a_chars, b_chars = collect(a), collect(b)
    n = length(b_chars)
    previous = collect(0:n)
    current = similar(previous)
    for (i, ca) in enumerate(a_chars)
        current[1] = i
        for j in 1:n
            cost = ca == b_chars[j] ? 0 : 1
            current[j + 1] = min(previous[j + 1] + 1, current[j] + 1, previous[j] + cost)
        end
        previous, current = current, previous
    end
    return previous[n + 1]
end

"""
    suggest(key, candidates)

Returns the candidate closest to `key` (ignoring case, spaces and
punctuation), or `nothing` if none is close enough.
"""
function suggest(key::AbstractString, candidates)
    normalized_key = _normalize(key)
    best = nothing
    best_distance = typemax(Int)
    for candidate in candidates
        distance = levenshtein(normalized_key, _normalize(candidate))
        if distance < best_distance
            best, best_distance = String(candidate), distance
        end
    end
    best === nothing && return nothing
    return best_distance <= max(2, length(normalized_key) ÷ 4) ? best : nothing
end

function unknown_key_message(key::AbstractString, candidates)
    suggestion = suggest(key, candidates)
    return suggestion === nothing ? "unknown key" :
           "unknown key — did you mean \"$suggestion\"?"
end

function format_errors(errors::Vector{InputError})
    io = IOBuffer()
    print(io, "Input errors (", length(errors), "):")
    for e in errors
        print(io, "\n  ", e.path, ": ", e.message)
    end
    return String(take!(io))
end

"""
    report!(ctx)

Logs all warnings, then aborts with every error at once if there are any.
"""
function report!(ctx::ParseContext)
    for w in ctx.errors
        w.severity == :warning && @warn "$(w.path): $(w.message)"
    end
    errors = filter(e -> e.severity == :error, ctx.errors)
    isempty(errors) || @abort format_errors(errors)
    return nothing
end

"""
    strict_mode(input; no_strict_flag = false)

Strict validation is on unless the command line flag `--no-strict` is given or
the input deck sets `Strict Validation: false`.
"""
function strict_mode(input::AbstractDict; no_strict_flag::Bool = false)
    no_strict_flag && return false
    value = get(input, "Strict Validation", true)
    value isa Bool || @abort "\"Strict Validation\" must be true or false, got $(repr(value))"
    return value
end
```

In `src/PeriLab.jl`, directly after the line `include("./IO/exceptions.jl")`, add:

```julia
include("./Support/Parameters/Spec/ParameterSpec.jl")
```

- [ ] **Step 5: Run test to verify it passes**

Run: `julia --project=. test/unit_tests/Support/Parameters/Spec/run_spec_tests.jl`
Expected: PASS (all `ut_errors.jl` testsets green).

---

### Task 2: Dependent values (`Constant`, `Table1D`, data files, `combine`)

**Files:**
- Create: `src/Support/Parameters/Spec/dependent.jl`
- Modify: `src/Support/Parameters/Spec/ParameterSpec.jl`
- Create: `test/unit_tests/Support/Parameters/Spec/ut_dependent.jl`
- Modify: `test/unit_tests/Support/Parameters/Spec/spec_tests.jl`

**Interfaces:**
- Consumes (Task 1): `ParseContext`, `add_error!`, `has_errors`.
- Produces:
  - `abstract type Dependent end`
  - `struct Constant <: Dependent; value::Float64; end`
  - `struct Table1D <: Dependent` with fields `field_name::String`, `x::Vector{Float64}`, `y::Vector{Float64}`, `spline::Spline1D`, `field::Base.RefValue{Vector{Float64}}`, `bound::Base.RefValue{Bool}`, `warn::Base.RefValue{Bool}`, `source::String`; constructor `Table1D(field_name, x, y, source)`
  - `value(c::Constant, iID::Int64)::Float64`, `value(t::Table1D, iID::Int64)::Float64` (throws `ArgumentError` if unbound)
  - `bind_table!(t::Table1D, field::Vector{Float64})`
  - `read_table(file::String, alias::String, path::String, ctx::ParseContext)::Union{Table1D,Nothing}`
  - `combine(f, a::Dependent, b::Dependent)::Dependent`

- [ ] **Step 1: Write the failing test**

`test/unit_tests/Support/Parameters/Spec/ut_dependent.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const PS = PeriLab.ParameterSpec

const UT_TABLE_DIR = mktempdir()
function ut_write(name, content)
    path = joinpath(UT_TABLE_DIR, name)
    write(path, content)
    return path
end

const UT_E_T = """
# Young's modulus and Poisson's ratio over temperature
header: Temperature Young's_Modulus Poisson's_Ratio
0.0 200.0 0.30
100.0 180.0 0.30
200.0 150.0 0.31
"""

ut_table(alias) = PS.read_table(ut_write("E_T_$(hash(alias)).txt", UT_E_T), alias, "p",
                                PS.ParseContext())

@testset "Constant" begin
    c = PS.Constant(5.0)
    @test PS.value(c, 1) == 5.0
    @test PS.value(c, 99) == 5.0
    @test (@inferred PS.value(c, 3)) == 5.0
end

@testset "read_table" begin
    file = ut_write("E_T.txt", UT_E_T)
    ctx = PS.ParseContext()
    t = PS.read_table(file, "Young's Modulus", "M.\"Young's Modulus\"", ctx)
    @test !PS.has_errors(ctx)
    @test t isa PS.Table1D
    @test t.field_name == "Temperature"
    @test t.x == [0.0, 100.0, 200.0]
    @test t.y == [200.0, 180.0, 150.0]
    @test t.source == file
    @test !t.bound[]
    nu = PS.read_table(file, "Poisson's Ratio", "p", ctx)
    @test nu.y == [0.30, 0.30, 0.31]
end

@testset "read_table tolerates CRLF, tabs and blank lines" begin
    file = ut_write("crlf.txt",
                    "# comment\r\nheader:\tTemperature Young's_Modulus\r\n\r\n0.0\t200.0\r\n100.0  180.0\r\n")
    ctx = PS.ParseContext()
    t = PS.read_table(file, "Young's Modulus", "p", ctx)
    @test !PS.has_errors(ctx)
    @test t.x == [0.0, 100.0]
    @test t.y == [200.0, 180.0]
end

@testset "read_table errors" begin
    cases = [("missing.txt", nothing, "data file"),
             ("no_header.txt", "0.0 1.0\n1.0 2.0\n", "line 1: expected 'header:"),
             ("no_column.txt", "header: Temperature Density\n0 1\n1 2\n",
              "has no column \"Young's_Modulus\""),
             ("bad_value.txt", "header: Temperature Young's_Modulus\n0 1\n1 abc\n",
              "line 3: non-numeric value"),
             ("bad_count.txt", "header: Temperature Young's_Modulus\n0 1\n1\n",
              "line 3: expected 2 values, got 1"),
             ("unsorted.txt", "header: Temperature Young's_Modulus\n10 1\n0 2\n",
              "must be strictly increasing"),
             ("one_row.txt", "header: Temperature Young's_Modulus\n0 1\n",
              "at least 2 data rows")]
    for (name, content, expected) in cases
        file = content === nothing ? joinpath(UT_TABLE_DIR, name) : ut_write(name, content)
        ctx = PS.ParseContext()
        @test PS.read_table(file, "Young's Modulus", "p", ctx) === nothing
        @test length(ctx.errors) == 1
        @test occursin(expected, ctx.errors[1].message)
        @test ctx.errors[1].path == "p"
    end
end

@testset "Table1D value and binding" begin
    t = ut_table("Young's Modulus")
    @test_throws ArgumentError PS.value(t, 1)
    temperature = [0.0, 100.0, 200.0, -50.0, 500.0]
    PS.bind_table!(t, temperature)
    @test t.bound[]
    @test PS.value(t, 1) ≈ 200.0
    @test PS.value(t, 2) ≈ 180.0
    @test PS.value(t, 3) ≈ 150.0
    @test t.warn[]
    @test PS.value(t, 4) ≈ 200.0          # below range: nearest boundary value
    @test !t.warn[]                        # warned once, never again
    @test PS.value(t, 5) ≈ 150.0          # above range: nearest boundary value
    temperature[1] = 100.0                 # bound by reference: sees field updates
    @test PS.value(t, 1) ≈ 180.0
    @test (@inferred PS.value(t, 2)) ≈ 180.0
end

@testset "combine" begin
    @test PS.combine(+, PS.Constant(1.0), PS.Constant(2.0)) == PS.Constant(3.0)
    E = ut_table("Young's Modulus")
    doubled = PS.combine(*, E, PS.Constant(2.0))
    @test doubled isa PS.Table1D
    @test doubled.x == E.x
    @test doubled.y == [400.0, 360.0, 300.0]
    @test doubled.field_name == "Temperature"
    @test !doubled.bound[]
    @test PS.combine(-, PS.Constant(1000.0), E).y == [800.0, 820.0, 850.0]
    nu = ut_table("Poisson's Ratio")
    G = PS.combine((e, n) -> e / (2 * (1 + n)), E, nu)
    @test G.x == [0.0, 100.0, 200.0]
    @test G.y ≈ [200.0 / 2.6, 180.0 / 2.6, 150.0 / 2.62]
    other = PS.read_table(ut_write("E_D.txt", "header: Damage Young's_Modulus\n0 1\n1 2\n"),
                          "Young's Modulus", "p", PS.ParseContext())
    @test_throws ArgumentError PS.combine(+, E, other)
end
```

Add `"ut_dependent.jl"` to the list in `spec_tests.jl`:

```julia
for file in ["ut_errors.jl", "ut_dependent.jl"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. test/unit_tests/Support/Parameters/Spec/run_spec_tests.jl`
Expected: FAIL in `ut_dependent.jl` with `UndefVarError: Constant not defined in PeriLab.ParameterSpec`.

- [ ] **Step 3: Write the implementation**

In `ParameterSpec.jl`, replace

```julia
using ..PeriLabExceptions: @abort

include("errors.jl")
```

with

```julia
using ..PeriLabExceptions: @abort
using Dierckx: Spline1D, evaluate

include("errors.jl")
include("dependent.jl")
```

`src/Support/Parameters/Spec/dependent.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

export Dependent, Constant, Table1D, value, combine

"""
    Dependent

A parameter that is either a constant (`Constant`) or depends on one node
field via a data table (`Table1D`). Read it with `value(d, iID)`.
"""
abstract type Dependent end

struct Constant <: Dependent
    value::Float64
end

"""
    Table1D

Values interpolated (spline) over one node field, read from a data file.
Must be bound to the field array (`bind_table!`) before `value` is called.
"""
struct Table1D <: Dependent
    field_name::String
    x::Vector{Float64}
    y::Vector{Float64}
    spline::Spline1D
    field::Base.RefValue{Vector{Float64}}
    bound::Base.RefValue{Bool}
    warn::Base.RefValue{Bool}
    source::String
end

function Table1D(field_name::AbstractString, x::Vector{Float64}, y::Vector{Float64},
                 source::AbstractString)
    k = min(3, length(x) - 1)
    return Table1D(String(field_name), x, y, Spline1D(x, y; k = k, bc = "nearest"),
                   Ref(Float64[]), Ref(false), Ref(true), String(source))
end

@inline value(c::Constant, ::Int64) = c.value

function value(t::Table1D, iID::Int64)
    t.bound[] ||
        throw(ArgumentError("dependent value from $(t.source) is not bound to field \"$(t.field_name)\"; call bind_dependents! first"))
    x = t.field[][iID]
    if t.warn[] && (x < t.x[1] || x > t.x[end])
        @warn "$(t.field_name) = $x is outside the data range [$(t.x[1]), $(t.x[end])] of $(t.source). Using the nearest boundary value."
        t.warn[] = false
    end
    return evaluate(t.spline, x)
end

function bind_table!(t::Table1D, field::Vector{Float64})
    t.field[] = field
    t.bound[] = true
    return t
end

"""
    read_table(file, alias, path, ctx)

Reads the data file for parameter `alias`. Format: optional `#` comment lines,
a line `header: <field> <column> ...`, then whitespace separated numbers. The
first column is the node field the value depends on; the column used is
`alias` with spaces replaced by underscores. Problems are added to `ctx` at
`path` and `nothing` is returned.
"""
function read_table(file::String, alias::String, path::String, ctx::ParseContext)
    if !isfile(file)
        add_error!(ctx, path, "data file \"$file\" not found")
        return nothing
    end
    header = String[]
    rows = Vector{Vector{Float64}}()
    for (line_number, raw_line) in enumerate(eachline(file))
        line = strip(raw_line)
        (isempty(line) || startswith(line, "#")) && continue
        if isempty(header)
            if startswith(line, "header:")
                header = String.(split(line)[2:end])
                continue
            end
            add_error!(ctx, path,
                       "$file line $line_number: expected 'header: <field> <column> ...' before the data")
            return nothing
        end
        parts = split(line)
        if length(parts) != length(header)
            add_error!(ctx, path,
                       "$file line $line_number: expected $(length(header)) values, got $(length(parts))")
            return nothing
        end
        row = tryparse.(Float64, parts)
        if any(isnothing, row)
            add_error!(ctx, path, "$file line $line_number: non-numeric value")
            return nothing
        end
        push!(rows, Float64.(row))
    end
    column = replace(alias, " " => "_")
    index = findfirst(==(column), header)
    if index === nothing || index == 1
        add_error!(ctx, path,
                   "$file has no column \"$column\" (header: $(join(header, " ")))")
        return nothing
    end
    if length(rows) < 2
        add_error!(ctx, path, "$file needs at least 2 data rows")
        return nothing
    end
    x = [row[1] for row in rows]
    y = [row[index] for row in rows]
    if !all(diff(x) .> 0)
        add_error!(ctx, path,
                   "$file: first column ($(header[1])) must be strictly increasing")
        return nothing
    end
    return Table1D(header[1], x, y, file)
end

"""
    combine(f, a, b)

Applies `f` pointwise to two dependent values, e.g. to derive a shear modulus
from Young's modulus and Poisson's ratio in `derive`. Tables must depend on the
same field. The result is unbound.
"""
combine(f, a::Constant, b::Constant) = Constant(f(a.value, b.value))
combine(f, a::Table1D, b::Constant) = Table1D(a.field_name, a.x, f.(a.y, b.value), a.source)
combine(f, a::Constant, b::Table1D) = Table1D(b.field_name, b.x, f.(a.value, b.y), b.source)
function combine(f, a::Table1D, b::Table1D)
    a.field_name == b.field_name ||
        throw(ArgumentError("cannot combine values depending on \"$(a.field_name)\" ($(a.source)) and \"$(b.field_name)\" ($(b.source))"))
    x = sort!(unique!(vcat(a.x, b.x)))
    return Table1D(a.field_name, x, f.(evaluate(a.spline, x), evaluate(b.spline, x)),
                   "$(a.source) + $(b.source)")
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. test/unit_tests/Support/Parameters/Spec/run_spec_tests.jl`
Expected: PASS.

---

### Task 3: Field declarations, value conversion, constraints

**Files:**
- Create: `src/Support/Parameters/Spec/field_spec.jl`
- Create: `src/Support/Parameters/Spec/convert.jl`
- Modify: `src/Support/Parameters/Spec/ParameterSpec.jl`
- Create: `test/unit_tests/Support/Parameters/Spec/ut_convert.jl`
- Modify: `test/unit_tests/Support/Parameters/Spec/spec_tests.jl`

**Interfaces:**
- Consumes: Task 1 (`ParseContext`, `add_error!`, `has_errors`, `_normalize`, `join_path`), Task 2 (`Dependent`, `Constant`, `Table1D`, `read_table`).
- Produces:
  - `struct NoDefault end`, `const NO_DEFAULT`
  - `struct FieldDecl` (`alias`, `required`, `default`, `min`, `max`, `allowed`, `quantity`, `description`)
  - `req(alias; min, max, allowed, quantity, description)::FieldDecl`, `opt(alias; default, min, max, allowed, quantity, description)::FieldDecl`
  - `struct FieldSpec` (`name::Symbol`, `type::Any`, `alias`, `required`, `default`, `min`, `max`, `allowed`, `quantity`, `description`), constructor `FieldSpec(name::Symbol, type, decl::FieldDecl; default = decl.default)`
  - `struct Failed end`, `const FAILED`
  - `convert_value(T, raw, path::String, ctx::ParseContext; alias::String = "")` → value of type `T`, or `FAILED` (error already recorded)
  - `enum_aliases(::Type{E}) where {E<:Enum}` → `Dict{String,E}` (modules extend it)
  - `check_constraints!(fs::FieldSpec, v, path::String, ctx::ParseContext)::Bool`
  - helpers `_describe(raw)::String`, `_fresh(x)`, `_nonnothing(T::Union)`, `_fmt(x::Real)::String`

- [ ] **Step 1: Write the failing test**

`test/unit_tests/Support/Parameters/Spec/ut_convert.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const PS = PeriLab.ParameterSpec

@enum UTHardening LinearHardening ExponentialHardening
@enum UTSymmetry PlaneStress PlaneStrain Full3D
PS.enum_aliases(::Type{UTSymmetry}) = Dict{String,UTSymmetry}("3D" => Full3D)

function ut_conv(T, raw; directory = "", alias = "")
    ctx = PS.ParseContext(directory = directory)
    v = PS.convert_value(T, raw, "k", ctx; alias = alias)
    return v, ctx
end

function ut_conv_error(T, raw)
    v, ctx = ut_conv(T, raw)
    @test v === PS.FAILED
    @test length(ctx.errors) == 1
    @test ctx.errors[1].path == "k"
    return ctx.errors[1].message
end

@testset "Float64" begin
    @test ut_conv(Float64, 2)[1] === 2.0
    @test ut_conv(Float64, 2.5)[1] === 2.5
    @test ut_conv_error(Float64, "abc") == "expected a number, got \"abc\""
    @test ut_conv_error(Float64, "E.txt") ==
          "expected a number, got file path \"E.txt\" — this parameter does not support dependent values"
    @test ut_conv_error(Float64, true) == "expected a number, got true"
    @test ut_conv_error(Float64, nothing) == "expected a number, got an empty value"
end

@testset "Int64" begin
    @test ut_conv(Int64, 3)[1] === 3
    @test ut_conv(Int64, 100.0)[1] === 100
    @test ut_conv_error(Int64, 2.5) == "expected an integer, got 2.5"
    @test ut_conv_error(Int64, "3") == "expected an integer, got \"3\""
end

@testset "Bool and String" begin
    @test ut_conv(Bool, true)[1] === true
    @test ut_conv_error(Bool, "yes") == "expected true or false, got \"yes\""
    @test ut_conv(String, "abc")[1] == "abc"
    @test ut_conv_error(String, 5) == "expected text, got 5"
end

@testset "Enum" begin
    @test ut_conv(UTHardening, "LinearHardening")[1] === LinearHardening
    @test ut_conv(UTHardening, "Linear Hardening")[1] === LinearHardening
    @test ut_conv(UTHardening, "exponential_hardening")[1] === ExponentialHardening
    @test ut_conv(UTHardening, ExponentialHardening)[1] === ExponentialHardening
    @test ut_conv(UTSymmetry, "3D")[1] === Full3D
    @test ut_conv(UTSymmetry, "plane stress")[1] === PlaneStress
    @test ut_conv_error(UTHardening, "Cubic") ==
          "\"Cubic\" is not one of: LinearHardening, ExponentialHardening"
    @test ut_conv_error(UTSymmetry, 3) ==
          "3 is not one of: PlaneStress, PlaneStrain, Full3D, 3D"
end

@testset "Vectors" begin
    v = ut_conv(Vector{Float64}, [1, 2.5])[1]
    @test v == [1.0, 2.5] && v isa Vector{Float64}
    @test ut_conv(Vector{Int64}, [1, 2])[1] == [1, 2]
    @test ut_conv(Vector{String}, ["a", "b"])[1] == ["a", "b"]
    v, ctx = ut_conv(Vector{Float64}, [1, "a"])
    @test v === PS.FAILED
    @test ctx.errors[1].path == "k[2]"
    @test ut_conv_error(Vector{Float64}, 3) == "expected a list, got 3"
end

@testset "Union{Nothing,T}" begin
    @test ut_conv(Union{Nothing,Float64}, nothing)[1] === nothing
    @test ut_conv(Union{Nothing,Float64}, 2)[1] === 2.0
    @test ut_conv_error(Union{Nothing,Float64}, "x") == "expected a number, got \"x\""
end

@testset "Dependent" begin
    @test ut_conv(PS.Dependent, 3)[1] === PS.Constant(3.0)
    dir = mktempdir()
    write(joinpath(dir, "E.txt"), "header: Temperature Young's_Modulus\n0 200\n100 180\n")
    t, ctx = ut_conv(PS.Dependent, "E.txt"; directory = dir, alias = "Young's Modulus")
    @test !PS.has_errors(ctx)
    @test t isa PS.Table1D
    @test t.y == [200.0, 180.0]
    v, ctx = ut_conv(PS.Dependent, "missing.txt"; directory = dir, alias = "Young's Modulus")
    @test v === PS.FAILED
    @test occursin("not found", ctx.errors[1].message)
    @test ut_conv_error(PS.Dependent, true) ==
          "expected a number or a data file path, got true"
end

@testset "req / opt / FieldSpec" begin
    d = PS.req("Horizon"; min = 0, quantity = :length, description = "radius")
    @test d.alias == "Horizon" && d.required && d.min === 0.0 && d.max === nothing
    @test d.quantity === :length && d.description == "radius"
    o = PS.opt("Steps"; default = 10, allowed = [10, 20])
    @test !o.required && o.default == 10 && o.allowed == Any[10, 20]
    @test PS.opt("X").default === PS.NO_DEFAULT
    fs = PS.FieldSpec(:steps, Int64, o)
    @test fs.name === :steps && fs.type === Int64 && fs.alias == "Steps" && fs.default == 10
end

@testset "check_constraints!" begin
    ctx = PS.ParseContext()
    fs = PS.FieldSpec(:h, Float64, PS.req("Horizon"; min = 0, max = 10))
    @test PS.check_constraints!(fs, 5.0, "h", ctx)
    @test !PS.check_constraints!(fs, -0.1, "h", ctx)
    @test ctx.errors[end].message == "-0.1 is below minimum 0"
    @test !PS.check_constraints!(fs, 11.0, "h", ctx)
    @test ctx.errors[end].message == "11.0 is above maximum 10"
    vec_fs = PS.FieldSpec(:v, Vector{Float64}, PS.req("V"; min = 0))
    @test !PS.check_constraints!(vec_fs, [1.0, -2.0], "v", ctx)
    @test ctx.errors[end].message == "-2.0 is below minimum 0"
    dir = mktempdir()
    write(joinpath(dir, "E.txt"), "header: Temperature Young's_Modulus\n0 200\n100 180\n")
    table = PS.read_table(joinpath(dir, "E.txt"), "Young's Modulus", "e", ctx)
    dep_fs = PS.FieldSpec(:e, PS.Dependent, PS.req("Young's Modulus"; max = 190))
    @test !PS.check_constraints!(dep_fs, table, "e", ctx)
    @test ctx.errors[end].message == "200.0 is above maximum 190"
    @test PS.check_constraints!(dep_fs, PS.Constant(100.0), "e", ctx)
    type_fs = PS.FieldSpec(:t, String, PS.req("Type"; allowed = ["Exodus", "CSV"]))
    @test PS.check_constraints!(type_fs, "CSV", "t", ctx)
    @test !PS.check_constraints!(type_fs, "VTK", "t", ctx)
    @test ctx.errors[end].message == "\"VTK\" is not one of: \"Exodus\", \"CSV\""
    union_fs = PS.FieldSpec(:u, Union{Nothing,Float64}, PS.opt("U"; default = nothing, min = 0))
    @test PS.check_constraints!(union_fs, nothing, "u", ctx)
end
```

Add `"ut_convert.jl"` to the list in `spec_tests.jl`:

```julia
for file in ["ut_errors.jl", "ut_dependent.jl", "ut_convert.jl"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. test/unit_tests/Support/Parameters/Spec/run_spec_tests.jl`
Expected: FAIL in `ut_convert.jl` with `UndefVarError: enum_aliases not defined in PeriLab.ParameterSpec`.

- [ ] **Step 3: Write the implementation**

In `ParameterSpec.jl`, replace

```julia
include("errors.jl")
include("dependent.jl")
```

with

```julia
include("errors.jl")
include("dependent.jl")
include("field_spec.jl")
include("convert.jl")
```

`src/Support/Parameters/Spec/field_spec.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

export req, opt

struct NoDefault end
const NO_DEFAULT = NoDefault()

"""
    FieldDecl

What a module author wrote for one field via `req(...)` or `opt(...)`.
"""
struct FieldDecl
    alias::String
    required::Bool
    default::Any
    min::Union{Nothing,Float64}
    max::Union{Nothing,Float64}
    allowed::Union{Nothing,Vector{Any}}
    quantity::Union{Nothing,Symbol}
    description::String
end

_bound(::Nothing) = nothing
_bound(x::Real) = Float64(x)
_allowed(::Nothing) = nothing
_allowed(values) = Vector{Any}(collect(values))

"""
    req(yaml_key; min, max, allowed, quantity, description)

A required parameter. `quantity` (e.g. `:stress`) is documentation only;
PeriLab has no fixed unit system.
"""
function req(alias::AbstractString; min = nothing, max = nothing, allowed = nothing,
             quantity = nothing, description::AbstractString = "")
    return FieldDecl(String(alias), true, NO_DEFAULT, _bound(min), _bound(max),
                     _allowed(allowed), quantity, String(description))
end

"""
    opt(yaml_key; default, min, max, allowed, quantity, description)

An optional parameter; `default` is required.
"""
function opt(alias::AbstractString; default = NO_DEFAULT, min = nothing, max = nothing,
             allowed = nothing, quantity = nothing, description::AbstractString = "")
    return FieldDecl(String(alias), false, default, _bound(min), _bound(max),
                     _allowed(allowed), quantity, String(description))
end

"""
    FieldSpec

A validated field declaration: struct field name, declared type, and the
`FieldDecl` metadata. `default` is already converted to `type`.
"""
struct FieldSpec
    name::Symbol
    type::Any
    alias::String
    required::Bool
    default::Any
    min::Union{Nothing,Float64}
    max::Union{Nothing,Float64}
    allowed::Union{Nothing,Vector{Any}}
    quantity::Union{Nothing,Symbol}
    description::String
end

function FieldSpec(name::Symbol, type, decl::FieldDecl; default = decl.default)
    return FieldSpec(name, type, decl.alias, decl.required, default, decl.min, decl.max,
                     decl.allowed, decl.quantity, decl.description)
end
```

`src/Support/Parameters/Spec/convert.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

"Returned by `convert_value` when conversion failed; the error is already in the context."
struct Failed end
const FAILED = Failed()

_fail(ctx::ParseContext, path::AbstractString, msg::AbstractString) = (add_error!(ctx, path,
                                                                                  msg);
                                                                       FAILED)

_describe(::Nothing) = "an empty value"
_describe(raw::AbstractString) = "\"$raw\""
_describe(raw) = repr(raw)

_fresh(x) = (x isa AbstractArray || x isa AbstractDict) ? copy(x) : x

_nonnothing(T::Union) = T.a === Nothing ? T.b : T.a

_fmt(x::Real) = isinteger(x) ? string(Int(x)) : string(x)

"""
    enum_aliases(::Type{E})

Extra YAML spellings for an `@enum`, for values that are not valid Julia
names. Extend it in the module defining the enum, e.g.
`ParameterSpec.enum_aliases(::Type{Symmetry}) = Dict("3D" => Full3D)`.
Without an alias, a YAML string matches an enum instance if both are equal
ignoring case, spaces and punctuation ("plane stress" matches `PlaneStress`).
"""
enum_aliases(::Type{E}) where {E<:Enum} = Dict{String,E}()

"""
    convert_value(T, raw, path, ctx; alias = "")

Converts a raw YAML value to the declared type `T`. On failure the error is
added to `ctx` and `FAILED` is returned. `alias` is the YAML key, needed to
find the column in a `Dependent` data file.
"""
function convert_value(T, raw, path::String, ctx::ParseContext; alias::String = "")
    raw isa T && return _fresh(raw)
    if T isa Union
        return convert_value(_nonnothing(T), raw, path, ctx; alias = alias)
    elseif T === Float64
        raw isa Real && !(raw isa Bool) && return Float64(raw)
        if raw isa AbstractString && endswith(lowercase(raw), ".txt")
            return _fail(ctx, path,
                         "expected a number, got file path \"$raw\" — this parameter does not support dependent values")
        end
        return _fail(ctx, path, "expected a number, got $(_describe(raw))")
    elseif T === Int64
        raw isa Integer && !(raw isa Bool) && return Int64(raw)
        raw isa AbstractFloat && isinteger(raw) && return Int64(raw)
        return _fail(ctx, path, "expected an integer, got $(_describe(raw))")
    elseif T === Bool
        return _fail(ctx, path, "expected true or false, got $(_describe(raw))")
    elseif T === String
        raw isa AbstractString && return String(raw)
        return _fail(ctx, path, "expected text, got $(_describe(raw))")
    elseif T === Dependent
        raw isa Real && !(raw isa Bool) && return Constant(Float64(raw))
        if raw isa AbstractString
            table = read_table(joinpath(ctx.directory, raw), alias, path, ctx)
            return table === nothing ? FAILED : table
        end
        return _fail(ctx, path, "expected a number or a data file path, got $(_describe(raw))")
    elseif T isa DataType && T <: Enum
        return _convert_enum(T, raw, path, ctx)
    elseif T isa DataType && T <: Vector
        return _convert_vector(T, raw, path, ctx)
    end
    throw(ArgumentError("convert_value: unsupported type $T"))
end

function _convert_enum(::Type{E}, raw, path::String, ctx::ParseContext) where {E<:Enum}
    aliases = enum_aliases(E)
    if raw isa AbstractString
        haskey(aliases, raw) && return aliases[raw]
        key = _normalize(raw)
        for instance in instances(E)
            _normalize(string(instance)) == key && return instance
        end
    end
    names = vcat([string(instance) for instance in instances(E)], collect(keys(aliases)))
    return _fail(ctx, path, "$(_describe(raw)) is not one of: $(join(names, ", "))")
end

function _convert_vector(::Type{Vector{S}}, raw, path::String, ctx::ParseContext) where {S}
    raw isa AbstractVector || return _fail(ctx, path, "expected a list, got $(_describe(raw))")
    out = Vector{S}(undef, length(raw))
    ok = true
    for (i, item) in enumerate(raw)
        v = convert_value(S, item, "$path[$i]", ctx)
        if v === FAILED
            ok = false
        else
            out[i] = v
        end
    end
    return ok ? out : FAILED
end

_numbers(v::Bool) = ()
_numbers(v::Real) = (v,)
_numbers(v::AbstractVector{<:Real}) = v
_numbers(v::Constant) = (v.value,)
_numbers(v::Table1D) = v.y
_numbers(v) = ()

"""
    check_constraints!(fs, v, path, ctx)

Checks `min` / `max` (on every number in `v`, including all table values) and
`allowed`. Adds the first violation to `ctx` and returns `false`.
"""
function check_constraints!(fs::FieldSpec, v, path::String, ctx::ParseContext)
    for x in _numbers(v)
        if fs.min !== nothing && x < fs.min
            add_error!(ctx, path, "$x is below minimum $(_fmt(fs.min))")
            return false
        end
        if fs.max !== nothing && x > fs.max
            add_error!(ctx, path, "$x is above maximum $(_fmt(fs.max))")
            return false
        end
    end
    if fs.allowed !== nothing && v !== nothing && !(v in fs.allowed)
        add_error!(ctx, path,
                   "$(_describe(v)) is not one of: $(join(_describe.(fs.allowed), ", "))")
        return false
    end
    return true
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. test/unit_tests/Support/Parameters/Spec/run_spec_tests.jl`
Expected: PASS.

---

### Task 4: `@params` macro, building structs from YAML dicts, unknown keys, `derive`

**Files:**
- Create: `src/Support/Parameters/Spec/params_macro.jl`
- Create: `src/Support/Parameters/Spec/build.jl`
- Modify: `src/Support/Parameters/Spec/convert.jl` (add nested-section branches)
- Modify: `src/Support/Parameters/Spec/ParameterSpec.jl`
- Create: `test/unit_tests/Support/Parameters/Spec/ut_params.jl`
- Modify: `test/unit_tests/Support/Parameters/Spec/spec_tests.jl`

**Interfaces:**
- Consumes: Tasks 1–3 (`ParamsDefinitionError`, `ParseContext`, `join_path`, `unknown_key_message`, `add_error!`, `add_warning!`, `has_errors`, `Dependent`, `FieldDecl`, `FieldSpec`, `NO_DEFAULT`, `convert_value`, `check_constraints!`, `FAILED`, `_describe`, `_fresh`, `_nonnothing`).
- Produces:
  - `@params struct Name ... end` — emits the struct (parametric `Name{T_field<:Dependent}` for each `Dependent` field), `parameter_spec(::Type{<:Name})::Vector{FieldSpec}`, `is_params(::Type{<:Name}) = true`
  - `is_params(::Type)::Bool` (fallback `false`), `parameter_spec(T)`, `derive(p)` (fallback returns `p`; modules extend it)
  - `build_spec(T, display::String, entries::Vector{Any})::Vector{FieldSpec}`, `supported_type(T)::Bool`, `_typename(T)::String`
  - `build(T, dict::AbstractDict, path::String, ctx::ParseContext; owner::String)` → instance or `nothing` (does **not** check unknown keys)
  - `check_unknown!(dict::AbstractDict, known, path::String, ctx::ParseContext)` (skips `"Globals"`)
  - `aliases(T)::Set{String}`
  - `parse_section(T, dict::AbstractDict, path::String, ctx::ParseContext)` → `derive(instance)` or `nothing`

- [ ] **Step 1: Write the failing test**

`test/unit_tests/Support/Parameters/Spec/ut_params.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const PS = PeriLab.ParameterSpec

@enum UTMode FastMode SafeMode

PS.@params struct UTInner
    radius::Float64 = req("Radius"; min = 0)
    allow_contact::Bool = opt("Allow Contact"; default = false)
end

"""Docstring for UTOuter."""
PS.@params struct UTOuter
    horizon::Float64 = req("Horizon"; min = 0, quantity = :length,
                           description = "Neighborhood radius")
    steps::Int64 = opt("Number of Steps"; default = 10, min = 1)
    mode::UTMode = opt("Mode"; default = SafeMode)
    note::Union{Nothing,String} = opt("Note"; default = nothing)
    weights::Vector{Float64} = opt("Weights"; default = [1.0, 2.0])
    filter::UTInner = req("Filter")
    sets::Dict{String,UTInner} = opt("Sets"; default = Dict{String,UTInner}())
end

PS.@params struct UTDependentMat
    youngs_modulus::Dependent = req("Young's Modulus"; min = 0)
    poissons_ratio::Float64 = req("Poisson's Ratio"; min = -1, max = 0.5)
end

PS.@params struct UTDerived
    a::Float64 = req("A")
    twice_a::Float64 = opt("Twice A"; default = 0.0)
end
PS.derive(p::UTDerived) = UTDerived(p.a, 2 * p.a)

ut_errors_by_path(ctx) = Dict(e.path => e.message for e in ctx.errors)

@testset "@params generates a plain struct and its spec" begin
    @test PS.is_params(UTOuter)
    @test !PS.is_params(Float64)
    @test fieldnames(UTOuter) ==
          (:horizon, :steps, :mode, :note, :weights, :filter, :sets)
    @test isconcretetype(UTOuter)
    spec = PS.parameter_spec(UTOuter)
    @test [fs.alias for fs in spec] ==
          ["Horizon", "Number of Steps", "Mode", "Note", "Weights", "Filter", "Sets"]
    @test spec[1].required && spec[1].min === 0.0 && spec[1].quantity === :length
    @test spec[1].description == "Neighborhood radius"
    @test spec[2].default === 10
    @test spec[3].default === SafeMode
    @test occursin("Docstring for UTOuter", string(@doc UTOuter))
end

@testset "Dependent fields make the struct parametric" begin
    @test !isconcretetype(UTDependentMat)
    @test fieldtype(UTDependentMat{PS.Constant}, :youngs_modulus) === PS.Constant
    @test PS.parameter_spec(UTDependentMat)[1].type === PS.Dependent
end

@testset "build: valid input" begin
    dict = Dict{String,Any}("Horizon" => 1,
                            "Filter" => Dict{String,Any}("Radius" => 0.5),
                            "Mode" => "Fast Mode",
                            "Sets" => Dict{String,Any}("left" => Dict{String,Any}("Radius" => 1.0,
                                                                                  "Allow Contact" => true)))
    ctx = PS.ParseContext()
    p = PS.parse_section(UTOuter, dict, "Disc", ctx)
    @test isempty(ctx.errors)
    @test p isa UTOuter
    @test p.horizon === 1.0
    @test p.steps === 10
    @test p.mode === FastMode
    @test p.note === nothing
    @test p.weights == [1.0, 2.0]
    @test p.filter == UTInner(0.5, false)
    @test p.sets["left"] == UTInner(1.0, true)
    p2 = PS.parse_section(UTOuter, dict, "Disc", ctx)
    @test p2.weights !== p.weights        # defaults are not shared between instances
end

@testset "build: collects all errors" begin
    dict = Dict{String,Any}("Horizon" => -0.1,
                            "Number of Steps" => 0,
                            "Filter" => Dict{String,Any}("Radius" => "big", "Radus" => 1),
                            "Poissons Ratio" => 0.3,
                            "Globals" => Dict{String,Any}("anything" => 1))
    ctx = PS.ParseContext()
    @test PS.parse_section(UTOuter, dict, "Disc", ctx) === nothing
    msgs = ut_errors_by_path(ctx)
    @test msgs["Disc.Horizon"] == "-0.1 is below minimum 0"
    @test msgs["Disc.\"Number of Steps\""] == "0 is below minimum 1"
    @test msgs["Disc.Filter.Radius"] == "expected a number, got \"big\""
    @test msgs["Disc.Filter.Radus"] == "unknown key — did you mean \"Radius\"?"
    @test msgs["Disc.\"Poissons Ratio\""] == "unknown key"
    @test length(ctx.errors) == 5         # "Globals" is never reported
end

@testset "missing and empty values" begin
    ctx = PS.ParseContext()
    @test PS.parse_section(UTOuter, Dict{String,Any}("Horizon" => nothing), "Disc", ctx) ===
          nothing
    msgs = ut_errors_by_path(ctx)
    @test msgs["Disc.Horizon"] == "expected a number, got an empty value"
    @test msgs["Disc.Filter"] == "missing (required by UTOuter)"
end

@testset "Int field accepts integral float" begin
    ctx = PS.ParseContext()
    p = PS.parse_section(UTOuter,
                         Dict{String,Any}("Horizon" => 1.0, "Number of Steps" => 100.0,
                                          "Filter" => Dict{String,Any}("Radius" => 1)),
                         "Disc", ctx)
    @test isempty(ctx.errors)
    @test p.steps === 100
    @test p.filter.radius === 1.0
end

@testset "case variants are not accepted silently" begin
    ctx = PS.ParseContext()
    @test PS.parse_section(UTInner, Dict{String,Any}("radius" => 1.0), "F", ctx) === nothing
    msgs = ut_errors_by_path(ctx)
    @test msgs["F.radius"] == "unknown key — did you mean \"Radius\"?"
    @test msgs["F.Radius"] == "missing (required by UTInner)"
end

@testset "non-strict mode downgrades unknown keys to warnings" begin
    ctx = PS.ParseContext(strict = false)
    p = PS.parse_section(UTInner, Dict{String,Any}("Radius" => 1.0, "Colour" => "red"), "F",
                         ctx)
    @test p == UTInner(1.0, false)
    @test !PS.has_errors(ctx)
    @test ctx.errors[1].severity == :warning
    @test ctx.errors[1].path == "F.Colour"
end

@testset "derive runs after building" begin
    ctx = PS.ParseContext()
    p = PS.parse_section(UTDerived, Dict{String,Any}("A" => 2.0), "D", ctx)
    @test p.twice_a == 4.0
end

@testset "Dependent fields via parse_section" begin
    dir = mktempdir()
    write(joinpath(dir, "E.txt"), "header: Temperature Young's_Modulus\n0 200\n100 180\n")
    ctx = PS.ParseContext(directory = dir)
    pc = PS.parse_section(UTDependentMat,
                          Dict{String,Any}("Young's Modulus" => 210.0,
                                           "Poisson's Ratio" => 0.3), "M", ctx)
    @test pc isa UTDependentMat{PS.Constant}
    @test isconcretetype(typeof(pc))
    pt = PS.parse_section(UTDependentMat,
                          Dict{String,Any}("Young's Modulus" => "E.txt",
                                           "Poisson's Ratio" => 0.3), "M", ctx)
    @test pt isa UTDependentMat{PS.Table1D}
    @test isempty(ctx.errors)
    @test PS.parse_section(UTInner, Dict{String,Any}("Radius" => "E.txt"), "F", ctx) ===
          nothing
    @test occursin("does not support dependent values", ctx.errors[end].message)
end

function ut_definition_error(ex)
    try
        Core.eval(@__MODULE__, ex)
    catch e
        return e
    end
    return nothing
end

@testset "definition-time errors" begin
    cases = [(:(PS.@params struct UTBad1
                    x::Float64
                end),
              "UTBad1: field `x` must be written as `x::Type = req(\"YAML key\"; ...)` or `x::Type = opt(\"YAML key\"; default = ...)`"),
             (:(PS.@params struct UTBad2
                    x::Float64 = 3.0
                end),
              "UTBad2: field `x` must be written as"),
             (:(PS.@params struct UTBad3
                    x::Matrix{Float64} = req("X")
                end),
              "UTBad3.x: unsupported field type Matrix{Float64}"),
             (:(PS.@params struct UTBad4
                    x::Float64 = opt("X")
                end),
              "UTBad4.x: opt(\"X\") needs a default"),
             (:(PS.@params struct UTBad5
                    x::Float64 = opt("X"; default = -1.0, min = 0)
                end),
              "UTBad5.x: default -1.0 is invalid: -1.0 is below minimum 0"),
             (:(PS.@params struct UTBad6
                    m::UTMode = opt("Mode"; default = "Turbo")
                end),
              "UTBad6.m: default \"Turbo\" is invalid: \"Turbo\" is not one of: FastMode, SafeMode"),
             (:(PS.@params struct UTBad7
                    a::Float64 = req("X")
                    b::Float64 = req("X")
                end),
              "UTBad7.b: YAML key \"X\" is already used by field `a`"),
             (:(PS.@params struct UTBad8
                    s::String = req("S"; min = 0)
                end),
              "UTBad8.s: min/max are only allowed on numeric fields"),
             (:(PS.@params mutable struct UTBad9
                    x::Float64 = req("X")
                end),
              "@params structs must be immutable"),
             (:(PS.@params struct UTBad10
                    inner::UTDependentMat = req("Inner")
                end),
              "UTBad10.inner: nested section type UTDependentMat contains Dependent fields; this is not supported"),
             (:(PS.@params x = 1),
              "@params must be applied to a struct definition")]
    for (ex, expected) in cases
        e = ut_definition_error(ex)
        @test e isa PS.ParamsDefinitionError
        @test e !== nothing && occursin(expected, e.msg)
    end
end
```

Add `"ut_params.jl"` to the list in `spec_tests.jl`:

```julia
for file in ["ut_errors.jl", "ut_dependent.jl", "ut_convert.jl", "ut_params.jl"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. test/unit_tests/Support/Parameters/Spec/run_spec_tests.jl`
Expected: FAIL in `ut_params.jl` with `UndefVarError: @params not defined in PeriLab.ParameterSpec`.

- [ ] **Step 3: Write the implementation**

In `ParameterSpec.jl`, replace

```julia
include("field_spec.jl")
include("convert.jl")
```

with

```julia
include("field_spec.jl")
include("convert.jl")
include("params_macro.jl")
include("build.jl")
```

`src/Support/Parameters/Spec/params_macro.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

export @params, derive

const _SELF = @__MODULE__

"`true` for structs declared with `@params`."
is_params(::Type) = false

"""
    parameter_spec(T) -> Vector{FieldSpec}

Field declarations of an `@params` struct, in field order.
"""
function parameter_spec end

"""
    derive(p) -> p′

Hook run once after a parameter struct is built from the input deck. Return a
struct of the same struct type (type parameters may differ) with derived or
normalized values filled in. The default returns `p` unchanged.
"""
derive(p) = p

_typename(T) = replace(string(T), r"(\w+\.)+" => "")

const _SCALAR_TYPES = (Float64, Int64, Bool, String)
const _VECTOR_TYPES = (Vector{Float64}, Vector{Int64}, Vector{String})

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

function _numeric_type(T)
    T isa Union && return _numeric_type(_nonnothing(T))
    return T in (Float64, Int64, Dependent, Vector{Float64}, Vector{Int64})
end

"""
    build_spec(T, display, entries)

Turns the `(field name, declared type, FieldDecl)` entries emitted by `@params`
into `FieldSpec`s, raising `ParamsDefinitionError` for invalid declarations.
"""
function build_spec(T, display::String, entries::Vector{Any})
    specs = FieldSpec[]
    used = Dict{String,Symbol}()
    for (fname, ftype, decl) in entries
        location = "$display.$fname"
        if ftype isa UnionAll && is_params(ftype)
            throw(ParamsDefinitionError("$location: nested section type $(_typename(ftype)) contains Dependent fields; this is not supported"))
        end
        supported_type(ftype) ||
            throw(ParamsDefinitionError("$location: unsupported field type $(_typename(ftype)). Supported: Float64, Int64, Bool, String, an @enum, Vector{Float64}, Vector{Int64}, Vector{String}, Dependent, Union{Nothing,T}, a nested @params struct, or Dict{String,<nested @params struct>}"))
        if haskey(used, decl.alias)
            throw(ParamsDefinitionError("$location: YAML key \"$(decl.alias)\" is already used by field `$(used[decl.alias])`"))
        end
        used[decl.alias] = fname
        if (decl.min !== nothing || decl.max !== nothing) && !_numeric_type(ftype)
            throw(ParamsDefinitionError("$location: min/max are only allowed on numeric fields"))
        end
        default = NO_DEFAULT
        if !decl.required
            if decl.default === NO_DEFAULT
                throw(ParamsDefinitionError("$location: opt(\"$(decl.alias)\") needs a default, e.g. opt(\"$(decl.alias)\"; default = ...). Use a Union{Nothing,T} field with default = nothing if the value may be absent"))
            end
            ctx = ParseContext()
            v = convert_value(ftype, decl.default, "", ctx; alias = decl.alias)
            v === FAILED || check_constraints!(FieldSpec(fname, ftype, decl), v, "", ctx)
            if has_errors(ctx)
                throw(ParamsDefinitionError("$location: default $(_describe(decl.default)) is invalid: $(ctx.errors[1].message)"))
            end
            default = v
        end
        push!(specs, FieldSpec(fname, ftype, decl; default = default))
    end
    return specs
end

_is_dependent_type(t) = t === :Dependent ||
                        (t isa Expr && t.head === :. && t.args[end] == QuoteNode(:Dependent))

_definition_error(msg) = :(throw($ParamsDefinitionError($msg)))

function _field_usage(name, fname)
    return "$name: field `$fname` must be written as `$fname::Type = req(\"YAML key\"; ...)` or `$fname::Type = opt(\"YAML key\"; default = ...)`"
end

"""
    @params struct Name
        field::Type = req("YAML key"; min, max, allowed, quantity, description)
        field::Type = opt("YAML key"; default, min, max, allowed, quantity, description)
    end

Declares a parameter struct. Generates the plain immutable struct (type
parameters are added for `Dependent` fields so that every instance is a
concrete type) and its `parameter_spec`. Invalid declarations raise a
`ParamsDefinitionError` naming the struct and field.
"""
macro params(structdef)
    if !(structdef isa Expr && structdef.head === :struct)
        return _definition_error("@params must be applied to a struct definition")
    end
    if structdef.args[1]
        return _definition_error("@params structs must be immutable: use `struct`, not `mutable struct`")
    end
    name = structdef.args[2]
    if !(name isa Symbol)
        return _definition_error("@params struct $(name): write a plain name without type parameters or supertype")
    end
    fields = Any[]
    typeparams = Any[]
    entries = Any[]
    for line in structdef.args[3].args
        if line isa LineNumberNode
            push!(fields, line)
            continue
        end
        line isa AbstractString && continue
        if !(line isa Expr && line.head === :(=) && line.args[1] isa Expr &&
             line.args[1].head === :(::) && length(line.args[1].args) == 2)
            fname = line isa Expr && line.head === :(::) ? line.args[1] : line
            return _definition_error(_field_usage(name, fname))
        end
        fname, ftype = line.args[1].args
        decl = line.args[2]
        if !(decl isa Expr && decl.head === :call && decl.args[1] in (:req, :opt))
            return _definition_error(_field_usage(name, fname))
        end
        call = Expr(:call, GlobalRef(_SELF, decl.args[1]), decl.args[2:end]...)
        if _is_dependent_type(ftype)
            typeparam = Symbol("T_", fname)
            push!(typeparams, Expr(:<:, typeparam, Dependent))
            push!(fields, Expr(:(::), fname, typeparam))
            push!(entries, Expr(:tuple, QuoteNode(fname), Dependent, call))
        else
            push!(fields, Expr(:(::), fname, ftype))
            push!(entries, Expr(:tuple, QuoteNode(fname), ftype, call))
        end
    end
    head = isempty(typeparams) ? name : Expr(:curly, name, typeparams...)
    specname = Symbol("__params_spec_", name)
    structexpr = Expr(:struct, false, head, Expr(:block, fields...))
    return esc(quote
                   Base.@__doc__ $structexpr
                   const $specname = $(_SELF).build_spec($name, $(string(name)),
                                                         Any[$(entries...)])
                   $(_SELF).parameter_spec(::Type{<:$name}) = $specname
                   $(_SELF).is_params(::Type{<:$name}) = true
                   nothing
               end)
end
```

`src/Support/Parameters/Spec/build.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

"YAML keys declared by an `@params` struct."
aliases(T) = Set{String}(fs.alias for fs in parameter_spec(T))

"""
    build(T, dict, path, ctx; owner)

Builds an instance of the `@params` struct `T` from a YAML dict. Every field
is checked, so all problems are collected in `ctx`; returns `nothing` if any
field failed. Does not check for unknown keys (see `check_unknown!`), so that
composite models can share one dict. `owner` names the struct or model in
"missing" messages.
"""
function build(T, dict::AbstractDict, path::String, ctx::ParseContext;
               owner::String = string(nameof(T)))
    values = Any[]
    failed = false
    for fs in parameter_spec(T)
        field_path = join_path(path, fs.alias)
        if !haskey(dict, fs.alias)
            if fs.required
                add_error!(ctx, field_path, "missing (required by $owner)")
                failed = true
            else
                push!(values, _fresh(fs.default))
            end
            continue
        end
        v = convert_value(fs.type, dict[fs.alias], field_path, ctx; alias = fs.alias)
        if v === FAILED || !check_constraints!(fs, v, field_path, ctx)
            failed = true
        else
            push!(values, v)
        end
    end
    return failed ? nothing : T(values...)
end

"""
    check_unknown!(dict, known, path, ctx)

Reports keys of `dict` that are not in `known`: errors in strict mode,
warnings otherwise. `Globals` is an explicit escape hatch and never reported.
"""
function check_unknown!(dict::AbstractDict, known, path::String, ctx::ParseContext)
    for k in keys(dict)
        key = string(k)
        (key in known || key == "Globals") && continue
        message = unknown_key_message(key, known)
        key_path = join_path(path, key)
        ctx.strict ? add_error!(ctx, key_path, message) : add_warning!(ctx, key_path, message)
    end
    return nothing
end

"""
    parse_section(T, dict, path, ctx)

Builds `T` from `dict`, reports unknown keys, and runs `derive`.
"""
function parse_section(T, dict::AbstractDict, path::String, ctx::ParseContext)
    p = build(T, dict, path, ctx)
    check_unknown!(dict, aliases(T), path, ctx)
    return p === nothing ? nothing : derive(p)
end
```

In `convert.jl`, inside `convert_value`, replace

```julia
    elseif T isa DataType && T <: Vector
        return _convert_vector(T, raw, path, ctx)
    end
    throw(ArgumentError("convert_value: unsupported type $T"))
end
```

with

```julia
    elseif T isa DataType && T <: Vector
        return _convert_vector(T, raw, path, ctx)
    elseif T isa DataType && T <: Dict && T.parameters[1] === String &&
           is_params(T.parameters[2])
        return _convert_named_sections(T, raw, path, ctx)
    elseif is_params(T)
        raw isa AbstractDict ||
            return _fail(ctx, path,
                         "expected a section of `key: value` entries, got $(_describe(raw))")
        section = parse_section(T, raw, path, ctx)
        return section === nothing ? FAILED : section
    end
    throw(ArgumentError("convert_value: unsupported type $T"))
end

function _convert_named_sections(::Type{Dict{String,V}}, raw, path::String,
                                 ctx::ParseContext) where {V}
    raw isa AbstractDict ||
        return _fail(ctx, path, "expected named entries, got $(_describe(raw))")
    out = Dict{String,V}()
    ok = true
    for (k, item) in raw
        v = convert_value(V, item, join_path(path, string(k)), ctx)
        if v === FAILED
            ok = false
        else
            out[string(k)] = v
        end
    end
    return ok ? out : FAILED
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. test/unit_tests/Support/Parameters/Spec/run_spec_tests.jl`
Expected: PASS. If the docstring test fails, check that `Base.@__doc__` wraps the struct expression in the macro output.

---

### Task 5: Model registry, `NoModel`, `Composite`, `parse_model`

**Files:**
- Create: `src/Support/Parameters/Spec/registry.jl`
- Create: `src/Support/Parameters/Spec/model.jl`
- Modify: `src/Support/Parameters/Spec/ParameterSpec.jl`
- Create: `test/unit_tests/Support/Parameters/Spec/ut_model.jl`
- Modify: `test/unit_tests/Support/Parameters/Spec/spec_tests.jl`

**Interfaces:**
- Consumes: Tasks 1–4 (`ParamsDefinitionError`, `ParseContext`, `add_error!`, `join_path`, `suggest`, `is_params`, `parameter_spec`, `build`, `check_unknown!`, `aliases`, `derive`, `_typename`).
- Produces:
  - `struct UnavailableModel; reason::String; end`
  - `register_model!(category::Symbol, name::AbstractString, T::Type)`; convenience `register_material`, `register_damage`, `register_thermal`, `register_additive`, `register_degradation`, `register_pre_calculation` (each `(name, T)`)
  - `register_unavailable!(category::Symbol, name::AbstractString; reason = "requires a license that is not available")`
  - `lookup_model(category, name)` → `Type`, `UnavailableModel` or `nothing`
  - `registered_names(category)::Vector{String}` (sorted, available models only)
  - `struct NoModel end`; `struct Composite{P<:Tuple}; parts::P; end`
  - `parse_model(category::Symbol, dict::Union{Nothing,AbstractDict}, path::String, ctx::ParseContext; name_key::String)` → `NoModel()`, a model struct, a `Composite`, or `nothing`

- [ ] **Step 1: Write the failing test**

`test/unit_tests/Support/Parameters/Spec/ut_model.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const PS = PeriLab.ParameterSpec

PS.@params struct UTElastic
    youngs_modulus::Dependent = req("Young's Modulus"; min = 0)
    poissons_ratio::Float64 = req("Poisson's Ratio")
end

PS.@params struct UTPlastic
    youngs_modulus::Dependent = req("Young's Modulus"; min = 0)
    yield_stress::Float64 = req("Yield Stress"; min = 0)
end

PS.@params struct UTConflicting
    youngs_modulus::Float64 = req("Young's Modulus")
end

PS.register_model!(:ut_material, "UT Elastic", UTElastic)
PS.register_model!(:ut_material, "UT Plastic", UTPlastic)
PS.register_model!(:ut_material, "UT Conflicting", UTConflicting)
PS.register_unavailable!(:ut_material, "UT Licensed")

const UT_KEY = "Material Model"
const UT_PATH = "Models.\"Material Models\".Steel"

function ut_parse(dict; strict = true)
    ctx = PS.ParseContext(strict = strict)
    return PS.parse_model(:ut_material, dict, UT_PATH, ctx; name_key = UT_KEY), ctx
end

@testset "registry" begin
    @test PS.lookup_model(:ut_material, "UT Elastic") === UTElastic
    @test PS.lookup_model(:ut_material, "nope") === nothing
    @test PS.lookup_model(:no_such_category, "UT Elastic") === nothing
    @test PS.lookup_model(:ut_material, "UT Licensed") isa PS.UnavailableModel
    @test PS.registered_names(:ut_material) == ["UT Conflicting", "UT Elastic", "UT Plastic"]
    PS.register_model!(:ut_material, "UT Elastic", UTElastic)      # same type again: fine
    e = try
        PS.register_model!(:ut_material, "UT Elastic", UTPlastic)
    catch err
        err
    end
    @test e isa PS.ParamsDefinitionError
    @test e.msg == "ut_material model \"UT Elastic\" is already registered by UTElastic"
    @test_throws PS.ParamsDefinitionError PS.register_model!(:ut_material, "X", Float64)
    PS.register_unavailable!(:ut_material, "UT Elastic")           # never hides a real model
    @test PS.lookup_model(:ut_material, "UT Elastic") === UTElastic
    PS.register_unavailable!(:ut_material, "UT Later")
    PS.register_model!(:ut_material, "UT Later", UTPlastic)        # real replaces stub
    @test PS.lookup_model(:ut_material, "UT Later") === UTPlastic
end

@testset "no model" begin
    m, ctx = ut_parse(nothing)
    @test m === PS.NoModel()
    @test isempty(ctx.errors)
end

@testset "single model" begin
    m, ctx = ut_parse(Dict{String,Any}(UT_KEY => "UT Elastic", "Young's Modulus" => 210.0,
                                       "Poisson's Ratio" => 0.3))
    @test isempty(ctx.errors)
    @test m isa UTElastic{PS.Constant}
end

@testset "composite model" begin
    m, ctx = ut_parse(Dict{String,Any}(UT_KEY => "UT Elastic + UT Plastic",
                                       "Young's Modulus" => 210.0,
                                       "Poisson's Ratio" => 0.3, "Yield Stress" => 5.0))
    @test isempty(ctx.errors)
    @test m isa PS.Composite{Tuple{UTElastic{PS.Constant},UTPlastic{PS.Constant}}}
    @test isconcretetype(typeof(m))
    @test m.parts[2].yield_stress == 5.0
    @test m.parts[1].youngs_modulus.value == m.parts[2].youngs_modulus.value == 210.0
end

@testset "composite: a key is unknown only if no part declares it" begin
    m, ctx = ut_parse(Dict{String,Any}(UT_KEY => "UT Elastic + UT Plastic",
                                       "Young's Modulus" => 210.0,
                                       "Poisson's Ratio" => 0.3, "Yeild Stress" => 5.0))
    @test m === nothing
    msgs = Dict(e.path => e.message for e in ctx.errors)
    @test msgs["$UT_PATH.\"Yield Stress\""] == "missing (required by UT Plastic)"
    @test msgs["$UT_PATH.\"Yeild Stress\""] == "unknown key — did you mean \"Yield Stress\"?"
    @test length(ctx.errors) == 2
end

@testset "model name errors" begin
    name_path = "$UT_PATH.\"Material Model\""
    cases = [("UT Elastc", "model \"UT Elastc\" not found — did you mean \"UT Elastic\"?"),
             ("Completely Different",
              "model \"Completely Different\" not found; it may require a licensed module"),
             ("UT Licensed",
              "model \"UT Licensed\" requires a license that is not available"),
             ("UT Elastic + ", "empty model name in \"UT Elastic + \""),
             (42, "expected a model name, got 42")]
    for (name, expected) in cases
        m, ctx = ut_parse(Dict{String,Any}(UT_KEY => name, "Young's Modulus" => 1.0,
                                           "Poisson's Ratio" => 0.3))
        @test m === nothing
        @test ctx.errors[1].path == name_path
        @test ctx.errors[1].message == expected
    end
    m, ctx = ut_parse(Dict{String,Any}("Young's Modulus" => 1.0))
    @test m === nothing
    @test ctx.errors[1].path == name_path
    @test ctx.errors[1].message == "missing (names the model to use)"
end

@testset "alias type conflict between combined models" begin
    m, ctx = ut_parse(Dict{String,Any}(UT_KEY => "UT Elastic + UT Conflicting",
                                       "Young's Modulus" => 1.0, "Poisson's Ratio" => 0.3))
    @test m === nothing
    @test ctx.errors[1].path == "$UT_PATH.\"Young's Modulus\""
    @test ctx.errors[1].message ==
          "declared as Dependent by \"UT Elastic\" but as Float64 by \"UT Conflicting\"; models combined with + must agree"
end
```

Add `"ut_model.jl"` to the list in `spec_tests.jl`:

```julia
for file in ["ut_errors.jl", "ut_dependent.jl", "ut_convert.jl", "ut_params.jl", "ut_model.jl"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. test/unit_tests/Support/Parameters/Spec/run_spec_tests.jl`
Expected: FAIL in `ut_model.jl` with `UndefVarError: register_model! not defined in PeriLab.ParameterSpec`.

- [ ] **Step 3: Write the implementation**

In `ParameterSpec.jl`, replace

```julia
include("params_macro.jl")
include("build.jl")
```

with

```julia
include("params_macro.jl")
include("build.jl")
include("registry.jl")
include("model.jl")
```

`src/Support/Parameters/Spec/registry.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

export register_model!, register_unavailable!, register_material, register_damage,
       register_thermal, register_additive, register_degradation, register_pre_calculation

"A model name that is known (e.g. from a license manifest) but not loaded."
struct UnavailableModel
    reason::String
end

# category => model name => @params struct type or UnavailableModel.
# Only mutated at runtime (module `__init__` functions, runtime-loaded modules),
# never during precompilation.
const REGISTRY = Dict{Symbol,Dict{String,Any}}()

"""
    register_model!(category, name, T)

Makes model `name` with parameter struct `T` available in `category`. Call it
from your module's `__init__()`.
"""
function register_model!(category::Symbol, name::AbstractString, T::Type)
    is_params(T) ||
        throw(ParamsDefinitionError("register_model!: $(_typename(T)) is not an @params struct"))
    models = get!(Dict{String,Any}, REGISTRY, category)
    existing = get(models, name, nothing)
    if existing isa Type && existing !== T
        throw(ParamsDefinitionError("$category model \"$name\" is already registered by $(_typename(existing))"))
    end
    models[String(name)] = T
    return nothing
end

"""
    register_unavailable!(category, name; reason)

Records a model name that exists but cannot be used (e.g. no license), so a
deck referencing it gets a specific message. Never hides a registered model.
"""
function register_unavailable!(category::Symbol, name::AbstractString;
                               reason::AbstractString = "requires a license that is not available")
    models = get!(Dict{String,Any}, REGISTRY, category)
    haskey(models, name) && models[name] isa Type && return nothing
    models[String(name)] = UnavailableModel(String(reason))
    return nothing
end

lookup_model(category::Symbol, name::AbstractString) = get(get(REGISTRY, category,
                                                                Dict{String,Any}()), name,
                                                            nothing)

function registered_names(category::Symbol)
    models = get(REGISTRY, category, Dict{String,Any}())
    return sort!([name for (name, entry) in models if entry isa Type])
end

register_material(name::AbstractString, T::Type) = register_model!(:material, name, T)
register_damage(name::AbstractString, T::Type) = register_model!(:damage, name, T)
register_thermal(name::AbstractString, T::Type) = register_model!(:thermal, name, T)
register_additive(name::AbstractString, T::Type) = register_model!(:additive, name, T)
register_degradation(name::AbstractString, T::Type) = register_model!(:degradation, name, T)
function register_pre_calculation(name::AbstractString, T::Type)
    register_model!(:pre_calculation, name, T)
end
```

`src/Support/Parameters/Spec/model.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

export NoModel, Composite

"Placeholder for a block that has no model of a category."
struct NoModel end

"Models combined with `+` in the input deck; each part keeps its own struct."
struct Composite{P<:Tuple}
    parts::P
end

function _check_alias_conflicts!(names, types, path::String, ctx::ParseContext)
    seen = Dict{String,Tuple{String,Any}}()
    ok = true
    for (name, T) in zip(names, types), fs in parameter_spec(T)
        if haskey(seen, fs.alias)
            other_name, other_type = seen[fs.alias]
            if other_type !== fs.type
                add_error!(ctx, join_path(path, fs.alias),
                           "declared as $(_typename(other_type)) by \"$other_name\" but as $(_typename(fs.type)) by \"$name\"; models combined with + must agree")
                ok = false
            end
        else
            seen[fs.alias] = (name, fs.type)
        end
    end
    return ok
end

"""
    parse_model(category, dict, path, ctx; name_key)

Builds the model(s) named by `dict[name_key]` (e.g. "Material Model"), which
may combine several registered models with `+`. All parts read their keys from
the same `dict`; a key is unknown only if no part declares it. Returns
`NoModel()` for `dict === nothing`, the single model struct, a `Composite`, or
`nothing` if there were errors.
"""
function parse_model(category::Symbol, dict::Union{Nothing,AbstractDict}, path::String,
                     ctx::ParseContext; name_key::String)
    dict === nothing && return NoModel()
    name_path = join_path(path, name_key)
    raw = get(dict, name_key, nothing)
    if !(raw isa AbstractString)
        add_error!(ctx, name_path,
                   raw === nothing ? "missing (names the model to use)" :
                   "expected a model name, got $(_describe(raw))")
        return nothing
    end
    names = String.(strip.(split(raw, "+")))
    if any(isempty, names)
        add_error!(ctx, name_path, "empty model name in \"$raw\"")
        return nothing
    end
    types = Any[]
    for name in names
        entry = lookup_model(category, name)
        if entry === nothing
            suggestion = suggest(name, registered_names(category))
            add_error!(ctx, name_path,
                       suggestion === nothing ?
                       "model \"$name\" not found; it may require a licensed module" :
                       "model \"$name\" not found — did you mean \"$suggestion\"?")
        elseif entry isa UnavailableModel
            add_error!(ctx, name_path, "model \"$name\" $(entry.reason)")
        else
            push!(types, entry)
        end
    end
    length(types) == length(names) || return nothing
    _check_alias_conflicts!(names, types, path, ctx) || return nothing
    parts = Any[]
    for (name, T) in zip(names, types)
        part = build(T, dict, path, ctx; owner = name)
        push!(parts, part === nothing ? nothing : derive(part))
    end
    known = Set{String}([name_key])
    for T in types
        union!(known, aliases(T))
    end
    check_unknown!(dict, known, path, ctx)
    any(isnothing, parts) && return nothing
    return length(parts) == 1 ? parts[1] : Composite(Tuple(parts))
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. test/unit_tests/Support/Parameters/Spec/run_spec_tests.jl`
Expected: PASS.

---

### Task 6: Binding dependent values, end-to-end module author flow, type stability, full regression

**Files:**
- Create: `src/Support/Parameters/Spec/bind.jl`
- Modify: `src/Support/Parameters/Spec/ParameterSpec.jl`
- Create: `test/unit_tests/Support/Parameters/Spec/ut_end_to_end.jl`
- Modify: `test/unit_tests/Support/Parameters/Spec/spec_tests.jl`

**Interfaces:**
- Consumes: Tasks 1–5 (`Table1D`, `bind_table!`, `Composite`, `is_params`, `parameter_spec`, `join_path`, `add_error!`, `parse_model`, `register_model!`, `combine`, `Constant`, `value`, `derive`, `@params`).
- Produces:
  - `bind_dependents!(x, lookup, path::String, ctx::ParseContext)` — walks `@params` structs, `Composite`s and `AbstractDict`s; binds every `Table1D` to `lookup(field_name)`, which must return a `Vector{Float64}` or `nothing`. Errors go to `ctx`. Phase 3 calls it with a lookup backed by `Data_Manager` (`NP1` fields).

- [ ] **Step 1: Write the failing test**

`test/unit_tests/Support/Parameters/Spec/ut_end_to_end.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const PS = PeriLab.ParameterSpec

# A module written the way a module author would, loaded at runtime like a
# licensed module (registration in __init__).
const UT_AUTHOR_MODULE = raw"""
module UTAuthorMaterial
using PeriLab.ParameterSpec: @params, register_model!, value, combine, Constant
import PeriLab.ParameterSpec: derive

@params struct UTLinearElastic
    youngs_modulus::Dependent = req("Young's Modulus"; min = 0, quantity = :stress)
    poissons_ratio::Float64 = req("Poisson's Ratio"; min = -1, max = 0.5)
    shear_modulus::Dependent = opt("Shear Modulus"; default = 0.0, quantity = :stress)
end

function derive(p::UTLinearElastic)
    G = combine((E, nu) -> E / (2 * (1 + nu)), p.youngs_modulus, Constant(p.poissons_ratio))
    return UTLinearElastic(p.youngs_modulus, p.poissons_ratio, G)
end

function sum_shear(p::UTLinearElastic, nodes)
    G = p.shear_modulus
    s = 0.0
    for iID in nodes
        s += value(G, iID)
    end
    return s
end

function __init__()
    register_model!(:ut_author, "UT Linear Elastic", UTLinearElastic)
end
end
"""
Base.include_string(@__MODULE__, UT_AUTHOR_MODULE, "UTAuthorMaterial.jl")

const UT_E2E_DIR = mktempdir()
write(joinpath(UT_E2E_DIR, "E_T.txt"),
      "header: Temperature Young's_Modulus\n0.0 200.0\n100.0 180.0\n")

function ut_author_model(youngs_modulus)
    ctx = PS.ParseContext(directory = UT_E2E_DIR)
    m = PS.parse_model(:ut_author,
                       Dict{String,Any}("Material Model" => "UT Linear Elastic",
                                        "Young's Modulus" => youngs_modulus,
                                        "Poisson's Ratio" => 0.3),
                       "Models.\"Material Models\".A", ctx; name_key = "Material Model")
    return m, ctx
end

@testset "runtime-loaded module registers in __init__" begin
    @test PS.lookup_model(:ut_author, "UT Linear Elastic") === UTAuthorMaterial.UTLinearElastic
end

@testset "constant parameters: derive, type stability, no allocations" begin
    m, ctx = ut_author_model(260.0)
    @test isempty(ctx.errors)
    @test m.shear_modulus isa PS.Constant
    @test m.shear_modulus.value ≈ 100.0
    @test isconcretetype(typeof(m))
    @test (@inferred UTAuthorMaterial.sum_shear(m, 1:2)) ≈ 200.0
    UTAuthorMaterial.sum_shear(m, 1:2)
    @test (@allocated UTAuthorMaterial.sum_shear(m, 1:2)) == 0
end

@testset "table parameters: derive, bind, evaluate" begin
    m, ctx = ut_author_model("E_T.txt")
    @test isempty(ctx.errors)
    @test m.youngs_modulus isa PS.Table1D
    @test m.shear_modulus isa PS.Table1D
    temperature = [0.0, 100.0]
    PS.bind_dependents!(m, name -> name == "Temperature" ? temperature : nothing,
                        "Models.\"Material Models\".A", ctx)
    @test isempty(ctx.errors)
    @test m.youngs_modulus.bound[] && m.shear_modulus.bound[]
    @test (@inferred UTAuthorMaterial.sum_shear(m, 1:2)) ≈ 200.0 / 2.6 + 180.0 / 2.6
end

@testset "binding errors" begin
    m, _ = ut_author_model("E_T.txt")
    ctx = PS.ParseContext()
    PS.bind_dependents!(m, name -> nothing, "B", ctx)
    @test ctx.errors[1].path == "B.\"Young's Modulus\""
    @test ctx.errors[1].message ==
          "field \"Temperature\" required by $(joinpath(UT_E2E_DIR, "E_T.txt")) does not exist"
    m2, _ = ut_author_model("E_T.txt")
    ctx2 = PS.ParseContext()
    PS.bind_dependents!(m2, name -> [1, 2], "B", ctx2)
    @test ctx2.errors[1].message ==
          "field \"Temperature\" must be a per-node Vector{Float64}, got Vector{Int64}"
end

@testset "binding walks composites and named entries" begin
    a, _ = ut_author_model("E_T.txt")
    b, _ = ut_author_model("E_T.txt")
    c, _ = ut_author_model("E_T.txt")
    temperature = [0.0, 100.0]
    lookup = name -> temperature
    ctx = PS.ParseContext()
    PS.bind_dependents!(PS.Composite((a, b)), lookup, "C", ctx)
    PS.bind_dependents!(Dict("block_1" => c), lookup, "D", ctx)
    @test isempty(ctx.errors)
    @test a.youngs_modulus.bound[] && b.youngs_modulus.bound[] && c.youngs_modulus.bound[]
end
```

Add `"ut_end_to_end.jl"` to the list in `spec_tests.jl`:

```julia
for file in ["ut_errors.jl", "ut_dependent.jl", "ut_convert.jl", "ut_params.jl", "ut_model.jl",
             "ut_end_to_end.jl"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. test/unit_tests/Support/Parameters/Spec/run_spec_tests.jl`
Expected: FAIL in `ut_end_to_end.jl` with `UndefVarError: bind_dependents! not defined in PeriLab.ParameterSpec` (the first two testsets pass).

- [ ] **Step 3: Write the implementation**

In `ParameterSpec.jl`, replace

```julia
include("registry.jl")
include("model.jl")
```

with

```julia
include("registry.jl")
include("model.jl")
include("bind.jl")
```

`src/Support/Parameters/Spec/bind.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

"""
    bind_dependents!(x, lookup, path, ctx)

Binds every `Table1D` inside `x` (an `@params` struct, a `Composite`, or a
dict of them) to its node field. `lookup(field_name)` returns the field's
`Vector{Float64}` or `nothing`. Call after the fields exist; problems are added
to `ctx`.
"""
function bind_dependents!(x, lookup, path::String, ctx::ParseContext)
    is_params(typeof(x)) || return nothing
    for fs in parameter_spec(typeof(x))
        bind_dependents!(getfield(x, fs.name), lookup, join_path(path, fs.alias), ctx)
    end
    return nothing
end

function bind_dependents!(t::Table1D, lookup, path::String, ctx::ParseContext)
    field = lookup(t.field_name)
    if field === nothing
        add_error!(ctx, path, "field \"$(t.field_name)\" required by $(t.source) does not exist")
    elseif !(field isa Vector{Float64})
        add_error!(ctx, path,
                   "field \"$(t.field_name)\" must be a per-node Vector{Float64}, got $(typeof(field))")
    else
        bind_table!(t, field)
    end
    return nothing
end

function bind_dependents!(c::Composite, lookup, path::String, ctx::ParseContext)
    for part in c.parts
        bind_dependents!(part, lookup, path, ctx)
    end
    return nothing
end

function bind_dependents!(d::AbstractDict, lookup, path::String, ctx::ParseContext)
    for (k, v) in d
        bind_dependents!(v, lookup, join_path(path, string(k)), ctx)
    end
    return nothing
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. test/unit_tests/Support/Parameters/Spec/run_spec_tests.jl`
Expected: PASS for all six spec test files.

- [ ] **Step 5: Run the full test suite (regression + Aqua)**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS, including the `Aqua` testset (no undefined exports, no piracy) and the new `ParameterSpec` testset under Support → Parameters. All pre-existing tests are unchanged and must still pass, since no existing code path uses `ParameterSpec` yet.

- [ ] **Step 6: Leave changes uncommitted**

Run: `git status --short`
Expected: the new files under `src/Support/Parameters/Spec/` and `test/unit_tests/Support/Parameters/Spec/`, plus modified `src/PeriLab.jl` and `test/runtests.jl`, all unstaged. Do not commit (user instruction).
