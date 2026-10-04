# Typed Input Parameters — Phase 3a: Material Models Declare and Validate Their Parameters

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Every material model declares its YAML parameters as an `@params` struct, and `Models."Material Models"` is parsed and validated against them. This replaces the legacy validator for materials. The runtime keeps reading the validated dict, which is the spec's phase 2–3 bridge.

**Architecture:**
- Shared base part:
  - The Material factory declares one `MaterialBaseParams` struct with the shared material vocabulary:
    - symmetry;
    - isotropic, anisotropic and C-matrix moduli;
    - State Factor ID, Accuracy Order, Zero Energy Control, Bond Associated, Linear Strain;
    - Flaw Function.
  - It registers that struct as the base part of the `:material` category.
- Parsing:
  - `ParameterSpec.parse_model` builds the base part from the same flat YAML block as the named model(s).
  - It returns `WithBase(base, model)`.
- Module declarations:
  - Each material module declares only its own extra keys (Yield Stress, UMAT File, …) and registers its struct in `__init__`.
  - Indexed keys such as `Property_1 … Property_N` are declared with a new `key_patterns` hook.
- Dependent values: optional dependent values (`Union{Nothing,Dependent}`) become supported, because anisotropic moduli may come from `.txt` tables.

**Tech Stack:** Julia 1.12, `PeriLab.ParameterSpec`, `PeriLab.InputDeck`.

**Spec:** `docs/superpowers/specs/2026-10-03-typed-input-parameters-design.md` (§2.1, 2.2, 2.4, 2.6, 2.7, 5)

## Decisions taken with the user for this phase

- The shared material keys go into **one base part declared in the Material factory**. It is added to every material model automatically, and modules do not repeat the moduli.
- **Validate first:** this plan covers declarations, registration, parsing and validation. BlockModels and the runtime switch to the structs follow in phase 3b. `derive` (completing K/E/G/ν) also moves in 3b.

## Global Constraints

- Full YAML backward compatibility. Every shipped deck parses in strict mode (`ut_golden_decks.jl`). A deck key that no model reads is removed from the deck, with a ledger note, because before this phase it was silently ignored.
- Unknown keys are errors (strict), with the `--no_strict` / `Strict Validation: false` opt-out.
- No fixed units; `quantity` is documentation only.
- Dependent values depend on exactly one node field; a `.txt` path given for a non-`Dependent` field is an input error.
- Registration runs at load time in `__init__`, never during precompilation.
- A licensed install registers its modules when loaded; a model name that is not registered is an input error ("not found; it may require a licensed module").
- Structs are immutable. The runtime still uses `Data_Manager.get_properties(block, "Material Model")` until phase 3b.
- Branch `feature/typed-input-parameters`; one commit per task; never merge. Read the full-suite result before committing a task that runs it.

## Review Focus

1. **A composite with a misspelled key, e.g. `Correspondence Elastic + Correspondence Plastic` with `Yeild Stress`**, reports both the missing `Yield Stress` and the unknown key with a suggestion. Pinned in Task 5 (`misspelled key in a composite`).
2. **An orthotropic deck reading `Young's Modulus X` from a `.txt` table** (`test_symmetry/symmetry_eng_dep.yaml`) parses to a `Table1D`. A table missing the column is an input error naming the file. Pinned in Task 1 (`optional Dependent from a file`) and by the golden test.
3. **A UMAT block with `Property_1 … Property_27`** passes. `Property_3: "abc"` is a type error on `Property_3`, and `Propery_3` is an unknown key. Pinned in Task 2 (`key patterns`) and Task 4 (`UMAT properties`).
4. **A block using a material name that is not registered (`"Bond based Elastic"`)** fails with a near-match suggestion, not a runtime abort. Pinned in Task 5 (`unknown material model`).
5. **A model that needs a key the base doesn't have (PD Solid Plastic without `Yield Stress`)** reports "missing (required by PD Solid Plastic)" at parse time. Pinned in Task 3 (`PD Solid Plastic requires Yield Stress`).

---

## Shared test commands

Create the runner once in the plan's workspace (`<workspace>` = output of `sdd-workspace PLAN_FILE`), as `<workspace>/unit.jl`:

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

- ParameterSpec tests: `julia --project=. test/unit_tests/Support/Parameters/Spec/run_spec_tests.jl 2>&1 | tail -15`
- Input tests (golden decks included): `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl 2>&1 | tail -15`
- Single file: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/<path>`
- Full suite (~30 min, background, read the tail): `julia --project=. -e 'using Pkg; Pkg.test()' > <workspace>/suite.log 2>&1`

`Logging.disable_logging(Logging.Warn)` hides `:warn` records, so assert only `:error` logs.

## Relative imports of `ParameterSpec`

A top-level package module is its own parent, so extra dots stay at the root. Use these depths:

| Module | Import |
|---|---|
| `PeriLab.Solver_Manager.Model_Factory.Material` (`Material_Factory.jl`) | `using .....ParameterSpec: …` |
| `….Material.<Model>` (bond-based, PD solid, Rigid, Material_template) | `using ......ParameterSpec: …` |
| `….Material.Correspondence.<Model>` (Correspondence_Elastic/Plastic/UMAT/VUMAT) | `using .......ParameterSpec: …` |

---

### Task 1: Optional `Dependent` fields

**Files:**
- Modify: `src/Support/Parameters/Spec/params_macro.jl` (`supported_type`, macro loop at the `_is_dependent_type(ftype)` branch, new `_is_optional_dependent_type`)
- Test: `test/unit_tests/Support/Parameters/Spec/ut_params.jl` (append)

**Interfaces:**
- Consumes: `@params`, `convert_value(Union{Nothing,Dependent}, …)` (already delegates to `Dependent`), `bind_dependents!` (already ignores `nothing`).
- Produces: fields `x::Union{Nothing,Dependent} = opt("Key"; default = nothing, …)`. The struct gets a type parameter `T_x <: Union{Nothing,Dependent}`, so every instance stays concrete. `parameter_spec(T)[i].type === Union{Nothing,Dependent}`.

- [ ] **Step 1: Write the failing tests** (append to `ut_params.jl`)

```julia
@testset "optional Dependent fields" begin
    ex = :(PS.@params struct UTOptionalDependent
               youngs_modulus_x::Union{Nothing,Dependent} = opt("Young's Modulus X";
                                                                default = nothing, min = 0)
               poissons_ratio::Float64 = req("Poisson's Ratio")
           end)
    @test ut_definition_error(ex) === nothing
    T = getfield(@__MODULE__, :UTOptionalDependent)
    @test PS.parameter_spec(T)[1].type === Union{Nothing,PS.Dependent}
    dir = mktempdir()
    write(joinpath(dir, "Ex.txt"), "header: Temperature Young's_Modulus_X\n0 200\n100 180\n")
    ctx = PS.ParseContext(directory = dir)
    absent = PS.parse_section(T, Dict{String,Any}("Poisson's Ratio" => 0.3), "M", ctx)
    constant = PS.parse_section(T, Dict{String,Any}("Young's Modulus X" => 210.0,
                                                    "Poisson's Ratio" => 0.3), "M", ctx)
    table = PS.parse_section(T, Dict{String,Any}("Young's Modulus X" => "Ex.txt",
                                                 "Poisson's Ratio" => 0.3), "M", ctx)
    @test isempty(ctx.errors)
    @test absent.youngs_modulus_x === nothing && isconcretetype(typeof(absent))
    @test constant.youngs_modulus_x isa PS.Constant && isconcretetype(typeof(constant))
    @test table.youngs_modulus_x isa PS.Table1D && isconcretetype(typeof(table))
end

@testset "optional Dependent from a file: missing column and min" begin
    T = getfield(@__MODULE__, :UTOptionalDependent)
    dir = mktempdir()
    write(joinpath(dir, "bad.txt"), "header: Temperature Other\n0 1\n1 2\n")
    write(joinpath(dir, "neg.txt"), "header: Temperature Young's_Modulus_X\n0 -1\n1 2\n")
    ctx = PS.ParseContext(directory = dir)
    PS.parse_section(T, Dict{String,Any}("Young's Modulus X" => "bad.txt",
                                         "Poisson's Ratio" => 0.3), "M", ctx)
    @test occursin("has no column \"Young's_Modulus_X\"", ctx.errors[end].message)
    PS.parse_section(T, Dict{String,Any}("Young's Modulus X" => "neg.txt",
                                         "Poisson's Ratio" => 0.3), "M", ctx)
    @test ctx.errors[end].message == "-1.0 is below minimum 0"
end
```

- [ ] **Step 2: Run to verify they fail**

Run: ParameterSpec tests
Expected: "optional Dependent fields" fails, because `ut_definition_error(ex)` returns a `ParamsDefinitionError` ("unsupported field type Union{Nothing, Dependent}"). The follow-up testset errors with `UndefVarError: UTOptionalDependent`.

- [ ] **Step 3: Implement** in `src/Support/Parameters/Spec/params_macro.jl`

In `supported_type`, replace

```julia
        return Nothing <: T && !(S isa Union) && S !== Dependent && supported_type(S)
```

with

```julia
        return Nothing <: T && !(S isa Union) && supported_type(S)
```

After `_is_dependent_type`, add

```julia
"`Union{Nothing,Dependent}` in a field declaration (either order)."
function _is_optional_dependent_type(t)
    (t isa Expr && t.head === :curly && t.args[1] === :Union && length(t.args) == 3) ||
        return false
    members = t.args[2:3]
    return any(==(:Nothing), members) && any(_is_dependent_type, members)
end
```

In the macro loop, replace

```julia
        if _is_dependent_type(ftype)
            typeparam = Symbol("T_", fname)
            push!(typeparams, Expr(:<:, typeparam, Dependent))
            push!(fields, Expr(:(::), fname, typeparam))
            push!(entries, Expr(:tuple, QuoteNode(fname), Dependent, call))
```

with

```julia
        if _is_dependent_type(ftype) || _is_optional_dependent_type(ftype)
            bound = _is_dependent_type(ftype) ? Dependent : Union{Nothing,Dependent}
            typeparam = Symbol("T_", fname)
            push!(typeparams, Expr(:<:, typeparam, bound))
            push!(fields, Expr(:(::), fname, typeparam))
            push!(entries, Expr(:tuple, QuoteNode(fname), bound, call))
```

- [ ] **Step 4: Run to verify they pass**

Run: ParameterSpec tests
Expected: all pass, including the existing definition-error tests (nested Dependent sections stay unsupported).

- [ ] **Step 5: Commit**

```bash
git add src/Support/Parameters/Spec/params_macro.jl test/unit_tests/Support/Parameters/Spec/ut_params.jl
git commit -m "ParameterSpec: optional Dependent fields"
```

---

### Task 2: Base part per category and indexed key patterns in `parse_model`

**Files:**
- Modify: `src/Support/Parameters/Spec/registry.jl` (`register_base!`, `base_model`, export)
- Modify: `src/Support/Parameters/Spec/model.jl` (`WithBase`, `_parse_model`, `_check_key_patterns!`)
- Modify: `src/Support/Parameters/Spec/params_macro.jl` (`key_patterns` hook, export)
- Test: `test/unit_tests/Support/Parameters/Spec/ut_model.jl` (append)

**Interfaces:**
- Consumes: `build`, `check!`, `derive`, `aliases`, `check_unknown!`, `convert_value`, `_check_alias_conflicts!`.
- Produces (exported from `PeriLab.ParameterSpec`):
  - `register_base!(category::Symbol, T::Type)`: one `@params` struct per category, read by every model of the category from the same YAML block. Registering a different type twice throws `ParamsDefinitionError`.
  - `base_model(category::Symbol)`: the registered type or `nothing`.
  - `struct WithBase{B,M}; base::B; model::M; end`: `parse_model` returns it when the category has a base. `model` is the single model struct or a `Composite`.
  - `key_patterns(::Type)::Vector{Pair{Regex,Any}}`: default empty. Extend it for indexed keys; matching keys are type-checked and are not unknown. Their values are not stored in the struct (phase 3b decides storage).

- [ ] **Step 1: Write the failing tests** (append to `ut_model.jl`)

```julia
PS.@params struct UTBase
    symmetry::Union{Nothing,String} = opt("Symmetry"; default = nothing)
    youngs_modulus::Union{Nothing,Float64} = opt("Young's Modulus"; default = nothing, min = 0)
end

PS.@params struct UTEmpty
end

PS.@params struct UTPatterned
    file::String = req("File")
end
PS.key_patterns(::Type{UTPatterned}) = [r"^Property_\d+$" => Float64]

PS.@params struct UTBaseConflict
    symmetry::Float64 = req("Symmetry")
end

PS.register_base!(:ut_based, UTBase)
PS.register_model!(:ut_based, "UT Empty", UTEmpty)
PS.register_model!(:ut_based, "UT Patterned", UTPatterned)
PS.register_model!(:ut_based, "UT Base Conflict", UTBaseConflict)

function ut_parse_based(dict; strict = true)
    ctx = PS.ParseContext(strict = strict)
    return PS.parse_model(:ut_based, dict, UT_PATH, ctx; name_key = UT_KEY), ctx
end

@testset "base part registration" begin
    @test PS.base_model(:ut_based) === UTBase
    @test PS.base_model(:ut_material) === nothing
    PS.register_base!(:ut_based, UTBase)                       # same type again: fine
    e = try
        PS.register_base!(:ut_based, UTPatterned)
    catch err
        err
    end
    @test e isa PS.ParamsDefinitionError
    @test e.msg == "ut_based base parameters are already registered by UTBase"
    @test_throws PS.ParamsDefinitionError PS.register_base!(:ut_other, Float64)
end

@testset "base part is read from the same block" begin
    m, ctx = ut_parse_based(Dict{String,Any}(UT_KEY => "UT Empty", "Symmetry" => "isotropic",
                                             "Young's Modulus" => 210.0))
    @test isempty(ctx.errors)
    @test m isa PS.WithBase{UTBase,UTEmpty}
    @test m.base.symmetry == "isotropic" && m.base.youngs_modulus == 210.0
    m, ctx = ut_parse_based(Dict{String,Any}(UT_KEY => "UT Empty + UT Patterned",
                                             "File" => "a.so"))
    @test isempty(ctx.errors)
    @test m isa PS.WithBase{UTBase,PS.Composite{Tuple{UTEmpty,UTPatterned}}}
    @test m.base.symmetry === nothing
end

@testset "base keys: errors, unknown keys and conflicts" begin
    m, ctx = ut_parse_based(Dict{String,Any}(UT_KEY => "UT Empty", "Young's Modulus" => -1.0,
                                             "Youngs Modulus" => 1.0))
    @test m === nothing
    msgs = Dict(e.path => e.message for e in ctx.errors)
    @test msgs["$UT_PATH.\"Young's Modulus\""] == "-1.0 is below minimum 0"
    @test msgs["$UT_PATH.\"Youngs Modulus\""] == "unknown key — did you mean \"Young's Modulus\"?"
    m, ctx = ut_parse_based(Dict{String,Any}(UT_KEY => "UT Base Conflict", "Symmetry" => 1.0))
    @test m === nothing
    @test ctx.errors[1].message ==
          "declared as Union{Nothing, String} by \"base parameters\" but as Float64 by \"UT Base Conflict\"; models combined with + must agree"
end

@testset "key patterns" begin
    m, ctx = ut_parse_based(Dict{String,Any}(UT_KEY => "UT Patterned", "File" => "a.so",
                                             "Property_1" => 1, "Property_27" => 2.5))
    @test isempty(ctx.errors)
    @test m isa PS.WithBase{UTBase,UTPatterned}
    m, ctx = ut_parse_based(Dict{String,Any}(UT_KEY => "UT Patterned", "File" => "a.so",
                                             "Property_3" => "abc", "Propery_4" => 1.0))
    @test m === nothing
    msgs = Dict(e.path => e.message for e in ctx.errors)
    @test msgs["$UT_PATH.Property_3"] == "expected a number, got \"abc\""
    @test startswith(msgs["$UT_PATH.Propery_4"], "unknown key")
    # a pattern of one model does not make the key known for another
    m, ctx = ut_parse_based(Dict{String,Any}(UT_KEY => "UT Empty", "Property_1" => 1.0))
    @test m === nothing
    @test startswith(ctx.errors[1].message, "unknown key")
end
```

- [ ] **Step 2: Run to verify they fail**

Run: ParameterSpec tests
Expected: errors such as `UndefVarError: key_patterns` / `register_base!` not defined in `PeriLab.ParameterSpec`. These happen at file load, so the whole `ut_model.jl` testset errors.

- [ ] **Step 3: Implement**

`src/Support/Parameters/Spec/params_macro.jl`: change the export line to `export @params, derive, check!, key_patterns`, and after the `check!` default add:

```julia
"""
    key_patterns(::Type{T})

Indexed YAML keys of the `@params` struct `T`, e.g. `Property_1 … Property_N`,
as `regex => type` pairs. Matching keys in a model block are type-checked and
are not reported as unknown; their values are not stored in the struct.
"""
key_patterns(::Type) = Pair{Regex,Any}[]
```

`src/Support/Parameters/Spec/registry.jl`: extend the export list with `register_base!, base_model` and append:

```julia
# category => @params struct read by every model of the category
const BASES = Dict{Symbol,Any}()

"""
    register_base!(category, T)

Declares `T` as the base part of `category`: every model of the category reads
`T`'s keys from its YAML block in addition to its own (e.g. the shared material
moduli). Call it from the factory's `__init__()`.
"""
function register_base!(category::Symbol, T::Type)
    is_params(T) ||
        throw(ParamsDefinitionError("register_base!: $(_typename(T)) is not an @params struct"))
    existing = get(BASES, category, nothing)
    if existing !== nothing && existing !== T
        throw(ParamsDefinitionError("$category base parameters are already registered by $(_typename(existing))"))
    end
    BASES[category] = T
    return nothing
end

base_model(category::Symbol) = get(BASES, category, nothing)
```

`src/Support/Parameters/Spec/model.jl`:
- Change the export line to `export NoModel, Composite, WithBase`.
- After `Composite`, add:

```julia
"A model together with the base part of its category (see `register_base!`)."
struct WithBase{B,M}
    base::B
    model::M
end

function _check_key_patterns!(dict::AbstractDict, known::Set{String}, types, path::String,
                              ctx::ParseContext)
    for k in keys(dict)
        key = string(k)
        key in known && continue
        for T in types
            pattern = findfirst(p -> occursin(first(p), key), key_patterns(T))
            pattern === nothing && continue
            convert_value(last(key_patterns(T)[pattern]), dict[k], join_path(path, key), ctx)
            push!(known, key)
            break
        end
    end
    return nothing
end
```

- In `_parse_model`, replace everything from `length(types) == length(names) || return nothing` to the end of the function with:

```julia
    length(types) == length(names) || return nothing
    base = base_model(category)
    all_names = base === nothing ? names : ["base parameters"; names]
    all_types = base === nothing ? types : Any[base; types]
    _check_alias_conflicts!(all_names, all_types, path, ctx) || return nothing
    base_part = nothing
    if base !== nothing
        base_part = build(base, dict, path, ctx; owner = "every $category model")
        base_part === nothing || check!(base_part, path, ctx)
        base_part = base_part === nothing ? nothing : derive(base_part)
    end
    parts = Any[]
    for (name, T) in zip(names, types)
        part = build(T, dict, path, ctx; owner = name)
        part === nothing || check!(part, path, ctx)
        push!(parts, part === nothing ? nothing : derive(part))
    end
    known = Set{String}([name_key])
    for T in all_types
        union!(known, aliases(T))
    end
    _check_key_patterns!(dict, known, all_types, path, ctx)
    check_unknown!(dict, known, path, ctx)
    any(isnothing, parts) && return nothing
    base !== nothing && base_part === nothing && return nothing
    model = length(parts) == 1 ? parts[1] : Composite(Tuple(parts))
    return base === nothing ? model : WithBase(base_part, model)
```

  Note: `check!` errors are recorded in `ctx` and abort via `report!`, but they do not make `build` return `nothing`. Tests therefore assert on `ctx.errors`.

- [ ] **Step 4: Run to verify they pass**

Run: ParameterSpec tests
Expected: all pass (existing `:ut_material` tests unchanged; it has no base).

- [ ] **Step 5: Commit**

```bash
git add src/Support/Parameters/Spec/registry.jl src/Support/Parameters/Spec/model.jl \
        src/Support/Parameters/Spec/params_macro.jl test/unit_tests/Support/Parameters/Spec/ut_model.jl
git commit -m "ParameterSpec: category base part and indexed key patterns"
```

---

### Task 3: Material base part and the non-correspondence material models

**Files:**
- Modify: `src/Models/Material/Material_Factory.jl` (`FlawFunctionParams`, `MaterialBaseParams`, `__init__`)
- Modify: `src/Models/Material/Material_Models/BondBased/Bondbased_Elastic.jl`, `1D_Bondbased_Elastic.jl`, `Unified_Bondbased_Elastic.jl`
- Modify: `src/Models/Material/Material_Models/Ordinary/PD_Solid_Elastic.jl`, `PD_Solid_Plastic.jl`
- Modify: `src/Models/Material/Material_Models/Rigid/Rigid.jl`
- Modify: `src/Models/Material/Material_Models/Material_template/material_template.jl`
- Create: `test/unit_tests/Support/Parameters/Input/ut_material_params.jl`
- Modify: `test/unit_tests/Support/Parameters/Input/input_tests.jl` (include the new file)

**Interfaces:**
- Consumes (Tasks 1–2): `Union{Nothing,Dependent}` fields, `register_base!`, `register_material`, `WithBase`, `parse_model`.
- Produces:
  - `Material.MaterialBaseParams` registered as `base_model(:material)`, plus `Material.FlawFunctionParams`.
  - Registered `:material` names and their structs (each in its module):
    - "Bond-based Elastic" → `BondbasedElasticParams`
    - "1D Bond-based Elastic" → `OneDBondbasedElasticParams` (`id1`, `id2`)
    - "Unified Bond-based Elastic" → `UnifiedBondbasedElasticParams`
    - "PD Solid Elastic" → `PDSolidElasticParams`
    - "PD Solid Plastic" → `PDSolidPlasticParams` (`yield_stress::Dependent`)
    - "Rigid" → `RigidParams`
    - "Material Template" → `MaterialTemplateParams`

- [ ] **Step 1: Write the failing tests**

Create `test/unit_tests/Support/Parameters/Input/ut_material_params.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

const MPS = PeriLab.ParameterSpec
const MAT = PeriLab.Solver_Manager.Model_Factory.Material

function ut_material(dict; directory = "")
    ctx = MPS.ParseContext(directory = directory)
    m = MPS.parse_model(:material, Dict{String,Any}(dict), "Models.\"Material Models\".M",
                        ctx; name_key = "Material Model")
    return m, ctx
end

@testset "material base part" begin
    @test MPS.base_model(:material) === MAT.MaterialBaseParams
    m, ctx = ut_material(Dict("Material Model" => "PD Solid Elastic", "Symmetry" => "isotropic",
                              "Bulk Modulus" => 2.5e3, "Shear Modulus" => 1.15e3,
                              "C11" => 1.0, "State Factor ID" => 2,
                              "Zero Energy Control" => "Global", "Bond Associated" => true,
                              "Flaw Function" => Dict("Active" => true, "Function" => "Pre-defined",
                                                      "Flaw Size" => 0.2, "Flaw Magnitude" => 0.5,
                                                      "Flaw Location X" => 1.0,
                                                      "Flaw Location Y" => 0.0)))
    @test isempty(ctx.errors)
    @test m isa MPS.WithBase
    @test m.base.bulk_modulus == 2.5e3 && m.base.c11 == 1.0 && m.base.c66 === nothing
    @test m.base.bond_associated && !m.base.linear_strain
    @test m.base.flaw_function.flaw_size == 0.2
end

@testset "non-correspondence material names" begin
    for (name, T) in [("Bond-based Elastic", MAT.Bondbased_Elastic.BondbasedElasticParams),
                      ("1D Bond-based Elastic",
                       MAT.OneD_Bond_Based_Elastic.OneDBondbasedElasticParams),
                      ("Unified Bond-based Elastic",
                       MAT.Unified_Bondbased_Elastic.UnifiedBondbasedElasticParams),
                      ("PD Solid Elastic", MAT.PD_Solid_Elastic.PDSolidElasticParams),
                      ("PD Solid Plastic", MAT.PD_Solid_Plastic.PDSolidPlasticParams),
                      ("Rigid", MAT.Rigid.RigidParams),
                      ("Material Template", MAT.Material_template.MaterialTemplateParams)]
        @test MPS.lookup_model(:material, name) === T
    end
    m, ctx = ut_material(Dict("Material Model" => "1D Bond-based Elastic",
                              "Young's Modulus" => 1.0, "Id1" => 1, "Id2" => 2))
    @test isempty(ctx.errors) && m.model.id2 === 2
    m, ctx = ut_material(Dict("Material Model" => "PD Solid Elastic + PD Solid Plastic",
                              "Bulk Modulus" => 1.0, "Shear Modulus" => 1.0,
                              "Yield Stress" => 5))
    @test isempty(ctx.errors)
    @test m.model isa MPS.Composite
    @test m.model.parts[2].yield_stress == MPS.Constant(5.0)
end

@testset "PD Solid Plastic requires Yield Stress" begin
    m, ctx = ut_material(Dict("Material Model" => "PD Solid Plastic", "Bulk Modulus" => 1.0))
    @test m === nothing
    @test only(ctx.errors).message == "missing (required by PD Solid Plastic)"
    @test only(ctx.errors).path == "Models.\"Material Models\".M.\"Yield Stress\""
end

@testset "material base constraints" begin
    m, ctx = ut_material(Dict("Material Model" => "Bond-based Elastic",
                              "Poisson's Ratio" => 0.7, "Shear Modulus" => -1.0,
                              "Flaw Function" => Dict("Active" => true, "Function" => "Gauss")))
    @test m === nothing
    msgs = Dict(e.path => e.message for e in ctx.errors)
    p = "Models.\"Material Models\".M"
    @test msgs["$p.\"Poisson's Ratio\""] == "0.7 is above maximum 0.5"
    @test msgs["$p.\"Shear Modulus\""] == "-1.0 is below minimum 0"
    @test msgs["$p.\"Flaw Function\".Function"] == "\"Gauss\" is not one of: \"Pre-defined\""
end
```

In `test/unit_tests/Support/Parameters/Input/input_tests.jl`, append `"ut_material_params.jl"` to the file list.

- [ ] **Step 2: Run to verify they fail**

Run: input tests
Expected: `ut_material_params.jl` fails: `base_model(:material) === nothing`, `UndefVarError: MaterialBaseParams`, and model names not found.

- [ ] **Step 3: Implement**

`src/Models/Material/Material_Factory.jl`: after `using ....ModuleLoader: find_module_files, create_module_specifics` add

```julia
using .....ParameterSpec: @params, Dependent, register_base!
```

Then, before `global module_list = …`, so that the base exists before the model files are included, add:

```julia
@params struct FlawFunctionParams
    active::Bool = req("Active")
    function_name::String = req("Function"; allowed = ["Pre-defined"])
    flaw_size::Union{Nothing,Float64} = opt("Flaw Size"; default = nothing, min = 0,
                                            quantity = :length)
    flaw_magnitude::Union{Nothing,Float64} = opt("Flaw Magnitude"; default = nothing, min = 0,
                                                 max = 1)
    flaw_location_x::Union{Nothing,Float64} = opt("Flaw Location X"; default = nothing)
    flaw_location_y::Union{Nothing,Float64} = opt("Flaw Location Y"; default = nothing)
    flaw_location_z::Union{Nothing,Float64} = opt("Flaw Location Z"; default = nothing)
end

"""
    MaterialBaseParams

Keys every material model may use (read from the same YAML block as the model's
own keys): symmetry, elastic constants, and options shared by the material
framework. Which elastic constants are needed depends on the model and the
symmetry; they are completed at initialisation.
"""
@params struct MaterialBaseParams
    symmetry::Union{Nothing,String} = opt("Symmetry"; default = nothing,
                                          description = "e.g. isotropic, isotropic plane strain, isotropic plane stress, orthotropic, anisotropic")
    youngs_modulus::Union{Nothing,Float64} = opt("Young's Modulus"; default = nothing, min = 0,
                                                 quantity = :stress)
    poissons_ratio::Union{Nothing,Float64} = opt("Poisson's Ratio"; default = nothing,
                                                 min = -1, max = 0.5)
    bulk_modulus::Union{Nothing,Float64} = opt("Bulk Modulus"; default = nothing, min = 0,
                                               quantity = :stress)
    shear_modulus::Union{Nothing,Float64} = opt("Shear Modulus"; default = nothing, min = 0,
                                                quantity = :stress)
    youngs_modulus_x::Union{Nothing,Dependent} = opt("Young's Modulus X"; default = nothing,
                                                     min = 0, quantity = :stress)
    youngs_modulus_y::Union{Nothing,Dependent} = opt("Young's Modulus Y"; default = nothing,
                                                     min = 0, quantity = :stress)
    youngs_modulus_z::Union{Nothing,Dependent} = opt("Young's Modulus Z"; default = nothing,
                                                     min = 0, quantity = :stress)
    poissons_ratio_xy::Union{Nothing,Dependent} = opt("Poisson's Ratio XY"; default = nothing)
    poissons_ratio_yz::Union{Nothing,Dependent} = opt("Poisson's Ratio YZ"; default = nothing)
    poissons_ratio_xz::Union{Nothing,Dependent} = opt("Poisson's Ratio XZ"; default = nothing)
    shear_modulus_xy::Union{Nothing,Dependent} = opt("Shear Modulus XY"; default = nothing,
                                                     min = 0, quantity = :stress)
    shear_modulus_yz::Union{Nothing,Dependent} = opt("Shear Modulus YZ"; default = nothing,
                                                     min = 0, quantity = :stress)
    shear_modulus_xz::Union{Nothing,Dependent} = opt("Shear Modulus XZ"; default = nothing,
                                                     min = 0, quantity = :stress)
    c11::Union{Nothing,Float64} = opt("C11"; default = nothing, quantity = :stress)
    c12::Union{Nothing,Float64} = opt("C12"; default = nothing, quantity = :stress)
    c13::Union{Nothing,Float64} = opt("C13"; default = nothing, quantity = :stress)
    c14::Union{Nothing,Float64} = opt("C14"; default = nothing, quantity = :stress)
    c15::Union{Nothing,Float64} = opt("C15"; default = nothing, quantity = :stress)
    c16::Union{Nothing,Float64} = opt("C16"; default = nothing, quantity = :stress)
    c22::Union{Nothing,Float64} = opt("C22"; default = nothing, quantity = :stress)
    c23::Union{Nothing,Float64} = opt("C23"; default = nothing, quantity = :stress)
    c24::Union{Nothing,Float64} = opt("C24"; default = nothing, quantity = :stress)
    c25::Union{Nothing,Float64} = opt("C25"; default = nothing, quantity = :stress)
    c26::Union{Nothing,Float64} = opt("C26"; default = nothing, quantity = :stress)
    c33::Union{Nothing,Float64} = opt("C33"; default = nothing, quantity = :stress)
    c34::Union{Nothing,Float64} = opt("C34"; default = nothing, quantity = :stress)
    c35::Union{Nothing,Float64} = opt("C35"; default = nothing, quantity = :stress)
    c36::Union{Nothing,Float64} = opt("C36"; default = nothing, quantity = :stress)
    c44::Union{Nothing,Float64} = opt("C44"; default = nothing, quantity = :stress)
    c45::Union{Nothing,Float64} = opt("C45"; default = nothing, quantity = :stress)
    c46::Union{Nothing,Float64} = opt("C46"; default = nothing, quantity = :stress)
    c55::Union{Nothing,Float64} = opt("C55"; default = nothing, quantity = :stress)
    c56::Union{Nothing,Float64} = opt("C56"; default = nothing, quantity = :stress)
    c66::Union{Nothing,Float64} = opt("C66"; default = nothing, quantity = :stress)
    state_factor_id::Union{Nothing,Int64} = opt("State Factor ID"; default = nothing, min = 1,
                                                description = "index of the state variable that scales the elastic constants")
    accuracy_order::Union{Nothing,Int64} = opt("Accuracy Order"; default = nothing, min = 1)
    zero_energy_control::Union{Nothing,String} = opt("Zero Energy Control";
                                                     default = nothing,
                                                     description = "zero energy control model, e.g. Global")
    bond_associated::Bool = opt("Bond Associated"; default = false)
    linear_strain::Bool = opt("Linear Strain"; default = false)
    flaw_function::Union{Nothing,FlawFunctionParams} = opt("Flaw Function"; default = nothing)
end

# registration runs at load time, never during precompilation
__init__() = register_base!(:material, MaterialBaseParams)
```

Each model module below gets, directly after its existing `using …Data_Manager` line:

```julia
using ......ParameterSpec: @params, register_material
```

(`PD_Solid_Plastic.jl` additionally imports `Dependent`.) Each module also gets, after its `material_name()` function, the struct and its `__init__`:

`Bondbased_Elastic.jl`:

```julia
"Parameters of Bond-based Elastic beyond the shared material keys (none)."
@params struct BondbasedElasticParams
end
__init__() = register_material("Bond-based Elastic", BondbasedElasticParams)
```

`1D_Bondbased_Elastic.jl`:

```julia
@params struct OneDBondbasedElasticParams
    id1::Int64 = req("Id1"; min = 1)
    id2::Int64 = req("Id2"; min = 1)
end
__init__() = register_material("1D Bond-based Elastic", OneDBondbasedElasticParams)
```

`Unified_Bondbased_Elastic.jl`:

```julia
"Parameters of Unified Bond-based Elastic beyond the shared material keys (none)."
@params struct UnifiedBondbasedElasticParams
end
__init__() = register_material("Unified Bond-based Elastic", UnifiedBondbasedElasticParams)
```

`PD_Solid_Elastic.jl`:

```julia
"Parameters of PD Solid Elastic beyond the shared material keys (none)."
@params struct PDSolidElasticParams
end
__init__() = register_material("PD Solid Elastic", PDSolidElasticParams)
```

`PD_Solid_Plastic.jl` (import line: `using ......ParameterSpec: @params, Dependent, register_material`):

```julia
@params struct PDSolidPlasticParams
    yield_stress::Dependent = req("Yield Stress"; min = 0, quantity = :stress)
end
__init__() = register_material("PD Solid Plastic", PDSolidPlasticParams)
```

`Rigid.jl`:

```julia
"Parameters of Rigid beyond the shared material keys (none)."
@params struct RigidParams
end
__init__() = register_material("Rigid", RigidParams)
```

`material_template.jl`:
- Struct and `__init__`:

```julia
"""
    MaterialTemplateParams

Declare the YAML keys your material needs beyond the shared material keys
(Symmetry, Young's Modulus, … are already available). Example:

    my_parameter::Float64 = req("My Parameter"; min = 0, description = "...")
"""
@params struct MaterialTemplateParams
end
__init__() = register_material("Material Template", MaterialTemplateParams)
```

- Also add one line to the template's module docstring or top comment:

```julia
# Declare your parameters with @params (see MaterialTemplateParams) and register them in __init__.
```

If a module already defines `__init__`, merge the registration into it. Run `grep -n '__init__' src/Models/Material -r` first; no existing module defines one.

- [ ] **Step 4: Run to verify they pass**

Run: input tests
Expected: `ut_material_params.jl` passes, and the rest stays green. Material Models are not parsed by `read_input` yet, so the golden test is unaffected.

- [ ] **Step 5: Commit**

```bash
git add src/Models/Material/Material_Factory.jl src/Models/Material/Material_Models/BondBased \
        src/Models/Material/Material_Models/Ordinary/PD_Solid_Elastic.jl \
        src/Models/Material/Material_Models/Ordinary/PD_Solid_Plastic.jl \
        src/Models/Material/Material_Models/Rigid/Rigid.jl \
        src/Models/Material/Material_Models/Material_template/material_template.jl \
        test/unit_tests/Support/Parameters/Input/ut_material_params.jl \
        test/unit_tests/Support/Parameters/Input/input_tests.jl
git commit -m "Material base part and parameter declarations of bond-based, PD solid and rigid models"
```

---

### Task 4: Correspondence material models

**Files:**
- Modify: `src/Models/Material/Material_Models/Correspondence/Correspondence_Elastic.jl`, `Correspondence_Plastic.jl`, `Correspondence_UMAT.jl`, `Correspondence_VUMAT.jl`
- Modify: `src/Models/Material/Material_Models/Material_template/correspondence_template.jl` (documentation only; this file is not loaded)
- Test: `test/unit_tests/Support/Parameters/Input/ut_material_params.jl` (append)

**Interfaces:**
- Consumes (Tasks 1–3): `register_material`, `key_patterns`, `Dependent`, base part.
- Produces (each module's struct registered under `:material`):
  - "Correspondence Elastic" → `CorrespondenceElasticParams`.
  - "Correspondence Plastic" → `CorrespondencePlasticParams` (`yield_stress::Dependent`).
  - "Correspondence UMAT" → `CorrespondenceUMATParams`, with fields:
    - `file::String`, `number_of_properties::Int64`;
    - `number_of_state_variables`, `predefined_field_names`, `umat_material_name`, `umat_name`, all `Union{Nothing,…}`;
    - key pattern `Property_N => Float64`.
  - "Correspondence VUMAT" → `CorrespondenceVUMATParams`: the same without `predefined_field_names`, with `vumat_material_name`, `vumat_name`.

- [ ] **Step 1: Write the failing tests** (append to `ut_material_params.jl`)

```julia
const CORR = MAT.Correspondence

@testset "correspondence material names" begin
    for (name, T) in [("Correspondence Elastic", CORR.Correspondence_Elastic.CorrespondenceElasticParams),
                      ("Correspondence Plastic", CORR.Correspondence_Plastic.CorrespondencePlasticParams),
                      ("Correspondence UMAT", CORR.Correspondence_UMAT.CorrespondenceUMATParams),
                      ("Correspondence VUMAT", CORR.Correspondence_VUMAT.CorrespondenceVUMATParams)]
        @test MPS.lookup_model(:material, name) === T
    end
    m, ctx = ut_material(Dict("Material Model" => "Correspondence Elastic + Correspondence Plastic",
                              "Symmetry" => "isotropic plane strain", "Bulk Modulus" => 1.0,
                              "Shear Modulus" => 1.0, "Yield Stress" => 2.0))
    @test isempty(ctx.errors)
    @test m.model.parts[2].yield_stress == MPS.Constant(2.0)
end

@testset "UMAT properties" begin
    props = Dict("Property_$i" => Float64(i) for i in 1:27)
    m, ctx = ut_material(merge(Dict{String,Any}("Material Model" => "Correspondence UMAT",
                                                "File" => "libusertest.so",
                                                "Number of Properties" => 27,
                                                "Number of State Variables" => 0,
                                                "UMAT Material Name" => "test"), props))
    @test isempty(ctx.errors)
    @test m.model.number_of_properties == 27 && m.model.umat_name === nothing
    m, ctx = ut_material(Dict("Material Model" => "Correspondence VUMAT", "File" => "a.so",
                              "Number of Properties" => 3, "Property_3" => "abc",
                              "Propery_2" => 1.0))
    @test m === nothing
    msgs = Dict(e.path => e.message for e in ctx.errors)
    p = "Models.\"Material Models\".M"
    @test msgs["$p.Property_3"] == "expected a number, got \"abc\""
    @test startswith(msgs["$p.Propery_2"], "unknown key")
    m, ctx = ut_material(Dict("Material Model" => "Correspondence UMAT", "File" => "a.so"))
    @test only(ctx.errors).message == "missing (required by Correspondence UMAT)"
end
```

- [ ] **Step 2: Run to verify they fail**

Run: input tests
Expected: `UndefVarError` for the four `…Params` types, and the model names are not found.

- [ ] **Step 3: Implement**

Each of the four correspondence modules gets, after its `using …Data_Manager` line,

```julia
using .......ParameterSpec: @params, register_material
```

(`Correspondence_Plastic.jl` adds `Dependent`; UMAT/VUMAT add `key_patterns` via `import .......ParameterSpec: key_patterns`.) After each module's `correspondence_name()` function add:

`Correspondence_Elastic.jl`:

```julia
"Parameters of Correspondence Elastic beyond the shared material keys (none)."
@params struct CorrespondenceElasticParams
end
__init__() = register_material("Correspondence Elastic", CorrespondenceElasticParams)
```

`Correspondence_Plastic.jl` (`using .......ParameterSpec: @params, Dependent, register_material`):

```julia
@params struct CorrespondencePlasticParams
    yield_stress::Dependent = req("Yield Stress"; min = 0, quantity = :stress)
end
__init__() = register_material("Correspondence Plastic", CorrespondencePlasticParams)
```

`Correspondence_UMAT.jl`:

```julia
@params struct CorrespondenceUMATParams
    file::String = req("File"; description = "UMAT library, relative to the input deck")
    number_of_properties::Int64 = req("Number of Properties"; min = 1,
                                      description = "number of Property_N values passed to the UMAT")
    number_of_state_variables::Union{Nothing,Int64} = opt("Number of State Variables";
                                                          default = nothing, min = 0)
    predefined_field_names::Union{Nothing,String} = opt("Predefined Field Names";
                                                        default = nothing)
    umat_material_name::Union{Nothing,String} = opt("UMAT Material Name"; default = nothing)
    umat_name::Union{Nothing,String} = opt("UMAT name"; default = nothing,
                                           description = "name of the UMAT routine, default UMAT")
end
key_patterns(::Type{CorrespondenceUMATParams}) = [r"^Property_\d+$" => Float64]
__init__() = register_material("Correspondence UMAT", CorrespondenceUMATParams)
```

`Correspondence_VUMAT.jl`:

```julia
@params struct CorrespondenceVUMATParams
    file::String = req("File"; description = "VUMAT library, relative to the input deck")
    number_of_properties::Int64 = req("Number of Properties"; min = 1,
                                      description = "number of Property_N values passed to the VUMAT")
    number_of_state_variables::Union{Nothing,Int64} = opt("Number of State Variables";
                                                          default = nothing, min = 0)
    vumat_material_name::Union{Nothing,String} = opt("VUMAT Material Name"; default = nothing)
    vumat_name::Union{Nothing,String} = opt("VUMAT name"; default = nothing,
                                            description = "name of the VUMAT routine, default VUMAT")
end
key_patterns(::Type{CorrespondenceVUMATParams}) = [r"^Property_\d+$" => Float64]
__init__() = register_material("Correspondence VUMAT", CorrespondenceVUMATParams)
```

`correspondence_template.jl`: below the `correspondence_name()` function, add this documentation block. The file is not loaded, so keep it a comment:

```julia
# Declare the YAML keys of your model beyond the shared material keys and register them:
#
#     using .......ParameterSpec: @params, register_material
#     @params struct CorrespondenceTemplateParams
#         my_parameter::Float64 = req("My Parameter"; min = 0)
#     end
#     __init__() = register_material("Correspondence Template", CorrespondenceTemplateParams)
```

- [ ] **Step 4: Run to verify they pass**

Run: input tests
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add src/Models/Material/Material_Models/Correspondence \
        src/Models/Material/Material_Models/Material_template/correspondence_template.jl \
        test/unit_tests/Support/Parameters/Input/ut_material_params.jl
git commit -m "Parameter declarations of the correspondence material models"
```

---

### Task 5: `read_input` validates `Material Models`

**Files:**
- Modify: `src/Support/Parameters/Input/input.jl` (`PeriLabInput.materials`, `read_input`)
- Modify: `src/Support/Parameters/parameter_handling.jl` (`validate_models` skips `Material Models`)
- Modify: decks:
  - `examples/Compact_Tension/CompactTension.yaml` (remove material `Density`)
  - `examples/DCB/DCBmodel_ba.yaml`, `examples/Dogbone/Dogbone_plastic_corr.yaml`, `examples/Dogbone/Dogbone_plastic_corr_static.yaml` (remove `Yield Stress` from the `Correspondence Elastic` materials)
  - `examples/Training/Input.yaml` (remove `my new parameter`)
- Modify tests:
  - `test/unit_tests/Support/Parameters/Input/ut_validate_yaml.jl`
  - `test/unit_tests/Support/Parameters/Input/ut_read_input.jl`
  - `test/unit_tests/Support/Parameters/ut_parameter_handling.jl` (material name `"a"` → registered name)
- Create: `test/unit_tests/Support/Parameters/Input/ut_material_models.jl` and add it to `input_tests.jl`

**Interfaces:**
- Consumes (Tasks 2–4): `parse_model(:material, …)` returns `WithBase` or `nothing`.
- Produces: `PeriLabInput.materials::Dict{String,Any}`, mapping material name to its `WithBase` struct (empty if there are no `Material Models`). The `PeriLabInput` constructor gains this argument after `contact`: `PeriLabInput(sections, contact, materials, models, globals)`.

- [ ] **Step 1: Write the failing tests**

Create `test/unit_tests/Support/Parameters/Input/ut_material_models.jl`:

```julia
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
    @test haskey(input.models, "Material Models")     # raw dict stays for the runtime bridge
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
```

Add `"ut_material_models.jl"` to `input_tests.jl`.

In `test/unit_tests/Support/Parameters/Input/ut_validate_yaml.jl`:
- In `ut_valid_params`, `"Material Model" => "a"` becomes `"Material Model" => "Bond-based Elastic"`.
- Rename `@testset "legacy model validation still applies"` to `@testset "material model names are validated"` (its body, `"Material Model" => 5`, still aborts).
- Append:

```julia
@testset "an unknown material key aborts with the typed error" begin
    params = ut_valid_params()
    params["PeriLab"]["Models"]["Material Models"]["mat_1"]["Youngs Modulus"] = 1.0
    @test_logs (:error, r"did you mean \"Young's Modulus\"") match_mode=:any @test_throws PeriLab.PeriLabError begin
        PeriLab.Parameter_Handling.validate_yaml(params)
    end
end
```

In `test/unit_tests/Support/Parameters/Input/ut_read_input.jl`, in `@testset "minimal deck"`, after `@test haskey(input.models, "Material Models")`, add:

```julia
    @test input.materials["Mat"] isa PS.WithBase
```

In `test/unit_tests/Support/Parameters/ut_parameter_handling.jl`, `@testset "ut_validate_yaml"`, last `params = …`: change `"Material Model" => "a"` to `"Material Model" => "Bond-based Elastic"`. Leave the earlier cases in that testset unchanged; they abort for other reasons.

- [ ] **Step 2: Run to verify they fail**

Run: input tests
Expected:
- `ut_material_models.jl` errors with `type PeriLabInput has no field materials`, or fails because a misspelled key produces no error.
- `ut_read_input.jl` "minimal deck" errors on `input.materials`.
- "an unknown material key aborts…" fails (legacy validation only warns).

- [ ] **Step 3: Implement**

`src/Support/Parameters/Input/input.jl`:
- Add the field to `PeriLabInput` after `contact`:

```julia
struct PeriLabInput
    sections::PeriLabSections
    contact::Union{Nothing,ContactInput}
    materials::Dict{String,Any}
    models::Dict{String,Any}
    globals::Dict{String,Any}
end
```

- Update its docstring sentence to: "`materials` holds the typed material models (`ParameterSpec.WithBase`), keyed by name; `models` stays the raw `Models` dict for the categories not yet migrated and for the runtime until phase 3b."
- Add after `_SPECIAL_KEYS`:

```julia
"Typed `Material Models`, keyed by name; errors go to `ctx`."
function parse_materials(models::AbstractDict, ctx::ParseContext)
    materials = Dict{String,Any}()
    raw = get(models, "Material Models", nothing)
    raw === nothing && return materials
    path = join_path("Models", "Material Models")
    if !(raw isa AbstractDict)
        add_error!(ctx, path,
                   "expected a section of `key: value` entries, got $(ParameterSpec._describe(raw))")
        return materials
    end
    for (name, entry) in raw
        entry_path = join_path(path, string(name))
        if !(entry isa AbstractDict)
            add_error!(ctx, entry_path,
                       "expected a section of `key: value` entries, got $(ParameterSpec._describe(entry))")
            continue
        end
        model = ParameterSpec.parse_model(:material,
                                          Dict{String,Any}(string(k) => v for (k, v) in entry),
                                          entry_path, ctx; name_key = "Material Model")
        model === nothing || (materials[string(name)] = model)
    end
    return materials
end
```

- In `read_input`, after the `models` checks (just before `globals = get(deck, "Globals", …)`), add

```julia
    materials = models isa AbstractDict ? parse_materials(models, ctx) : Dict{String,Any}()
```

- Change the constructor call to

```julia
    input = PeriLabInput(sections, contact, materials,
                         Dict{String,Any}(string(k) => v for (k, v) in models),
                         globals isa AbstractDict ?
                         Dict{String,Any}(string(k) => v for (k, v) in globals) :
                         Dict{String,Any}())
```

`src/Support/Parameters/parameter_handling.jl`, `validate_models`: after `models isa Dict || return true`, add

```julia
    # Material Models are validated by the typed declarations (InputDeck.parse_materials)
    models = Dict{Any,Any}(k => v for (k, v) in models if k != "Material Models")
```

and update its docstring's first sentence to "Legacy validation of the `Models` section, except `Material Models` (typed since phase 3a)".

Then run `grep -rn 'PeriLabInput(' src test` and update any other constructor call. The plan expects only `input.jl`.

Decks: remove exactly these lines (find them with `grep -n`):
- `examples/Compact_Tension/CompactTension.yaml`: the `Density:` line inside `Models: Material Models: Aluminium` (not the one under `Blocks`).
- `examples/DCB/DCBmodel_ba.yaml`, `examples/Dogbone/Dogbone_plastic_corr.yaml`, `examples/Dogbone/Dogbone_plastic_corr_static.yaml`: the `Yield Stress:` line of the `Correspondence Elastic` material.
- `examples/Training/Input.yaml`: the `my new parameter: 20` line.

Before deleting, check that no `Damage Models`/other block of the same deck relies on that line (it is inside the material block, so it was read by no one). Ledger one `Ruling:` line listing these removals.

- [ ] **Step 4: Run to verify they pass**

Run: input tests
Expected: all pass, including `ut_golden_decks.jl`, where every shipped deck's materials now parse in strict mode. If a golden deck fails on a material key, do not widen a declaration without checking the code reads that key. Ledger what you found and either declare it (if read) or remove it from the deck (if not read).

Then run the legacy parameter test: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Support/Parameters/ut_parameter_handling`
Expected: pass.

- [ ] **Step 5: Full suite, then commit**

Full suite in background; read the tail. Expected: 0 failures, 0 errors. If an in-memory test deck uses an unregistered material name and is validated (`validate_yaml`/`validate_input`/`read_input`), change it to a registered name with matching keys. Ledger it.

```bash
git add src/Support/Parameters/Input/input.jl src/Support/Parameters/parameter_handling.jl \
        examples/Compact_Tension/CompactTension.yaml examples/DCB/DCBmodel_ba.yaml \
        examples/Dogbone/Dogbone_plastic_corr.yaml examples/Dogbone/Dogbone_plastic_corr_static.yaml \
        examples/Training/Input.yaml test/unit_tests/Support/Parameters
git commit -m "Material Models are validated against the typed declarations"
```
