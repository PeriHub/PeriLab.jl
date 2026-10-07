# Typed Input Parameters — Phase 3f: Additive, Degradation, Thermal, Pre Calculation

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** The four remaining model categories (Additive, Degradation, Thermal incl. HETVAL, Pre Calculation) are validated into typed structs when the input is read. They are stored per block as typed values and read through struct fields. Afterwards no model category is stored as a properties dict. Surface Correction is already typed since phase 2.

**Architecture:**
- **Pattern.** Every category follows phase 3d (Damage):
  - `@params` structs registered per model (`register_additive`, `register_degradation`, `register_thermal`, `register_pre_calculation`);
  - a typed `input.<category>` dict;
  - one typed value per block in `Data_Manager`;
  - dispatch through `parentmodule(typeof(model))` instead of name lookups.
- **Storage.**
  - The new categories share one generic slot, `Data_Manager.set_block_model(category, block, model)` / `get_block_model(category, block)`.
  - Material and Damage keep their own slots. Phase 4 folds all slots into `BlockModels`.
- **Thermal.**
  - A shared base part holds `Thermal Conductivity`, which the critical time step reads for every thermal model.
  - `+` composites are kept.
  - `Thermal Flow`'s `Type` becomes an enum (spec §2.7).
  - `Thermal Expansion Coefficient` is a number or a list. `ParameterSpec` gains that field type.
- **Additive.**
  - The only additive models are licensed. They are loaded once at startup, before the input is read (spec §2.5), so they can register.
  - The golden-deck test treats "requires a license" errors as expected.
- **Pre Calculation.**
  - The sections are on/off switches: `Pre Calculation Global` and `Pre Calculation Models`.
  - Each block stores its active pre-calculation names, in run order.
  - The pre-calculation modules take no parameters any more.
- **No property dicts.** `read_properties` stops building property dicts altogether. The `Model_Factory` helpers become typed-only.

**Tech Stack:** Julia 1.12, `PeriLab.ParameterSpec`, `PeriLab.InputDeck`, `PeriLab.Data_Manager`.

**Spec:** `docs/superpowers/specs/2026-10-03-typed-input-parameters-design.md` (§2.2, §2.3, §2.5 license-gated modules, §2.7 module interface and enums, §3 Use, §5 phase 3)

## Global Constraints

- **Numerical results unchanged.** The full suite stays green, including:
  - `test_thermal_flow`, `test_thermal_expansion`, `test_thermal_decomp`, `test_Hetval`, `test_Multistep`;
  - `test_Newmark` (thermal expansion) and the matrix-based thermal decks;
  - every deck that runs pre-calculations (correspondence, bond-associated, `test_calculation`, `test_PD_Solid_Elastic`, `test_Surface_Correction`).
- **Golden decks.** Every deck under `test/` and `examples/` still parses in strict mode (`ut_golden_decks.jl`).
  - The only accepted error is "requires a license that is not available" for licensed model names.
  - A deck key that no code reads for the named model is removed from the deck, with a `Ruling:` per deck.
  - The allowlist is never extended.
- **Licensed modules** (the additive model "Simple") are maintained in-house. They must adopt the typed module interface of the Additive template (Task 2). Nothing in this repo can test them.
- **Read fields directly; functions only for real logic.**
- **Imports.** Factories (`Thermal`, `Additive`, `Degradation`, `Pre_Calculation`) import `ParameterSpec` as `.....ParameterSpec`, like `Damage`. Model modules inside a factory use `......ParameterSpec`.
- **`ut_MPI.jl`** runs standalone under mpiexec without `helper.jl`.
- **`/dev/shm` is 64 MB.** Delete stale `/dev/shm/mpich_shm_*` files before a suite run when no MPI process runs.
- **Git.** Branch `feature/typed-input-parameters`; one commit per task; never merge. Read the full-suite result before committing.

## Review Focus

1. **Composite thermal models.** `Thermal Flow + Heat Transfer`, other orders, and `Thermal Flow+Heat Transfer` without spaces run every part, in deck order. Pinned in Task 5 (`thermal composite runs every part`).
2. **`Environmental Temperature`** given as an expression string (evaluated with `t`) still works and equals the numeric form. Pinned in Task 5 (`environmental temperature expression`).
3. **A deck naming an unavailable licensed additive model** gets a clear input error, not a crash at init. Pinned in Task 2 (`additive input errors`).
4. **Pre Calculation switches.**
   - A block-specific `Pre Calculation Model` replaces `Pre Calculation Global`.
   - `Bond Associated Deformation Gradient: false` (56 decks) is accepted.
   - `Bond Associated Deformation Gradient: true` is an error that names the replacement.
   - An unknown switch gets a suggestion.

   Pinned in Task 6 (`pre-calculation switches`) and Task 7 (`read_properties stores pre-calculations`).
5. **Thermal critical time step.** It reads `Thermal Conductivity` of the block's thermal model, whatever the model. A block without a thermal model gives a mechanical-only time step. Pinned in Task 5 (`thermal critical time step`).

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
- Input tests: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Support/Parameters/Input/ut_category_models unit_tests/Support/Parameters/Input/ut_golden_decks`
- Full suite (~30 min, background, read the tail): `julia --project=. -e 'using Pkg; Pkg.test()' > <workspace>/suite.log 2>&1`

---

### Task 1: Shared machinery — number-or-list fields, generic block slot, block reference helper

**Files:**
- Modify: `src/Support/Parameters/Spec/params_macro.jl`, `src/Support/Parameters/Spec/convert.jl` (field type `Union{Float64,Vector{Float64}}`)
- Modify: `src/Core/Data_manager.jl`, `src/Core/Data_manager/data_manager_checkpoint.jl` (slot `Block Models`)
- Modify: `src/Models/Model_Factory.jl` (`typed_block_model`, `block_typed_model`)
- Modify: `src/Support/Parameters/Input/input.jl` (`parse_single_models`; `parse_damages` uses it)
- Modify: `test/helper.jl` (`typed_model`)
- Test: `test/unit_tests/Support/Parameters/Spec/ut_convert.jl` (append). Create `test/unit_tests/Models/ut_block_models.jl` and include it in `test/runtests.jl` directly after `ut_Model_Factory.jl`.

**Interfaces:**
- Produces:
  - `@params` field type `Union{Float64,Vector{Float64}}`: a YAML number becomes a `Float64`, a list becomes a `Vector{Float64}`, anything else is the error "expected a number or a list of numbers, got …".
  - `Data_Manager.set_block_model(category::String, block::Int64, model)`.
  - `Data_Manager.get_block_model(category::String, block::Int64)`, which returns `nothing` if the block has none.
  - `Model_Factory.typed_block_model(block::Int64, name::String)`.
  - `Model_Factory.block_typed_model(input, block_name::String, field::Symbol, category::String, models::AbstractDict)`. It returns the parsed model or `nothing`, and aborts with the legacy messages:
    - "<category> is defined in blocks, but no <category>s definition block exists"
    - "<category> model with name <name> is defined in blocks, but missing in the <category>s definition."
  - `InputDeck.parse_single_models(models, section, category, name_key, ctx)`: like `parse_models`, plus the error "<category> models cannot be combined with +" (e.g. "damage models cannot be combined with +").
  - Test helper `typed_model(category::Symbol, raw; name_key::String)`: the parsed model, or `WithBase` for categories with a base part.

- [ ] **Step 1: Write the failing tests**

Append to `test/unit_tests/Support/Parameters/Spec/ut_convert.jl`:

```julia
@params struct UTNumberOrList
    alpha::Union{Float64,Vector{Float64}} = req("Alpha"; min = 0)
end

@testset "number or list fields" begin
    ctx = PS.ParseContext()
    @test PS.parse_section(UTNumberOrList, Dict("Alpha" => 2), "t", ctx).alpha === 2.0
    @test PS.parse_section(UTNumberOrList, Dict("Alpha" => [1, 2.5]), "t", ctx).alpha ==
          [1.0, 2.5]
    @test isempty(ctx.errors)
    PS.parse_section(UTNumberOrList, Dict("Alpha" => "x"), "t", ctx)
    @test ctx.errors[end].message == "expected a number or a list of numbers, got \"x\""
    PS.parse_section(UTNumberOrList, Dict("Alpha" => [1.0, -1.0]), "t", ctx)
    @test occursin("must be ≥ 0", ctx.errors[end].message)
end
```

Before writing it, check how `ut_convert.jl` names the `ParameterSpec` module (`PS`) and how it imports `@params`, and match that. Also check the min-violation message text in `convert.jl` (`grep -n 'must be' src/Support/Parameters/Spec/*.jl`) and use the real wording in the last assertion.

Create `test/unit_tests/Models/ut_block_models.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

const UT_MF = PeriLab.Solver_Manager.Model_Factory

@testset "generic block model slot" begin
    PeriLab.Data_Manager.initialize_data()
    @test PeriLab.Data_Manager.get_block_model("Thermal Model", 1) === nothing
    PeriLab.Data_Manager.set_block_model("Thermal Model", 1, :x)
    @test PeriLab.Data_Manager.get_block_model("Thermal Model", 1) === :x
    @test PeriLab.Data_Manager.get_block_model("Thermal Model", 2) === nothing
    @test PeriLab.Data_Manager.get_block_model("Additive Model", 1) === nothing
    @test "Block Models" in PeriLab.Data_Manager.CHECKPOINT_EXCLUDED_KEYS
end

@testset "block_typed_model" begin
    blocks = Dict("block_1" => Dict("Block ID" => 1, "Density" => 1.0, "Horizon" => 1.0,
                                    "Damage Model" => "Dam"),
                  "block_2" => Dict("Block ID" => 2, "Density" => 1.0, "Horizon" => 1.0))
    models = Dict("Damage Models" => Dict("Dam" => Dict("Damage Model" => "Critical Stretch",
                                                        "Critical Value" => 0.1)))
    input = typed_input(Dict("Blocks" => blocks, "Models" => models))
    d = UT_MF.block_typed_model(input, "block_1", :damage_model, "Damage Model",
                                input.damages)
    @test d === input.damages["Dam"]
    @test UT_MF.block_typed_model(input, "block_2", :damage_model, "Damage Model",
                                  input.damages) === nothing
    @test_logs (:error,
                "Damage Model model with name Dam is defined in blocks, but missing in the Damage Models definition.") @test_throws PeriLab.PeriLabError UT_MF.block_typed_model(input,
                                                                                                                                                                              "block_1",
                                                                                                                                                                              :damage_model,
                                                                                                                                                                              "Damage Model",
                                                                                                                                                                              Dict{String,Any}())
end

@testset "typed_model helper" begin
    m = typed_model(:damage, Dict("Damage Model" => "Critical Stretch", "Critical Value" => 0.1);
                    name_key = "Damage Model")
    @test m isa PeriLab.ParameterSpec.WithBase
end
```

In `test/runtests.jl`, include `unit_tests/Models/ut_block_models.jl` in the `Models` testset, directly after the line that includes `ut_Model_Factory.jl`.

- [ ] **Step 2: Run to verify they fail**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Support/Parameters/Spec/ut_convert unit_tests/Models/ut_block_models`
Expected:
- The `@params` definition fails with a `ParamsDefinitionError`, because the type is not supported.
- `UndefVarError`s for `get_block_model`, `block_typed_model` and `typed_model`.

- [ ] **Step 3: Implement**

`params_macro.jl`:
- In `supported_type(T)`, add as the first line: `T === Union{Float64,Vector{Float64}} && return true`.
- In `_numeric_type(T)`, add before the `_nonnothing` line: `T === Union{Float64,Vector{Float64}} && return true`.

`convert.jl`, in `convert_value`, add after `raw isa T && return _fresh(raw)`:

```julia
    if T === Union{Float64,Vector{Float64}}
        raw isa Real && !(raw isa Bool) && return Float64(raw)
        raw isa AbstractVector && return _convert_vector(Vector{Float64}, raw, path, ctx)
        return _fail(ctx, path, "expected a number or a list of numbers, got $(_describe(raw))")
    end
```

The min/max check runs on the result through `_numbers`, which already handles numbers and vectors.

`Data_manager.jl`:
- Next to `data["Block Damages"]`, add `data["Block Models"] = Dict{String,Dict{Int64,Any}}()`.
- Export and add below `get_block_damage`:

```julia
"""
	set_block_model(category, block, model)

Stores the typed model of `category` (e.g. "Thermal Model") of a block.
"""
function set_block_model(category::String, block::Int64, model)
    get!(Dict{Int64,Any}, data["Block Models"], category)[block] = model
end

"""
	get_block_model(category, block)

The typed model of `category` of a block, or `nothing`.
"""
function get_block_model(category::String, block::Int64)
    return get(get(data["Block Models"], category, Dict{Int64,Any}()), block, nothing)
end
```

`data_manager_checkpoint.jl`: add `"Block Models"` to `CHECKPOINT_EXCLUDED_KEYS`, and change the comment to `# typed models: rebuilt from the input by read_properties`.

`Model_Factory.jl`:
- `typed_block_model` becomes:

```julia
function typed_block_model(block::Int64, name::String)
    name == "Material Model" && return Data_Manager.get_block_material(block)
    name == "Damage Model" && return Data_Manager.get_block_damage(block)
    return Data_Manager.get_block_model(name, block)
end
```

- Add before `read_properties`:

```julia
"""
    block_typed_model(input, block_name, field, category, models)

The parsed model a block names in `field` (e.g. `:damage_model`), or `nothing`.
Aborts if the name or the category's models section is missing.
"""
function block_typed_model(input::PeriLabInput, block_name::String, field::Symbol,
                           category::String, models::AbstractDict)
    block_params = get(input.sections.blocks, block_name, nothing)
    name = block_params === nothing ? nothing : getfield(block_params, field)
    name === nothing && return nothing
    if !haskey(input.models, category * "s")
        @abort "$category is defined in blocks, but no $(category)s definition block exists"
    end
    if !haskey(models, name)
        @abort "$category model with name $name is defined in blocks, but missing in the $(category)s definition."
    end
    return models[name]
end
```

- In `read_properties`, the material loop and the damage loop use it:

```julia
    if material_model
        dof = Data_Manager.get_dof()
        for (block_name, block) in zip(block_name_list, block_id_list)
            wb = block_typed_model(input, block_name, :material_model, "Material Model",
                                   input.materials)
            wb === nothing && continue
            material_name = input.sections.blocks[block_name].material_model
            model_name = String(input.models["Material Models"][material_name]["Material Model"])
            material = Material.block_material(wb, model_name, dof)
            Material.check_material_symmetry(material, dof)
            Data_Manager.set_block_material(block, material)
        end
    end
    for (block_name, block) in zip(block_name_list, block_id_list)
        d = block_typed_model(input, block_name, :damage_model, "Damage Model", input.damages)
        d === nothing || Data_Manager.set_block_damage(block, Damage.block_damage(d))
    end
```

`input.jl`:
- Add `parse_single_models` below `parse_models`. `parse_damages` becomes a one-line call of it.

```julia
"Typed models of `section` that cannot be combined with `+`."
function parse_single_models(models::AbstractDict, section::String, category::Symbol,
                             name_key::String, ctx::ParseContext)
    parsed = parse_models(models, section, category, name_key, ctx)
    for (name, m) in parsed
        part = m isa ParameterSpec.WithBase ? m.model : m
        part isa ParameterSpec.Composite || continue
        add_error!(ctx, join_path(join_path(join_path("Models", section), name), name_key),
                   "$category models cannot be combined with +")
    end
    return parsed
end

parse_damages(models::AbstractDict, ctx::ParseContext) = parse_single_models(models,
                                                                             "Damage Models",
                                                                             :damage,
                                                                             "Damage Model",
                                                                             ctx)
```

`test/helper.jl`, append:

```julia
"""
    typed_model(category, raw; name_key)

The model of `category` parsed from a raw model block the way the input reader
does (`WithBase` for categories with a base part); aborts on input errors.
"""
function typed_model(category::Symbol, raw::AbstractDict; name_key::String)
    ctx = PeriLab.ParameterSpec.ParseContext()
    model = PeriLab.ParameterSpec.parse_model(category,
                                              Dict{String,Any}(string(k) => v
                                                               for (k, v) in raw),
                                              "test", ctx; name_key = name_key)
    PeriLab.ParameterSpec.report!(ctx)
    return model
end
```

- [ ] **Step 4: Run to verify they pass**

Run the Step 2 command, plus `unit_tests/Models/ut_Model_Factory`, `unit_tests/Support/Parameters/Input/ut_damage_models` and `unit_tests/Models/Damage/ut_block_damage`.
Expected: all pass. The existing material and damage `read_properties` tests keep their messages.

- [ ] **Step 5: Full suite, then commit**

```bash
git add src test
git commit -m "Shared typed-model machinery: number-or-list fields, generic block model slot, block reference helper"
```

---

### Task 2: Additive on typed structs; licensed additive models load before the input is read

**Files:**
- Modify: `src/Models/Additive/Additive_Factory.jl`
- Modify: `src/Models/Additive/Additive_template/additive_template.jl`
- Modify: `src/PeriLab.jl` (load licensed models before `IO.initialize_data`)
- Modify: `src/Support/Parameters/Input/input.jl` (`PeriLabInput.additives`)
- Modify: `src/Models/Model_Factory.jl` (`TYPED_CATEGORIES`, `read_properties`)
- Modify: `test/unit_tests/Support/Parameters/Input/ut_golden_decks.jl`
- Test: create `test/unit_tests/Support/Parameters/Input/ut_category_models.jl` and add it to `input_tests.jl` after `"ut_damage_models.jl"`. Append to `test/unit_tests/Models/ut_block_models.jl`.

**Interfaces:**
- Consumes (Task 1): `parse_single_models`, `set_block_model` / `get_block_model`, `block_typed_model`, `typed_block_model`.
- Produces:
  - `PeriLabInput.additives::Dict{String,Any}`, placed after `damages`.
  - `Additive.load_licensed_models()`: includes and registers the licensed additive modules once; returns their module list.
  - Module interface of an additive model:
    - its `@params` struct registered with `register_additive(name, T)`;
    - `init_model(nodes::AbstractVector{Int64}, p::T, block::Int64)`;
    - `compute_model(nodes::AbstractVector{Int64}, p::T, block::Int64, time::Float64, dt::Float64)`;
    - `fields_for_local_synchronization(model::String)`.
  - `Additive_template.AdditiveTemplateParams` (unregistered).

- [ ] **Step 1: Write the failing tests**

Create `test/unit_tests/Support/Parameters/Input/ut_category_models.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

const CMID = PeriLab.InputDeck

function ut_category_deck(models; block = Dict{String,Any}())
    return Dict{String,Any}("Discretization" => Dict{String,Any}("Type" => "Text File",
                                                                 "Input Mesh File" => "m.txt"),
                            "Blocks" => Dict{String,Any}("block_1" => merge(Dict{String,Any}("Block ID" => 1,
                                                                                             "Density" => 1.0,
                                                                                             "Horizon" => 1.0),
                                                                            block)),
                            "Models" => Dict{String,Any}(models),
                            "Solver" => Dict{String,Any}("Initial Time" => 0.0,
                                                         "Final Time" => 1.0,
                                                         "Verlet" => Dict{String,Any}()))
end
ut_category_messages(ctx) = Dict(e.path => e.message for e in ctx.errors)

@testset "additive input errors" begin
    _, ctx = CMID.read_input(ut_category_deck(Dict("Additive Models" => Dict("Add" => Dict("Additive Model" => "Simple",
                                                                                             "Print Temperature" => 100.0)))))
    msg = ut_category_messages(ctx)["Models.\"Additive Models\".Add.\"Additive Model\""]
    @test startswith(msg, "model \"Simple\"")
    input, ctx = CMID.read_input(ut_category_deck(Dict{String,Any}()))
    @test isempty(input.additives)
end
```

`Simple` is not registered without a license. The message is either "not found; it may require a licensed module …" or, when registered unavailable, "requires a license …". Both start with `model "Simple"`.

Append to `test/unit_tests/Models/ut_block_models.jl`:

```julia
@testset "additive dispatch reads the block model" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(3)
    AT = UT_MF.Additive.Additive_template
    p = AT.AdditiveTemplateParams()
    PeriLab.Data_Manager.set_block_model("Additive Model", 1, p)
    @test UT_MF.has_block_model(1, "Additive Model")
    @test UT_MF.block_model_parameters(1, "Additive Model") === p
    UT_MF.Additive.init_model([1, 2, 3], 1)
    UT_MF.Additive.compute_model([1, 2, 3], p, 1, 0.0, 1.0)
    UT_MF.Additive.fields_for_local_synchronization("Additive Model", 1)
    @test length(methods(AT.compute_model)) == 1
end

@testset "licensed additive models load without a license" begin
    withenv("LICENSE_SERVER_URL" => nothing, "PERIHUB_LICENSE_KEY" => nothing,
            "LICENSED_MODULES_CONFIG" => nothing, "LICENSED_MODULES_DIR" => nothing) do
        @test isempty(UT_MF.Additive.load_licensed_models())
    end
end
```

If the template module is not reachable as `Additive.Additive_template`, check `grep -n 'local_user_modules\|include' src/Models/Additive/Additive_Factory.jl`. If the template is not included at all, include it the same way the other factories include theirs (`find_module_files(@__DIR__, "additive_name")`), and record a `Ruling:`.

- [ ] **Step 2: Run to verify they fail**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Support/Parameters/Input/ut_category_models unit_tests/Models/ut_block_models`
Expected:
- `type PeriLabInput has no field additives`;
- `UndefVarError: AdditiveTemplateParams`;
- `UndefVarError: load_licensed_models`.

- [ ] **Step 3: Implement**

`Additive_Factory.jl`:
- Keep `local_user_modules` and their `include`.
- Replace `compute_model`, `init_model` and `fields_for_local_synchronization` (docstrings updated) with:

```julia
global licensed = Any[]

"""
    load_licensed_models()

Includes the licensed additive modules (if a license is configured) so that they
register their parameter structs before the input deck is read. Runs once.
"""
function load_licensed_models()
    isempty(licensed) || return licensed
    global licensed = licensed_modules(@__MODULE__, "Additive")
    return licensed
end

model_module(p) = parentmodule(typeof(p))

function compute_model(nodes::AbstractVector{Int64}, p, block::Int64, time::Float64,
                       dt::Float64)
    # invokelatest: licensed models are defined at runtime
    return Base.@invokelatest model_module(p).compute_model(nodes, p, block, time, dt)
end

function init_model(nodes::AbstractVector{Int64}, block::Int64)
    p = Data_Manager.get_block_model("Additive Model", block)
    Base.@invokelatest model_module(p).init_model(nodes, p, block)
end

function fields_for_local_synchronization(model, block)
    p = Data_Manager.get_block_model("Additive Model", block)
    return Base.@invokelatest model_module(p).fields_for_local_synchronization(model)
end
```

Drop `create_module_specifics` from the `ModuleLoader` import.

`additive_template.jl`, after its `using` lines:

```julia
using ......ParameterSpec: @params, register_additive

"""
    AdditiveTemplateParams

Declare the YAML keys your additive model needs. Register it under your model name
by uncommenting `__init__` (the template stays unregistered so that a copy never
collides with it). Licensed modules register the same way; they are loaded before
the input deck is read (`Additive.load_licensed_models`).
"""
@params struct AdditiveTemplateParams
end

# __init__() = register_additive("Additive Template", AdditiveTemplateParams)
```

Use the template's `Data_Manager` dot count plus one. Its `init_model` / `compute_model` get the signatures from the Interfaces block, replacing the `Dict` argument with `p::AdditiveTemplateParams`. Update their docstrings and the `@info` lines that mention the parameter dict, as in the damage template.

`src/PeriLab.jl`: directly before the line `@timeit "IO.initialize_data" params, input,`, add:

```julia
            # licensed models register their parameters before the input deck is read
            Solver_Manager.Model_Factory.Additive.load_licensed_models()
```

`input.jl`:
- Add the field `additives::Dict{String,Any}` after `damages` and mention it in the docstring.
- In `read_input`, add `additives = models isa AbstractDict ? parse_single_models(models, "Additive Models", :additive, "Additive Model", ctx) : Dict{String,Any}()` and pass it to the constructor after `damages`.

`Model_Factory.jl`:
- `const TYPED_CATEGORIES = ("Material Model", "Damage Model", "Additive Model")`.
- In `read_properties`, after the damage loop:

```julia
    for (block_name, block) in zip(block_name_list, block_id_list)
        a = block_typed_model(input, block_name, :additive_model, "Additive Model",
                              input.additives)
        a === nothing || Data_Manager.set_block_model("Additive Model", block, a)
    end
```

`ut_golden_decks.jl`:
- At the top, after `UT_DECK_ALLOWLIST`, add:

```julia
# licensed models: without a license a deck naming them cannot be validated
PeriLab.ParameterSpec.register_unavailable!(:additive, "Simple")
ut_license_error(e) = endswith(e.message, "requires a license that is not available")
```

- Change the error filter to `errors = filter(e -> e.severity == :error && !ut_license_error(e), ctx.errors)`.

- [ ] **Step 4: Run to verify they pass**

Run the Step 2 command, plus `unit_tests/Support/Parameters/Input/ut_golden_decks`.
Expected: all pass.

- [ ] **Step 5: Full suite, then commit**

```bash
git add src test
git commit -m "Additive models on typed structs; licensed additive models load before the input is read"
```

---

### Task 3: Degradation on typed structs

**Files:**
- Modify: `src/Models/Degradation/Degradation_Factory.jl`, `thermal_decomposition.jl`, `BondBased_Corrosion.jl`, `Degradation_template/degradation_template.jl`
- Modify: `src/Support/Parameters/Input/input.jl` (`PeriLabInput.degradations`)
- Modify: `src/Models/Model_Factory.jl` (`TYPED_CATEGORIES`, `read_properties`)
- Test: append to `ut_category_models.jl` and `ut_block_models.jl`

**Interfaces:**
- Consumes (Task 1): as Task 2.
- Produces:
  - `PeriLabInput.degradations::Dict{String,Any}`, placed after `additives`.
  - `Thermal_Decomposition.ThermalDecompositionParams(decomposition_temperature::Float64)`, registered as "Thermal Decomposition".
  - `Bondbased_Corrosion.BondbasedCorrosionParams()`, registered as "Bond-based Corrosion".
  - `Degradation_template.DegradationTemplateParams` (unregistered).
  - Model interface: `init_model(nodes, p, block)`, `compute_model(nodes, p, block, time, dt)`, `fields_for_local_synchronization(model)`.

- [ ] **Step 1: Write the failing tests**

Append to `ut_category_models.jl`:

```julia
@testset "degradation models are parsed into typed structs" begin
    input, ctx = CMID.read_input(ut_category_deck(Dict("Degradation Models" => Dict("Deg" => Dict("Degradation Model" => "Thermal Decomposition",
                                                                                                    "Decomposition Temperature" => 50)))))
    @test isempty(ctx.errors)
    @test input.degradations["Deg"].decomposition_temperature == 50.0
    _, ctx = CMID.read_input(ut_category_deck(Dict("Degradation Models" => Dict("Deg" => Dict("Degradation Model" => "Thermal Decomposition")))))
    @test ut_category_messages(ctx)["Models.\"Degradation Models\".Deg.\"Decomposition Temperature\""] ==
          "missing (required by Thermal Decomposition)"
end
```

Append to `ut_block_models.jl`:

```julia
@testset "thermal decomposition on the typed interface" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(2)
    PeriLab.Data_Manager.set_dof(2)
    nn = PeriLab.Data_Manager.create_constant_node_scalar_field("Number of Neighbors", Int64)
    nn .= 1
    nlist = PeriLab.Data_Manager.create_constant_bond_scalar_state("Neighborhoodlist", Int64)
    nlist[1] = [2]
    nlist[2] = [1]
    PeriLab.Data_Manager.create_bond_scalar_state("Bond Damage", Float64; default_value = 1)
    PeriLab.Data_Manager.create_constant_node_scalar_field("Active", Bool;
                                                           default_value = true)
    _, temperature = PeriLab.Data_Manager.create_node_scalar_field("Temperature", Float64)
    temperature .= [10.0, 60.0]
    UT_MF.Degradation.init_fields()
    p = typed_model(:degradation, Dict("Degradation Model" => "Thermal Decomposition",
                                       "Decomposition Temperature" => 50);
                    name_key = "Degradation Model")
    PeriLab.Data_Manager.set_block_model("Degradation Model", 1, p)
    UT_MF.Degradation.init_model([1, 2], 1)
    UT_MF.Degradation.compute_model([1, 2], p, 1, 0.0, 1.0)
    @test PeriLab.Data_Manager.get_field("Active") == [true, false]
    bd = PeriLab.Data_Manager.get_bond_damage("NP1")
    @test bd[2][1] == 0.0 && bd[1][1] == 0.0
    @test length(methods(UT_MF.Degradation.Thermal_Decomposition.compute_model)) == 1
end
```

- [ ] **Step 2: Run to verify they fail**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Support/Parameters/Input/ut_category_models unit_tests/Models/ut_block_models`
Expected:
- `type PeriLabInput has no field degradations`;
- `model "Thermal Decomposition" not found` (from `typed_model`).

- [ ] **Step 3: Implement**

Each model file gets its struct after the `export` lines:

```julia
# thermal_decomposition.jl
using ......ParameterSpec: @params, register_degradation
@params struct ThermalDecompositionParams
    decomposition_temperature::Float64 = req("Decomposition Temperature";
                                             quantity = :temperature)
end
__init__() = register_degradation("Thermal Decomposition", ThermalDecompositionParams)

# BondBased_Corrosion.jl
using ......ParameterSpec: @params, register_degradation
@params struct BondbasedCorrosionParams
end
__init__() = register_degradation("Bond-based Corrosion", BondbasedCorrosionParams)
```

If `quantity = :temperature` is not a known quantity, drop it: `grep -n 'quantity' src/Support/Parameters/Spec/field_spec.jl`.

In the models:
- The dict `init_model(nodes, degradation_parameter::Dict, block)` becomes `init_model(nodes::AbstractVector{Int64}, p::<Struct>, block::Int64)`.
- The dict `compute_model(nodes, degradation_parameter::Dict, block, time, dt)` becomes `compute_model(nodes::AbstractVector{Int64}, p::<Struct>, block::Int64, time::Float64, dt::Float64)`.
- `degradation_parameter["Decomposition Temperature"]` becomes `p.decomposition_temperature`.
- Update the docstrings: `p` is the model parameters.

`Degradation_Factory.jl`:
- Replace `compute_model`, `init_model` and `fields_for_local_synchronization` with the Additive versions from Task 2, without `@invokelatest` and with `"Degradation Model"` as the slot name.
- Drop `create_module_specifics` from the import.

`degradation_template.jl`: add `DegradationTemplateParams`, unregistered, with the commented `__init__`. Use the typed signatures and update the docstrings, as in Task 2.

`input.jl`: add the field `degradations` after `additives`, filled by `parse_single_models(models, "Degradation Models", :degradation, "Degradation Model", ctx)`.

`Model_Factory.jl`:
- Add `"Degradation Model"` to `TYPED_CATEGORIES`.
- In `read_properties`, store the block model the same way as Additive (`:degradation_model`, `input.degradations`).

- [ ] **Step 4: Run to verify they pass**

Run the Step 2 command, plus `unit_tests/Support/Parameters/Input/ut_golden_decks` and `unit_tests/Models/Degradation/ut_Degradation_Factory`.
Expected: all pass.

- [ ] **Step 5: Full suite, then commit**

The suite includes `test_thermal_decomp`.

```bash
git add src test
git commit -m "Degradation models on typed structs"
```

---

### Task 4: Thermal parameter structs and `input.thermals`

**Files:**
- Modify: `src/Models/Thermal/Thermal_Factory.jl` (`ThermalBaseParams`, `__init__`)
- Modify: `thermal_flow.jl`, `thermal_expansion.jl`, `heat_transfer.jl`, `hetval.jl` (structs and registration only)
- Modify: `src/Support/Parameters/Input/input.jl` (`PeriLabInput.thermals`)
- Modify (decks; keys no code reads for the named model — one `Ruling:` per deck):
  - `test/fullscale_tests/test_Additive/additive_2d_heat.yaml`: remove `Type`, `Specific Heat Capacity` from the thermal model.
  - `test/fullscale_tests/test_Hetval/hetval.yaml`: remove `Type`, `Thermal Expansion Coefficient`.
  - `test/fullscale_tests/test_Multistep/multistep.yaml`: remove `Type`.
  - `test/fullscale_tests/test_Newmark/Newmark_thermal_expansion.yaml`, `test/fullscale_tests/test_linear_static_matrix_based_solver/linear_static_matrix_based_thermal_expansion.yaml` and `test/fullscale_tests/test_matrix_based_additive/matrix_based_only_thermal.yaml`: remove `Heat Transfer Coefficient`, `Thermal Conductivity Print Bed`, `Environmental Temperature`.
- Test: append to `ut_category_models.jl`

**Interfaces:**
- Consumes (Task 1): `parse_models`, the number-or-list field type.
- Produces:
  - `Thermal.ThermalBaseParams(thermal_conductivity::Union{Nothing,Float64})`, registered with `register_base!(:thermal, …)`.
  - `Thermal_Flow.ThermalFlowParams` with fields:
    - `type::ThermalFlowType`, where `@enum ThermalFlowType BondBased Correspondence`;
    - `print_bed_temperature::Union{Nothing,Float64}`;
    - `thermal_conductivity_print_bed::Union{Nothing,Float64}`;
    - `print_bed_z_coordinate::Float64`.
  - `Thermal_Expansion.ThermalExpansionParams` with fields:
    - `thermal_expansion_coefficient::Union{Float64,Vector{Float64}}`;
    - `reference_temperature::Union{Nothing,Float64}`.
  - `Heat_Transfer.HeatTransferParams` with fields:
    - `heat_transfer_coefficient::Float64`;
    - `environmental_temperature::Union{Float64,String}`;
    - `allow_surface_change::Bool`.
  - `HETVAL.HETVALParams` with fields:
    - `file::String`;
    - `number_of_state_variables::Int64`;
    - `hetval_material_name::Union{Nothing,String}`;
    - `hetval_name::String`;
    - `predefined_field_names::Union{Nothing,String}`;
    - plus `key_patterns` `Property_N`.
  - `PeriLabInput.thermals::Dict{String,Any}` (`WithBase`, composites allowed), placed after `degradations`.

- [ ] **Step 1: Write the failing tests**

Append to `ut_category_models.jl`:

```julia
function ut_thermal(entry)
    return CMID.read_input(ut_category_deck(Dict("Thermal Models" => Dict("Th" => Dict{String,Any}(entry)))))
end
const TH_PATH = "Models.\"Thermal Models\".Th"

@testset "thermal models are parsed into typed structs" begin
    input, ctx = ut_thermal(Dict("Thermal Model" => "Thermal Flow + Heat Transfer",
                                 "Thermal Conductivity" => 2.0, "Type" => "Correspondence",
                                 "Heat Transfer Coefficient" => 1.0,
                                 "Environmental Temperature" => "20+t"))
    @test isempty(ctx.errors)
    th = input.thermals["Th"]
    @test th.base.thermal_conductivity == 2.0
    flow, transfer = th.model.parts
    @test string(flow.type) == "Correspondence"
    @test transfer.environmental_temperature == "20+t"
    @test transfer.allow_surface_change
    input, _ = ut_thermal(Dict("Thermal Model" => "Thermal Flow+Heat Transfer",
                               "Heat Transfer Coefficient" => 1.0,
                               "Environmental Temperature" => 30))
    flow, transfer = input.thermals["Th"].model.parts
    @test string(flow.type) == "BondBased"                     # default
    @test transfer.environmental_temperature === 30.0
    @test input.thermals["Th"].base.thermal_conductivity === nothing
    input, _ = ut_thermal(Dict("Thermal Model" => "Thermal Expansion",
                               "Thermal Expansion Coefficient" => [1.0, 2.0]))
    @test input.thermals["Th"].model.thermal_expansion_coefficient == [1.0, 2.0]
    @test input.thermals["Th"].model.reference_temperature === nothing
    input, _ = ut_thermal(Dict("Thermal Model" => "HETVAL", "File" => "h.so",
                               "Property_1" => 2.0))
    @test input.thermals["Th"].model.hetval_name == "HETVAL"
    @test input.thermals["Th"].model.number_of_state_variables == 1
end

@testset "thermal input errors" begin
    _, ctx = ut_thermal(Dict("Thermal Model" => "Thermal Flow", "Type" => "Bond"))
    @test haskey(ut_category_messages(ctx), "$TH_PATH.Type")
    _, ctx = ut_thermal(Dict("Thermal Model" => "Thermal Flow", "Thermal Conductivity" => 1.0,
                             "Print Bed Temperature" => 300.0))
    @test ut_category_messages(ctx)["$TH_PATH.\"Thermal Conductivity Print Bed\""] ==
          "required when Print Bed Temperature is given"
    _, ctx = ut_thermal(Dict("Thermal Model" => "HETVAL", "File" => "h.so",
                             "HETVAL Material Name" => repeat("x", 81)))
    @test ut_category_messages(ctx)["$TH_PATH.\"HETVAL Material Name\""] ==
          "at most 80 characters (Fortran)"
    _, ctx = ut_thermal(Dict("Thermal Model" => "Heat Transfer",
                             "Heat Transfer Coefficient" => 1.0,
                             "Environmental Temperature" => 30, "Type" => "Bond based"))
    @test startswith(ut_category_messages(ctx)["$TH_PATH.Type"], "unknown key")
end
```

- [ ] **Step 2: Run to verify they fail**

Run the input tests command.
Expected: `type PeriLabInput has no field thermals`, plus errors that the thermal models are not found.

- [ ] **Step 3: Implement**

`Thermal_Factory.jl`, after the `using` lines and before `global module_list`:

```julia
using .....ParameterSpec: @params, register_base!

"""
    ThermalBaseParams

Keys every thermal model may use: `Thermal Conductivity` (read by the critical
time step for any thermal model, and by Thermal Flow).
"""
@params struct ThermalBaseParams
    thermal_conductivity::Union{Nothing,Float64} = opt("Thermal Conductivity";
                                                       default = nothing, min = 0)
end

# registration runs at load time, never during precompilation
__init__() = register_base!(:thermal, ThermalBaseParams)
```

`thermal_flow.jl`, after its `export` lines:

```julia
using ......ParameterSpec: @params, register_thermal, ParseContext, add_error!, join_path
import ......ParameterSpec: check!

@enum ThermalFlowType BondBased Correspondence

@params struct ThermalFlowParams
    type::ThermalFlowType = opt("Type"; default = BondBased,
                                description = "Bond based or Correspondence")
    print_bed_temperature::Union{Nothing,Float64} = opt("Print Bed Temperature";
                                                        default = nothing)
    thermal_conductivity_print_bed::Union{Nothing,Float64} = opt("Thermal Conductivity Print Bed";
                                                                 default = nothing,
                                                                 min = 0)
    print_bed_z_coordinate::Float64 = opt("Print Bed Z Coordinate"; default = 0.0)
end

function check!(p::ThermalFlowParams, path::String, ctx::ParseContext)
    if p.print_bed_temperature !== nothing && p.thermal_conductivity_print_bed === nothing
        add_error!(ctx, join_path(path, "Thermal Conductivity Print Bed"),
                   "required when Print Bed Temperature is given")
    end
    return nothing
end

__init__() = register_thermal("Thermal Flow", ThermalFlowParams)
```

`thermal_expansion.jl`:

```julia
using ......ParameterSpec: @params, register_thermal
@params struct ThermalExpansionParams
    thermal_expansion_coefficient::Union{Float64,Vector{Float64}} = req("Thermal Expansion Coefficient";
                                                                        description = "a number, or one value per direction")
    reference_temperature::Union{Nothing,Float64} = opt("Reference Temperature";
                                                        default = nothing,
                                                        description = "0 if missing (with a warning)")
end
__init__() = register_thermal("Thermal Expansion", ThermalExpansionParams)
```

`heat_transfer.jl`:

```julia
using ......ParameterSpec: @params, register_thermal
@params struct HeatTransferParams
    heat_transfer_coefficient::Float64 = req("Heat Transfer Coefficient"; min = 0)
    environmental_temperature::Union{Float64,String} = req("Environmental Temperature";
                                                           description = "a number, or an expression in the time t")
    allow_surface_change::Bool = opt("Allow Surface Change"; default = true)
end
__init__() = register_thermal("Heat Transfer", HeatTransferParams)
```

`hetval.jl`:

```julia
using ......ParameterSpec: @params, register_thermal, ParseContext, add_error!, join_path
import ......ParameterSpec: check!, key_patterns
@params struct HETVALParams
    file::String = req("File"; description = "HETVAL library, relative to the input deck")
    number_of_state_variables::Int64 = opt("Number of State Variables"; default = 1, min = 1)
    hetval_material_name::Union{Nothing,String} = opt("HETVAL Material Name";
                                                      default = nothing)
    hetval_name::String = opt("HETVAL name"; default = "HETVAL")
    predefined_field_names::Union{Nothing,String} = opt("Predefined Field Names";
                                                        default = nothing,
                                                        description = "space separated node field names")
end
key_patterns(::Type{HETVALParams}) = [r"^Property_\d+$" => Float64]
function check!(p::HETVALParams, path::String, ctx::ParseContext)
    if p.hetval_material_name !== nothing && length(p.hetval_material_name) > 80
        add_error!(ctx, join_path(path, "HETVAL Material Name"),
                   "at most 80 characters (Fortran)")
    end
    return nothing
end
__init__() = register_thermal("HETVAL", HETVALParams)
```

The model modules sit one level below the factory. The factory uses `.....Data_Manager`, the models use `.....Data_Manager` too (re-exported), so check one existing `ParameterSpec` import at the same depth before choosing the dots. `Damage.Critical_Stretch` uses `......ParameterSpec`; use the same.

`input.jl`: add the field `thermals` after `degradations`, filled by `parse_models(models, "Thermal Models", :thermal, "Thermal Model", ctx)`.

Decks: remove the keys listed under Files. Each of them is ignored by the named models today (see Global Constraints), so the results do not change. Ledger one `Ruling:` per deck.

- [ ] **Step 4: Run to verify they pass**

Run the input tests command.
Expected:
- the new tests pass;
- `ut_golden_decks` passes, so every thermal deck parses.

If a golden deck fails on a thermal key that a thermal model reads, add the key to that model's struct and ledger a `Ruling:`.

- [ ] **Step 5: Full suite, then commit**

```bash
git add src test
git commit -m "Thermal models are validated into typed structs (input.thermals)"
```

---

### Task 5: Thermal dispatch, models, heat capacity and critical time step on the typed values

**Files:**
- Modify: `src/Models/Thermal/Thermal_Factory.jl` (`init_model`, `compute_model`, `fields_for_local_synchronization`)
- Modify: `thermal_flow.jl`, `thermal_expansion.jl`, `heat_transfer.jl`, `hetval.jl`, `Thermal_template/thermal_template.jl`
- Modify: `src/Models/Model_Factory.jl`:
  - `TYPED_CATEGORIES`;
  - `read_properties`;
  - `set_heat_capacity(input, block_nodes, heat_capacity)`;
  - thermal part of `compute_crititical_time_step`.
- Test: `test/unit_tests/Models/Thermal/ut_Thermal_Factory.jl`, `ut_Thermal_flow.jl`, `ut_Thermal_expansion.jl`, `ut_Heat_transfer.jl`, `ut_HETVAL.jl`; append to `ut_block_models.jl`

**Interfaces:**
- Consumes (Task 4): the thermal structs and `input.thermals`.
- Produces:
  - Thermal model interface:
    - `init_model(nodes::AbstractVector{Int64}, p::<Params>, thermal, block::Int64)`;
    - `compute_model(nodes::AbstractVector{Int64}, p::<Params>, thermal, block::Int64, time::Float64, dt::Float64)`;
    - `fields_for_local_synchronization(model::String)`.

    `thermal` is the block's `WithBase`; `thermal.base.thermal_conductivity` is the shared key.
  - `Thermal.compute_model(nodes, thermal::WithBase, block, time, dt)` runs every part of `thermal.model`, in deck order.
  - `Model_Factory.set_heat_capacity(input::PeriLabInput, block_nodes, heat_capacity)`. It aborts with "Specific Heat Capacity of <block name> is not defined" (the legacy message).

- [ ] **Step 1: Write the failing tests**

Append to `ut_block_models.jl`:

```julia
function ut_thermal_model(raw)
    return typed_model(:thermal, raw; name_key = "Thermal Model")
end

@testset "thermal composite runs every part" begin
    PeriLab.Data_Manager.initialize_data()
    th = ut_thermal_model(Dict("Thermal Model" => "Thermal Flow + Heat Transfer",
                               "Thermal Conductivity" => 1.0,
                               "Heat Transfer Coefficient" => 1.0,
                               "Environmental Temperature" => 30))
    parts = UT_MF.Thermal.model_parts(th.model)
    @test nameof.(typeof.(parts)) == (:ThermalFlowParams, :HeatTransferParams)
    th2 = ut_thermal_model(Dict("Thermal Model" => "Heat Transfer+Thermal Flow",
                                "Thermal Conductivity" => 1.0,
                                "Heat Transfer Coefficient" => 1.0,
                                "Environmental Temperature" => 30))
    @test nameof.(typeof.(UT_MF.Thermal.model_parts(th2.model))) ==
          (:HeatTransferParams, :ThermalFlowParams)
    single = ut_thermal_model(Dict("Thermal Model" => "Thermal Expansion",
                                   "Thermal Expansion Coefficient" => 1.0))
    @test UT_MF.Thermal.model_parts(single.model) == (single.model,)
end

@testset "environmental temperature expression" begin
    HT = UT_MF.Thermal.Heat_Transfer
    @test HT.environmental_temperature(30.0, 5.0) == 30.0
    @test HT.environmental_temperature("20+t", 5.0) == 25.0
end

@testset "thermal critical time step" begin
    PeriLab.Data_Manager.initialize_data()
    th = ut_thermal_model(Dict("Thermal Model" => "Heat Transfer",
                               "Thermal Conductivity" => 0.12,
                               "Heat Transfer Coefficient" => 1.0,
                               "Environmental Temperature" => 30))
    PeriLab.Data_Manager.set_block_model("Thermal Model", 1, th)
    @test UT_MF.block_thermal_conductivity(1) == 0.12
    @test UT_MF.block_thermal_conductivity(2) === nothing
end

@testset "heat capacity from the typed blocks" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(3)
    input = typed_input(Dict("Blocks" => Dict("block_1" => Dict("Block ID" => 1,
                                                                "Density" => 1.0,
                                                                "Horizon" => 1.0,
                                                                "Specific Heat Capacity" => 5.0),
                                              "block_2" => Dict("Block ID" => 2,
                                                                "Density" => 1.0,
                                                                "Horizon" => 1.0))))
    hc = zeros(3)
    UT_MF.set_heat_capacity(input, Dict(1 => [1, 2]), hc)
    @test hc == [5.0, 5.0, 0.0]
    @test_logs (:error,
                "Specific Heat Capacity of block_2 is not defined") @test_throws PeriLab.PeriLabError UT_MF.set_heat_capacity(input,
                                                                                                                                Dict(2 => [3]),
                                                                                                                                hc)
end
```

Rewrite the dict calls in the thermal unit tests. Every `Dict(...)` parameter becomes a typed model built with `typed_model(:thermal, Dict("Thermal Model" => "<model>", <same keys>); name_key = "Thermal Model")`. The call passes `(nodes, th.model, th, 1)` for `init_model` and `(nodes, th.model, th, 1, time, dt)` for `compute_model`. In detail:

- `ut_Thermal_Factory.jl`, `@testset "init_model"`: replace the two `set_properties` calls and the "Missing" abort with:

```julia
    th = typed_model(:thermal, Dict("Thermal Model" => "Heat Transfer",
                                    "Heat Transfer Coefficient" => 1.0,
                                    "Environmental Temperature" => 30);
                     name_key = "Thermal Model")
    PeriLab.Data_Manager.set_block_model("Thermal Model", 1, th)
    PeriLab.Solver_Manager.Model_Factory.Thermal.init_model([1], 1)
    @test PeriLab.Data_Manager.has_key("Bond Norm")
```

  An unknown thermal model name is a parse error now ("thermal input errors" covers the input side).

- `ut_Thermal_flow.jl`, `@testset "ut_init_model"`: replace the body from the first `init_model` call to the end with:

```julia
    TF = PeriLab.Solver_Manager.Model_Factory.Thermal.Thermal_Flow
    flow(raw) = typed_model(:thermal, merge(Dict{String,Any}("Thermal Model" => "Thermal Flow"), raw);
                            name_key = "Thermal Model")
    for type in ("Bond based", "Correspondence")
        th = flow(Dict("Type" => type, "Thermal Conductivity" => 100))
        TF.init_model(Vector{Int64}(1:3), th.model, th, 1)
        th = flow(Dict("Type" => type))
        @test_logs (:error,
                    "Thermal Conductivity not defined.") @test_throws PeriLab.PeriLabError TF.init_model(Vector{Int64}(1:3),
                                                                                                         th.model,
                                                                                                         th,
                                                                                                         1)
    end
    coordinates = PeriLab.Data_Manager.create_constant_node_vector_field("Coordinates",
                                                                         Float64, 3)
    coordinates .= [0 0 1; 1 0 1; 1 1 1]
    th = flow(Dict("Type" => "Correspondence", "Thermal Conductivity" => 100,
                   "Print Bed Temperature" => 2, "Thermal Conductivity Print Bed" => 1,
                   "Print Bed Z Coordinate" => -1))
    PeriLab.Data_Manager.set_dof(3)
    TF.init_model(Vector{Int64}(1:3), th.model, th, 1)
    @test TF.print_bed_active(th.model, 3)
    PeriLab.Data_Manager.set_dof(2)
    TF.init_model(Vector{Int64}(1:3), th.model, th, 1)    # warns: warnings are off in the suite
    @test !TF.print_bed_active(th.model, 2)
```

- `ut_Thermal_expansion.jl`: replace the `thermal_parameter = Dict(...)` and the `compute_model` call with:

```julia
    th = typed_model(:thermal, Dict("Thermal Model" => "Thermal Expansion",
                                    "Thermal Expansion Coefficient" => 1.0,
                                    "Reference Temperature" => 0.0);
                     name_key = "Thermal Model")

    temperature_NP1 .= 0
    PeriLab.Solver_Manager.Model_Factory.Thermal.Thermal_Expansion.compute_model(nodes,
                                                                                 th.model,
                                                                                 th, 1,
                                                                                 1.0, 1.0)
```

- `ut_Heat_transfer.jl`, `@testset "ut_compute_thermal_model"`: both `compute_model` calls pass a typed model instead of the dict:

```julia
    th = typed_model(:thermal, Dict("Thermal Model" => "Heat Transfer",
                                    "Heat Transfer Coefficient" => 1,
                                    "Environmental Temperature" => 1.2,
                                    "Allow Surface Change" => false);
                     name_key = "Thermal Model")
    PeriLab.Solver_Manager.Model_Factory.Thermal.Heat_Transfer.compute_model(Vector{Int64}(1:10),
                                                                             th.model, th, 1,
                                                                             1.0, 1.0)
```

- `ut_HETVAL.jl`, `@testset "init exceptions"`:
  - Remove the "HETVAL file is not defined." block (missing `File` is an input error now).
  - Build each remaining dict as `typed_model(:thermal, merge(Dict{String,Any}("Thermal Model" => "HETVAL"), <dict without "Number of Properties">); name_key = "Thermal Model")`.
  - Call `HETVAL.init_model(Vector{Int64}(1:nodes), th.model, th, 1)`.
  - Replace the assertions on the mutated dict (`mat_dict["HETVAL name"] == "test_sub"`, `!haskey`/`haskey(mat_dict, "HETVAL name")`, `== "HETVAL"`) with `@test th.model.hetval_name == "test_sub"` for the first model and `@test th.model.hetval_name == "HETVAL"` for the second.

- [ ] **Step 2: Run to verify they fail**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Models/ut_block_models unit_tests/Models/Thermal/ut_Thermal_Factory unit_tests/Models/Thermal/ut_Thermal_flow unit_tests/Models/Thermal/ut_Thermal_expansion unit_tests/Models/Thermal/ut_Heat_transfer unit_tests/Models/Thermal/ut_HETVAL`
Expected:
- `UndefVarError`s: `model_parts`, `environmental_temperature`, `block_thermal_conductivity`, `print_bed_active`;
- `MethodError`s for the typed model functions and for `set_heat_capacity(::PeriLabInput, …)`.

- [ ] **Step 3: Implement**

**Thermal_Factory.jl.**
- Add `WithBase, Composite` to the `ParameterSpec` import.
- Drop `create_module_specifics`.
- Replace `compute_model`, `init_model` and `fields_for_local_synchronization` (docstrings updated) with:

```julia
model_parts(model::Composite) = model.parts
model_parts(model) = (model,)

function compute_model(nodes::AbstractVector{Int64}, thermal::WithBase, block::Int64,
                       time::Float64, dt::Float64)
    for part in model_parts(thermal.model)
        @timeit "$(nameof(parentmodule(typeof(part))))" parentmodule(typeof(part)).compute_model(nodes,
                                                                                                  part,
                                                                                                  thermal,
                                                                                                  block,
                                                                                                  time,
                                                                                                  dt)
    end
end

function init_model(nodes::AbstractVector{Int64}, block::Int64)
    thermal = Data_Manager.get_block_model("Thermal Model", block)
    for part in model_parts(thermal.model)
        parentmodule(typeof(part)).init_model(nodes, part, thermal, block)
    end
end

function fields_for_local_synchronization(model, block)
    thermal = Data_Manager.get_block_model("Thermal Model", block)
    for part in model_parts(thermal.model)
        parentmodule(typeof(part)).fields_for_local_synchronization(model)
    end
end
```

**The models.** Replace the dict methods with typed ones. The bodies stay the same apart from these reads:

| File | Old | New |
|---|---|---|
| thermal_flow init | the `Type` check with warning and default | delete (validated; default `BondBased`) |
| thermal_flow init | `haskey(thermal_parameter, "Print Bed Temperature")` / `delete!` | `if p.print_bed_temperature !== nothing` … in 2D only `@warn "Print bed temperature can only be defined for 3D problems. Its deactivated."` (no mutation) |
| thermal_flow init | `get(thermal_parameter, "Print Bed Z Coordinate", 0.0)` | `p.print_bed_z_coordinate` |
| thermal_flow init | `!haskey(thermal_parameter, "Thermal Conductivity")` | `thermal.base.thermal_conductivity === nothing` |
| thermal_flow compute | `thermal_parameter["Thermal Conductivity"]` | `thermal.base.thermal_conductivity` |
| thermal_flow compute | `haskey(…, "Print Bed Temperature")` | `print_bed_active(p, dof)` |
| thermal_flow compute | `thermal_parameter["Print Bed Temperature"]` / `["Thermal Conductivity Print Bed"]` | `p.print_bed_temperature` / `p.thermal_conductivity_print_bed` |
| thermal_flow compute | `thermal_parameter["Type"] == "Bond based"` / `"Correspondence"` | `p.type == BondBased` / `p.type == Correspondence` |
| thermal_expansion init | `!haskey(thermal_parameter, "Reference Temperature")` | `p.reference_temperature === nothing` |
| thermal_expansion compute | `thermal_parameter["Thermal Expansion Coefficient"]` | `p.thermal_expansion_coefficient` |
| thermal_expansion compute | `get(thermal_parameter, "Reference Temperature", 0.0)` | `something(p.reference_temperature, 0.0)` |
| heat_transfer compute | `thermal_parameter["Heat Transfer Coefficient"]` | `p.heat_transfer_coefficient` |
| heat_transfer compute | the `isa String` / `eval` block | `Tenv::Float64 = environmental_temperature(p.environmental_temperature, time)` |
| heat_transfer compute | `get(thermal_parameter, "Allow Surface Change", true)` | `p.allow_surface_change` |

Add to `thermal_flow.jl`:

```julia
# the print bed only exists in 3D
print_bed_active(p::ThermalFlowParams, dof::Int64) = p.print_bed_temperature !== nothing &&
                                                     dof == 3
```

Add to `heat_transfer.jl`, with the legacy evaluation of an expression in the time `t`:

```julia
environmental_temperature(T::Float64, time::Float64) = T
function environmental_temperature(expression::String, time::Float64)
    global t = time
    return Float64(eval(Meta.parse(expression)))
end
```

`hetval.jl` `init_model(nodes, p::HETVALParams, thermal, block)`:
- Drop the missing-file abort; keep the file-path join and the existence check with `p.file`.
- `num_state_vars = p.number_of_state_variables`.
- The material name:

```julia
    if p.hetval_material_name === nothing
        @warn "No HETVAL Material Name is defined. Please check if you use it as method to check different material in your HETVAL."
    end
    global hetval_cmname = malloc_cstring(something(p.hetval_material_name, ""))
```

  Legacy set `hetval_cmname` only when the name was missing, which left it undefined when a name was given. Ledger this as a `Ruling:` (bug fix).
- `field_names = p.predefined_field_names === nothing ? ["Volume"] : split(p.predefined_field_names, " ")`.
- The 80-character abort and the `HETVAL name` default are gone: they are validated.

`thermal_template.jl`:
- Add `ThermalTemplateParams`, unregistered, with the commented `__init__ = register_thermal(...)`.
- Use the typed signatures and update the docstrings, as in Task 2.

Afterwards run `grep -rn '::Dict\|thermal_parameter' src/Models/Thermal`. Expected: no matches.

**Model_Factory.jl.**
- Add `"Thermal Model"` to `TYPED_CATEGORIES`.
- In `read_properties`, store the block model as in Task 2 (`:thermal_model`, `input.thermals`).
- Add:

```julia
# Thermal Conductivity of a block's thermal model, or `nothing`
function block_thermal_conductivity(block::Int64)
    thermal = Data_Manager.get_block_model("Thermal Model", block)
    return thermal === nothing ? nothing : thermal.base.thermal_conductivity
end
```

- In `compute_crititical_time_step`, `lambda = block_thermal_conductivity(iblock)`.
- `set_heat_capacity` becomes:

```julia
function set_heat_capacity(input::PeriLabInput, block_nodes::Dict,
                           heat_capacity::NodeScalarField{Float64})
    for block in eachindex(block_nodes)
        name, params = InputDeck.block_by_id(input.sections.blocks, block)
        if params.specific_heat_capacity === nothing
            @abort "Specific Heat Capacity of $name is not defined"
        end
        heat_capacity[block_nodes[block]] .= params.specific_heat_capacity
    end
    return heat_capacity
end
```

  The call in `init_models` passes `input`. Remove `get_heat_capacity` from the `Parameter_Handling` import. Check how `Model_Factory` reaches `InputDeck` (`using ...InputDeck: PeriLabInput, ContactInput`) and add `block_by_id` to that import if needed.

- [ ] **Step 4: Run to verify they pass**

Run the Step 2 command, plus `unit_tests/Support/Parameters/Input/ut_category_models`.
Expected: all pass.

- [ ] **Step 5: Full suite, then commit**

The suite includes `test_thermal_flow`, `test_thermal_expansion`, `test_Hetval` and `test_Multistep`.

```bash
git add src test
git commit -m "Thermal dispatch, models, heat capacity and critical time step read the typed thermal model"
```

---

### Task 6: Pre-calculation switches validated (`input.pre_calculation_global`, `input.pre_calculations`)

**Files:**
- Modify: `src/Models/Pre_calculation/Pre_Calculation_Factory.jl` (registration helper), and each pre-calculation module plus the template (empty `@params` struct, registered)
- Modify: `src/Support/Parameters/Input/input.jl`
- Test: append to `ut_category_models.jl`

**Interfaces:**
- Produces:
  - Every pre-calculation module registers an empty `@params` struct under its `pre_calculation_name()`. The template stays unregistered. Registered names:
    - "Deformed Bond Geometry", "Shape Tensor", "Deformation Gradient";
    - "Bond Associated Correspondence", "Axis Symmetric".
  - `PeriLabInput.pre_calculation_global::Union{Nothing,Dict{String,Bool}}` and `PeriLabInput.pre_calculations::Dict{String,Dict{String,Bool}}`, placed after `thermals`.
  - `InputDeck.parse_pre_calculation_switches(raw, path, ctx)::Union{Nothing,Dict{String,Bool}}`.

- [ ] **Step 1: Write the failing tests**

Append to `ut_category_models.jl`:

```julia
@testset "pre-calculation switches" begin
    global_switches = Dict("Deformed Bond Geometry" => true, "Shape Tensor" => false,
                           "Bond Associated Deformation Gradient" => false)
    input, ctx = CMID.read_input(ut_category_deck(Dict("Pre Calculation Global" => global_switches,
                                                       "Pre Calculation Models" => Dict("Pre" => Dict("Deformation Gradient" => true)))))
    @test isempty(ctx.errors)
    @test input.pre_calculation_global ==
          Dict("Deformed Bond Geometry" => true, "Shape Tensor" => false)
    @test input.pre_calculations["Pre"] == Dict("Deformation Gradient" => true)
    _, ctx = CMID.read_input(ut_category_deck(Dict("Pre Calculation Global" => Dict("Bond Associated Deformation Gradient" => true))))
    @test ut_category_messages(ctx)["Models.\"Pre Calculation Global\".\"Bond Associated Deformation Gradient\""] ==
          "no longer supported; use \"Bond Associated Correspondence\""
    _, ctx = CMID.read_input(ut_category_deck(Dict("Pre Calculation Global" => Dict("Shape Tenser" => true))))
    @test ut_category_messages(ctx)["Models.\"Pre Calculation Global\".\"Shape Tenser\""] ==
          "unknown key — did you mean \"Shape Tensor\"?"
    input, _ = CMID.read_input(ut_category_deck(Dict{String,Any}()))
    @test input.pre_calculation_global === nothing
    @test isempty(input.pre_calculations)
end
```

- [ ] **Step 2: Run to verify they fail**

Run the input tests command.
Expected: `type PeriLabInput has no field pre_calculation_global`.

- [ ] **Step 3: Implement**

In each pre-calculation module, after its `export` lines, add the following. Use one more dot than the module's own `Data_Manager` import; the modules differ (`grep -n '^using .*Data_Manager' src/Models/Pre_calculation/*.jl`).

```julia
using .......ParameterSpec: @params, register_pre_calculation
"Switch only: this pre-calculation has no parameters."
@params struct ShapeTensorParams
end
__init__() = register_pre_calculation("Shape Tensor", ShapeTensorParams)
```

The struct names are `DeformedBondGeometryParams`, `ShapeTensorParams`, `DeformationGradientParams`, `BondAssociatedCorrespondenceParams` and `AxisSymmetricParams`. Each module registers its own `pre_calculation_name()` string. The template gets `PreCalculationTemplateParams` with a commented `__init__`.

`input.jl`:

```julia
# old switch names: accepted when off, an error naming the replacement when on
const DEPRECATED_PRE_CALCULATIONS = Dict("Bond Associated Deformation Gradient" => "Bond Associated Correspondence")

"""
    parse_pre_calculation_switches(raw, path, ctx)

Pre-calculation switches (`name: true/false`); names are the registered
pre-calculation models.
"""
function parse_pre_calculation_switches(raw, path::String, ctx::ParseContext)
    switches = ParameterSpec.convert_value(Dict{String,Bool}, raw, path, ctx)
    switches === ParameterSpec.FAILED && return nothing
    for (name, replacement) in DEPRECATED_PRE_CALCULATIONS
        haskey(switches, name) || continue
        switches[name] &&
            add_error!(ctx, join_path(path, name), "no longer supported; use \"$replacement\"")
        delete!(switches, name)
    end
    ParameterSpec.check_unknown!(switches,
                                 Set(ParameterSpec.registered_names(:pre_calculation)),
                                 path, ctx)
    return switches
end
```

- Add the fields `pre_calculation_global::Union{Nothing,Dict{String,Bool}}` and `pre_calculations::Dict{String,Dict{String,Bool}}` after `thermals`, and mention them in the docstring.
- In `read_input`:

```julia
    pre_global = models isa AbstractDict && haskey(models, "Pre Calculation Global") ?
                 parse_pre_calculation_switches(models["Pre Calculation Global"],
                                                join_path("Models", "Pre Calculation Global"),
                                                ctx) : nothing
    pre_calculations = Dict{String,Dict{String,Bool}}()
    raw_pre = models isa AbstractDict ? get(models, "Pre Calculation Models", nothing) : nothing
    if raw_pre isa AbstractDict
        for (name, entry) in raw_pre
            s = parse_pre_calculation_switches(entry,
                                               join_path(join_path("Models",
                                                                   "Pre Calculation Models"),
                                                         string(name)), ctx)
            s === nothing || (pre_calculations[string(name)] = s)
        end
    elseif raw_pre !== nothing
        add_error!(ctx, join_path("Models", "Pre Calculation Models"),
                   "expected a section of `key: value` entries, got $(ParameterSpec._describe(raw_pre))")
    end
```

  Pass both to the constructor.

If `check_unknown!` does not suggest across `known` the same way for switch dicts, compare its message with the test and adjust only the test's expected text to the real format.

- [ ] **Step 4: Run to verify they pass**

Run the input tests command.
Expected: all pass, including the golden decks (56 decks carry `Bond Associated Deformation Gradient: false`).

- [ ] **Step 5: Full suite, then commit**

```bash
git add src test
git commit -m "Pre-calculation switches are validated (input.pre_calculation_global, input.pre_calculations)"
```

---

### Task 7: Pre-calculation runtime on the typed switches

**Files:**
- Modify: `src/Models/Pre_calculation/Pre_Calculation_Factory.jl` (`init_model`, `compute_model`, `fields_for_local_synchronization`, `check_dependencies`, `active_pre_calculations`, `order_pre_calculations`)
- Modify: every pre-calculation module (`init_model(nodes, block)`, `compute(nodes, block)`) and the template
- Modify: `src/Models/Model_Factory.jl` (`TYPED_CATEGORIES`, `read_properties`)
- Test: `test/unit_tests/Models/Material/ut_block_material.jl` (`pre-calculation dependencies`); append to `ut_block_models.jl`

**Interfaces:**
- Consumes (Task 6): `input.pre_calculation_global`, `input.pre_calculations`.
- Produces:
  - Block slot `"Pre Calculation Model"`: a `Vector{String}` of active names, in run order. It is stored only when non-empty.
  - `Pre_Calculation.order_pre_calculations(names)::Vector{String}`:
    - names first in `Data_Manager.get_pre_calculation_order()`;
    - then the remaining names, sorted.
  - `Pre_Calculation.active_pre_calculations(switches::AbstractDict{String,Bool})::Vector{String}`.
  - `Pre_Calculation.compute_model(nodes, names::Vector{String}, block, time, dt)`.
  - Module interface: `init_model(nodes::AbstractVector{Int64}, block::Int64)`, `compute(nodes::AbstractVector{Int64}, block::Int64)`, `fields_for_local_synchronization(model::String)`.

- [ ] **Step 1: Write the failing tests**

In `ut_block_material.jl`, `@testset "pre-calculation dependencies"`:
- `ut_dependencies` returns `PeriLab.Data_Manager.get_block_model("Pre Calculation Model", 1)`, and drops the `init_properties` call.
- The assertions become:

```julia
    p = ut_dependencies(merge(corr, Dict{String,Any}("Bond Associated" => true)))
    @test p == ["Deformed Bond Geometry", "Bond Associated Correspondence"]
    p = ut_dependencies(corr)
    @test p == ["Deformed Bond Geometry", "Shape Tensor", "Deformation Gradient"]
    p = ut_dependencies(Dict{String,Any}("Material Model" => "PD Solid Elastic",
                                         "Bulk Modulus" => 1.0, "Shear Modulus" => 1.0))
    @test p == ["Deformed Bond Geometry"]
```

Append to `ut_block_models.jl`:

```julia
@testset "read_properties stores pre-calculations" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(2)
    PeriLab.Data_Manager.set_dof(2)
    PeriLab.Data_Manager.set_block_name_list(["block_1", "block_2", "block_3"])
    PeriLab.Data_Manager.set_block_id_list([1, 2, 3])
    blocks = Dict("block_1" => Dict("Block ID" => 1, "Density" => 1.0, "Horizon" => 1.0),
                  "block_2" => Dict("Block ID" => 2, "Density" => 1.0, "Horizon" => 1.0,
                                    "Pre Calculation Model" => "Pre"),
                  "block_3" => Dict("Block ID" => 3, "Density" => 1.0, "Horizon" => 1.0,
                                    "Pre Calculation Model" => "Off"))
    models = Dict("Pre Calculation Global" => Dict("Shape Tensor" => true,
                                                   "Deformed Bond Geometry" => true),
                  "Pre Calculation Models" => Dict("Pre" => Dict("Deformation Gradient" => true),
                                                   "Off" => Dict("Shape Tensor" => false)))
    input = typed_input(Dict("Blocks" => blocks, "Models" => models))
    UT_MF.read_properties(Dict{String,Any}("Blocks" => blocks, "Models" => input.models),
                          input, false)
    get(b) = PeriLab.Data_Manager.get_block_model("Pre Calculation Model", b)
    @test get(1) == ["Deformed Bond Geometry", "Shape Tensor"]
    @test get(2) == ["Deformation Gradient"]          # replaces the global switches
    @test get(3) === nothing                          # nothing active
    @test UT_MF.has_block_model(1, "Pre Calculation Model")
    @test !UT_MF.has_block_model(3, "Pre Calculation Model")
    @test UT_MF.Pre_Calculation.order_pre_calculations(["Axis Symmetric", "Shape Tensor",
                                                        "Deformed Bond Geometry"]) ==
          ["Deformed Bond Geometry", "Shape Tensor", "Axis Symmetric"]
end
```

- [ ] **Step 2: Run to verify they fail**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Models/ut_block_models unit_tests/Models/Material/ut_block_material`
Expected:
- the dependencies test fails: `get_block_model` returns `nothing`;
- `read_properties stores pre-calculations` fails on `get(1)`;
- `UndefVarError: order_pre_calculations`.

- [ ] **Step 3: Implement**

`Pre_Calculation_Factory.jl`:
- Remove `create_module_specifics` and `DataStructures`. Add `using .....ParameterSpec`.
- Replace `init_model`, `compute_model`, `fields_for_local_synchronization` and `check_dependencies` (docstrings updated) with:

```julia
pre_calculation_module(name::String) = parentmodule(ParameterSpec.lookup_model(:pre_calculation,
                                                                               name))

"Active names of `switches` (on/off per pre-calculation), in run order."
active_pre_calculations(switches::AbstractDict{String,Bool}) = order_pre_calculations([name
                                                                                       for (name, on) in switches
                                                                                       if on])

"`names` in run order: the fixed order first, then the others sorted."
function order_pre_calculations(names)
    order = Data_Manager.get_pre_calculation_order()
    return vcat([n for n in order if n in names], sort!([n for n in names if !(n in order)]))
end

function init_model(nodes::AbstractVector{Int64}, block::Int64)
    for name in Data_Manager.get_block_model("Pre Calculation Model", block)
        pre_calculation_module(name).init_model(nodes, block)
    end
end

function compute_model(nodes::AbstractVector{Int64}, names::Vector{String}, block::Int64,
                       time::Float64, dt::Float64)
    for name in names
        @timeit "compute $name" pre_calculation_module(name).compute(nodes, block)
    end
end

function fields_for_local_synchronization(model, block)
    for name in Data_Manager.get_block_model("Pre Calculation Model", block)
        pre_calculation_module(name).fields_for_local_synchronization(model)
    end
end

function check_dependencies(block_nodes::Dict{Int64,Vector{Int64}})
    for block_id in eachindex(block_nodes)
        material = Data_Manager.get_block_material(block_id)
        material === nothing && continue
        current = Data_Manager.get_block_model("Pre Calculation Model", block_id)
        names = Set{String}(current === nothing ? String[] : current)
        push!(names, "Deformed Bond Geometry")
        if material.correspondence
            if material.base.bond_associated
                push!(names, "Bond Associated Correspondence")
            else
                push!(names, "Shape Tensor", "Deformation Gradient")
            end
        end
        "Deformation Gradient" in names && push!(names, "Shape Tensor")
        Data_Manager.set_block_model("Pre Calculation Model", block_id,
                                     order_pre_calculations(collect(names)))
    end
end
```

Ledger two `Ruling:`s:
- Legacy `check_dependencies` set a bond-associated block's switches by merging the last block's switches (a loop over all blocks writing into `block_id`). The typed version keeps the block's own switches. The result only differs when blocks have different switch sets.
- Legacy dropped names that are not in the fixed order (`Axis Symmetric`) for material blocks. They now run after the ordered ones. No deck uses `Axis Symmetric`.

Modules:
- `init_model(nodes, parameter, block)` becomes `init_model(nodes::AbstractVector{Int64}, block::Int64)`.
- `compute(nodes, parameter, block)` becomes `compute(nodes::AbstractVector{Int64}, block::Int64)`.
- The bodies never read `parameter`.
- Update the docstrings: remove the `parameter` argument lines.
- `axissymmetric.jl`'s `init(nodes, Pre_calculation_parameter::Dict)` loses the dict argument if it is unused (`grep -n 'Pre_calculation_parameter' src/Models/Pre_calculation/axissymmetric.jl`).
- The template gets the same signatures.

Afterwards run `grep -rn 'parameter::\|::Dict\|OrderedDict' src/Models/Pre_calculation`. Expected: no matches.

`Model_Factory.jl`:
- Add `"Pre Calculation Model"` to `TYPED_CATEGORIES`.
- In `read_properties`, add:

```julia
    for (block_name, block) in zip(block_name_list, block_id_list)
        switches = block_typed_model(input, block_name, :pre_calculation_model,
                                     "Pre Calculation Model", input.pre_calculations)
        switches === nothing && (switches = input.pre_calculation_global)
        switches === nothing && continue
        names = Pre_Calculation.active_pre_calculations(switches)
        isempty(names) ||
            Data_Manager.set_block_model("Pre Calculation Model", block, names)
    end
```

- [ ] **Step 4: Run to verify they pass**

Run the Step 2 command, plus `unit_tests/Models/Pre_calculation/ut_pre_bond_associated_correspondence` and `unit_tests/Models/ut_Model_Factory`.
Expected: all pass.

- [ ] **Step 5: Full suite, then commit**

The pre-calculation, correspondence and bond-associated decks pin the numbers.

```bash
git add src test
git commit -m "Pre-calculations run from the typed per-block switch list"
```

---

### Task 8: No property dicts for any model category

**Files:**
- Modify: `src/Models/Model_Factory.jl`:
  - `read_properties` no longer calls `get_block_model_definition` or `init_properties`;
  - delete `get_block_model_definition`;
  - `has_block_model` / `block_model_parameters` become typed-only;
  - drop the `get_model_parameter` import.
- Test: `test/unit_tests/Models/ut_Model_Factory.jl` (delete `ut_get_block_model_definition`; adapt `ut_read_properties`)

**Interfaces:**
- Consumes: all `TYPED_CATEGORIES` (Tasks 2, 3, 5, 7).
- Produces:
  - `has_block_model(block, name) = typed_block_model(block, name) !== nothing`.
  - `block_model_parameters(block, name) = typed_block_model(block, name)`.
  - `TYPED_CATEGORIES` is removed: every category is typed.

- [ ] **Step 1: Write the failing test**

In `ut_Model_Factory.jl`:
- Delete `@testset "ut_get_block_model_definition"`; the function is removed.
- Replace the body of `@testset "ut_read_properties"` with:

```julia
@testset "ut_read_properties" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_block_name_list(["block_1"])
    PeriLab.Data_Manager.set_block_id_list([1])
    input = typed_input(Dict("Blocks" => Dict("block_1" => Dict("Block ID" => 1,
                                                                "Density" => 1.0,
                                                                "Horizon" => 1.0,
                                                                "Thermal Model" => "therm")),
                             "Models" => Dict("Thermal Models" => Dict("therm" => Dict("Thermal Model" => "Heat Transfer",
                                                                                       "Heat Transfer Coefficient" => 1.0,
                                                                                       "Environmental Temperature" => 30)))))
    params = Dict{String,Any}("Blocks" => Dict{String,Any}(), "Models" => input.models)
    PeriLab.Solver_Manager.Model_Factory.read_properties(params, input, false)
    @test PeriLab.Data_Manager.get_block_model("Thermal Model", 1) === input.thermals["therm"]
    @test !haskey(PeriLab.Data_Manager.data["properties"], 1)     # no property dicts
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Models/ut_Model_Factory`
Expected: `ut_read_properties` fails on the last assertion: `init_properties` still creates the dicts. Also, `params["Blocks"]` is empty, so the legacy `get_block_model_definition` is the only reader of it.

- [ ] **Step 3: Implement**

`Model_Factory.jl`:
- In `read_properties`, delete `prop_keys = Data_Manager.init_properties()`, the `directory = …` line and the `get_block_model_definition(...)` call.
- Delete `get_block_model_definition` and its docstring.
- `has_block_model` / `block_model_parameters` become the one-liners from the Interfaces block. Delete `TYPED_CATEGORIES` and its comment.
- Remove `get_model_parameter` from the `Parameter_Handling` import. Drop the import line if nothing is left.

Run `grep -rn 'get_properties\|get_property\|check_property\|set_properties\|init_properties' src --include=*.jl | grep -v 'src/Core/Data_manager.jl\|src/Support/Parameters/'`.
Expected:
- Only `Matrix_linear_static.jl`'s commented line.
- `Data_Manager`'s own definitions and `Parameter_Handling` stay for phase 4.

If `init_models` or a solver still calls `Data_Manager.get_properties` for any model, switch it to `block_model_parameters` and ledger a `Ruling:`.

- [ ] **Step 4: Run to verify it passes**

Run the Step 2 command, plus `unit_tests/Models/ut_block_models`, `unit_tests/Models/Material/ut_block_material` and `unit_tests/Models/Damage/ut_block_damage`.
Expected: all pass.

- [ ] **Step 5: Full suite, then commit**

```bash
git add src test
git commit -m "No property dicts for any model category; Model_Factory dispatch is typed-only"
```
