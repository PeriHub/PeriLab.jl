# Typed Input Parameters — Phase 3d: Damage Models

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** The damage category is typed end to end. `Damage Models` are validated into `@params` structs when the input is read. Every block gets a typed `BlockDamage`. The damage factory, the three damage models, the template and local damping read struct fields. The damage properties dict is no longer stored.

**Architecture:**
- **Validation (as Material in 3a).**
  - A shared base part `DamageBaseParams` is declared in the damage factory and registered with `register_base!(:damage, …)`. It holds `Critical Value` (a `Dependent`), `Interblock Damage`, `Anisotropic Damage` and `Local Damping`.
  - Each model registers its own struct for its own keys: `Only Tension` (the default differs per model) and `Thickness`.
  - `read_input` parses `Models."Damage Models"` into `input.damages`. Combining damage models with `+` is an input error.
- **Runtime (as Material in 3b/3c).**
  - `BlockDamage{B,M}(base, model, tables)` is stored per block in `Data_Manager` (`set_block_damage` / `get_block_damage`).
  - The factory dispatches on `parentmodule(typeof(damage.model))`. The model functions take `(nodes, p::<Model>Params, damage::BlockDamage, block, …)`.
  - Hot loops get the concrete `Critical Value` type through a function barrier.
  - `Dependent` tables are bound to their node field before every compute, by a binder shared with Material and moved to `Data_Manager`.
- **Removal.** `read_properties` stops storing the damage dict. All dict methods of the damage modules are deleted. The template moves to the typed interface.

**Tech Stack:** Julia 1.12, `PeriLab.ParameterSpec`, `PeriLab.InputDeck`, `PeriLab.Data_Manager`.

**Spec:** `docs/superpowers/specs/2026-10-03-typed-input-parameters-design.md` (§2.2 registry, §2.3 InputReader, §3 Use, §5 phase 3 "then Damage … templates and licensed modules of that category are updated in the same phase")

## Global Constraints

- **Numerical results unchanged.** The full suite stays green. It includes the damage decks:
  - `test_critical_stretch` (incl. interface and tension/aniso);
  - `test_critical_energy` (incl. interface);
  - `test_aniso_damage`;
  - `test_DCB` (incl. temperature-dependent Critical Value and local damping);
  - `test_Dogbone`, `test_Newmark`, `test_bond_value_export`, `test_3D_aniso_material`;
  - the matrix-based decks.
- **Golden decks.** Every deck under `test/` and `examples/` still parses in strict mode (`ut_golden_decks.jl`). No deck is added to the allowlist to hide a validator gap.
- **Other categories keep the dict path.** `Data_Manager.get_properties` / `check_property` stay for thermal, additive, degradation and pre-calculation (later phases).
- **Read fields directly; functions only for real logic.**
- **Per-category slots, not `BlockModels` yet.** Like `Block Materials`, the damage gets its own slot `Block Damages`. Folding the slots into the spec's `BlockModels` struct is phase 4 work, once every category is typed.
- **Licensed damage modules** are maintained in-house. They follow the typed interface of the template in Task 4; nothing in this repo references them.
- **`ut_MPI.jl`** runs standalone under mpiexec without `helper.jl`.
- **`/dev/shm` is 64 MB.** If MPI tests die with "Bus error (signal 7)", check that no MPI process runs, then delete the stale `/dev/shm/mpich_shm_*` files.
- **Git.** Branch `feature/typed-input-parameters`; one commit per task; never merge. Read the full-suite result before committing a task that runs it.

## Review Focus

1. **Interblock Damage** builds the same critical-value matrix as before:
   - the matrix is filled with `Critical Value`, then the `Interblock Critical Value i_j` entries are set for the block's slice;
   - block ids beyond the block count are ignored.

   Pinned in Task 2 (`typed interface critical values`).
2. **A field-dependent Critical Value** (data file over Temperature) is bound to the current field before every compute and gives the spline value. Pinned in Task 2 (`block damage binds its tables`) and by `test_DCB` (temperature-dependent deck).
3. **Anisotropic damage.**
   - `Critical Energy Anisotropic` without `Anisotropic Damage` aborts with a clear message, not a `KeyError` in the time loop.
   - In 3D, `Anisotropic Damage` without `Critical Value Z` uses `Critical Value Y`.

   Pinned in Task 2 (`typed anisotropic critical values`) and Task 3 (`anisotropic energy needs Anisotropic Damage`).
4. **A block that names an undefined damage model** aborts with the legacy message. So does a deck without a `Damage Models` section. Pinned in Task 2 (`read_properties builds block damages`).
5. **Local Damping.**
   - Missing keys are input errors at parse time.
   - A block without damage is skipped by the local-damping compute.

   Pinned in Task 1 (`local damping keys are required`) and Task 3 (`local damping reads the block damage`).

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
- Input tests: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Support/Parameters/Input/ut_damage_models unit_tests/Support/Parameters/Input/ut_golden_decks`
- Full suite (~30 min, background, read the tail): `julia --project=. -e 'using Pkg; Pkg.test()' > <workspace>/suite.log 2>&1`

---

### Task 1: Damage parameter structs and `input.damages`

**Files:**
- Modify: `src/Models/Damage/Damage_Factory.jl` (base structs, `check!`, `__init__`)
- Modify: `src/Models/Damage/Critical_Stretch.jl`, `Energy_release.jl`, `Energy_release_aniso.jl` (model structs + registration)
- Modify: `src/Support/Parameters/Input/input.jl` (`PeriLabInput.damages`, generic `parse_models`, `parse_damages`)
- Create: `test/unit_tests/Support/Parameters/Input/ut_damage_models.jl`
- Modify: `test/unit_tests/Support/Parameters/Input/input_tests.jl` (add the file to the list)

**Interfaces:**
- Produces:
  - `Damage.DamageBaseParams` with fields:
    - `critical_value::Dependent`;
    - `interblock_damage::Union{Nothing,Dict{String,Float64}}`;
    - `anisotropic_damage::Union{Nothing,AnisotropicDamageParams}`;
    - `local_damping::Union{Nothing,LocalDampingParams}`.
  - `Damage.AnisotropicDamageParams(critical_value_x::Float64, critical_value_y::Float64, critical_value_z::Union{Nothing,Float64})`.
  - `Damage.LocalDampingParams(representative_youngs_modulus::Float64, damping_coefficient::Float64)`.
  - Model structs:
    - `Critical_Stretch.CriticalStretchParams(only_tension::Bool = true)`;
    - `Critical_Energy.CriticalEnergyParams(only_tension::Bool = true, thickness::Float64 = 1.0)`;
    - `Critical_Energy_Aniso.CriticalEnergyAnisotropicParams(only_tension::Bool = false, thickness::Float64 = 1.0)`.

    Registered as `"Critical Stretch"`, `"Critical Energy"` and `"Critical Energy Anisotropic"` (the existing `damage_name()` strings). Check the module names with `grep -n '^module' src/Models/Damage/*.jl` and use the real ones.
  - `PeriLabInput.damages::Dict{String,Any}` (name ⇒ `WithBase`), placed directly after `materials`.
  - `InputDeck.parse_models(models, section, category, name_key, ctx)`. `parse_materials(models, ctx)` stays as a one-line call of it.

- [ ] **Step 1: Write the failing tests**

Create `test/unit_tests/Support/Parameters/Input/ut_damage_models.jl`:

```julia
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
```

Add `"ut_damage_models.jl"` to the file list in `input_tests.jl`, after `"ut_material_models.jl"`.

- [ ] **Step 2: Run to verify they fail**

Run the input tests command.
Expected:
- `ut_damage_models` fails with `type PeriLabInput has no field damages`;
- the error-message tests fail because no damage model is registered: `model "Critical Stretch" not found`.

- [ ] **Step 3: Implement**

`Damage_Factory.jl`: add after the `using` lines, before `global module_list`:

```julia
using .....ParameterSpec: @params, Dependent, Constant, register_base!, ParseContext,
                          add_error!, join_path
import .....ParameterSpec: check!

@params struct AnisotropicDamageParams
    critical_value_x::Float64 = req("Critical Value X")
    critical_value_y::Float64 = req("Critical Value Y")
    critical_value_z::Union{Nothing,Float64} = opt("Critical Value Z"; default = nothing,
                                                   description = "defaults to Critical Value Y")
end

@params struct LocalDampingParams
    representative_youngs_modulus::Float64 = req("Representative Young's modulus"; min = 0,
                                                 quantity = :stress)
    damping_coefficient::Float64 = req("Damping coefficient"; min = 0)
end

"""
    DamageBaseParams

Keys every damage model may use (read from the same YAML block as the model's
own keys).
"""
@params struct DamageBaseParams
    critical_value::Dependent = req("Critical Value"; min = 0)
    interblock_damage::Union{Nothing,Dict{String,Float64}} = opt("Interblock Damage";
                                                                 default = nothing,
                                                                 description = "Interblock Critical Value <block>_<block> entries")
    anisotropic_damage::Union{Nothing,AnisotropicDamageParams} = opt("Anisotropic Damage";
                                                                     default = nothing)
    local_damping::Union{Nothing,LocalDampingParams} = opt("Local Damping"; default = nothing)
end

const INTERBLOCK_KEY = r"^Interblock Critical Value \d+_\d+$"

function check!(p::DamageBaseParams, path::String, ctx::ParseContext)
    p.interblock_damage === nothing && return nothing
    for name in keys(p.interblock_damage)
        occursin(INTERBLOCK_KEY, name) ||
            add_error!(ctx, join_path(join_path(path, "Interblock Damage"), name),
                       "unknown key — expected \"Interblock Critical Value <block>_<block>\"")
    end
    p.critical_value isa Constant ||
        add_error!(ctx, join_path(path, "Critical Value"),
                   "must be a number when Interblock Damage is used")
    return nothing
end

# registration runs at load time, never during precompilation
__init__() = register_base!(:damage, DamageBaseParams)
```

Each model file adds its struct after its `using` lines. The model modules sit one level deeper than the factory, so they use `......ParameterSpec`: one more dot than their `.....Data_Manager`.

```julia
# Critical_Stretch.jl
using ......ParameterSpec: @params, register_damage
@params struct CriticalStretchParams
    only_tension::Bool = opt("Only Tension"; default = true)
end
__init__() = register_damage("Critical Stretch", CriticalStretchParams)

# Energy_release.jl
using ......ParameterSpec: @params, register_damage
@params struct CriticalEnergyParams
    only_tension::Bool = opt("Only Tension"; default = true)
    thickness::Float64 = opt("Thickness"; default = 1.0, min = 0, quantity = :length)
end
__init__() = register_damage("Critical Energy", CriticalEnergyParams)

# Energy_release_aniso.jl
using ......ParameterSpec: @params, register_damage
@params struct CriticalEnergyAnisotropicParams
    only_tension::Bool = opt("Only Tension"; default = false)
    thickness::Float64 = opt("Thickness"; default = 1.0, min = 0, quantity = :length)
end
__init__() = register_damage("Critical Energy Anisotropic", CriticalEnergyAnisotropicParams)
```

If a module already has an `__init__`, append the call to it instead of adding a second one.

`input.jl`:
- Add the field `damages::Dict{String,Any}` to `PeriLabInput` after `materials`, and update the docstring: "`damages` holds the typed damage models".
- Replace `parse_materials` with a generic reader and two one-liners:

```julia
"Typed models of `section` (e.g. \"Material Models\"), keyed by name; errors go to `ctx`."
function parse_models(models::AbstractDict, section::String, category::Symbol,
                      name_key::String, ctx::ParseContext)
    parsed = Dict{String,Any}()
    raw = get(models, section, nothing)
    raw === nothing && return parsed
    path = join_path("Models", section)
    if !(raw isa AbstractDict)
        add_error!(ctx, path,
                   "expected a section of `key: value` entries, got $(ParameterSpec._describe(raw))")
        return parsed
    end
    for (name, entry) in raw
        entry_path = join_path(path, string(name))
        if !(entry isa AbstractDict)
            add_error!(ctx, entry_path,
                       "expected a section of `key: value` entries, got $(ParameterSpec._describe(entry))")
            continue
        end
        model = ParameterSpec.parse_model(category,
                                          Dict{String,Any}(string(k) => v for (k, v) in entry),
                                          entry_path, ctx; name_key = name_key)
        model === nothing || (parsed[string(name)] = model)
    end
    return parsed
end

parse_materials(models::AbstractDict, ctx::ParseContext) = parse_models(models,
                                                                        "Material Models",
                                                                        :material,
                                                                        "Material Model", ctx)

"Typed `Damage Models`; combining damage models with `+` is an error."
function parse_damages(models::AbstractDict, ctx::ParseContext)
    damages = parse_models(models, "Damage Models", :damage, "Damage Model", ctx)
    for (name, d) in damages
        d.model isa ParameterSpec.Composite || continue
        add_error!(ctx,
                   join_path(join_path(join_path("Models", "Damage Models"), name),
                             "Damage Model"),
                   "damage models cannot be combined with +")
    end
    return damages
end
```

- In `read_input`:
  - next to `materials = …`, add `damages = models isa AbstractDict ? parse_damages(models, ctx) : Dict{String,Any}()`;
  - pass `damages` to `PeriLabInput(sections, contact, materials, damages, …)`.

- [ ] **Step 4: Run to verify they pass**

Run the input tests command.
Expected:
- `ut_damage_models` passes.
- `ut_golden_decks` passes: every shipped damage deck parses.

If a golden deck fails on a damage key that the legacy code read, add that key to the struct whose code reads it. A ledgered example is "Critical Value Z", which `init_aniso_crit_values` reads. If the key was never read by any damage code, the deck is wrong: report it with a `Ruling:`. Never allowlist a deck.

- [ ] **Step 5: Full suite, then commit**

Full suite in the background; read the tail. Expected: 0 failures, 0 errors.

```bash
git add src/Models/Damage src/Support/Parameters/Input/input.jl test/unit_tests/Support/Parameters/Input
git commit -m "Damage models are validated into typed structs (input.damages)"
```

---

### Task 2: `BlockDamage` runtime struct, shared table binder, typed critical-value setup

**Files:**
- Modify: `src/Core/Data_manager.jl`:
  - slot `Block Damages`;
  - `set_block_damage` / `get_block_damage`;
  - `bind_dependent_tables!` (moved from `Material_Factory`).
- Modify: `src/Core/Data_manager/data_manager_checkpoint.jl` (exclude `Block Damages`)
- Modify: `src/Support/Parameters/Spec/dependent.jl`: `dependent_tables`, moved from `Material_Factory._tables` and exported.
- Modify: `src/Models/Material/Material_Factory.jl` (use `dependent_tables` and `Data_Manager.bind_dependent_tables!`; delete `_tables` and `_dependent_field`)
- Modify: `src/Models/Damage/Damage_Factory.jl`:
  - `BlockDamage`, `block_damage`;
  - typed `init_interface_crit_values(damage::BlockDamage, block)`;
  - typed `init_aniso_crit_values(aniso::AnisotropicDamageParams, block, dof)`.
- Modify: `src/Models/Model_Factory.jl` (`read_properties` builds and stores block damages; abort on an undefined name)
- Test: create `test/unit_tests/Models/Damage/ut_block_damage.jl` and include it in `runtests.jl` after `ut_Energy_release.jl`. Modify `test/unit_tests/Models/ut_Model_Factory.jl`.

**Interfaces:**
- Consumes (Task 1): `input.damages`, `DamageBaseParams`, `AnisotropicDamageParams`.
- Produces:
  - `Damage.BlockDamage{B,M}` with fields `base::B`, `model::M`, `tables::Vector{Table1D}`.
  - `Damage.block_damage(wb::WithBase)::BlockDamage`.
  - `Data_Manager.set_block_damage(block::Int64, damage)` and `Data_Manager.get_block_damage(block::Int64)` (returns `nothing` if the block has none).
  - `Data_Manager.bind_dependent_tables!(tables::Vector{Table1D})`. It aborts with `Field "<f>" required by <source> does not exist or is not a per-node Vector{Float64}.` (the message is unchanged from Material).
  - `ParameterSpec.dependent_tables(x)::Vector{Table1D}`. `x` is an `@params` struct or a `Composite`.
  - `Damage.init_interface_crit_values(damage::BlockDamage, block_id::Int64)`.
  - `Damage.init_aniso_crit_values(aniso::AnisotropicDamageParams, block_id::Int64, dof::Int64)`.
  - The dict methods of both stay until Task 3.

- [ ] **Step 1: Write the failing tests**

Create `test/unit_tests/Models/Damage/ut_block_damage.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

const BDAM = PeriLab.Solver_Manager.Model_Factory.Damage

function ut_typed_damage(raw)
    ctx = PeriLab.ParameterSpec.ParseContext()
    wb = PeriLab.ParameterSpec.parse_model(:damage,
                                           Dict{String,Any}(string(k) => v for (k, v) in raw),
                                           "test", ctx; name_key = "Damage Model")
    PeriLab.ParameterSpec.report!(ctx)
    return BDAM.block_damage(wb)
end

@testset "block damage binds its tables" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(2)
    file = joinpath(mktempdir(), "crit.txt")
    write(file, "header: Temperature Critical_Value\n0.0 1.0\n100.0 2.0\n")
    d = ut_typed_damage(Dict("Damage Model" => "Critical Energy", "Critical Value" => file))
    @test d isa BDAM.BlockDamage
    @test length(d.tables) == 1
    @test_logs (:error,
                "Field \"Temperature\" required by $(d.tables[1].source) does not exist or is not a per-node Vector{Float64}.") @test_throws PeriLab.PeriLabError PeriLab.Data_Manager.bind_dependent_tables!(d.tables)
    _, temperature = PeriLab.Data_Manager.create_node_scalar_field("Temperature", Float64)
    temperature .= [0.0, 100.0]
    PeriLab.Data_Manager.bind_dependent_tables!(d.tables)
    @test PeriLab.ParameterSpec.value(d.base.critical_value, 1) ≈ 1.0
    @test PeriLab.ParameterSpec.value(d.base.critical_value, 2) ≈ 2.0
    @test isempty(ut_typed_damage(Dict("Damage Model" => "Critical Stretch",
                                       "Critical Value" => 0.1)).tables)
end

@testset "block damage slot" begin
    PeriLab.Data_Manager.initialize_data()
    @test PeriLab.Data_Manager.get_block_damage(1) === nothing
    d = ut_typed_damage(Dict("Damage Model" => "Critical Stretch", "Critical Value" => 0.1))
    PeriLab.Data_Manager.set_block_damage(1, d)
    @test PeriLab.Data_Manager.get_block_damage(1) === d
    @test "Block Damages" in PeriLab.Data_Manager.CHECKPOINT_EXCLUDED_KEYS
end

@testset "typed interface critical values" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_block_id_list([2, 3, 1])
    d = ut_typed_damage(Dict("Damage Model" => "Critical Stretch", "Critical Value" => 1.0,
                             "Interblock Damage" => Dict("Interblock Critical Value 1_2" => 0.2,
                                                         "Interblock Critical Value 2_3" => 0.3,
                                                         "Interblock Critical Value 2_1" => 0.4,
                                                         "Interblock Critical Value 4_1" => 0.9)))
    BDAM.init_interface_crit_values(d, 1)
    m = PeriLab.Data_Manager.get_crit_values_matrix()
    @test size(m) == (3, 3, 3)
    @test m[1, 2, 1] == 0.2 && m[2, 3, 1] == 0.3 && m[2, 1, 1] == 0.4
    @test m[1, 1, 1] == 1.0 && m[3, 3, 2] == 1.0          # filled with Critical Value
    plain = ut_typed_damage(Dict("Damage Model" => "Critical Stretch", "Critical Value" => 1.0))
    PeriLab.Data_Manager.initialize_data()
    BDAM.init_interface_crit_values(plain, 1)
    @test PeriLab.Data_Manager.get_crit_values_matrix() == fill(-1, (1, 1, 1))
end

@testset "typed anisotropic critical values" begin
    PeriLab.Data_Manager.initialize_data()
    d = ut_typed_damage(Dict("Damage Model" => "Critical Energy Anisotropic",
                             "Critical Value" => 1.0,
                             "Anisotropic Damage" => Dict("Critical Value X" => 1.0,
                                                          "Critical Value Y" => 2.0)))
    BDAM.init_aniso_crit_values(d.base.anisotropic_damage, 1, 2)
    BDAM.init_aniso_crit_values(d.base.anisotropic_damage, 2, 3)
    z = ut_typed_damage(Dict("Damage Model" => "Critical Energy Anisotropic",
                             "Critical Value" => 1.0,
                             "Anisotropic Damage" => Dict("Critical Value X" => 1.0,
                                                          "Critical Value Y" => 2.0,
                                                          "Critical Value Z" => 3.0)))
    BDAM.init_aniso_crit_values(z.base.anisotropic_damage, 3, 3)
    aniso = PeriLab.Data_Manager.get_aniso_crit_values()
    @test aniso[1] == [1.0, 2.0]
    @test aniso[2] == [1.0, 2.0, 2.0]          # Z defaults to Y
    @test aniso[3] == [1.0, 2.0, 3.0]
end
```

Append to `test/unit_tests/Models/ut_Model_Factory.jl`:

```julia
@testset "read_properties builds block damages" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(3)
    PeriLab.Data_Manager.set_dof(2)
    PeriLab.Data_Manager.set_block_name_list(["block_1", "block_2"])
    PeriLab.Data_Manager.set_block_id_list([1, 2])
    blocks = Dict("block_1" => Dict("Block ID" => 1, "Density" => 1.0, "Horizon" => 1.0,
                                    "Damage Model" => "Dam"),
                  "block_2" => Dict("Block ID" => 2, "Density" => 1.0, "Horizon" => 1.0))
    models = Dict("Damage Models" => Dict("Dam" => Dict("Damage Model" => "Critical Stretch",
                                                        "Critical Value" => 0.1)))
    input = typed_input(Dict("Blocks" => blocks, "Models" => models))
    params = Dict{String,Any}("Blocks" => blocks, "Models" => input.models)
    PeriLab.Solver_Manager.Model_Factory.read_properties(params, input, false)
    d = PeriLab.Data_Manager.get_block_damage(1)
    @test d isa PeriLab.Solver_Manager.Model_Factory.Damage.BlockDamage
    @test PeriLab.Data_Manager.get_block_damage(2) === nothing

    blocks["block_1"]["Damage Model"] = "Dmg"
    input = typed_input(Dict("Blocks" => blocks, "Models" => models))
    params = Dict{String,Any}("Blocks" => blocks, "Models" => input.models)
    @test_logs (:error,
                "Damage Model model with name Dmg is defined in blocks, but missing in the Damage Models definition.") match_mode=:any @test_throws PeriLab.PeriLabError PeriLab.Solver_Manager.Model_Factory.read_properties(params,
                                                                                                                                                                                                                                input,
                                                                                                                                                                                                                                false)
    input = typed_input(Dict("Blocks" => blocks, "Models" => Dict()))
    params = Dict{String,Any}("Blocks" => blocks, "Models" => input.models)
    @test_logs (:error,
                "Damage Model is defined in blocks, but no Damage Models definition block exists") match_mode=:any @test_throws PeriLab.PeriLabError PeriLab.Solver_Manager.Model_Factory.read_properties(params,
                                                                                                                                                                                                                 input,
                                                                                                                                                                                                                 false)
end
```

In `ut_read_properties`, remove `"Damage Model" => "a"` from `block_3` and the two damage assertions (`get_property(3, "Damage Model", "value")`). The test passes an empty typed input, so that undefined damage name would now abort. The new testset above covers the damage path.

- [ ] **Step 2: Run to verify they fail**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Models/Damage/ut_block_damage unit_tests/Models/ut_Model_Factory`
Expected:
- `UndefVarError: block_damage`, `bind_dependent_tables!`, `get_block_damage`;
- `read_properties builds block damages` fails: no block damage is stored, and the legacy abort comes from `get_model_parameter` with the same text for the first case.

Watch that the second and third cases now fail or pass for the right reason, then keep going.

- [ ] **Step 3: Implement**

`dependent.jl`: add, and export `dependent_tables`:

```julia
"""
    dependent_tables(x)

The `Table1D` parameters of an `@params` struct or of every part of a `Composite`.
"""
function dependent_tables(x)
    found = Table1D[]
    for part in (x isa Composite ? x.parts : (x,))
        part === nothing && continue
        for fs in parameter_spec(typeof(part))
            v = getfield(part, fs.name)
            v isa Table1D && push!(found, v)
        end
    end
    return found
end
```

If `Composite` or `parameter_spec` is defined in a file that `ParameterSpec.jl` includes after `dependent.jl`, put the function at the end of `model.jl` instead. It must live in a file that sees both.

`Data_manager.jl`:
- Add `using ..ParameterSpec: Table1D, bind_table!`. Use the same number of dots as the existing `using ...PeriLabExceptions` line.
- Export `set_block_damage`, `get_block_damage` and `bind_dependent_tables!`.
- Next to `data["Block Materials"]`, add `data["Block Damages"] = Dict{Int64,Any}()`.
- Below `get_block_material`, add:

```julia
"""
	set_block_damage(block, damage)

Stores the typed damage model (`BlockDamage`) of a block.
"""
function set_block_damage(block::Int64, damage)
    data["Block Damages"][block] = damage
end

"""
	get_block_damage(block)

The typed damage model of a block, or `nothing`.
"""
function get_block_damage(block::Int64)
    return get(data["Block Damages"], block, nothing)
end

# node field a dependent value reads (the NP1 state if the field has one)
function _dependent_field(name::String)
    has_key(name * "NP1") && return get_field(name, "NP1")
    has_key(name) && return get_field(name)
    return nothing
end

"""
	bind_dependent_tables!(tables)

Binds every table to its node field (call before each compute: the N/NP1 switch
replaces the arrays). Aborts if a field is missing.
"""
function bind_dependent_tables!(tables::Vector{Table1D})
    for table in tables
        field = _dependent_field(table.field_name)
        if !(field isa Vector{Float64})
            @abort "Field \"$(table.field_name)\" required by $(table.source) does not exist or is not a per-node Vector{Float64}."
        end
        bind_table!(table, field)
    end
    return tables
end
```

`data_manager_checkpoint.jl`: add `"Block Damages"` to `CHECKPOINT_EXCLUDED_KEYS`, under the material comment, as `# typed materials and damages: rebuilt from the input by read_properties`.

`Material_Factory.jl`:
- Delete `_tables` and `_dependent_field`.
- In `block_material`, use `vcat(dependent_tables(wb.base), dependent_tables(wb.model))`.
- `bind_material!(material)` becomes:

```julia
function bind_material!(material::BlockMaterial)
    Data_Manager.bind_dependent_tables!(material.tables)
    return material
end
```

- Add `dependent_tables` to the `ParameterSpec` import and drop `bind_table!` from it if it is no longer used there.

`Damage_Factory.jl`: extend the `ParameterSpec` import with `WithBase, Table1D, dependent_tables, value`, and add:

```julia
"""
    BlockDamage

The typed damage model of a block: the shared base part, the model's own part,
and the field-dependent tables to bind before each compute.
"""
struct BlockDamage{B,M}
    base::B
    model::M
    tables::Vector{Table1D}
end

block_damage(wb::WithBase) = BlockDamage(wb.base, wb.model,
                                         vcat(dependent_tables(wb.base),
                                              dependent_tables(wb.model)))

function init_interface_crit_values(damage::BlockDamage, block_id::Int64)
    interblock = damage.base.interblock_damage
    interblock === nothing && return
    max_block_id = maximum(Data_Manager.get_block_id_list())
    inter_critical_value = Data_Manager.get_crit_values_matrix()
    if inter_critical_value == fill(-1, (1, 1, 1))
        inter_critical_value = fill(value(damage.base.critical_value, 1),
                                    (max_block_id, max_block_id, max_block_id))
    end
    for block_iId in 1:max_block_id, block_jId in 1:max_block_id
        name = "Interblock Critical Value $(block_iId)_$block_jId"
        haskey(interblock, name) &&
            (inter_critical_value[block_iId, block_jId, block_id] = interblock[name])
    end
    Data_Manager.set_crit_values_matrix(inter_critical_value)
end

function init_aniso_crit_values(aniso::AnisotropicDamageParams, block_id::Int64,
                                dof::Int64)
    aniso_crit::Dict{Int64,Any} = Data_Manager.get_aniso_crit_values()
    aniso_crit[block_id] = dof == 2 ?
                           [aniso.critical_value_x, aniso.critical_value_y] :
                           [aniso.critical_value_x, aniso.critical_value_y,
                            something(aniso.critical_value_z, aniso.critical_value_y)]
    Data_Manager.set_aniso_crit_values(aniso_crit)
end
```

`Model_Factory.jl`, in `read_properties`, after the `if material_model … end` block:

```julia
    for (block_name, block) in zip(block_name_list, block_id_list)
        block_params = get(input.sections.blocks, block_name, nothing)
        damage_name = block_params === nothing ? nothing : block_params.damage_model
        damage_name === nothing && continue
        if !haskey(input.models, "Damage Models")
            @abort "Damage Model is defined in blocks, but no Damage Models definition block exists"
        end
        if !haskey(input.damages, damage_name)
            @abort "Damage Model model with name $damage_name is defined in blocks, but missing in the Damage Models definition."
        end
        Data_Manager.set_block_damage(block, Damage.block_damage(input.damages[damage_name]))
    end
```

The legacy `get_model_parameter` call for "Damage Model" in `get_block_model_definition` still runs before this loop and aborts with the same messages. That is fine until Task 3 removes it.

- [ ] **Step 4: Run to verify they pass**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Models/Damage/ut_block_damage unit_tests/Models/ut_Model_Factory unit_tests/Models/Material/ut_block_material unit_tests/Models/Damage/ut_Damage_Factory`
Expected: all pass. The material binding tests keep their message.

- [ ] **Step 5: Full suite, then commit**

Full suite in the background; read the tail. Expected: 0 failures, 0 errors.

```bash
git add src/Core src/Support/Parameters/Spec src/Models test/unit_tests/Models test/runtests.jl
git commit -m "Typed BlockDamage per block; table binding shared with Material"
```

---

### Task 3: Damage dispatch and models on the typed interface; the damage dict is no longer stored

**Files:**
- Modify: `src/Models/Damage/Damage_Factory.jl` (`init_model`, `compute_model`, `fields_for_local_synchronization`; delete the dict `init_interface_crit_values` / `init_aniso_crit_values`)
- Modify: `src/Models/Damage/Critical_Stretch.jl`, `Energy_release.jl`, `Energy_release_aniso.jl` (typed `init_model` / `compute_model`; dict methods deleted)
- Modify: `src/Models/Model_Factory.jl`:
  - `has_block_model`, `block_model_parameters`;
  - local damping init and compute;
  - `get_block_model_definition` skips "Damage Model".
- Modify: `src/Models/Material/Material_Basis.jl` (`init_local_damping_due_to_damage`, `local_damping_due_to_damage` read `LocalDampingParams` fields)
- Test: `test/unit_tests/Models/Damage/ut_block_damage.jl` (append), `ut_Damage_Factory.jl`, `test/unit_tests/Models/Material/ut_material_basis.jl`, `test/unit_tests/Models/ut_Model_Factory.jl`

**Interfaces:**
- Consumes (Task 2): `BlockDamage`, `get_block_damage`, `bind_dependent_tables!`, typed crit-value init.
- Produces:
  - Every damage model module provides:
    - `init_model(nodes::AbstractVector{Int64}, p::<Model>Params, damage, block::Int64)`;
    - `compute_model(nodes::AbstractVector{Int64}, p::<Model>Params, damage, block::Int64, time::Float64, dt::Float64)`;
    - `fields_for_local_synchronization(model::String)`.
  - `Damage.compute_model(nodes, damage::BlockDamage, block, time, dt)`.
  - `Material_Basis.init_local_damping_due_to_damage(nodes, symmetry::String, local_damping)` and `local_damping_due_to_damage(nodes, local_damping, dt)`. `local_damping` has the fields `representative_youngs_modulus` and `damping_coefficient`.

- [ ] **Step 1: Write the failing tests**

Append to `ut_block_damage.jl`:

```julia
function ut_stretch_setup(deformed)
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(2)
    PeriLab.Data_Manager.set_dof(2)
    nn = PeriLab.Data_Manager.create_constant_node_scalar_field("Number of Neighbors", Int64)
    nn .= 1
    nlist = PeriLab.Data_Manager.create_constant_bond_scalar_state("Neighborhoodlist", Int64)
    nlist[1] = [2]
    nlist[2] = [1]
    PeriLab.Data_Manager.create_constant_node_scalar_field("Block_Id", Int64) .= 1
    PeriLab.Data_Manager.create_constant_node_scalar_field("Update", Bool)
    PeriLab.Data_Manager.create_bond_scalar_state("Bond Damage", Float64; default_value = 1)
    len = PeriLab.Data_Manager.create_constant_bond_scalar_state("Bond Length", Float64)
    len[1] .= 1.0
    len[2] .= 1.0
    _, dlen = PeriLab.Data_Manager.create_bond_scalar_state("Deformed Bond Length", Float64)
    dlen[1] .= deformed[1]
    dlen[2] .= deformed[2]
    return nothing
end

function ut_stretch_damage(raw; deformed = (1.2, 0.9))
    ut_stretch_setup(deformed)
    d = ut_typed_damage(raw)
    PeriLab.Data_Manager.set_block_damage(1, d)
    BDAM.Critical_Stretch.init_model([1, 2], d.model, d, 1)
    BDAM.Critical_Stretch.compute_model([1, 2], d.model, d, 1, 0.0, 1.0)
    bd = PeriLab.Data_Manager.get_bond_damage("NP1")
    return (bd[1][1], bd[2][1])
end

@testset "critical stretch on the typed interface" begin
    # stretches: node 1 +0.2, node 2 -0.1
    @test ut_stretch_damage(Dict("Damage Model" => "Critical Stretch",
                                 "Critical Value" => 0.15)) == (0.0, 1.0)
    @test ut_stretch_damage(Dict("Damage Model" => "Critical Stretch",
                                 "Critical Value" => 0.05, "Only Tension" => false)) ==
          (0.0, 0.0)
    @test ut_stretch_damage(Dict("Damage Model" => "Critical Stretch",
                                 "Critical Value" => 0.05)) == (0.0, 1.0)
end

@testset "anisotropic energy needs Anisotropic Damage" begin
    ut_stretch_setup((1.0, 1.0))
    PeriLab.Data_Manager.create_constant_node_scalar_field("Horizon", Float64) .= 1.0
    d = ut_typed_damage(Dict("Damage Model" => "Critical Energy Anisotropic",
                             "Critical Value" => 1.0))
    @test_logs (:error,
                "Critical Energy Anisotropic requires Anisotropic Damage.") @test_throws PeriLab.PeriLabError BDAM.Critical_Energy_Aniso.init_model([1, 2],
                                                                                                                                                       d.model,
                                                                                                                                                       d,
                                                                                                                                                       1)
end

@testset "damage dispatch reads the block damage" begin
    ut_stretch_setup((1.2, 0.9))
    PeriLab.Data_Manager.create_constant_node_scalar_field("Volume", Float64) .= 1.0
    PeriLab.Data_Manager.create_node_scalar_field("Damage", Float64)
    d = ut_typed_damage(Dict("Damage Model" => "Critical Stretch", "Critical Value" => 0.15))
    PeriLab.Data_Manager.set_block_damage(1, d)
    BDAM.init_model([1, 2], 1)
    BDAM.compute_model([1, 2], d, 1, 0.0, 1.0)
    @test PeriLab.Data_Manager.get_damage("NP1") == [1.0, 0.0]
    @test length(methods(BDAM.Critical_Stretch.compute_model)) == 1
end
```

Use the real model module names from Task 1 (`grep -n '^module' src/Models/Damage/*.jl`) for `Critical_Stretch`, `Critical_Energy_Aniso`.

`ut_material_basis.jl`: replace `@testset "ut_init_local_damping_due_to_damage"` with:

```julia
@testset "ut_init_local_damping_due_to_damage" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(3)
    PeriLab.Data_Manager.set_dof(2)
    PeriLab.Data_Manager.create_constant_node_scalar_field("Horizon", Float64) .= 1.0
    local_damping = (representative_youngs_modulus = 70.0e9, damping_coefficient = 0.5)
    PeriLab.Solver_Manager.Material_Basis.init_local_damping_due_to_damage(collect(1:2),
                                                                           "plane strain",
                                                                           local_damping)
    @test PeriLab.Data_Manager.has_key("Bond Based Constant")
end
```

The missing-key aborts are input errors now (Task 1, `local damping keys are required`).

Append to `ut_Model_Factory.jl`:

```julia
@testset "local damping reads the block damage" begin
    PeriLab.Data_Manager.initialize_data()
    MF = PeriLab.Solver_Manager.Model_Factory
    @test !MF.has_block_model(1, "Damage Model")
    d = MF.Damage.block_damage(typed_input(Dict("Blocks" => Dict("block_1" => Dict("Block ID" => 1,
                                                                                    "Density" => 1.0,
                                                                                    "Horizon" => 1.0,
                                                                                    "Damage Model" => "Dam")),
                                                "Models" => Dict("Damage Models" => Dict("Dam" => Dict("Damage Model" => "Critical Stretch",
                                                                                                       "Critical Value" => 0.1,
                                                                                                       "Local Damping" => Dict("Representative Young's modulus" => 1.0,
                                                                                                                               "Damping coefficient" => 0.5))))).damages["Dam"])
    PeriLab.Data_Manager.set_block_damage(1, d)
    @test MF.has_block_model(1, "Damage Model")
    @test MF.block_model_parameters(1, "Damage Model") === d
    @test MF.block_local_damping(1).damping_coefficient == 0.5
    @test MF.block_local_damping(2) === nothing
end
```

In `ut_get_block_model_definition`, replace the `Damage Model` assertion with `@test isempty(PeriLab.Data_Manager.get_properties(3, "Damage Model"))`.

In `ut_Damage_Factory.jl`:
- delete `@testset "ut_Damage_factory_exceptions"` (an unknown name is a parse error: Task 1, `damage input errors`);
- delete `@testset "ut_init_interface_crit_values"` (typed version: Task 2, `typed interface critical values`).

- [ ] **Step 2: Run to verify they fail**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Models/Damage/ut_block_damage unit_tests/Models/Damage/ut_Damage_Factory unit_tests/Models/Material/ut_material_basis unit_tests/Models/ut_Model_Factory`
Expected:
- `MethodError`s for the typed `init_model` / `compute_model` of the models and the factory;
- `UndefVarError: block_local_damping`;
- `has_block_model` false after `set_block_damage`;
- local damping init fails reading `["Local Damping"]` of a NamedTuple.

- [ ] **Step 3: Implement**

**Damage_Factory.jl.**
- Delete the dict `init_interface_crit_values(damage_parameter::Dict, …)` and `init_aniso_crit_values(damage_parameter::Dict, …)`.
- Remove `create_module_specifics` from the `ModuleLoader` import: the module list is still used to include the files, but the name lookup is gone.
- Replace `compute_model`, `fields_for_local_synchronization` and `init_model` with:

```julia
model_module(damage::BlockDamage) = parentmodule(typeof(damage.model))

"""
    compute_model(nodes, damage, block, time, dt)

Binds the damage's dependent values, computes the block's damage model and the
damage index.
"""
function compute_model(nodes::AbstractVector{Int64}, damage::BlockDamage, block::Int64,
                       time::Float64, dt::Float64)
    Data_Manager.bind_dependent_tables!(damage.tables)
    model_module(damage).compute_model(nodes, damage.model, damage, block, time, dt)
    if isnothing(Data_Manager.get_filtered_nlist())
        @timeit "compute index" return damage_index(nodes)
    end
    @timeit "compute index" return damage_index(nodes, Data_Manager.get_filtered_nlist())
end

function fields_for_local_synchronization(model, block)
    damage = Data_Manager.get_block_damage(block)
    model_module(damage).fields_for_local_synchronization(model)
end

function init_model(nodes::AbstractVector{Int64}, block::Int64)
    damage = Data_Manager.get_block_damage(block)
    model_module(damage).init_model(nodes, damage.model, damage, block)
    model_module(damage).fields_for_local_synchronization("Damage Model")
    init_interface_crit_values(damage, block)
    damage.base.anisotropic_damage === nothing ||
        init_aniso_crit_values(damage.base.anisotropic_damage, block, Data_Manager.get_dof())
end
```

Keep the existing docstrings, updated to the new arguments.

**The models.** In each file, replace the dict `init_model` / `compute_model` with typed methods. A thin entry function passes the concrete `Critical Value` to a kernel, which is the old body with the parameter reads replaced:

| File | Old read | New |
|---|---|---|
| all | `damage_parameter["Critical Value"]` / `critical_value_fn(iID)` / `get_dependent_value(...)` / `is_dependent("Critical Value", …)` + `interpol_data(dependent_field[iID], …)` | `value(critical_value, iID)` (`critical_value` = kernel argument) |
| all | `get(damage_parameter, "Only Tension", …)` | `p.only_tension` |
| all | `haskey(damage_parameter, "Interblock Damage")` / `Data_Manager.haskey(...)` | `damage.base.interblock_damage !== nothing` |
| energy (both) init | `get(damage_parameter, "Thickness", 1)` | `p.thickness` |
| aniso init | `haskey(damage_parameter, "Anisotropic Damage")` | `damage.base.anisotropic_damage !== nothing` |
| aniso compute | the `is_dependent(param_name, damage_parameter)` block for `Interblock Critical Value i_j` | delete (keep `critical_energy_value = inter_critical_energy[…]`) |

Shape, shown for Critical Stretch; do the same in the two energy files:

```julia
function init_model(nodes::AbstractVector{Int64}, p::CriticalStretchParams, damage,
                    block::Int64)
    Data_Manager.create_constant_bond_scalar_state("Bond Stretch", Float64)
end

function compute_model(nodes::AbstractVector{Int64}, p::CriticalStretchParams, damage,
                       block::Int64, time::Float64, dt::Float64)
    return _critical_stretch!(nodes, damage.base.critical_value, p.only_tension,
                              damage.base.interblock_damage !== nothing, block)
end

# function barrier: `critical_value` has a concrete type (Constant or Table1D) here
function _critical_stretch!(nodes::AbstractVector{Int64}, critical_value, tension::Bool,
                            inter_block_damage::Bool, block::Int64)
    # old body; `crit_stretch = critical_field ? critical_stretch[iID] :
    #   inter_block_damage ? inter_critical_stretch[…] : value(critical_value, iID)`
end
```

- In the Critical Stretch body, rename the field variable so it never holds both a field and a number: `critical_stretch = Data_Manager.get_field("Critical_Value")` only inside `if critical_field`.
- In the Critical Energy Anisotropic `init_model`, add as the first lines: `damage.base.anisotropic_damage === nothing && @abort "Critical Energy Anisotropic requires Anisotropic Damage."`. Keep the existing rotation check (`Anisotropic damage requires Angles field`).
- Add `value` to each model's `ParameterSpec` import and `@abort` where it is now used.
- Drop imports that became unused: `get_dependent_value`, `is_dependent`, `interpol_data`. Grep each file before removing.

Afterwards run `grep -rn '::Dict\|damage_parameter\|get_dependent_value\|is_dependent' src/Models/Damage --include=*.jl | grep -v Damage_template`. Expected: no matches.

**Model_Factory.jl.**

```julia
# typed categories (Data_Manager slots); the other categories still use the property dicts
function typed_block_model(block::Int64, name::String)
    name == "Material Model" && return Data_Manager.get_block_material(block)
    name == "Damage Model" && return Data_Manager.get_block_damage(block)
    return nothing
end
const TYPED_CATEGORIES = ("Material Model", "Damage Model")
has_block_model(block::Int64, name::String) = name in TYPED_CATEGORIES ?
                                              typed_block_model(block, name) !== nothing :
                                              Data_Manager.check_property(block, name)
block_model_parameters(block::Int64, name::String) = name in TYPED_CATEGORIES ?
                                                     typed_block_model(block, name) :
                                                     Data_Manager.get_properties(block, name)
# local damping of a block's damage model, or `nothing`
function block_local_damping(block::Int64)
    damage = Data_Manager.get_block_damage(block)
    return damage === nothing ? nothing : damage.base.local_damping
end
```

- In `init_models`, replace the `haskey(block_model_parameters(…), "Local Damping")` branch with:

```julia
if active_model_name == "Damage Model" && block_local_damping(block) !== nothing
    Material.init_local_damping(block_nodes[block], local_damping_symmetry(block),
                                block_local_damping(block))
end
```

- In `compute_models`, the local-damping loop becomes `local_damping = block_local_damping(block)`, then `local_damping === nothing && continue`. Pass `local_damping` to `Material.compute_local_damping(active_nodes, local_damping, dt)`.
- In `get_block_model_definition`, replace `model == "Material Model" && continue` with `model in TYPED_CATEGORIES && continue   # typed: input.materials / input.damages`.

**Material_Basis.jl.**
- `init_local_damping_due_to_damage(nodes, symmetry::String, local_damping)`:
  - drop the two missing-key aborts;
  - log `@info "Local damping is active with damping coefficient $(local_damping.damping_coefficient)"`;
  - keep the rest.
- `local_damping_due_to_damage`:
  - `local_damping = params.damping_coefficient`;
  - `E = params.representative_youngs_modulus`.
- In `Material_Factory.jl`, rename the argument of `init_local_damping(nodes, symmetry, damage_parameter)` to `local_damping`.

- [ ] **Step 4: Run to verify they pass**

Run the Step 2 command, plus `unit_tests/Support/Parameters/Input/ut_damage_models`.
Expected: all pass.

- [ ] **Step 5: Full suite, then commit**

Full suite in the background; read the tail. Expected: 0 failures, 0 errors. The damage fullscale decks pin the numbers of all three models, interblock, anisotropic, temperature-dependent and local damping.

```bash
git add src/Models test/unit_tests/Models
git commit -m "Damage dispatch and models read the block damage; the damage dict is no longer stored"
```

---

### Task 4: Damage template on the typed interface

**Files:**
- Modify: `src/Models/Damage/Damage_template/damage_template.jl`
- Modify: `test/unit_tests/Models/ut_templates.jl` (`ut_damage_template`; the file is commented out in `runtests.jl`, run it standalone)

**Interfaces:**
- Consumes (Task 3): the typed model interface.
- Produces: `Damage_template.DamageTemplateParams` (unregistered) and the typed `init_model` / `compute_model` of the template.

- [ ] **Step 1: Write the failing test**

In `ut_templates.jl`, replace the body of `@testset "ut_damage_template"` with:

```julia
@testset "ut_damage_template" begin
    DT = PeriLab.Solver_Manager.Model_Factory.Damage.Damage_template
    @test DT.damage_name() == "Damage Template"
    p = DT.DamageTemplateParams()
    damage = PeriLab.Solver_Manager.Model_Factory.Damage.BlockDamage(nothing, p,
                                                                      PeriLab.ParameterSpec.Table1D[])
    DT.init_model(Vector{Int64}(1:3), p, damage, 1)
    DT.compute_model(Vector{Int64}(1:3), p, damage, 1, 0.0, 0.0)
    @test length(methods(DT.compute_model)) == 1
    DT.fields_for_local_synchronization("")
end
```

Adjust the module path if the template module is reachable under another name (`grep -rn 'Damage_template' test/unit_tests/Models/ut_templates.jl`).

- [ ] **Step 2: Run to verify it fails**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Models/ut_templates`
Expected: `ut_damage_template` errors with `UndefVarError: DamageTemplateParams`. The file's other testsets error as before (the template modules are not loaded there); ledger that.

- [ ] **Step 3: Implement**

In `damage_template.jl`:
- After the `using` lines, add:

```julia
using ......ParameterSpec: @params, register_damage

"""
    DamageTemplateParams

Declare the YAML keys your damage model needs beyond the shared damage keys
(Critical Value, Interblock Damage, Anisotropic Damage, Local Damping are in
`damage.base`). Register it under your model name by uncommenting `__init__` (the
template stays unregistered so that a copy never collides with it).
"""
@params struct DamageTemplateParams
end

# __init__() = register_damage("Damage Template", DamageTemplateParams)
```

  Use the same dot count as the template's `Data_Manager` import plus one.
- `compute_model(nodes, damage_parameter::Dict, block, time, dt)` becomes `compute_model(nodes::AbstractVector{Int64}, p::DamageTemplateParams, damage, block::Int64, time::Float64, dt::Float64)`. Update the docstring argument lines:
  - `p`: the model parameters;
  - `damage::BlockDamage`: base and model parameters.

  The two `@info` lines that mention `damage_parameter` become:
  - "Fill the compute_model(nodes, p, damage, block, time, dt) function."
  - "The Data_Manager, p and damage hold all you need to solve your problem on material level."
- `init_model(nodes, damage_parameter::Dict, block)` becomes `init_model(nodes::AbstractVector{Int64}, p::DamageTemplateParams, damage, block::Int64)`, with the docstring updated the same way.
- Fix the docstring typo `c ompute_damage(...)` to `compute_model(nodes, p, damage, block, time, dt)`.

Afterwards run `grep -rn '::Dict\|damage_parameter' src/Models/Damage`. Expected: no matches.

- [ ] **Step 4: Run to verify it passes**

Run the Step 2 command.
Expected: `ut_damage_template` passes. The other testsets of the file error as before.

- [ ] **Step 5: Full suite, then commit**

Full suite in the background; read the tail. Expected: 0 failures, 0 errors.

```bash
git add src/Models/Damage test/unit_tests/Models/ut_templates.jl
git commit -m "Damage template on the typed interface"
```
