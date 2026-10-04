# Typed Input Parameters — Phase 2b-5: Contact, FEM/Coupling, Surface Correction

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Contact, FEM/coupling and surface correction read the typed `ContactInput`, `FEMParams` and `SurfaceCorrectionParams` from `PeriLabInput` instead of the params `Dict`. Two pre-existing `eval_bc!` crashes found in the 2b-4 review are fixed as well.

**Architecture:**
- Code-level defaults that used to be written into the params dict at init time become declared defaults in the `@params` structs:
  - Penalty: `Contact Stiffness` 1e8, `Friction Coefficient` 0, `Symmetry` "3D".
  - Arlequin: `PD Weight` 0.5, `Kappa` 1.0.
- Checks that need only the section itself become `check!` input errors. The checks that need the block list stay at runtime.
- The typed section is stored once in a `Data_Manager` slot (`Contact Properties`, `FEM Parameters`, `Surface Correction`). Compute code reads its fields through a function barrier.

**Tech Stack:** Julia 1.12, `PeriLab.ParameterSpec` (`@params`, `check!`), `PeriLab.InputDeck`.

**Spec:** `docs/superpowers/specs/2026-10-03-typed-input-parameters-design.md`

## Global Constraints

- Full YAML backward compatibility: every shipped deck still parses in strict mode (golden test `ut_golden_decks.jl`).
- Unknown keys are errors (strict), with the `--no_strict` / `Strict Validation: false` opt-out.
- No fixed units; `quantity` is documentation only.
- Read fields directly; functions only for real logic (fallbacks, parsing, list building).
- Structs are immutable; per-node state stays in `Data_Manager` fields.
- Numerical results unchanged: the full suite (`Pkg.test()`) must stay green after every task.
- Branch `feature/typed-input-parameters`; one commit per task; never merge.
- Read the full-suite result before committing a task that runs it.

## Review Focus

1. **A contact deck without `Globals`** gets search frequency 1 and surface-only contact nodes, as before. Pinned in Task 2 (`contact defaults`).
2. **A group with its own `Global Search Frequency`** overrides `Globals`, and a group without one inherits it. Pinned in Task 2 (`contact_search_frequency`).
3. **An FEM `Degree` string such as `"2 1 1"`** gives one degree per direction. A malformed string like `"2 a"` is an input error, not a parse crash. Pinned in Task 2 (`FEM degree`).
4. **A coupled FEM deck without `PD Weight`/`Kappa`** uses 0.5 / 1.0 both at init and in every `compute_coupling` call. The old code relied on init mutating a dict that compute later re-read. Pinned in Task 4 (`coupling defaults reach compute`).
5. **A deck without `Surface Correction`** runs no correction. A deck with `Update: true` recomputes the correction each step. Pinned in Task 5.

Known pre-existing behaviour left unchanged (out of scope, report only): `compute_contact_model` computes `n = get_search_step(cg) + 1` but never stores it, so the global search runs every step regardless of `Global Search Frequency`.

---

## Shared test commands

Single unit-test file (run from repo root):

```bash
cd test && julia --project=.. -e 'using Test, Logging, MPI; import PeriLab; Logging.disable_logging(Logging.Warn); include("helper.jl"); MPI.Init(); @testset "t" begin include("unit_tests/<PATH>.jl") end' 2>&1 | tail -30
```

`Logging.disable_logging(Logging.Warn)` also hides `:warn` records, so never assert `@test_logs (:warn, …)` in these files. `:error` records are still captured.

Input tests:

```bash
julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl 2>&1 | tail -15
```

Full suite (~30 min; run in background, redirect to a file, read its tail):

```bash
julia --project=. -e 'using Pkg; Pkg.test()' > <workspace>/suite.log 2>&1
```

## Relative imports of `InputDeck`

`PeriLab.InputDeck` lives at the package root. A top-level package module is its own parent, so extra dots stay at the root. Use these depths:

| File | Module path | Import |
|---|---|---|
| `src/Models/Model_Factory.jl` | PeriLab.Solver_Manager.Model_Factory | `using ...InputDeck: …` |
| `src/Models/Contact/Contact_Factory.jl` | …Model_Factory.Contact | `using .....InputDeck: …` |
| `src/Models/Contact/Contact_search.jl` | …Contact.Contact_Search | `using .....InputDeck: …` |
| `src/Models/Contact/Penalty_model.jl` | …Contact.Penalty_Model | `using .....InputDeck: …` |
| `src/Models/Surface_correction/Surface_correction.jl` | …Model_Factory.Surface_Correction | `using .....InputDeck: …` |
| `src/FEM/FEM_Factory.jl` | PeriLab.Solver_Manager.FEM | `using ...InputDeck: …` |
| `src/FEM/FEM_basis.jl` (module FEM_Basis inside FEM) | …FEM.FEM_Basis | (no import needed) |
| `src/FEM/Coupling/Coupling_Factory.jl` | …FEM.Coupling | `using ....InputDeck: …` |
| `src/FEM/Coupling/Arlequin_coupling.jl` | …Coupling.Arlequin_Coupling | `using ......InputDeck: …` |

---

### Task 1: `eval_bc!` fixes and the z-in-2D regression test

**Files:**
- Modify: `src/Core/BC_manager.jl` (String method of `eval_bc!`, around lines 312 and 346)
- Test: `test/unit_tests/Core/ut_BC_manager.jl` (append three testsets)

**Interfaces:**
- Consumes: `init_BCs(::Dict{String,BoundaryConditionParams})`, `eval_bc!(field_values, bc, coordinates, time, step_time, dof, initial, name = "BC_1", neumann = false)`. These are unchanged.
- Produces: nothing new. `eval_bc!` now returns the compiled function on every path that compiled one.

- [ ] **Step 1: Write the failing tests** (append to `test/unit_tests/Core/ut_BC_manager.jl`)

```julia
@testset "skipped Initial condition still returns the compiled function" begin
    field = zeros(3)
    coordinates = [0.0 0.0 0.0; 1.0 0.0 0.0; 2.0 0.0 0.0]
    bc = PeriLab.Solver_Manager.Boundary_Conditions.eval_bc!(field, "10*sin(t)", coordinates,
                                                             1.0, 1.0, 3, true)
    @test bc isa Function
    @test field == zeros(3)
    # a second apply must not re-clean and re-parse the expression
    bc = PeriLab.Solver_Manager.Boundary_Conditions.eval_bc!(field, bc, coordinates,
                                                             2.0, 1.0, 3, true)
    @test bc isa Function
end

@testset "z in a 2D expression is an input error" begin
    field = zeros(2)
    coordinates = [0.0 0.0; 1.0 0.0]
    @test_logs (:error, "z is not valid in a 2D problem.") @test_throws PeriLab.PeriLabError begin
        PeriLab.Solver_Manager.Boundary_Conditions.eval_bc!(field, "10*z", coordinates,
                                                            0.0, 0.0, 2, false)
    end
    # names that merely contain a z (e.g. a function) are not the coordinate z
    @test PeriLab.Solver_Manager.Boundary_Conditions.eval_bc!(field, "10*t", coordinates,
                                                              1.0, 1.0, 2, false) isa Function
end

@testset "a z condition in 2D skips only itself" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(3)
    PeriLab.Data_Manager.set_dof(2)
    PeriLab.Data_Manager.set_nset("Nset_1", [1, 2])
    PeriLab.Data_Manager.set_glob_to_loc(Dict(1 => 1, 2 => 2, 3 => 3))
    PeriLab.Data_Manager.create_node_vector_field("Displacements", Float64, 2)
    raw = Dict{String,Any}("BC_$i" => Dict("Variable" => "Displacements",
                                           "Node Set" => "Nset_1",
                                           "Coordinate" => "x", "Value" => 0.0)
                           for i in 1:6)
    raw["BC_z"] = Dict("Variable" => "Displacements", "Node Set" => "Nset_1",
                       "Coordinate" => "z", "Value" => 0.0)
    params = typed_section(Dict{String,PeriLab.InputDeck.BoundaryConditionParams}, raw)
    # the old `break` dropped every condition iterated after the z condition;
    # the test only discriminates if some condition comes after it
    @test findfirst(==("BC_z"), collect(keys(params))) < length(params)
    bcs = PeriLab.Solver_Manager.Boundary_Conditions.init_BCs(params)
    @test sort(collect(keys(bcs))) == sort(["BC_$i" for i in 1:6])
end
```

- [ ] **Step 2: Run them to verify the first two fail**

Run: the single-file command with `<PATH>` = `Core/ut_BC_manager`
Expected:
- "skipped Initial condition…": FAILS at `@test bc isa Function`, because a `String` is returned.
- "z in a 2D expression…": FAILS, because an `UndefVarError: z` is thrown instead of `PeriLabError`.
- "a z condition in 2D skips only itself": PASSES. The `continue` fix landed in 2b-4, so this is a regression pin and cannot be watched failing. Ledger a `Ruling:` line saying so.

- [ ] **Step 3: Implement**

In `src/Core/BC_manager.jl`, String method of `eval_bc!`, replace

```julia
    if dof < 2 && "z" in bc
        @abort "z is not valid in a 2D problem."
        return nothing
    end
```

with

```julia
    if dof < 3 && occursin(r"\bz\b", bc)
        @abort "z is not valid in a 2D problem."
        return nothing
    end
```

In the same method, replace

```julia
    if isnothing(value) || (initial && time != 0.0)
        return bc
    end
```

with

```julia
    if isnothing(value) || (initial && time != 0.0)
        return bc_out
    end
```

Leave the Function method's `return bc` unchanged; there `bc` already is the function.

- [ ] **Step 4: Run to verify they pass**

Run: the single-file command with `<PATH>` = `Core/ut_BC_manager`
Expected: all testsets pass (previous 92 + new).

- [ ] **Step 5: Commit**

```bash
git add src/Core/BC_manager.jl test/unit_tests/Core/ut_BC_manager.jl
git commit -m "BC expressions: z check in 2D, compiled function kept on skipped Initial"
```

---

### Task 2: Typed input: defaults, checks and helpers for contact, FEM and surface correction

**Files:**
- Modify: `src/Support/Parameters/Input/contact.jl`
- Modify: `src/Support/Parameters/Input/blocks.jl` (`FEMCouplingParams`, `FEMParams`, `SurfaceCorrectionParams`)
- Modify: `src/Support/Parameters/Input/InputDeck.jl` (exports)
- Create: `test/unit_tests/Support/Parameters/Input/ut_contact_fem_typed.jl`
- Modify: `test/unit_tests/Support/Parameters/Input/input_tests.jl` (include the new file)
- Modify: `test/unit_tests/Support/Parameters/Input/ut_outputs_conditions_contact.jl` (friction default assertion)
- Modify: `test/helper.jl` (add `typed_contact`)

**Interfaces:**
- Consumes: `@params`, `opt`, `check!`, `add_error!`, `join_path`, `parse_contact`. All of these exist already.
- Produces (exported from `PeriLab.InputDeck`):
  - `ContactInput`, `ContactModelParams`, `ContactGroupParams`, `ContactGlobalsParams` (types; field names as today).
  - `ContactModelParams.contact_stiffness::Float64` (default 1e8), `.friction_coefficient::Float64` (default 0.0), `.symmetry::String` (default "3D").
  - `contact_blocks(contact::ContactInput)::Vector{Int64}`: sorted, unique master and slave block ids.
  - `contact_search_frequency(group::ContactGroupParams, globals::ContactGlobalsParams)::Int64`.
  - `FEMParams`, `FEMCouplingParams`; `FEMCouplingParams.pd_weight::Float64` (default 0.5) and `.kappa::Float64` (default 1.0).
  - `fem_degree(fem::FEMParams)::Vector{Int64}`: the degrees as written, one entry or one per direction.
  - `SurfaceCorrectionParams`; `type` allowed only `"Volume Correction"`.
  - Test helper `typed_contact(raw::AbstractDict)::ContactInput`.

- [ ] **Step 1: Add the test helper** (append to `test/helper.jl`)

```julia
"""
    typed_contact(raw)

The typed `Contact` section from its raw dict; aborts on input errors.
"""
function typed_contact(raw::AbstractDict)
    ctx = PeriLab.ParameterSpec.ParseContext()
    value = PeriLab.InputDeck.parse_contact(Dict{String,Any}(string(k) => v for (k, v) in raw),
                                            "Contact", ctx)
    PeriLab.ParameterSpec.report!(ctx)
    return value
end
```

- [ ] **Step 2: Write the failing tests**

Create `test/unit_tests/Support/Parameters/Input/ut_contact_fem_typed.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

const CF = PeriLab.InputDeck

cf_group(master, slave; extra...) = Dict{String,Any}("Master Block ID" => master,
                                                     "Slave Block ID" => slave,
                                                     "Search Radius" => 0.01,
                                                     (string(k) => v for (k, v) in extra)...)

function cf_contact(raw)
    ctx = PeriLab.ParameterSpec.ParseContext()
    return CF.parse_contact(raw, "Contact", ctx), ctx
end

@testset "contact defaults" begin
    c, ctx = cf_contact(Dict{String,Any}("C" => Dict{String,Any}("Type" => "Penalty Contact",
                                                                 "Contact Radius" => 0.005,
                                                                 "Contact Groups" => Dict{String,Any}("g" => cf_group(2, 1)))))
    @test isempty(ctx.errors)
    @test c.globals.global_search_frequency === 1
    @test c.globals.only_surface_contact_nodes
    m = c.models["C"]
    @test m.contact_stiffness === 1e8
    @test m.friction_coefficient === 0.0
    @test m.symmetry == "3D"
end

@testset "contact_search_frequency" begin
    c, ctx = cf_contact(Dict{String,Any}("Globals" => Dict{String,Any}("Global Search Frequency" => 3),
                                         "C" => Dict{String,Any}("Type" => "Penalty Contact",
                                                                 "Contact Radius" => 0.005,
                                                                 "Contact Groups" => Dict{String,Any}("own" => cf_group(2, 1; var"Global Search Frequency" = 5),
                                                                                                      "inherit" => cf_group(3, 4)))))
    @test isempty(ctx.errors)
    groups = c.models["C"].contact_groups
    @test CF.contact_search_frequency(groups["own"], c.globals) === 5
    @test CF.contact_search_frequency(groups["inherit"], c.globals) === 3
end

@testset "contact_blocks" begin
    c, ctx = cf_contact(Dict{String,Any}("C" => Dict{String,Any}("Type" => "Penalty Contact",
                                                                 "Contact Radius" => 0.005,
                                                                 "Contact Groups" => Dict{String,Any}("a" => cf_group(1, 2),
                                                                                                      "c" => cf_group(3, 5),
                                                                                                      "q" => cf_group(8, 5))),
                                         "D" => Dict{String,Any}("Type" => "Penalty Contact",
                                                                 "Contact Radius" => 0.005,
                                                                 "Contact Groups" => Dict{String,Any}("ba" => cf_group(2, 1)))))
    @test isempty(ctx.errors)
    @test CF.contact_blocks(c) == [1, 2, 3, 5, 8]
end

@testset "contact group checks" begin
    c, ctx = cf_contact(Dict{String,Any}("C" => Dict{String,Any}("Type" => "Penalty Contact",
                                                                 "Contact Radius" => 0.005,
                                                                 "Contact Groups" => Dict{String,Any}("self" => cf_group(1, 1),
                                                                                                      "flat" => Dict{String,Any}("Master Block ID" => 1,
                                                                                                                                 "Slave Block ID" => 2,
                                                                                                                                 "Search Radius" => 0.0)))))
    @test c === nothing
    msgs = Dict(e.path => e.message for e in ctx.errors)
    @test msgs["Contact.C.\"Contact Groups\".self"] ==
          "Master Block ID and Slave Block ID are equal; self contact is not implemented"
    @test msgs["Contact.C.\"Contact Groups\".flat.\"Search Radius\""] ==
          "must be greater than zero"
end

@testset "FEM degree" begin
    fem(degree) = Dict{String,Any}("Element Type" => "Lagrange", "Degree" => degree,
                                   "Material Model" => "Elastic")
    @test CF.fem_degree(typed_section(CF.FEMParams, fem(2))) == [2]
    @test CF.fem_degree(typed_section(CF.FEMParams, fem("2 1 1"))) == [2, 1, 1]
    for bad in ("2 a", "", "0 1")
        ctx = PeriLab.ParameterSpec.ParseContext()
        PeriLab.ParameterSpec.convert_value(CF.FEMParams, fem(bad), "FEM", ctx)
        @test only(ctx.errors).path == "FEM.Degree"
    end
    ctx = PeriLab.ParameterSpec.ParseContext()
    PeriLab.ParameterSpec.convert_value(CF.FEMParams, fem(0), "FEM", ctx)
    @test only(ctx.errors).path == "FEM.Degree"
end

@testset "coupling defaults" begin
    f = typed_section(CF.FEMParams,
                      Dict{String,Any}("Element Type" => "Lagrange", "Degree" => 1,
                                       "Material Model" => "Elastic",
                                       "Coupling" => Dict{String,Any}("Coupling Type" => "Arlequin")))
    @test f.coupling.pd_weight === 0.5
    @test f.coupling.kappa === 1.0
    @test f.coupling.coupling_block === nothing
end

@testset "surface correction type" begin
    ctx = PeriLab.ParameterSpec.ParseContext()
    PeriLab.ParameterSpec.convert_value(CF.SurfaceCorrectionParams,
                                        Dict{String,Any}("Type" => "Area Correction"),
                                        "Surface Correction", ctx)
    @test only(ctx.errors).path == "\"Surface Correction\".Type"
end
```

Add `include("ut_contact_fem_typed.jl")` to `test/unit_tests/Support/Parameters/Input/input_tests.jl`, after the line that includes `ut_bc_typed.jl`.

In `test/unit_tests/Support/Parameters/Input/ut_outputs_conditions_contact.jl`, change

```julia
    @test m.contact_stiffness === 1e8 && m.friction_coefficient === nothing
```

to

```julia
    @test m.contact_stiffness === 1e8 && m.friction_coefficient === 0.0
```

- [ ] **Step 3: Run to verify they fail**

Run: the input-tests command
Expected: failures and errors in the new file:
- `UndefVarError: contact_search_frequency` / `contact_blocks` / `fem_degree`.
- `friction_coefficient === 0.0` fails (it is `nothing`).
- `pd_weight === 0.5` fails.
- The checks produce no errors.

If an error path string differs only in quoting from the expected one, print `ctx.errors` once and match the format the existing tests in `ut_outputs_conditions_contact.jl` use (keys with spaces are quoted). Ledger the corrected expectation.

- [ ] **Step 4: Implement**

`src/Support/Parameters/Input/contact.jl`: replace `ContactModelParams` with

```julia
@params struct ContactModelParams
    type::String = req("Type")
    contact_radius::Float64 = req("Contact Radius"; min = 0, quantity = :length)
    contact_stiffness::Float64 = opt("Contact Stiffness"; default = 1e8, min = 0)
    friction_coefficient::Float64 = opt("Friction Coefficient"; default = 0.0, min = 0)
    symmetry::String = opt("Symmetry"; default = "3D",
                           description = "\"plane stress\", \"plane strain\" or 3D (anything else)")
    contact_groups::Dict{String,ContactGroupParams} = req("Contact Groups")
end
```

Append to `contact.jl`:

```julia
function check!(g::ContactGroupParams, path::String, ctx::ParseContext)
    if g.master_block_id == g.slave_block_id
        add_error!(ctx, path,
                   "Master Block ID and Slave Block ID are equal; self contact is not implemented")
    end
    g.search_radius > 0 ||
        add_error!(ctx, join_path(path, "Search Radius"), "must be greater than zero")
    return nothing
end

"Sorted ids of all blocks that take part in a contact group."
function contact_blocks(contact::ContactInput)
    ids = Int64[]
    for model in values(contact.models), group in values(model.contact_groups)
        push!(ids, group.master_block_id, group.slave_block_id)
    end
    return sort!(unique!(ids))
end

"Global search frequency of a contact group; falls back to `Globals`."
contact_search_frequency(group::ContactGroupParams, globals::ContactGlobalsParams) = something(group.global_search_frequency,
                                                                                               globals.global_search_frequency)
```

`src/Support/Parameters/Input/blocks.jl`: replace `FEMCouplingParams` and `SurfaceCorrectionParams` with

```julia
@params struct FEMCouplingParams
    coupling_type::String = req("Coupling Type")
    pd_weight::Float64 = opt("PD Weight"; default = 0.5)
    kappa::Float64 = opt("Kappa"; default = 1.0)
    coupling_block::Union{Nothing,Int64} = opt("Coupling Block"; default = nothing)
end
```

```julia
@params struct SurfaceCorrectionParams
    type::String = req("Type"; allowed = ["Volume Correction"])
    update::Bool = opt("Update"; default = false)
end
```

and append after `FEMParams`:

```julia
function check!(f::FEMParams, path::String, ctx::ParseContext)
    degrees = f.degree isa Int64 ? [f.degree] :
              [tryparse(Int64, part) for part in split(f.degree)]
    if isempty(degrees) || any(d -> d === nothing || d < 1, degrees)
        add_error!(ctx, join_path(path, "Degree"),
                   "expected a positive integer or positive integers separated by spaces, got \"$(f.degree)\"")
    end
    return nothing
end

"Polynomial degrees of the FEM elements as written: one value, or one per direction."
fem_degree(f::FEMParams) = f.degree isa Int64 ? [f.degree] : parse.(Int64, split(f.degree))
```

`src/Support/Parameters/Input/InputDeck.jl`: extend the export list with

```julia
       ContactInput, ContactModelParams, ContactGroupParams, ContactGlobalsParams,
       contact_blocks, contact_search_frequency, FEMParams, FEMCouplingParams, fem_degree,
       SurfaceCorrectionParams
```

- [ ] **Step 5: Run to verify they pass**

Run: the input-tests command
Expected: all pass, including `ut_golden_decks.jl` (every shipped deck still parses in strict mode). If a golden deck fails on a new check, stop and ledger it. Every shipped deck has search radius > 0, distinct master/slave, valid degrees and `Volume Correction`, so a failure means the check is wrong.

- [ ] **Step 6: Commit**

```bash
git add src/Support/Parameters/Input/contact.jl src/Support/Parameters/Input/blocks.jl \
        src/Support/Parameters/Input/InputDeck.jl test/helper.jl \
        test/unit_tests/Support/Parameters/Input/ut_contact_fem_typed.jl \
        test/unit_tests/Support/Parameters/Input/input_tests.jl \
        test/unit_tests/Support/Parameters/Input/ut_outputs_conditions_contact.jl
git commit -m "Typed contact/FEM/surface correction: declared defaults, checks, helpers"
```

---

### Task 3: Contact consumes `ContactInput`

**Files:**
- Modify: `src/Core/Data_manager.jl:194` (`data["Contact Properties"] = nothing`)
- Modify: `src/Core/Data_manager/data_manager_contact.jl:166` (`set_contact_properties`)
- Modify: `src/Models/Model_Factory.jl` (`init_models` signature, both `check_contact` methods)
- Modify: `src/Core/Solver/Solver_manager.jl:126` (pass `input` to `init_models`)
- Modify: `src/Models/Contact/Contact_Factory.jl`
- Modify: `src/Models/Contact/Contact_search.jl`
- Modify: `src/Models/Contact/Penalty_model.jl`
- Test: `test/unit_tests/Models/Contact/ut_Contact_Factory.jl`, `test/unit_tests/Models/Contact/ut_Penalty_model.jl`

**Interfaces:**
- Consumes (Task 2): `ContactInput`, `ContactModelParams`, `ContactGroupParams`, `contact_blocks`, `contact_search_frequency`, `typed_contact`.
- Produces:
  - `Data_Manager.set_contact_properties(contact)`: untyped; stores `ContactInput` or `nothing`. `get_contact_properties()` returns it; the default is `nothing`.
  - `Model_Factory.init_models(params::Dict, input::PeriLabInput, block_nodes, solver_options, synchronise_field)`.
  - `Model_Factory.check_contact(::Nothing)`, `check_contact(::ContactInput)`, `check_contact(::Nothing, time, dt)`, `check_contact(::ContactInput, time, dt)`.
  - `Contact.init_contact_model(contact::ContactInput)`, `Contact.compute_contact_model(contact::ContactInput, time::Float64, dt::Float64)`, `Contact.check_valid_contact_model(contact::ContactInput, block_ids)` (untyped: `Data_Manager.get_all_blocks()` starts as `Vector{Any}`).
  - `Contact_Search.init_contact_search(cm::String)`, `compute_contact_pairs(cg::String, group::ContactGroupParams, contact_radius::Float64)`, `global_contact_search(group::ContactGroupParams)`, `local_contact_search(contact_radius::Float64, master_nodes, slave_nodes)`.
  - Contact model module interface (`Penalty_Model`): `init_contact_model(params::ContactModelParams)`, `compute_contact_model(cg, params::ContactModelParams, compute_master_force_density, compute_slave_force_density)`.

- [ ] **Step 1: Write the failing tests**

Replace the body of `test/unit_tests/Models/Contact/ut_Contact_Factory.jl` from `@testset "ut_check_valid_contact_model"` through the end of `@testset "ut_get_all_contact_blocks"` with:

```julia
ut_cf_group(master, slave) = Dict{String,Any}("Master Block ID" => master,
                                              "Slave Block ID" => slave,
                                              "Search Radius" => 0.01)
ut_cf_model(groups) = Dict{String,Any}("Type" => "Penalty Contact", "Contact Radius" => 0.005,
                                       "Contact Groups" => groups)

@testset "ut_check_valid_contact_model" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_dof(2)
    PeriLab.Data_Manager.set_num_controller(4)
    block_id = PeriLab.Data_Manager.create_constant_node_scalar_field("Block_Id", Int64)
    block_id .= 1
    contact = typed_contact(Dict("cm" => ut_cf_model(Dict("cg" => ut_cf_group(1, 2)))))
    @test_logs (:error,
                "Block defintion in slave does not exist.") @test_throws PeriLab.PeriLabError begin
        PeriLab.Solver_Manager.Model_Factory.Contact.check_valid_contact_model(contact,
                                                                               block_id)
    end
    contact = typed_contact(Dict("cm" => ut_cf_model(Dict("cg" => ut_cf_group(2, 1)))))
    @test_logs (:error,
                "Block defintion in master does not exist.") @test_throws PeriLab.PeriLabError begin
        PeriLab.Solver_Manager.Model_Factory.Contact.check_valid_contact_model(contact,
                                                                               block_id)
    end
    block_id[2] = 2
    contact = typed_contact(Dict("cm" => ut_cf_model(Dict("cg" => ut_cf_group(1, 2))),
                                 "cm2" => ut_cf_model(Dict("cg" => ut_cf_group(2, 1)))))
    @test_logs (:error,
                "Master and Slave should be defined in an inverse way, e.g. Master = 1, Slave = 2 in model 1 and Master = 2, Slave = 1 in model 2.") @test_throws PeriLab.PeriLabError begin
        PeriLab.Solver_Manager.Model_Factory.Contact.check_valid_contact_model(contact,
                                                                               block_id)
    end
    contact = typed_contact(Dict("cm" => ut_cf_model(Dict("cg" => ut_cf_group(1, 2)))))
    @test isnothing(PeriLab.Solver_Manager.Model_Factory.Contact.check_valid_contact_model(contact,
                                                                                          block_id))
end
```

(The old "needs a Slave", "master and slave are equal", "needs a Search Radius" and "Search Radius > 0" cases are now input errors, tested in Task 2 `contact group checks` and by `req`. `get_all_contact_blocks` is replaced by `contact_blocks`, tested in Task 2.)

Replace the body of `test/unit_tests/Models/Contact/ut_Penalty_model.jl` (the `contact_initialize_data` testset) with:

```julia
@testset "contact_initialize_data" begin
    PeriLab.Data_Manager.initialize_data()
    contact = typed_contact(Dict("C" => Dict{String,Any}("Type" => "Penalty Contact",
                                                         "Contact Radius" => 0.005,
                                                         "Contact Groups" => Dict{String,Any}("g" => Dict{String,Any}("Master Block ID" => 2,
                                                                                                                      "Slave Block ID" => 1,
                                                                                                                      "Search Radius" => 0.01)))))
    params = contact.models["C"]
    @test isnothing(PeriLab.Solver_Manager.Model_Factory.Contact.Penalty_Model.init_contact_model(params))
    @test params.contact_stiffness == 1e8
end
```

- [ ] **Step 2: Run them to verify they fail**

Run: the single-file command with `<PATH>` = `Models/Contact/ut_Contact_Factory`, then with `Models/Contact/ut_Penalty_model`
Expected: `MethodError` for `check_valid_contact_model(::ContactInput, …)`, and for Penalty either a `MethodError` from `haskey(::ContactModelParams, …)` or a `setindex!` error.

- [ ] **Step 3: Implement**

`src/Core/Data_manager.jl:194`: `data["Contact Properties"] = Dict()` becomes `data["Contact Properties"] = nothing`.

`src/Core/Data_manager/data_manager_contact.jl:166`:

```julia
function set_contact_properties(contact)
    data["Contact Properties"] = contact
end
```

`src/Models/Model_Factory.jl`:
- Add `using ...InputDeck: PeriLabInput, ContactInput` after the `using ...Helpers: …` block.
- Change the `init_models` signature (and its docstring line) to

```julia
function init_models(params::Dict,
                     input::PeriLabInput,
                     block_nodes::Dict{Int64,Vector{Int64}},
                     solver_options::Dict,
                     synchronise_field)
```

- Inside it, replace `check_contact(params)` with `check_contact(input.contact)`.
- Replace both `check_contact` methods with

```julia
check_contact(::Nothing) = nothing
check_contact(contact::ContactInput) = Contact.init_contact_model(contact)
check_contact(::Nothing, time::Float64, dt::Float64) = nothing
check_contact(contact::ContactInput, time::Float64, dt::Float64) = Contact.compute_contact_model(contact,
                                                                                                time,
                                                                                                dt)
```

`src/Core/Solver/Solver_manager.jl:126`:

```julia
    @timeit "init_models" init_models(params,
                                      input,
                                      block_nodes,
                                      solver_options,
                                      synchronise_field)
```

`src/Models/Contact/Contact_Factory.jl`:
- Add `using .....InputDeck: ContactInput, contact_blocks, contact_search_frequency` after the `using .....ModuleLoader` line.
- Replace `init_contact_model(params)` from its first line through the `"Set contact models"` loop with:

```julia
function init_contact_model(contact::ContactInput)
    @info "Init Contact Model"
    Data_Manager.set_contact_properties(contact)

    check_valid_contact_model(contact, Data_Manager.get_all_blocks())
    contact_block_list = contact_blocks(contact)
    # get all the contact block surface global ids and reduce the exchange positions to this points.
    # all following functions deal with the local contact position id.
    only_surface = contact.globals.only_surface_contact_nodes

    global_contact_ids = identify_contact_block_nodes(contact_block_list,
                                                      only_surface)
    contact_nodes = Data_Manager.create_constant_node_scalar_field("Contact Nodes", Int64)
    block_list = Data_Manager.get_all_blocks()

    mapping = contact_block_ids(global_contact_ids, block_list, contact_block_list)

    Data_Manager.set_contact_block_ids(mapping)
    # identify all surface which have no neighboring nodes
    if only_surface
        free_surfaces = identify_free_contact_surfaces(contact_block_list)
    end
    points = Data_Manager.get_all_positions()

    block_nodes = get_block_nodes(block_list, length(block_list)) # all ids

    for model in values(contact.models)
        for (cg, group) in pairs(model.contact_groups)
            Data_Manager.set_search_step(cg, 0)
            slave_id = group.slave_block_id
            master_id = group.master_block_id
            @info "Contact pair Master block $master_id - Slave block $slave_id"
            if !only_surface
                Data_Manager.set_free_contact_nodes(master_id, block_nodes[master_id])
                Data_Manager.set_free_contact_nodes(slave_id, block_nodes[slave_id])
                continue
            end

            compute_and_set_free_contact_nodes(points,
                                               block_nodes[slave_id],
                                               free_surfaces[slave_id],
                                               global_contact_ids,
                                               slave_id)
            compute_and_set_free_contact_nodes(points,
                                               block_nodes[master_id],
                                               free_surfaces[master_id],
                                               global_contact_ids,
                                               master_id)
        end
    end
```

  Keep the following lines unchanged: the block_list reduction, the shared volume/horizon fields, `create_local_contact_id_mapping` and the two `synchronize_contact_points` calls. Then replace the `"Set contact models"` loop with:

```julia
    @info "Set contact models"
    for (cm, model) in pairs(contact.models)
        init_contact_search(cm)

        mod = create_module_specifics(model.type,
                                      module_list,
                                      @__MODULE__,
                                      "contact_model_name")
        if isnothing(mod)
            @abort "No contact model of type " * model.type * " exists."
            return
        end
        Data_Manager.set_model_module(model.type, mod)
        mod.init_contact_model(model)
    end

    @info "Finish Init Contact Model"
end
```

- Delete `get_all_contact_blocks` (and its docstring if any).
- In `compute_contact_model`, change the signature to `function compute_contact_model(contact::ContactInput, time::Float64, dt::Float64)` and replace the `@timeit "Contact search"` block with:

```julia
    @timeit "Contact search" begin
        for model in values(contact.models)
            mod = Data_Manager.get_model_module(model.type)
            for (cg, group) in pairs(model.contact_groups)
                n = Data_Manager.get_search_step(cg) + 1
                Data_Manager.set_contact_dict(cg, Dict())

                @timeit "compute_contact_pairs" compute_contact_pairs(cg, group,
                                                                      model.contact_radius)
                @timeit "compute_contact_model" mod.compute_contact_model(cg,
                                                                          model,
                                                                          compute_master_force_density,
                                                                          compute_slave_force_density)
                if n == contact_search_frequency(group, contact.globals)
                    Data_Manager.set_search_step(cg, 0)
                end
            end
        end
    end
```

- Replace `check_valid_contact_model` with:

```julia
function check_valid_contact_model(contact::ContactInput, block_ids)
    # an inverse pair (1-2 in one group, 2-1 in another) is not allowed
    check_dict = Dict{Int64,Int64}()
    for model in values(contact.models), group in values(model.contact_groups)
        master = group.master_block_id
        slave = group.slave_block_id
        if !(master in block_ids)
            @abort "Block defintion in master does not exist."
            return
        end
        if !(slave in block_ids)
            @abort "Block defintion in slave does not exist."
            return
        end
        check_dict[master] = slave
        if get(check_dict, slave, nothing) == master
            @abort "Master and Slave should be defined in an inverse way, e.g. Master = 1, Slave = 2 in model 1 and Master = 2, Slave = 1 in model 2."
            return
        end
    end
    return nothing
end
```

`src/Models/Contact/Contact_search.jl`:
- Add `using .....InputDeck: ContactGroupParams` after the `using ....Helpers` block.
- Then:

```julia
function init_contact_search(cm::String)
    Data_Manager.set_search_step(cm, 0)
end

function global_contact_search(group::ContactGroupParams)
    all_positions = Data_Manager.get_all_positions()
    #-------------
    dof = Data_Manager.get_dof()
    # ids are exchange vector ids (''all position'' ids)
    master_nodes = Data_Manager.get_free_contact_nodes(group.master_block_id)
    slave_nodes = Data_Manager.get_free_contact_nodes(group.slave_block_id)

    near_points,
    no_pairs = find_potential_contact_pairs(dof,
                                            all_positions[master_nodes, :],
                                            all_positions[slave_nodes, :],
                                            group.search_radius)
```

  The rest of `global_contact_search` stays unchanged.
- `compute_contact_pairs(cg::String, contact_params::Dict)` becomes `compute_contact_pairs(cg::String, group::ContactGroupParams, contact_radius::Float64)`. Inside it, `global_contact_search(contact_params)` becomes `global_contact_search(group)`, and the `local_contact_search(contact_params, …)` call becomes `local_contact_search(contact_radius, Data_Manager.get_global_search_master_nodes(cg), Data_Manager.get_global_search_slave_nodes(cg))`.
- `local_contact_search(contact_params, master_nodes, slave_nodes)` becomes `local_contact_search(contact_radius::Float64, master_nodes, slave_nodes)`, and `contact_params["Contact Radius"]` inside it becomes `contact_radius`.

`src/Models/Contact/Penalty_model.jl`:
- Add `using .....InputDeck: ContactModelParams` after the Helpers import.
- Replace `init_contact_model(params)` with:

```julia
function init_contact_model(params::ContactModelParams)
    @info "Contact Stiffness $(params.contact_stiffness)"
    return nothing
end
```

- In `compute_contact_model`:
  - Change the signature to `function compute_contact_model(cg, params::ContactModelParams, compute_master_force_density::Function, compute_slave_force_density::Function)`.
  - `params["Contact Stiffness"]` becomes `params.contact_stiffness`.
  - `params["Contact Radius"]` becomes `params.contact_radius`.
  - The two `params["Symmetry"] == …` comparisons become `params.symmetry == …`.
  - `params["Friction Coefficient"]` becomes `params.friction_coefficient`.

Then run `grep -rn 'get_all_contact_blocks\|"Contact Properties"\|check_contact(params' src test` and fix any leftover caller. Expected: only `data_manager_contact.jl`'s getter/setter and the `Data_manager.jl` init line mention `"Contact Properties"`.

- [ ] **Step 4: Run to verify they pass**

Run: the single-file command with `<PATH>` = `Models/Contact/ut_Contact_Factory`, then `Models/Contact/ut_Penalty_model`
Expected: both pass.

Then run the contact fullscale tests:

```bash
cd test && julia --project=.. -e 'using Test, Logging, MPI; import PeriLab; Logging.disable_logging(Logging.Warn); include("helper.jl"); MPI.Init(); @testset "t" begin include("fullscale_tests/test_Penalty_Contact/test_Penalty_Contact.jl") end' 2>&1 | tail -15
```

Expected: pass (numerical results unchanged).

- [ ] **Step 5: Full suite, then commit**

Run the full suite in the background and read the tail. Expected: 0 failures, 0 errors.

```bash
git add src/Core/Data_manager.jl src/Core/Data_manager/data_manager_contact.jl \
        src/Models/Model_Factory.jl src/Core/Solver/Solver_manager.jl \
        src/Models/Contact/Contact_Factory.jl src/Models/Contact/Contact_search.jl \
        src/Models/Contact/Penalty_model.jl \
        test/unit_tests/Models/Contact/ut_Contact_Factory.jl \
        test/unit_tests/Models/Contact/ut_Penalty_model.jl
git commit -m "Contact reads typed ContactInput"
```

---

### Task 4: FEM and coupling consume `FEMParams`

**Files:**
- Modify: `src/Core/Data_manager.jl` (new `FEM Parameters` slot plus `set_fem_params`/`get_fem_params`, next to `fem_active`)
- Modify: `src/FEM/FEM_Factory.jl` (`init_FEM`, delete `valid_models`, `eval_FEM`)
- Modify: `src/FEM/FEM_basis.jl` (`compute_FEM` drops `params`; `get_polynomial_degree(degree::Vector{Int64}, dof)`)
- Modify: `src/FEM/Coupling/Coupling_Factory.jl`, `src/FEM/Coupling/Arlequin_coupling.jl`
- Modify: `src/FEM/Element_formulation/Lagrange_element.jl`, `src/FEM/FEM_template/FEM_template.jl` (`init_element` takes untyped `element_params`)
- Modify: `src/Models/Model_Factory.jl` (eval_FEM and compute_coupling calls)
- Modify: `src/Core/Solver/Solver_manager.jl:158-162`
- Test: `test/unit_tests/FEM/ut_FEM_Factory.jl`, `test/unit_tests/FEM/ut_FEM_routines.jl`, new testset in `test/unit_tests/FEM/Coupling/ut_Arlequin.jl`

**Interfaces:**
- Consumes (Task 2): `FEMParams`, `FEMCouplingParams`, `fem_degree`.
- Produces:
  - `Data_Manager.set_fem_params(fem)` / `Data_Manager.get_fem_params()`: `FEMParams` or `nothing` (default).
  - `FEM.init_FEM(fem::FEMParams, material_models::AbstractDict)`; `FEM.init_FEM(::Nothing, ::AbstractDict)` aborts with `"Invalid FEM parameters"`.
  - `FEM.eval_FEM(elements::AbstractVector{Int64}, time::Float64, dt::Float64)`.
  - `FEM_Basis.compute_FEM(elements, compute_stresses!::Function, time, dt)`.
  - `FEM_Basis.get_polynomial_degree(degree::Vector{Int64}, dof::Int64)::Vector{Int64}`.
  - `Coupling.init_coupling(nodes, fem::FEMParams)`, `Coupling.compute_coupling(fem::FEMParams)`.
  - Coupling module interface (`Arlequin_Coupling`): `init_coupling_model(nodes, fem::FEMParams)`, `compute_coupling(coupling::FEMCouplingParams)`.

- [ ] **Step 1: Write the failing tests**

In `test/unit_tests/FEM/ut_FEM_Factory.jl`:
- Delete the whole `@testset "ut_valid_models"` testset.
- Add after `PeriLab.Data_Manager.initialize_data()` at the top:

```julia
ut_fem(material) = typed_section(PeriLab.InputDeck.FEMParams,
                                 Dict{String,Any}("Degree" => 1, "Element Type" => "Lagrange",
                                                  "Material Model" => material))
```

- In `ut_init_FEM`, replace the two early abort checks with

```julia
    @test_logs (:error, "Invalid FEM parameters") @test_throws PeriLab.PeriLabError begin
        PeriLab.Solver_Manager.FEM.init_FEM(nothing, Dict{String,Any}())
    end
    @test_logs (:error,
                "The FEM material model b is not defined") @test_throws PeriLab.PeriLabError begin
        PeriLab.Solver_Manager.FEM.init_FEM(ut_fem("b"), Dict{String,Any}("a" => "a"))
    end
```

- Replace each later `params = Dict{String,Any}("FEM" => …, "Models" => Dict(… "Material Models" => M))` with `material_models = M` (keep `M` verbatim), and each `PeriLab.Solver_Manager.FEM.init_FEM(params)` with `PeriLab.Solver_Manager.FEM.init_FEM(ut_fem("Elastic Model"), material_models)`.
- After the successful `init_FEM` in `ut_init_FEM` add `@test PeriLab.Data_Manager.get_fem_params().element_type == "Lagrange"`.
- Replace both `eval_FEM(elements, PeriLab.Data_Manager.get_properties(1, "FEM"), 0.0, 1.0e-6)` calls with `eval_FEM(elements, 0.0, 1.0e-6)`.

In `test/unit_tests/FEM/ut_FEM_routines.jl`:
- Replace both `get_polynomial_degree(params["FEM"]["FE_1"], dof)` calls with `get_polynomial_degree([1], dof)`.
- Replace the body of `@testset "ut_get_polynomial_degree"` with

```julia
@testset "ut_get_polynomial_degree" begin
    gpd = PeriLab.Solver_Manager.FEM.FEM_Basis.get_polynomial_degree
    @test gpd([1], 2) == [1, 1]
    @test gpd([1], 3) == [1, 1, 1]
    @test gpd([2], 2) == [2, 2]
    @test gpd([2, 1, 1], 3) == [2, 1, 1]
    @test_logs (:error,
                "Degree must be defined with length one or number of dof.") @test_throws PeriLab.PeriLabError begin
        gpd([2, 3, 1], 2)
    end
end
```

Append to `test/unit_tests/FEM/Coupling/ut_Arlequin.jl`:

```julia
@testset "coupling defaults reach compute" begin
    fem = typed_section(PeriLab.InputDeck.FEMParams,
                        Dict{String,Any}("Degree" => 1, "Element Type" => "Lagrange",
                                         "Material Model" => "m",
                                         "Coupling" => Dict{String,Any}("Coupling Type" => "Arlequin")))
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_fem_params(fem)
    # compute reads the same immutable struct init saw; no init-time mutation needed
    @test PeriLab.Data_Manager.get_fem_params().coupling.pd_weight === 0.5
    @test PeriLab.Data_Manager.get_fem_params().coupling.kappa === 1.0
    @test hasmethod(PeriLab.Solver_Manager.FEM.Coupling.Arlequin_Coupling.compute_coupling,
                    Tuple{PeriLab.InputDeck.FEMCouplingParams})
end
```

- [ ] **Step 2: Run to verify they fail**

Run: the single-file command for `FEM/ut_FEM_Factory`, `FEM/ut_FEM_routines`, `FEM/Coupling/ut_Arlequin`
Expected:
- `MethodError` for `init_FEM(::Nothing, …)` / `init_FEM(::FEMParams, …)`, `get_polynomial_degree(::Vector{Int64}, …)` and `eval_FEM(…, ::Float64, ::Float64)`.
- `UndefVarError`/`MethodError` for `set_fem_params`.
- The `hasmethod` test is false.

- [ ] **Step 3: Implement**

`src/Core/Data_manager.jl`:
- In `initialize_data`, after `data["coupling_dict"] = Dict{Int64,Int64}()`, add `data["FEM Parameters"] = nothing`.
- After the `fem_active()` function add:

```julia
"""
	set_fem_params(fem)

Stores the typed `FEM` section (`FEMParams`, or `nothing` without FEM).
"""
function set_fem_params(fem)
    data["FEM Parameters"] = fem
end

"""
	get_fem_params()

The typed `FEM` section stored by `set_fem_params`.
"""
function get_fem_params()
    return data["FEM Parameters"]
end
```

- Add `export set_fem_params` and `export get_fem_params` next to `export fem_active`.

`src/FEM/FEM_basis.jl`:
- Replace `get_polynomial_degree(params::Dict{String,Any}, dof::Int64)` with

```julia
function get_polynomial_degree(degree::Vector{Int64}, dof::Int64)
    if length(degree) == 1
        return fill(degree[1], dof)
    elseif length(degree) == dof
        return copy(degree)
    end
    @abort "Degree must be defined with length one or number of dof."
end
```

- Change the `compute_FEM` signature to

```julia
function compute_FEM(elements::AbstractVector{Int64},
                     compute_stresses!::Function,
                     time::Float64,
                     dt::Float64)
```

  The body never used `params`.

`src/FEM/FEM_Factory.jl`:
- Add `using ...InputDeck: FEMParams, fem_degree` after `using ...Helpers: …`.
- Replace the beginning of `init_FEM` up to (but not including) the `"Material Gradient"` check with:

```julia
function init_FEM(::Nothing, material_models::AbstractDict)
    @abort "Invalid FEM parameters"
    return
end

function init_FEM(fem::FEMParams, material_models::AbstractDict)
    Data_Manager.set_fem_params(fem)
    if !haskey(material_models, fem.material_model)
        @abort "The FEM material model $(fem.material_model) is not defined"
        return
    end
```

- In the rest of `init_FEM`:
  - `p = get_polynomial_degree(params, dof)` becomes `p = get_polynomial_degree(fem_degree(fem), dof)`.
  - Both `params["Element Type"]` become `fem.element_type`.
  - The `init_element` argument tuple `(elements, params, p)` becomes `(elements, fem, p)`.
- Delete `valid_models`. `Material Model` is required by `FEMParams`, and other keys are unknown-key input errors.
- Replace `eval_FEM` with:

```julia
function eval_FEM(elements::AbstractVector{Int64},
                  time::Float64,
                  dt::Float64)
    return compute_FEM(elements,
                       compute_stresses!,
                       time,
                       dt)
end
```

`src/FEM/Element_formulation/Lagrange_element.jl` and `src/FEM/FEM_template/FEM_template.jl`: in `init_element`, `element_params::Dict` becomes `element_params`. In both docstrings, `element_params::Dict` becomes `element_params::FEMParams`.

`src/FEM/Coupling/Coupling_Factory.jl`:
- Add `using ....InputDeck: FEMParams` after the Data_Manager import.
- Replace `init_coupling` and `compute_coupling` with:

```julia
function init_coupling(nodes, fem::FEMParams)
    Data_Manager.create_constant_node_scalar_field("PD Nodes", Int64)
    fem.coupling === nothing && return
    coupling_model = fem.coupling.coupling_type

    mod = create_module_specifics(coupling_model, module_list,
                                  @__MODULE__, "coupling_name")
    if isnothing(mod)
        @abort "No coupling model of name " * coupling_model * " exists."
        return
    end
    Data_Manager.set_model_module(coupling_model, mod)

    ###TODO nodes and blocks
    mod.init_coupling_model(nodes, fem)
end

function compute_coupling(fem::FEMParams)
    fem.coupling === nothing && return
    mod = Data_Manager.get_model_module(fem.coupling.coupling_type)
    return mod.compute_coupling(fem.coupling)
end
```

`src/FEM/Coupling/Arlequin_coupling.jl`:
- Add `using ......InputDeck: FEMParams, FEMCouplingParams, fem_degree` after the Data_Manager import.
- In `init_coupling_model`:
  - The signature becomes `function init_coupling_model(nodes, fem::FEMParams)`.
  - `p = get_polynomial_degree(fe_params, dof)` becomes `p = get_polynomial_degree(fem_degree(fem), dof)`.
  - Delete the `if !haskey(fe_params["Coupling"], "PD Weight") … end` block and the `if !haskey(fe_params["Coupling"], "Kappa") … end` block.
  - Replace the `Coupling Block` branch with

```julia
    if fem.coupling.coupling_block !== nothing
        @info "Pre defined coupling region is block $(fem.coupling.coupling_block)"
        pd_nodes = get_coupling_zone(pd_nodes,
                                     Data_Manager.get_field("Block_Id"),
                                     fem.coupling.coupling_block)
    end
```

  - Replace `kappa = fe_params["Coupling"]["Kappa"]` and `weight_coefficient = fe_params["Coupling"]["PD Weight"]` with `kappa = fem.coupling.kappa` and `weight_coefficient = fem.coupling.pd_weight`.
- In `compute_coupling`: the signature becomes `function compute_coupling(coupling::FEMCouplingParams)`, and `weight_coefficient = fem_params["Coupling"]["PD Weight"]` becomes `weight_coefficient = coupling.pd_weight`.

`src/Models/Model_Factory.jl`:
- `FEM.eval_FEM(Vector{Int64}(1:nelements), Data_Manager.get_properties(1, "FEM"), time, dt)` becomes `FEM.eval_FEM(Vector{Int64}(1:nelements), time, dt)`.
- `FEM.Coupling.compute_coupling(Data_Manager.get_properties(1, "FEM"))` becomes `FEM.Coupling.compute_coupling(Data_Manager.get_fem_params())`.

`src/Core/Solver/Solver_manager.jl:158-162`:

```julia
    if Data_Manager.fem_active()
        @timeit "init_FEM" FEM.init_FEM(input.sections.fem,
                                        get(input.models, "Material Models",
                                            Dict{String,Any}()))
        @timeit "init_coupling" FEM.Coupling.init_coupling(1:Data_Manager.get_nnodes(),
                                                           input.sections.fem)
    end
```

Then run `grep -rn 'get_properties([^)]*"FEM"\|valid_models\|get_polynomial_degree(params\|set_properties("FEM"' src test`. Expected: no matches.

- [ ] **Step 4: Run to verify they pass**

Run: the single-file command for `FEM/ut_FEM_Factory`, `FEM/ut_FEM_routines`, `FEM/Coupling/ut_Arlequin`
Expected: all pass.

Then run the FEM fullscale tests. Find the runner file names with `ls test/fullscale_tests/test_FEM test/fullscale_tests/test_FEM_Coupling`, and include each `test_*.jl` there via the single-file pattern with `fullscale_tests/...`.
Expected: pass.

- [ ] **Step 5: Full suite, then commit**

Full suite in the background; read the tail. Expected: 0 failures, 0 errors.

```bash
git add src/Core/Data_manager.jl src/FEM src/Models/Model_Factory.jl \
        src/Core/Solver/Solver_manager.jl test/unit_tests/FEM
git commit -m "FEM and coupling read typed FEMParams"
```

---

### Task 5: Surface correction consumes `SurfaceCorrectionParams`

**Files:**
- Modify: `src/Core/Data_manager.jl` (new `Surface Correction` slot plus `set_surface_correction`/`get_surface_correction`)
- Modify: `src/Models/Surface_correction/Surface_correction.jl`
- Modify: `src/Models/Model_Factory.jl` (`init_surface_correction` call)
- Test: `test/unit_tests/Surface_correction/ut_Surface_correction.jl`

**Interfaces:**
- Consumes (Task 2): `SurfaceCorrectionParams`. Consumes (Task 3): the `init_models(params, input, …)` signature.
- Produces:
  - `Data_Manager.set_surface_correction(sc)` / `Data_Manager.get_surface_correction()`: `SurfaceCorrectionParams` or `nothing` (default).
  - `Surface_Correction.init_surface_correction(sc::Union{Nothing,SurfaceCorrectionParams}, local_synch, synchronise_field)`.
  - `compute_surface_correction(nodes, local_synch, synchronise_field)`: signature unchanged.

(Phase 3 moves this into the per-block `BlockModels`; this task only removes the dict.)

- [ ] **Step 1: Write the failing test**

Replace the body of `test/unit_tests/Surface_correction/ut_Surface_correction.jl`'s testset with:

```julia
@testset "ut_init_surface_correction" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_block_id_list([1])
    PeriLab.Data_Manager.init_properties()
    PeriLab.Data_Manager.set_dof(3)
    PeriLab.Data_Manager.set_num_controller(4)
    block_iD = PeriLab.Data_Manager.create_constant_node_scalar_field("Block_Id", Int64)
    block_iD .= 1
    mod_struct = PeriLab.Solver_Manager.Model_Factory

    # no section: nothing is stored and the compute step is a no-op
    @test isnothing(mod_struct.init_surface_correction(nothing, "local_synch",
                                                       "synchronise_field"))
    @test isnothing(PeriLab.Data_Manager.get_surface_correction())
    @test isnothing(mod_struct.compute_surface_correction([1, 2], "local_synch",
                                                          "synchronise_field"))

    sc = typed_section(PeriLab.InputDeck.SurfaceCorrectionParams,
                       Dict{String,Any}("Type" => "Volume Correction", "Update" => true))
    PeriLab.Data_Manager.set_surface_correction(sc)
    @test PeriLab.Data_Manager.get_surface_correction().update
end
```

- [ ] **Step 2: Run to verify it fails**

Run: the single-file command with `<PATH>` = `Surface_correction/ut_Surface_correction`
Expected: `MethodError` for `init_surface_correction(::Nothing, …)` or `UndefVarError` for `get_surface_correction`.

- [ ] **Step 3: Implement**

`src/Core/Data_manager.jl`:
- In `initialize_data`, next to `data["FEM Parameters"] = nothing`, add `data["Surface Correction"] = nothing`.
- Next to `set_fem_params` add:

```julia
"""
	set_surface_correction(sc)

Stores the typed `Surface Correction` section (or `nothing`).
"""
function set_surface_correction(sc)
    data["Surface Correction"] = sc
end

"""
	get_surface_correction()

The typed `Surface Correction` section stored by `set_surface_correction`.
"""
function get_surface_correction()
    return data["Surface Correction"]
end
```

- Add `export set_surface_correction` and `export get_surface_correction`.

`src/Models/Surface_correction/Surface_correction.jl`:
- Add `using .....InputDeck: SurfaceCorrectionParams` after the PeriLabExceptions import.
- Replace `compute_surface_correction` and `init_surface_correction` with:

```julia
function compute_surface_correction(nodes,
                                    local_synch,
                                    synchronise_field)
    sc = Data_Manager.get_surface_correction()
    sc === nothing && return
    # SurfaceCorrectionParams allows only "Volume Correction"
    if sc.update
        volumen_correction(nodes, local_synch, synchronise_field)
    end
    return compute_surface_volume_correction(nodes)
end
```

```julia
function init_surface_correction(sc::Union{Nothing,SurfaceCorrectionParams},
                                 local_synch,
                                 synchronise_field)
    Data_Manager.set_surface_correction(sc)
    sc === nothing && return
    return init_volumen_correction(local_synch, synchronise_field)
end
```

`src/Models/Model_Factory.jl`: `init_surface_correction(params, local_synch, synchronise_field)` becomes `init_surface_correction(input.sections.surface_correction, local_synch, synchronise_field)`.

Then run `grep -rn '"Surface Correction"' src`. Expected: the `local_fields_to_synch` key in `Data_manager.jl`, the new slot, the `set_local_synch`/`local_synch` calls in `Surface_correction.jl` and `parameter_handling.jl` (the legacy schema; phase 4 deletes it). No `get_properties`/`set_properties` with `"Surface Correction"` remain.

- [ ] **Step 4: Run to verify it passes**

Run: the single-file command with `<PATH>` = `Surface_correction/ut_Surface_correction`
Expected: pass.

Then run the surface correction fullscale test. Find the runner with `ls test/fullscale_tests/test_Surface_Correction` and include it via the single-file pattern.
Expected: pass.

- [ ] **Step 5: Full suite, then commit**

Full suite in the background; read the tail. Expected: 0 failures, 0 errors (the DCB surface-correction deck with `Update: True` is covered).

```bash
git add src/Core/Data_manager.jl src/Models/Surface_correction/Surface_correction.jl \
        src/Models/Model_Factory.jl test/unit_tests/Surface_correction/ut_Surface_correction.jl
git commit -m "Surface correction reads typed SurfaceCorrectionParams"
```
