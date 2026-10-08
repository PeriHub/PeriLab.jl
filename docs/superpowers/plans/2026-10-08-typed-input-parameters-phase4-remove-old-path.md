# Typed Input Parameters — Phase 4: Remove the Old Path

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** The legacy input path is gone. This covers:
- the module `Parameter_Handling` (all `parameter_handling*.jl` files);
- `validate_yaml` / `validate_models` and the dict schema files (`Schema.jl`, `PeriLabInputSchema.jl`, `apply_defaults.jl`, `example_usage.jl`);
- the old dependent-value helpers in `Helpers.jl`;
- the `Data_Manager` properties API;
- the raw deck dict that is still threaded through IO and the solver.

The per-category block slots fold into one `BlockModels` value per block.

**Architecture:**
- **Typed helpers move out first.** A few typed helpers currently live in `parameter_handling_*.jl`:
  - node sets, topology file and header go to IO;
  - output and compute helpers go to InputDeck.

  `validate_input` becomes a small IO function without the legacy model validation.
- **The raw deck stops travelling.** `IO.initialize_data`, `init_data`, `Solver_Manager.init`, `read_properties` and `init_models` lose their `params` argument. `contact_basis` asks `input.contact`.
- **Then the legacy code is deleted, file by file**, with its tests. Tests that compared typed against dict results keep only their typed assertions.
- **Finally `BlockModels`** (spec §3) replaces the `Block Materials`, `Block Damages` and generic `Block Models` slots. Code reads its fields directly (`Data_Manager.get_block_models(block).material`).

**Tech Stack:** Julia 1.12, `PeriLab.InputDeck`, `PeriLab.IO`, `PeriLab.Data_Manager`.

**Spec:** `docs/superpowers/specs/2026-10-03-typed-input-parameters-design.md`:
- §5 phase 4: "delete `parameter_handling*.jl`, `validate_yaml`, `get_properties` & co., `find_data_files`, old dependent-value helpers, `Schema.jl`, `PeriLabInputSchema.jl`, `example_usage.jl`, `apply_defaults.jl`";
- §3 Use (`BlockModels`);
- "Bridge during phases 2–3 … Removed in phase 4".

## Global Constraints

- **Numerical results unchanged.** The full suite stays green. Its count drops only by deleted legacy tests; every deletion is listed in a `Ruling:`.
- **Golden decks still parse in strict mode** (`ut_golden_decks.jl`).
- **Read fields directly; functions only for real logic.**
- **No new behaviour.** This phase only removes code and moves code. Any behaviour change a reviewer could notice gets a `Ruling:`.
- **Surface Correction** stays the global typed section `input.sections.surface_correction`. It is not part of `BlockModels`, because it is not per block. Record that as a `Ruling:` against the spec's `S` type parameter.
- **`PeriLabInput.models`** (the raw `Models` dict) stays. Material reads its model name from it and `block_typed_model` checks section presence. Removing it is out of scope.
- **The test-only `legacy_material_oracle.jl`** keeps working. It receives verbatim copies of any helper it imports that this phase deletes.
- **Memory.** The machine has ~8 GB. Before every full-suite run, check for orphaned MPI ranks (`PeriLab.run`, `hydra_pmi_proxy`, `mpiexec` in `/proc/*/cmdline`) and kill them with `kill -9`. Also delete stale `/dev/shm/mpich_shm_*` files.
- **Git.** Branch `feature/typed-input-parameters`; one commit per task; never merge. Read the full-suite result before committing. Add files by path, never `git add test` or `git add src` wholesale: the user keeps untracked files in `test/`.

## Review Focus

1. **A missing or non-YAML input file** still aborts with "`<file>` can not be found. Make sure the file exist and is readable." / "Not a supported filetype `<file>`". Pinned in Task 1 (`read_input_deck errors`).
2. **A deck with `Contact`** still sets the contact positions and blocks on every rank; one without does not. Pinned in Task 2 (`contact basis follows input.contact`).
3. **Node sets** given as id lists, ranges, coordinate expressions, `All` or a node-set file give the same sets as before the move. Pinned in Task 1 (`ut_node_sets_typed` keeps its file and expression cases against fixed expected sets).
4. **A block without any model** behaves as before:
   - it is skipped by every model category;
   - it is a 3D block for local damping;
   - it gets no pre-calculation.

   Pinned in Task 6 (`empty block models`).
5. **Pre-calculations added by `check_dependencies`** survive in the block's models and run in order. Pinned in Task 6 (`check_dependencies updates block models`).

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

### Task 1: Move the typed helpers out of `Parameter_Handling`; `validate_input` without legacy model validation

**Files:**
- Create: `src/IO/node_sets.jl` with `get_header`, `external_topology_file` and `read_node_sets`, moved verbatim from `src/Support/Parameters/parameter_handling_mesh.jl`. Include it in `src/IO/IO.jl` before `mesh_data.jl`.
- Modify: `src/Support/Parameters/Input/outputs.jl`. Move here, verbatim:
  - `check_for_duplicates`, `output_filenames`, `output_frequencies` and `output_fieldnames` from `parameter_handling_output.jl`;
  - `compute_names` and `active_computes` from `parameter_handling_computes.jl`.

  Export them from `InputDeck.jl` if it lists exports.
- Modify: `src/IO/read_inputdeck.jl`:
  - `validate_input(params; directory, no_strict)` becomes an IO function: `read_input` + `report!`, without `validate_models`;
  - delete `read_input_file`;
  - `read_input_deck` keeps its aborts.
- Modify: `src/IO/IO.jl` and `src/IO/mesh_data.jl` (imports)
- Delete from `Parameter_Handling` the moved functions, plus `validate_yaml` and `validate_input`. Delete their export lines too.
- Tests:
  - `test/unit_tests/Support/Parameters/Input/ut_node_sets_typed.jl` and `ut_outputs_typed.jl`: point `PH.read_node_sets` etc. to the new homes;
  - `ut_validate_yaml.jl`: keep the `validate_input` / `read_input_deck` testsets, pointed to `PeriLab.IO`;
  - `test/unit_tests/IO/ut_read_inputdeck.jl`: `read_input_file` tests become `read_input_deck` tests.

**Interfaces:**
- Produces:
  - `IO.get_header(filename)`;
  - `IO.external_topology_file(d::DiscretizationParams, path)`;
  - `IO.read_node_sets(d::DiscretizationParams, path, mesh_df)`;
  - `InputDeck.check_for_duplicates`, `output_filenames`, `output_frequencies`, `output_fieldnames`, `compute_names`, `active_computes`;
  - `IO.validate_input(params; directory = "", no_strict = false) -> (deck, input)`;
  - `IO.read_input_deck(filename; directory, no_strict) -> (deck, input)` (unchanged).

- [ ] **Step 1: Point the tests to the new homes (they fail)**

In `ut_node_sets_typed.jl`, use `PeriLab.IO.read_node_sets` and `PeriLab.IO.external_topology_file` instead of `PH.read_node_sets` / `PH.external_topology_file`. Keep the `PH.get_node_sets` reference comparison until Task 3. In that file's comparison, replace the reference by fixed expected sets now, so it survives Task 3. Print each `PH.get_node_sets` result once and write it down as a literal `Dict`.

In `ut_outputs_typed.jl`, use `PeriLab.InputDeck.output_filenames`, `output_frequencies`, `output_fieldnames`, `compute_names` and `active_computes`.

In `ut_validate_yaml.jl`, `PeriLab.Parameter_Handling.validate_input` becomes `PeriLab.IO.validate_input`.

In `ut_read_inputdeck.jl`, the testsets calling `PeriLab.IO.read_input_file`:
- the missing-file and wrong-type aborts call `PeriLab.IO.read_input_deck` with the same expected messages;
- the success case asserts on the first element of `read_input_deck`'s result, the deck dict.

Add to `ut_read_inputdeck.jl`:

```julia
@testset "read_input_deck errors" begin
    @test_logs (:error,
                "missing.yaml can not be found. Make sure the file exist and is readable.") @test_throws PeriLab.PeriLabError PeriLab.IO.read_input_deck("missing.yaml")
    file = joinpath(mktempdir(), "deck.txt")
    write(file, "x")
    @test_logs (:error, "Not a supported filetype $file") @test_throws PeriLab.PeriLabError PeriLab.IO.read_input_deck(file)
end
```

Check the exact abort texts in `read_input_deck` (`src/IO/read_inputdeck.jl`) and copy them.

- [ ] **Step 2: Run to verify they fail**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Support/Parameters/Input/ut_node_sets_typed unit_tests/Support/Parameters/Input/ut_outputs_typed unit_tests/Support/Parameters/Input/ut_validate_yaml unit_tests/IO/ut_read_inputdeck`
Expected: `UndefVarError`s for the new homes (`read_node_sets not defined in PeriLab.IO`, `output_filenames not defined in PeriLab.InputDeck`, `validate_input not defined in PeriLab.IO`). "read_input_deck errors" passes already: it pins the existing aborts.

- [ ] **Step 3: Move**

- Cut the listed functions with their docstrings and paste them into the new files.
- `node_sets.jl` needs `CSV`, `DataFrames`, `Exodus` and `@abort`; IO already uses them. Check its `using` lines and add what is missing.
- `read_inputdeck.jl`:

```julia
using ..ParameterSpec: report!, strict_mode
using ..InputDeck: read_input as read_typed_input

"""
    validate_input(params; directory = "", no_strict = false) -> (deck, input)

Validates a loaded input deck against the typed input declarations. Reports
every problem at once and aborts if there is an error. Returns the
`params["PeriLab"]` dict and the typed `PeriLabInput`.
"""
function validate_input(params::Dict; directory::AbstractString = "", no_strict::Bool = false)
    if !haskey(params, "PeriLab") || !(params["PeriLab"] isa AbstractDict) ||
       length(params["PeriLab"]) < 2
        @abort "Yaml file is not valid."
        return
    end
    deck = params["PeriLab"]
    input, ctx = read_typed_input(deck, directory;
                                  strict = strict_mode(deck; no_strict_flag = no_strict))
    report!(ctx)
    return deck, input
end
```

  Use the dot counts of the file's existing `using ...Parameter_Handling` line for `ParameterSpec` and `InputDeck`. IO's own `read_input(filename)` keeps its name, hence the `as` import.
- Delete `read_input_file` and its export.
- In `IO.jl` and `mesh_data.jl`, replace the `Parameter_Handling` imports by the new homes.

`validate_models` is no longer called anywhere. Its legacy "Key not known" warnings for model keys are replaced by the typed strict checks of phase 3. Record a `Ruling:`.

- [ ] **Step 4: Run to verify they pass**

Run the Step 2 command, plus `unit_tests/Support/Parameters/Input/ut_golden_decks`.
Expected: all pass.

- [ ] **Step 5: Full suite, then commit**

```bash
git add src/IO src/Support/Parameters test/unit_tests/Support/Parameters/Input test/unit_tests/IO
git commit -m "Typed helpers leave Parameter_Handling; validate_input without the legacy model validation"
```

---

### Task 2: The raw deck no longer travels through IO and the solver

**Files:**
- Modify: `src/IO/IO.jl` (`initialize_data` returns `(input, steps)`)
- Modify: `src/IO/mesh_data.jl`:
  - `init_data(input, path, comm)` returns nothing;
  - `contact_basis(input::PeriLabInput, mesh, comm, rank)`.
- Modify: `src/PeriLab.jl` (call sites; drop `using .Parameter_Handling: get_initial_time`)
- Modify: `src/Core/Solver/Solver_manager.jl` (`init(input, step_id)`)
- Modify: `src/Models/Model_Factory.jl` (`read_properties(input, material_model)`, `init_models(input, block_nodes, solver_options, synchronise_field)`)
- Tests: every caller listed by `grep -rn 'read_properties(\|init_models(\|Solver_Manager.init(\|IO.initialize_data(\|init_data(' test --include=*.jl`

**Interfaces:**
- Produces:
  - `IO.initialize_data(filename, filedirectory, comm; no_strict) -> (input, steps)`;
  - `IO.init_data(input, path, comm)`;
  - `Solver_Manager.init(input, step_id)`;
  - `Model_Factory.read_properties(input, material_model::Bool)`;
  - `Model_Factory.init_models(input, block_nodes, solver_options, synchronise_field)`.

- [ ] **Step 1: Write the failing tests**

Change every test call of `read_properties(params, input, flag)` to `read_properties(input, flag)`, and drop the now unused `params = …` lines next to them. Add to `test/unit_tests/IO/ut_IO.jl` (or the file that tests `mesh_data.jl`; check with `grep -rln 'contact_basis\|init_data' test/unit_tests/IO`):

```julia
@testset "contact basis follows input.contact" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_dof(2)
    comm = MPI.COMM_WORLD
    mesh = PeriLab.IO.DataFrame(x = [0.0, 1.0], y = [0.0, 0.0], block_id = [1, 2])
    no_contact = typed_input(Dict())
    PeriLab.IO.contact_basis(no_contact, mesh, comm, 1)
    @test PeriLab.Data_Manager.get_all_positions() === nothing
    with_contact = typed_input(Dict("Contact" => Dict("Contact Group" => Dict("Contact Model" => "Penalty Contact",
                                                                               "Master Block ID" => 1,
                                                                               "Slave Block ID" => 2,
                                                                               "Search Radius" => 1.0,
                                                                               "Contact Radius" => 0.5))))
    PeriLab.IO.contact_basis(with_contact, mesh, comm, 1)
    @test PeriLab.Data_Manager.get_all_positions() == [0.0 0.0; 1.0 0.0]
    @test PeriLab.Data_Manager.get_all_blocks() == [1, 2]
end
```

- Check the getter names for the positions and blocks the function sets (`grep -n 'set_all_positions\|get_all_positions\|set_all_blocks' src/Core/Data_manager*.jl src/Core/Data_manager/*.jl`) and their default before `contact_basis` runs.
- Build `with_contact` from an existing typed contact test deck (`grep -rn 'typed_contact\|Contact Model' test/unit_tests | head`), so the keys are valid.
- If `DataFrame` is not reachable as `PeriLab.IO.DataFrame`, use `using DataFrames` in the test file as the other IO tests do.

- [ ] **Step 2: Run to verify they fail**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/IO/ut_IO unit_tests/Models/ut_Model_Factory unit_tests/Models/ut_block_models`
Expected: `MethodError`s for `contact_basis(::PeriLabInput, …)` and `read_properties(::PeriLabInput, ::Bool)`.

- [ ] **Step 3: Implement**

- `contact_basis`: the first lines become `input.contact === nothing && return`. Change the signature as listed.
- `init_data`: drop `params`. It passes `input` to `contact_basis` and returns nothing. In `IO.initialize_data`:
  - `deck, input = read_input_deck(...)` becomes `_, input = read_input_deck(...)`;
  - call `init_data(input, filedirectory, comm)`;
  - `return input, steps`.

  Update the docstrings.
- `PeriLab.jl`:
  - `params, input, steps = IO.initialize_data(...)` becomes `input, steps = IO.initialize_data(...)`;
  - `Solver_Manager.init(params, input, step_id)` becomes `Solver_Manager.init(input, step_id)`;
  - remove `using .Parameter_Handling: get_initial_time` (unused).
- `Solver_manager.jl`: `init(input, step_id)`; its calls become `read_properties(input, …)` and `init_models(input, …)`. Update the docstring.
- `Model_Factory.jl`: drop the `params` argument of `read_properties` and `init_models`, and their docstring lines.

Run `grep -rn 'params' src/IO/mesh_data.jl src/IO/IO.jl src/Core/Solver/Solver_manager.jl src/PeriLab.jl | grep -v 'solver_params\|step_solver_params\|block_params\|params::.*Params\|PARAMS'`.
Expected: no remaining raw-deck uses.

- [ ] **Step 4: Run to verify they pass**

Run the Step 2 command.
Expected: all pass.

- [ ] **Step 5: Full suite, then commit**

```bash
git add src/IO src/PeriLab.jl src/Core/Solver/Solver_manager.jl src/Models/Model_Factory.jl test/unit_tests/IO test/unit_tests/Models
git commit -m "The raw input deck no longer travels through IO and the solver"
```

---

### Task 3: Delete `Parameter_Handling` and the dict schema files

**Files:**
- Delete:
  - `src/Support/Parameters/parameter_handling.jl` and every `parameter_handling_*.jl`;
  - `src/Support/Parameters/Schema.jl`, `PeriLabInputSchema.jl`, `apply_defaults.jl`, `example_usage.jl`.
- Modify: `src/PeriLab.jl` (remove `include("./Support/Parameters/parameter_handling.jl")` and any `Parameter_Handling` import)
- Delete tests: `test/unit_tests/Support/Parameters/ut_parameter_handling.jl` and its include in `runtests.jl`. Its `ut_check_for_duplicates` testset moves first into `ut_outputs_typed.jl`, calling `PeriLab.InputDeck.check_for_duplicates`.
- Modify tests:
  - `ut_solver_typed.jl`, `ut_blocks_mesh_typed.jl`, `ut_outputs_typed.jl`, `ut_bc_typed.jl`, `ut_node_sets_typed.jl`: delete every assertion against `PH.*`. A testset with only `PH` comparisons is deleted, e.g. the "… equal the Dict getters on every shipped deck" testsets; the golden-deck test covers parsing every deck.
  - `ut_validate_yaml.jl`: delete the testsets that call `validate_yaml`, and "legacy model validation still applies to the other categories".
  - Remove now unused `const PH = …` lines and test data files only `ut_parameter_handling.jl` used (`test/unit_tests/Support/Parameters/test_data_file.txt`, if nothing else reads it: `grep -rn test_data_file test`).

**Interfaces:**
- Consumes (Tasks 1–2): nothing in `src` or in the kept tests references `Parameter_Handling`.

- [ ] **Step 1: Prove nothing depends on the module any more**

Run `grep -rn 'Parameter_Handling\|PH\.\|validate_yaml\|validate_models\|Schema\b\|apply_defaults' src test --include=*.jl | grep -v '^src/Support/Parameters/parameter_handling\|^src/Support/Parameters/Schema.jl\|^src/Support/Parameters/PeriLabInputSchema.jl\|^src/Support/Parameters/apply_defaults.jl\|^src/Support/Parameters/example_usage.jl'`.
Expected:
- only the test references this task deletes or rewrites;
- `src/PeriLab.jl`'s include line;
- no other `src` hit.

- [ ] **Step 2: Delete and rewrite**

Delete the files. Remove the include and imports. Rewrite the tests as listed. Ledger one `Ruling:` that lists every deleted testset and why each is covered (legacy-only function removed, or covered by its typed counterpart / the golden decks).

- [ ] **Step 3: Run to verify**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Support/Parameters/Input/input_tests`
Expected: all pass. The file includes every input test, the golden decks among them.

- [ ] **Step 4: Full suite, then commit**

```bash
git add -A src/Support/Parameters src/PeriLab.jl test/unit_tests/Support/Parameters test/runtests.jl
git commit -m "Delete Parameter_Handling and the dict schema files"
```

`git add -A` is limited to those paths, so it records the deletions without touching the user's untracked files elsewhere in `test/`.

---

### Task 4: Delete the old dependent-value helpers

**Files:**
- Modify: `src/Support/Helpers.jl`: delete these, with their exports and docstrings:
  - `is_dependent`, `interpolation`, `interpol_data`;
  - `AbstractDependentValue`, `ConstantValue`, `InterpolatedValue` and their call methods;
  - `get_dependent_value`, `get_dependent_value_with_ID`.
- Modify: `src/Models/Material/Material_Basis.jl`: delete the third `apply_pointwise_E(nodes, bond_force, dependent_field)` method (it references an undefined `damage_parameter` and nothing calls it). Remove `interpol_data` from the `Helpers` import.
- Modify: `test/unit_tests/Models/Material/legacy_material_oracle.jl`: replace the import of `get_dependent_value_with_ID` by verbatim copies of the deleted chain it needs (`get_dependent_value_with_ID`, `get_dependent_value`, `is_dependent`, `AbstractDependentValue`, `ConstantValue`, `InterpolatedValue` with their call methods, `interpol_data`). Copy them before deleting.
- Modify: `test/unit_tests/Support/ut_helpers.jl`: delete the testsets of the deleted helpers.

- [ ] **Step 1: Copy the chain into the oracle (oracle tests still pass)**

Extract the functions from `Helpers.jl` verbatim into the oracle module. They need `Dierckx: evaluate` for `interpol_data`; add `using Dierckx: evaluate` if the oracle module does not see it. Drop the oracle's `using PeriLab.Helpers: get_dependent_value_with_ID`.

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Models/Material/ut_block_material`
Expected: all pass, now through the oracle's own copies.

- [ ] **Step 2: Delete**

Delete the helpers, the `apply_pointwise_E` method and the `ut_helpers.jl` testsets.

Run `grep -rnw 'is_dependent\|interpolation\|interpol_data\|ConstantValue\|InterpolatedValue\|AbstractDependentValue\|get_dependent_value\|get_dependent_value_with_ID' src test --include=*.jl | grep -v legacy_material_oracle`.
Expected: no matches. Docstring mentions of the word "interpolation" in `Lagrange_element.jl` / `Gcode_Mesh.jl` are prose; make sure they are not function calls.

Ledger the deleted `ut_helpers.jl` testsets in a `Ruling:`.

- [ ] **Step 3: Run to verify**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Support/ut_helpers unit_tests/Models/Material/ut_block_material unit_tests/Models/Material/ut_material_basis`
Expected: all pass.

- [ ] **Step 4: Full suite, then commit**

```bash
git add src/Support/Helpers.jl src/Models/Material/Material_Basis.jl test/unit_tests/Models/Material/legacy_material_oracle.jl test/unit_tests/Support/ut_helpers.jl
git commit -m "Delete the old dependent-value helpers"
```

---

### Task 5: Delete the `Data_Manager` properties API

**Files:**
- Modify: `src/Core/Data_manager.jl`:
  - delete `get_properties`, `get_property`, `set_property` (both methods), `set_properties`, `check_property` and `init_properties`, with their exports and docstrings;
  - delete `data["properties"]` in `initialize_data`.
- Modify: the commented line in `src/Core/Solver/Matrix_linear_static.jl` that mentions `get_property`: delete it.
- Modify tests:
  - `test/unit_tests/Core/ut_data_manager.jl`: delete the testsets / assertions on the deleted functions;
  - remove the `PeriLab.Data_Manager.init_properties()` lines in `ut_Surface_correction.jl`, `ut_FEM_Factory.jl`, `ut_block_material.jl` and `ut_Material_Factory.jl`;
  - `ut_Model_Factory.jl`: the assertions `!haskey(PeriLab.Data_Manager.data["properties"], 1)` are deleted (there is no properties store).

- [ ] **Step 1: Prove no `src` user remains**

Run `grep -rn 'get_properties\|get_property\|set_property\|set_properties\|check_property\|init_properties\|"properties"' src --include=*.jl | grep -v 'src/Core/Data_manager.jl'`.
Expected: only the `Matrix_linear_static.jl` comment.

- [ ] **Step 2: Delete**

Delete the functions, the store and the test lines as listed. Ledger the deleted `ut_data_manager.jl` testsets in a `Ruling:`.

- [ ] **Step 3: Run to verify**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Core/ut_data_manager unit_tests/Models/ut_Model_Factory unit_tests/Surface_correction/ut_Surface_correction unit_tests/FEM/ut_FEM_Factory unit_tests/Models/Material/ut_Material_Factory unit_tests/Models/Material/ut_block_material`
Expected: all pass.

- [ ] **Step 4: Full suite, then commit**

```bash
git add src/Core src/Core/Solver/Matrix_linear_static.jl test/unit_tests/Core test/unit_tests/Models test/unit_tests/Surface_correction test/unit_tests/FEM
git commit -m "Delete the Data_Manager properties API"
```

---

### Task 6: `BlockModels` — one typed value per block

**Files:**
- Modify: `src/Core/Data_manager.jl`:
  - add `BlockModels` and `set_block_models` / `get_block_models`;
  - delete `set_block_material`, `get_block_material`, `set_block_damage`, `get_block_damage`, `set_block_model` and `get_block_model`;
  - delete the slots `Block Materials` and `Block Damages`;
  - `Block Models` becomes `Dict{Int64,BlockModels}`.
- Modify: `src/Core/Data_manager/data_manager_checkpoint.jl` (excluded keys: only `"Block Models"`)
- Modify: every `src` caller (`grep -rn 'get_block_material\|set_block_material\|get_block_damage\|set_block_damage\|get_block_model\b\|set_block_model\b' src --include=*.jl`):
  - `Model_Factory.jl` (`read_properties` builds one `BlockModels` per block; `typed_block_model`; `block_local_damping`; `local_damping_symmetry`; `block_thermal_conductivity`);
  - `Material_Factory.jl`, `Damage_Factory.jl`, `Thermal_Factory.jl`, `Additive_Factory.jl`, `Degradation_Factory.jl`, `Pre_Calculation_Factory.jl` (`check_dependencies`, `init_model`, `fields_for_local_synchronization`);
  - `compute_field_values.jl`, `Matrix_Verlet.jl`, `Matrix_linear_static.jl`, `Newmark.jl`, `Correspondence_matrix_based.jl`.
- Modify tests: every caller in `test/` (same grep).

**Interfaces:**
- Produces, in `Data_Manager`:

```julia
"""
    BlockModels

The typed models of one block, `nothing` where the block has none of a category.
`pre_calculation` lists the active pre-calculations in run order (empty: none).
Surface correction is global (`input.sections.surface_correction`), not per block.
"""
struct BlockModels{M,D,T,A,G}
    material::M
    damage::D
    thermal::T
    additive::A
    degradation::G
    pre_calculation::Vector{String}
end

function BlockModels(; material = nothing, damage = nothing, thermal = nothing,
                     additive = nothing, degradation = nothing,
                     pre_calculation::Vector{String} = String[])
    return BlockModels(material, damage, thermal, additive, degradation, pre_calculation)
end

# a block without models; never mutated
const NO_BLOCK_MODELS = BlockModels()

set_block_models(block::Int64, models::BlockModels) = (data["Block Models"][block] = models)
"The typed models of a block (`NO_BLOCK_MODELS` if it has none)."
get_block_models(block::Int64) = get(data["Block Models"], block, NO_BLOCK_MODELS)
```

  Export `BlockModels`, `set_block_models` and `get_block_models`.

- [ ] **Step 1: Write the failing tests**

Create `test/unit_tests/Core/ut_block_models_store.jl` and include it in `runtests.jl` next to `ut_data_manager.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

@testset "empty block models" begin
    PeriLab.Data_Manager.initialize_data()
    m = PeriLab.Data_Manager.get_block_models(1)
    @test m.material === nothing && m.damage === nothing && m.thermal === nothing
    @test m.additive === nothing && m.degradation === nothing
    @test isempty(m.pre_calculation)
    MF = PeriLab.Solver_Manager.Model_Factory
    for name in ("Material Model", "Damage Model", "Thermal Model", "Additive Model",
                 "Degradation Model", "Pre Calculation Model")
        @test !MF.has_block_model(1, name)
    end
    @test MF.local_damping_symmetry(1) == "3D"
    @test MF.block_local_damping(1) === nothing
    @test MF.block_thermal_conductivity(1) === nothing
end

@testset "check_dependencies updates block models" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(2)
    PeriLab.Data_Manager.set_dof(3)
    m = typed_block_material(Dict("Material Model" => "Correspondence Elastic",
                                  "Symmetry" => "isotropic", "Bulk Modulus" => 1.0,
                                  "Shear Modulus" => 1.0))
    th = typed_model(:thermal, Dict("Thermal Model" => "Heat Transfer",
                                    "Heat Transfer Coefficient" => 1.0,
                                    "Environmental Temperature" => 30);
                     name_key = "Thermal Model")
    PeriLab.Data_Manager.set_block_models(1,
                                          PeriLab.Data_Manager.BlockModels(material = m,
                                                                           thermal = th,
                                                                           pre_calculation = ["Shape Tensor"]))
    PeriLab.Solver_Manager.Model_Factory.Pre_Calculation.check_dependencies(Dict(1 => [1, 2]))
    b = PeriLab.Data_Manager.get_block_models(1)
    @test b.pre_calculation == ["Deformed Bond Geometry", "Shape Tensor", "Deformation Gradient"]
    @test b.material === m && b.thermal === th             # other parts kept
end
```

Rewrite the existing tests that use the deleted accessors:
- `set_block_material(b, m)` becomes `set_block_models(b, BlockModels(material = m))`;
- `get_block_material(b)` becomes `get_block_models(b).material`;
- the same for damage, and `get_block_model("Thermal Model", b)` becomes `get_block_models(b).thermal`, etc.;
- the "generic block model slot" testset in `ut_block_models.jl` is replaced by "empty block models" above (ledger it).

A test that sets two categories on the same block builds one `BlockModels` with both.

- [ ] **Step 2: Run to verify they fail**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Core/ut_block_models_store unit_tests/Models/ut_block_models`
Expected: `UndefVarError: get_block_models` / `BlockModels`.

- [ ] **Step 3: Implement**

`Data_Manager`: add the Interfaces code, remove the six old accessors and the two old slots, and set `data["Block Models"] = Dict{Int64,BlockModels}()`.

`Model_Factory.jl`:

```julia
const _BLOCK_MODEL_FIELDS = Dict("Material Model" => :material, "Damage Model" => :damage,
                                 "Thermal Model" => :thermal, "Additive Model" => :additive,
                                 "Degradation Model" => :degradation)

# the typed model of category `name` of a block, or `nothing`
function typed_block_model(block::Int64, name::String)
    models = Data_Manager.get_block_models(block)
    if name == "Pre Calculation Model"
        return isempty(models.pre_calculation) ? nothing : models.pre_calculation
    end
    return getfield(models, _BLOCK_MODEL_FIELDS[name])
end
```

`read_properties` collects the block's parts in locals: the material (as before, when `material_model`), the damage, additive, degradation, thermal and the pre-calculation names, each from `block_typed_model`. It then calls `Data_Manager.set_block_models(block, BlockModels(material = …, damage = …, thermal = …, additive = …, degradation = …, pre_calculation = names))` once per block, in one loop over the blocks.

`check_dependencies`:
- reads `models = Data_Manager.get_block_models(block_id)` and starts from `Set{String}(models.pre_calculation)`;
- stores `Data_Manager.set_block_models(block_id, BlockModels(models.material, models.damage, models.thermal, models.additive, models.degradation, order_pre_calculations(collect(names))))`.

Every other caller reads the field directly:

| Old | New |
|---|---|
| `Data_Manager.get_block_material(b)` | `Data_Manager.get_block_models(b).material` |
| `Data_Manager.get_block_damage(b)` | `Data_Manager.get_block_models(b).damage` |
| `Data_Manager.get_block_model("Thermal Model", b)` | `Data_Manager.get_block_models(b).thermal` |
| `Data_Manager.get_block_model("Additive Model", b)` | `Data_Manager.get_block_models(b).additive` |
| `Data_Manager.get_block_model("Degradation Model", b)` | `Data_Manager.get_block_models(b).degradation` |
| `Data_Manager.get_block_model("Pre Calculation Model", b)` | `Data_Manager.get_block_models(b).pre_calculation` |

Afterwards run `grep -rn 'get_block_material\|set_block_material\|get_block_damage\|set_block_damage\|get_block_model\b\|set_block_model\b\|Block Materials\|Block Damages' src test --include=*.jl`.
Expected: no matches.

- [ ] **Step 4: Run to verify they pass**

Run the Step 2 command, plus:
- `unit_tests/Models/ut_Model_Factory`, `unit_tests/Models/Material/ut_block_material`, `unit_tests/Models/Damage/ut_block_damage`;
- `unit_tests/Models/Thermal/ut_Thermal_Factory`, `unit_tests/Core/ut_data_manager`.

Expected: all pass.

- [ ] **Step 5: Full suite, then commit**

```bash
git add src test/unit_tests/Core test/unit_tests/Models test/runtests.jl
git commit -m "BlockModels: one typed value per block replaces the per-category slots"
```

Here `git add src` is fine: the user's untracked files are under `test/`.
