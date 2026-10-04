# Typed Input Parameters — Phase 2b-2 (Mesh, Discretization and Blocks Consume Typed Input) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Block setup, orientations, the block summary, the influence function, mesh import (including the format importers), neighbourhood search, surface extrusion, external topology, node sets and bond filters read the typed `BlockParams` / `DiscretizationParams` instead of the params `Dict`.

**Architecture:** Same rule as 2b-1: read fields directly; functions only for real logic. New logic helpers live next to the types in `InputDeck` (`block_by_id`, `block_angles`, `block_names_and_ids`, `mesh_scaling`, `gcode_block_ids`) or, where they need file / DataFrame / Exodus access, in `parameter_handling_mesh.jl` (`read_node_sets`, `external_topology_file`). Each consumer takes the smallest typed piece it needs (a `Dict{String,BlockParams}`, a `SurfaceExtrusionParams`, a `Bool`, a `String`); the mesh importers — an extension point — receive the whole `PeriLabInput`. The old `Dict` getters stay (deleted in phase 4) and are the reference in an equivalence test over every shipped deck.

**Tech Stack:** Julia 1.12, existing dependencies only.

**Spec:** `docs/superpowers/specs/2026-10-03-typed-input-parameters-design.md` (§3 step 4, phase 2 of §5). Previous slice: `docs/superpowers/plans/2026-10-03-typed-input-parameters-phase2b1-solver.md` (typed input already reaches `Solver_Manager.init` and `IO.initialize_data`). Remaining slices after this one: 2b-3 outputs and compute classes, 2b-4 boundary conditions, 2b-5 contact, FEM / coupling, surface correction.

## Global Constraints

- Julia `1.12`; no new packages.
- Branch `feature/typed-input-parameters`; one commit per task, message ending with the session's `Co-Authored-By` line; commit identity `-c user.name="Jan-Timo Hesse" -c user.email="jthesse@gmx.de"`. Never merge.
- Simulation results must not change: every fullscale test must still pass.
- The `Dict` getters in `parameter_handling_blocks.jl` / `parameter_handling_mesh.jl` stay unchanged (phase 4 deletes them); no typed methods are added to them.
- No thin wrappers around single fields.
- Relative imports: write the `InputDeck` import with the same number of dots as the file's existing `PeriLabExceptions` import (if the file has none, its `Data_Manager` import). Do not copy the dot count of `Parameter_Handling` imports — some of those resolve through re-exports.
- Every new source file starts with the SPDX header (`# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>` / `#` / `# SPDX-License-Identifier: BSD-3-Clause`).
- Test commands (from the repository root):
  - Input tests: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl`
  - Single existing test file, as `test/runtests.jl` runs it (cwd `test/`, warnings disabled, helper loaded):
    `cd test && julia --project=.. -e 'using Test, Logging, MPI; import PeriLab; Logging.disable_logging(Logging.Warn); include("helper.jl"); MPI.Init(); @testset "t" begin include("unit_tests/<path>.jl") end'`
  - Full suite (≈35 min, background, output to a file): `julia --project=. -e 'using Pkg; Pkg.test()'`

## Review Focus

1. A deck whose blocks have `Angle X` but not `Angle Y` / `Angle Z` must work in 2D and abort with "Angle Y of <block> is not defined" in 3D, as before — test in Task 1 ("block angles").
2. Gcode `Blocks` keys are block IDs written as YAML integers (`2: "z > 1"`); the importer must still assign integer block IDs — test in Task 1 ("gcode block ids") and Task 3 (Gcode importer uses `gcode_block_ids`).
3. A bond filter missing a field its type needs (e.g. a `Disk` without `Radius`) must be reported during input validation, not as a `KeyError` / `MethodError` deep in neighbourhood search — test in Task 1 ("bond filter required fields").
4. Node sets given as an id, an id list, a range (`"1:3"`), a coordinate expression (`"x > 0.5"`), `"All"` or a file must give the same nodes as before — test in Task 3 ("read_node_sets equals get_node_sets").
5. A block ID present in the deck but missing from the mesh must still only warn (serial) and be skipped — covered by the equivalence test in Task 1 (`block_names_and_ids` vs `get_block_names_and_ids` with the deck's IDs) plus a unit test ("block names and ids").

Known pre-existing bug, preserved: in `load_and_evaluate_mesh`, `params["Discretization"]["Type"] == "Text File"` (a comparison, intended as an assignment) means an existing `.txt` mesh is never reused when `Gcode.Overwrite Mesh` is false. Task 3 drops the no-op line without changing behaviour; fixing the intent is a separate decision for the maintainers.

---

## File Structure

| File | Change |
|---|---|
| `src/Support/Parameters/Input/blocks.jl` | `block_by_id`, `block_angles`, `block_names_and_ids` |
| `src/Support/Parameters/Input/discretization.jl` | `mesh_scaling`, `gcode_block_ids`; `check!` for `GcodeParams`, `BondFilterParams` |
| `src/Support/Parameters/Input/InputDeck.jl` | exports |
| `test/helper.jl` | `typed_section`, `typed_input` test helpers |
| `src/Core/Solver/Solver_manager.jl` | block setup reads `BlockParams`; influence function gets the expression |
| `src/Core/Influence_function.jl` | `init_influence_function(nodes, expr_str)` |
| `src/IO/IO.jl` | `init_orientations(blocks)`, `show_block_summary(solver_options, blocks, …)`, `initialize_data` passes `input` to `init_data` |
| `src/IO/mesh_data.jl` | `set_angles(blocks)`, `init_data(params, input, …)`, `load_and_evaluate_mesh(input, …)`, neighbours, extrusion, topology |
| `src/Support/Parameters/parameter_handling_mesh.jl`, `parameter_handling.jl` | `read_node_sets`, `external_topology_file` |
| `src/IO/Mesh_Import/*.jl` | `read_mesh(input::PeriLabInput, …)` |
| `src/IO/bond_filter.jl`, `src/IO/Bond_Filter/*.jl` | filters take `BondFilterParams` |
| `src/PeriLab.jl` | `run` passes `input.sections.blocks` |
| tests: `Input/ut_blocks_mesh_typed.jl` (new), `ut_IO.jl`, `ut_MPI.jl`, `ut_Influence_function.jl`, `ut_mesh_data.jl`, `ut_bond_filter.jl`, `Input/input_tests.jl` | |

---

### Task 1: Block and discretization helpers next to the types, proven equivalent on every shipped deck

**Files:**
- Modify: `src/Support/Parameters/Input/blocks.jl`, `discretization.jl`, `InputDeck.jl`
- Modify: `test/helper.jl`
- Create: `test/unit_tests/Support/Parameters/Input/ut_blocks_mesh_typed.jl`
- Modify: `test/unit_tests/Support/Parameters/Input/input_tests.jl`

**Interfaces:**
- Consumes: phase 2a `BlockParams`, `DiscretizationParams`, `GcodeParams`, `BondFilterParams`, `check!`, `read_input`; `Dict` getters (reference only).
- Produces (exported from `InputDeck`):
  - `block_by_id(blocks::Dict{String,BlockParams}, block_id::Integer)::Tuple{String,BlockParams}` — aborts `"Block with ID $block_id is not defined"`.
  - `block_angles(name::AbstractString, block::BlockParams, dof::Int64)` — `nothing` if `Angle X` not given; `Float64` (`Angle X`) for `dof == 2`; `Vector{Float64}` `[X, Y, Z]` for `dof == 3` (aborts `"Angle Y of $name is not defined"` / Z); `nothing` for other `dof`.
  - `block_names_and_ids(blocks::Dict{String,BlockParams}, mesh_block_ids::AbstractVector{<:Integer}, mpi::Bool)::Tuple{Vector{String},Vector{Int64}}`.
  - `mesh_scaling(d::DiscretizationParams)::Vector{Float64}` (`[x, y, z]`, 1.0 where not given).
  - `gcode_block_ids(g::GcodeParams)::Union{Nothing,Dict{Int64,String}}`.
  - `check!(::GcodeParams)`: every `Blocks` key is an integer. `check!(::BondFilterParams)`: fields required by type `"Disk"` (Center X/Y/Z, Normal Z, Radius) and `"Rectangular_Plane"` (Lower Left Corner X/Y, Bottom Unit Vector X/Y, Bottom Length, Side Length); other types unchecked (filter modules are an extension point).
  - Exports also `BlockParams`, `DiscretizationParams`, `GcodeParams`, `BondFilterParams`, `SurfaceExtrusionParams`, `ExternalTopologyParams`.
- `test/helper.jl`: `typed_section(T, raw)` (converts and validates a section, aborting on errors) and `typed_input(sections::AbstractDict)::PeriLabInput` (a minimal valid deck with `sections` merged in).

- [ ] **Step 1: Write the failing tests**

`test/unit_tests/Support/Parameters/Input/ut_blocks_mesh_typed.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const PS = PeriLab.ParameterSpec
const ID = PeriLab.InputDeck
const PH = PeriLab.Parameter_Handling

function ut_bm(T, raw)
    ctx = PS.ParseContext()
    value = PS.convert_value(T, raw, "X", ctx)
    return value, ctx
end

ut_bm_blocks(raw) = first(ut_bm(Dict{String,ID.BlockParams}, raw))

const UT_BM_BLOCKS = ut_bm_blocks(Dict{String,Any}("left" => Dict{String,Any}("Block ID" => 1,
                                                                              "Density" => 2.0,
                                                                              "Horizon" => 0.5,
                                                                              "Angle X" => 30.0),
                                                   "right" => Dict{String,Any}("Block ID" => 3,
                                                                               "Density" => 4.0,
                                                                               "Horizon" => 0.7,
                                                                               "Angle X" => 10.0,
                                                                               "Angle Y" => 20.0,
                                                                               "Angle Z" => 0.0)))

@testset "block by id" begin
    name, block = ID.block_by_id(UT_BM_BLOCKS, 3)
    @test name == "right" && block.density === 4.0
    @test_throws PeriLab.PeriLabError ID.block_by_id(UT_BM_BLOCKS, 2)
end

@testset "block angles" begin
    @test ID.block_angles(ID.block_by_id(UT_BM_BLOCKS, 1)..., 2) === 30.0
    @test_throws PeriLab.PeriLabError ID.block_angles(ID.block_by_id(UT_BM_BLOCKS, 1)..., 3)
    @test ID.block_angles(ID.block_by_id(UT_BM_BLOCKS, 3)..., 3) == [10.0, 20.0, 0.0]
    no_angles = ut_bm_blocks(Dict{String,Any}("b" => Dict{String,Any}("Block ID" => 1,
                                                                      "Density" => 1.0,
                                                                      "Horizon" => 1.0)))
    @test ID.block_angles(ID.block_by_id(no_angles, 1)..., 3) === nothing
end

@testset "block names and ids" begin
    names, ids = ID.block_names_and_ids(UT_BM_BLOCKS, [1, 3], true)
    @test names == ["left", "right"] && ids == [1, 3]
    names, ids = ID.block_names_and_ids(UT_BM_BLOCKS, [3], true)   # block 1 not in mesh
    @test names == ["right"] && ids == [3]
end

@testset "mesh scaling" begin
    d, ctx = ut_bm(ID.DiscretizationParams,
                   Dict{String,Any}("Type" => "Text File", "Input Mesh File" => "m.txt",
                                    "Horizon Mesh Scaling Y" => 2))
    @test isempty(ctx.errors)
    @test ID.mesh_scaling(d) == [1.0, 2.0, 1.0]
end

@testset "gcode block ids" begin
    gcode(blocks) = Dict{String,Any}("Overwrite Mesh" => true, "Sampling" => 1.0,
                                     "Width" => 1.0, "Height" => 1.0, "Blocks" => blocks)
    g, ctx = ut_bm(ID.GcodeParams, gcode(Dict{Any,Any}(2 => "z > 1", 1 => "z <= 1")))
    @test isempty(ctx.errors)
    @test ID.gcode_block_ids(g) == Dict(2 => "z > 1", 1 => "z <= 1")
    g, ctx = ut_bm(ID.GcodeParams, gcode(Dict{Any,Any}("top" => "z > 1")))
    @test ctx.errors[1].path == "X.Blocks.top"
    @test ctx.errors[1].message == "expected a block id (integer) as key"
    g, _ = ut_bm(ID.GcodeParams, Dict{String,Any}("Overwrite Mesh" => true, "Sampling" => 1.0,
                                                  "Width" => 1.0, "Height" => 1.0))
    @test ID.gcode_block_ids(g) === nothing
end

@testset "bond filter required fields" begin
    _, ctx = ut_bm(ID.BondFilterParams,
                   Dict{String,Any}("Type" => "Disk", "Normal X" => 0.0, "Normal Y" => 0.0,
                                    "Normal Z" => 1.0, "Center X" => 0.0, "Center Y" => 0.0,
                                    "Center Z" => 0.0))
    @test ctx.errors[1].path == "X"
    @test ctx.errors[1].message == "\"Disk\" bond filter requires: Radius"
    _, ctx = ut_bm(ID.BondFilterParams,
                   Dict{String,Any}("Type" => "Rectangular_Plane", "Normal X" => 0.0,
                                    "Normal Y" => 1.0))
    @test ctx.errors[1].message ==
          "\"Rectangular_Plane\" bond filter requires: Lower Left Corner X, Lower Left Corner Y, Bottom Unit Vector X, Bottom Unit Vector Y, Bottom Length, Side Length"
    _, ctx = ut_bm(ID.BondFilterParams,
                   Dict{String,Any}("Type" => "My_Filter", "Normal X" => 0.0, "Normal Y" => 1.0))
    @test isempty(ctx.errors)
end

const UT_BM_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", "..", ".."))
const UT_BM_SKIP = ("test/fullscale_tests/test_PD_Solid_Elastic/strain_xx_external.yaml",)

function ut_bm_decks()
    decks = String[]
    for root in ("examples", "test"), (dir, _, files) in walkdir(joinpath(UT_BM_ROOT, root))
        append!(decks, [joinpath(dir, f) for f in files if endswith(f, ".yaml")])
    end
    return sort!(decks)
end

function ut_bm_compare_all_decks()
    for file in ut_bm_decks()
        relpath(file, UT_BM_ROOT) in UT_BM_SKIP && continue
        raw = PeriLab.IO.read_input(file)
        (raw isa AbstractDict && haskey(raw, "PeriLab")) || continue
        deck = raw["PeriLab"]
        input, _ = ID.read_input(deck, dirname(file))
        blocks = input.sections.blocks
        ids = Int64[b["Block ID"] for b in values(deck["Blocks"])]
        @test ID.block_names_and_ids(blocks, ids, true) ==
              PH.get_block_names_and_ids(deck, ids, true)
        @test ID.mesh_scaling(input.sections.discretization) == PH.get_mesh_scaling(deck)
        for id in unique(ids)
            name, block = ID.block_by_id(blocks, id)
            @test block.density == PH.get_density(deck, id)
            @test block.horizon == PH.get_horizon(deck, id)
            @test something(block.fem, false) == PH.get_fem_block(deck, id)
            if block.specific_heat_capacity !== nothing
                @test block.specific_heat_capacity == PH.get_heat_capacity(deck, id)
            end
            for dof in (2, 3)
                if block.angle_x === nothing
                    @test PH.get_angles(deck, id, dof) === nothing
                elseif dof == 2 || (block.angle_y !== nothing && block.angle_z !== nothing)
                    @test ID.block_angles(name, block, dof) == PH.get_angles(deck, id, dof)
                end
            end
        end
    end
end

@testset "block and mesh helpers equal the Dict getters on every shipped deck" begin
    ut_bm_compare_all_decks()
end
```

Update `input_tests.jl`:

```julia
for file in ["ut_mesh_blocks.jl", "ut_solver.jl", "ut_outputs_conditions_contact.jl",
             "ut_read_input.jl", "ut_golden_decks.jl", "ut_validate_yaml.jl",
             "ut_solver_typed.jl", "ut_model_reduction_params.jl", "ut_blocks_mesh_typed.jl"]
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl`
Expected: FAIL in `ut_blocks_mesh_typed.jl` with `UndefVarError: block_by_id not defined in PeriLab.InputDeck`.

- [ ] **Step 3: Write the implementation**

Append to `src/Support/Parameters/Input/blocks.jl`:

```julia
"""
    block_by_id(blocks, block_id) -> (name, block)

The block with `Block ID` `block_id`.
"""
function block_by_id(blocks::Dict{String,BlockParams}, block_id::Integer)
    for (name, block) in blocks
        block.block_id == block_id && return name, block
    end
    @abort "Block with ID $block_id is not defined"
end

"""
    block_angles(name, block, dof)

Block rotation: `nothing` if `Angle X` is not given, `Angle X` in 2D, and
`[Angle X, Angle Y, Angle Z]` in 3D (all three are required then).
"""
function block_angles(name::AbstractString, block::BlockParams, dof::Int64)
    block.angle_x === nothing && return nothing
    dof == 2 && return block.angle_x
    dof == 3 || return nothing
    for (key, value) in (("Angle Y", block.angle_y), ("Angle Z", block.angle_z))
        value === nothing && @abort "$key of $name is not defined"
    end
    return [block.angle_x, block.angle_y, block.angle_z]
end

"""
    block_names_and_ids(blocks, mesh_block_ids, mpi) -> (names, ids)

Names and IDs of the blocks present in the mesh, ordered by ID. Blocks of the
input deck that are missing in the mesh are skipped (with a warning unless
running under MPI).
"""
function block_names_and_ids(blocks::Dict{String,BlockParams},
                             mesh_block_ids::AbstractVector{<:Integer}, mpi::Bool)
    names = String[]
    ids = Int64[]
    for id in 1:maximum(block.block_id for block in values(blocks))
        if !(id in mesh_block_ids)
            mpi || @warn "Block with ID $id is not defined in the provided mesh"
            continue
        end
        for (name, block) in blocks
            if block.block_id == id
                push!(names, name)
                push!(ids, id)
            end
        end
    end
    return names, ids
end
```

Append to `src/Support/Parameters/Input/discretization.jl`:

```julia
"Horizon scaling per direction, 1.0 where not given."
function mesh_scaling(d::DiscretizationParams)
    return [something(d.horizon_mesh_scaling_x, 1.0), something(d.horizon_mesh_scaling_y, 1.0),
            something(d.horizon_mesh_scaling_z, 1.0)]
end

function check!(g::GcodeParams, path::String, ctx::ParseContext)
    g.blocks === nothing && return nothing
    for key in keys(g.blocks)
        tryparse(Int64, key) === nothing &&
            add_error!(ctx, join_path(join_path(path, "Blocks"), key),
                       "expected a block id (integer) as key")
    end
    return nothing
end

"Gcode block assignment: block id => condition, or `nothing`."
function gcode_block_ids(g::GcodeParams)
    g.blocks === nothing && return nothing
    return Dict{Int64,String}(parse(Int64, key) => condition for (key, condition) in g.blocks)
end

# Fields each built-in bond filter needs (Z components only in 3D, not checked here).
const _BOND_FILTER_REQUIRED = Dict("Disk" => (("Center X", :center_x), ("Center Y", :center_y),
                                              ("Center Z", :center_z), ("Normal Z", :normal_z),
                                              ("Radius", :radius)),
                                   "Rectangular_Plane" => (("Lower Left Corner X",
                                                            :lower_left_corner_x),
                                                           ("Lower Left Corner Y",
                                                            :lower_left_corner_y),
                                                           ("Bottom Unit Vector X",
                                                            :bottom_unit_vector_x),
                                                           ("Bottom Unit Vector Y",
                                                            :bottom_unit_vector_y),
                                                           ("Bottom Length", :bottom_length),
                                                           ("Side Length", :side_length)))

function check!(f::BondFilterParams, path::String, ctx::ParseContext)
    required = get(_BOND_FILTER_REQUIRED, f.type, ())
    missing = [key for (key, field) in required if getfield(f, field) === nothing]
    isempty(missing) ||
        add_error!(ctx, path, "\"$(f.type)\" bond filter requires: $(join(missing, ", "))")
    return nothing
end
```

In `src/Support/Parameters/Input/InputDeck.jl`, extend the `export` statement with:

```julia
       BlockParams, DiscretizationParams, GcodeParams, BondFilterParams,
       SurfaceExtrusionParams, ExternalTopologyParams, block_by_id, block_angles,
       block_names_and_ids, mesh_scaling, gcode_block_ids
```

(append these names after `solver_step,` — keep one `export` statement).

Append to `test/helper.jl`:

```julia
"""
    typed_section(T, raw)

Converts and validates `raw` as `T` (an `@params` struct or `Dict{String,<@params>}`)
the way the input reader does; aborts on input errors.
"""
function typed_section(T, raw)
    ctx = PeriLab.ParameterSpec.ParseContext()
    value = PeriLab.ParameterSpec.convert_value(T, raw, "test", ctx)
    PeriLab.ParameterSpec.report!(ctx)
    return value
end

"""
    typed_input(sections)

A typed `PeriLabInput` from a minimal valid deck with `sections` merged in.
"""
function typed_input(sections::AbstractDict)
    deck = Dict{String,Any}("Discretization" => Dict{String,Any}("Type" => "Text File",
                                                                 "Input Mesh File" => "mesh.txt"),
                            "Blocks" => Dict{String,Any}("block_1" => Dict{String,Any}("Block ID" => 1,
                                                                                       "Density" => 1.0,
                                                                                       "Horizon" => 1.0)),
                            "Models" => Dict{String,Any}(),
                            "Solver" => Dict{String,Any}("Initial Time" => 0.0,
                                                         "Final Time" => 1.0,
                                                         "Verlet" => Dict{String,Any}()))
    merge!(deck, Dict{String,Any}(string(k) => v for (k, v) in sections))
    input, ctx = PeriLab.InputDeck.read_input(deck)
    PeriLab.ParameterSpec.report!(ctx)
    return input
end
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl`
Expected: PASS — including the golden deck test, which now also applies the Gcode and bond filter checks to every shipped deck. If a shipped deck fails a new check, the check is wrong unless the deck could not run (record a ruling either way).

- [ ] **Step 5: Commit**

```bash
git add src/Support/Parameters/Input test/helper.jl test/unit_tests/Support/Parameters/Input
git -c user.name="Jan-Timo Hesse" -c user.email="jthesse@gmx.de" commit -m "Typed block and discretization helpers; Gcode and bond filter checks"
```

---

### Task 2: Block consumers and the influence function read typed input

**Files:**
- Modify: `src/Core/Solver/Solver_manager.jl`
- Modify: `src/Core/Influence_function.jl`
- Modify: `src/IO/IO.jl` (`init_orientations`, `show_block_summary`)
- Modify: `src/IO/mesh_data.jl` (`set_angles(params::Dict)` only)
- Modify: `src/PeriLab.jl` (`run`)
- Modify: `test/unit_tests/IO/ut_IO.jl`, `test/unit_tests/MPI_communication/ut_MPI.jl`, `test/unit_tests/Core/ut_Influence_function.jl`

**Interfaces:**
- Consumes: Task 1 `block_by_id`, `block_angles`, `block_names_and_ids`, `BlockParams`; `typed_section` (tests).
- Produces:
  - `Solver_Manager.set_density(blocks::Dict{String,BlockParams}, block_nodes, density)`, `set_horizon(blocks, …)`, `set_fem_block(blocks, …)`, `set_angles(blocks::Dict{String,BlockParams}, block_nodes::Dict)`.
  - `Influence_Function.init_influence_function(nodes::AbstractVector{Int64}, expr_str::Union{Nothing,String})`.
  - `IO.init_orientations(blocks::Dict{String,BlockParams})`; `IO.set_angles(blocks::Dict{String,BlockParams})` (mesh_data).
  - `IO.show_block_summary(solver_options::Dict, blocks::Dict{String,BlockParams}, log_file::String, silent::Bool, comm::MPI.Comm)`.

- [ ] **Step 1: Update the existing tests to the new signatures (they become the failing tests)**

`test/unit_tests/Core/ut_Influence_function.jl` — every call passes the expression instead of a dict:

```bash
perl -0pi -e 's/(init_influence_function\([^,]+,\s*)params\)/$1params["Influence Function"])/g; s/(init_influence_function\([^,]+,\s*)Dict\(\)\)/$1nothing)/g' test/unit_tests/Core/ut_Influence_function.jl
grep -n 'init_influence_function(' -A1 test/unit_tests/Core/ut_Influence_function.jl | grep -c 'params\["Influence Function"\])\|nothing)'
```

Expected count: 9.

`test/unit_tests/IO/ut_IO.jl` — replace

```julia
    PeriLab.IO.init_orientations(Dict("Blocks" => Dict("block_1" => Dict("Block ID" => 1),
                                                       "block_2" => Dict("Block ID" => 2))))
```

with

```julia
    PeriLab.IO.init_orientations(typed_section(Dict{String,PeriLab.InputDeck.BlockParams},
                                               Dict("block_1" => Dict("Block ID" => 1,
                                                                      "Density" => 1.0,
                                                                      "Horizon" => 1.0),
                                                    "block_2" => Dict("Block ID" => 2,
                                                                      "Density" => 1.0,
                                                                      "Horizon" => 1.0))))
```

replace

```julia
    PeriLab.IO.init_orientations(Dict("Blocks" => Dict("block_1" => Dict("Block ID" => 1))))
```

with

```julia
    PeriLab.IO.init_orientations(typed_section(Dict{String,PeriLab.InputDeck.BlockParams},
                                               Dict("block_1" => Dict("Block ID" => 1,
                                                                      "Density" => 1.0,
                                                                      "Horizon" => 1.0))))
```

and replace the `params = Dict("Blocks" => Dict("block_1" => Dict("Material Models" => true, …)))` assignment before `PeriLab.IO.show_block_summary(solver_options,` (the whole statement through its closing `)))`) with

```julia
    params = typed_section(Dict{String,PeriLab.InputDeck.BlockParams},
                           Dict("block_1" => Dict("Block ID" => 1, "Density" => 1.0,
                                                  "Horizon" => 1.0, "Material Model" => "m",
                                                  "Damage Model" => "d",
                                                  "Additive Model" => "a",
                                                  "Thermal Model" => "t",
                                                  "Degradation Model" => "g"),
                                "block_2" => Dict("Block ID" => 2, "Density" => 1.0,
                                                  "Horizon" => 1.0, "Material Model" => "m")))
```

`test/unit_tests/MPI_communication/ut_MPI.jl` — replace

```julia
    params = Dict("Blocks" => Dict("block_1" => Dict("Block ID" => 1,
                                                     "Material Model" => "Test 1"),
                                   "block_2" => Dict("Block ID" => 2,
                                                     "Material Model" => "Test 2")))
```

with

```julia
    params = typed_section(Dict{String,PeriLab.InputDeck.BlockParams},
                           Dict("block_1" => Dict("Block ID" => 1, "Density" => 1.0,
                                                  "Horizon" => 1.0,
                                                  "Material Model" => "Test 1"),
                                "block_2" => Dict("Block ID" => 2, "Density" => 1.0,
                                                  "Horizon" => 1.0,
                                                  "Material Model" => "Test 2")))
```

- [ ] **Step 2: Run the updated tests to verify they fail**

Run each with the single-file command from Global Constraints for `unit_tests/Core/ut_Influence_function.jl` and `unit_tests/IO/ut_IO.jl`.
Expected: FAIL with `MethodError: no method matching init_influence_function(::Vector{Int64}, ::String)` and `MethodError: no method matching init_orientations(::Dict{String, PeriLab.InputDeck.BlockParams})`.

- [ ] **Step 3: Write the implementation**

`src/Core/Influence_function.jl` — replace

```julia
function init_influence_function(nodes::AbstractVector{Int64},
                                 params::Dict)
    if !haskey(params, "Influence Function")
        return
    end
```

with

```julia
function init_influence_function(nodes::AbstractVector{Int64},
                                 expr_str::Union{Nothing,String})
    expr_str === nothing && return
```

delete the line `    expr_str = params["Influence Function"]`, and in the docstring replace `init_influence_function(nodes::AbstractVector{Int64}, params::Dict)` with `init_influence_function(nodes::AbstractVector{Int64}, expr_str::Union{Nothing,String})` and the `params` argument description with ``- `expr_str`: The `Discretization` `Influence Function` expression, or `nothing` ``.

`src/Core/Solver/Solver_manager.jl`:
- delete the whole `using ..Parameter_Handling:` statement (its remaining names — `get_density`, `get_horizon`, `get_fem_block`, `get_angles`, `get_block_names_and_ids` — are all replaced);
- extend `using ..InputDeck: PeriLabInput, SolverParams, solver_name, model_options, solver_step` with `, BlockParams, block_by_id, block_angles, block_names_and_ids`;
- in `init`, replace

```julia
    block_name_list,
    block_id_list = get_block_names_and_ids(params, block_ids,
                                            Data_Manager.get_mpi_active())
```

with

```julia
    blocks = input.sections.blocks
    block_name_list,
    block_id_list = block_names_and_ids(blocks, block_ids, Data_Manager.get_mpi_active())
```

and replace `set_fem_block(params,`, `set_density(params,`, `set_horizon(params,`, `set_angles(params,` with `set_fem_block(blocks,`, `set_density(blocks,`, `set_horizon(blocks,`, `set_angles(blocks,`; replace

```julia
        Influence_Function.init_influence_function(block_nodes[iblock],
                                                   params["Discretization"])
```

with

```julia
        Influence_Function.init_influence_function(block_nodes[iblock],
                                                   input.sections.discretization.influence_function)
```

- replace the four setter definitions:

```julia
function set_density(blocks::Dict{String,BlockParams}, block_nodes::Dict,
                     density::NodeScalarField{Float64})
    for block in eachindex(block_nodes)
        density[block_nodes[block]] .= block_by_id(blocks, block)[2].density
    end
    return density
end
```

```julia
function set_fem_block(blocks::Dict{String,BlockParams}, block_nodes::Dict,
                       fem_block::Vector{Bool})
    for block in eachindex(block_nodes)
        fem_block[block_nodes[block]] .= something(block_by_id(blocks, block)[2].fem, false)
    end
    return fem_block
end
```

```julia
function set_horizon(blocks::Dict{String,BlockParams}, block_nodes::Dict,
                     horizon::NodeScalarField{Float64})
    for block in eachindex(block_nodes)
        horizon[block_nodes[block]] .= block_by_id(blocks, block)[2].horizon
    end
    return horizon
end
```

and in `function set_angles(params::Dict, block_nodes::Dict)` change the signature to `function set_angles(blocks::Dict{String,BlockParams}, block_nodes::Dict)` and replace both `get_angles(params, block, dof)` with `block_angles(block_by_id(blocks, block)..., dof)`. Update the four docstrings: first argument `blocks::Dict{String,BlockParams}`: "The blocks of the input deck".

`src/IO/mesh_data.jl` — in `function set_angles(params::Dict)` (the one without `block_nodes`) change the signature to `function set_angles(blocks::Dict{String,BlockParams})`, replace `get_angles(params, block, dof)` with `block_angles(block_by_id(blocks, block)..., dof)` and `get_angles(params, block_ids[iID], dof)` with `block_angles(block_by_id(blocks, block_ids[iID])..., dof)`; remove `get_angles` from the file's `using ..Parameter_Handling:` list and add `using ..InputDeck: BlockParams, block_by_id, block_angles` (mesh_data is included into `IO`, which uses `..` for PeriLab-level modules).

`src/IO/IO.jl`:
- remove `get_fem_block` from the `using ..Parameter_Handling:` list and extend `using ..InputDeck: solver_steps` with `, BlockParams`;
- `function init_orientations(params::Dict)` → `function init_orientations(blocks::Dict{String,BlockParams})`, and its first line `set_angles(params)` → `set_angles(blocks)`;
- in `show_block_summary`, change the second argument `params::Dict,` to `blocks::Dict{String,BlockParams},`, add before the function:

```julia
# Block summary columns that name a model of the block
const _BLOCK_MODEL_COLUMNS = Dict("Material" => :material_model, "Damage" => :damage_model,
                                  "Thermal" => :thermal_model, "Additive" => :additive_model,
                                  "Degradation" => :degradation_model)

function _block_summary_cell(block::BlockParams, column::String)
    column == "Density" && return @sprintf("%.3e", block.density)
    column == "Horizon" && return @sprintf("%.3e", block.horizon)
    field = get(_BLOCK_MODEL_COLUMNS, column, nothing)
    field === nothing && return ""
    return something(getfield(block, field), "")
end
```

and replace, inside the row loop,

```julia
            elseif name == "PD/FEM"
                fem_block = get_fem_block(params, id)
                if fem_block
                    push!(row, "FEM")
                else
                    push!(row, "PD")
                end
                # elseif !(name in solver_options["Models"])
                #     push!(row, "")
            elseif haskey(params["Blocks"][block_name_list[id]], name * " Model")
                push!(row, params["Blocks"][block_name_list[id]][name*" Model"])
            elseif haskey(params["Blocks"][block_name_list[id]], name)
                push!(row, @sprintf("%.3e", (params["Blocks"][block_name_list[id]][name])))
            else
                push!(row, "")
            end
```

with

```julia
            elseif name == "PD/FEM"
                fem = something(blocks[block_name_list[id]].fem, false)
                push!(row, fem ? "FEM" : "PD")
            else
                push!(row, _block_summary_cell(blocks[block_name_list[id]], name))
            end
```

(The old code looked the FEM flag up by list index `id` as if it were a block ID; reading it from the block itself is what the row shows. Update the docstring argument to `blocks`.)

`src/PeriLab.jl`, in `run`: `IO.init_orientations(params)` → `IO.init_orientations(input.sections.blocks)`, and in `IO.show_block_summary(solver_options,` the next argument `params,` → `input.sections.blocks,`.

- [ ] **Step 4: Run tests to verify they pass**

Run the single-file command for `unit_tests/Core/ut_Influence_function.jl`, `unit_tests/IO/ut_IO.jl` and (with 2 ranks not required — it runs serially in `runtests.jl`) `unit_tests/MPI_communication/ut_MPI.jl`.
Expected: PASS.
Run: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl`
Expected: PASS.
Run (background, ≈35 min): `julia --project=. -e 'using Pkg; Pkg.test()' > /tmp/full_2b2_t2.log 2>&1; tail -5 /tmp/full_2b2_t2.log`
Expected: `Testing PeriLab tests passed`.

- [ ] **Step 5: Commit**

```bash
git add src test
git -c user.name="Jan-Timo Hesse" -c user.email="jthesse@gmx.de" commit -m "Block setup, orientations, block summary and influence function read typed input"
```

---

### Task 3: Mesh pipeline reads typed input (importers, node sets, neighbours, extrusion, topology)

**Files:**
- Modify: `src/Support/Parameters/parameter_handling_mesh.jl`, `src/Support/Parameters/parameter_handling.jl`
- Modify: `src/IO/mesh_data.jl`, `src/IO/IO.jl` (`initialize_data`)
- Modify: `src/IO/Mesh_Import/Mesh_Import.jl`, `Text_Mesh.jl`, `Exodus_Mesh.jl`, `Gmsh_Mesh.jl`, `Abaqus_Mesh.jl`, `Gcode_Mesh.jl`
- Modify: `test/unit_tests/IO/ut_mesh_data.jl`
- Create: `test/unit_tests/Support/Parameters/Input/ut_node_sets_typed.jl` (+ `input_tests.jl`)

**Interfaces:**
- Consumes: Task 1 helpers; Task 2; 2b-1 `IO.initialize_data` has `input`.
- Produces:
  - `Parameter_Handling.read_node_sets(d::DiscretizationParams, path::String, mesh_df::DataFrame)::Dict{String,Vector{Int64}}` (same results as `get_node_sets`).
  - `Parameter_Handling.external_topology_file(d::DiscretizationParams, path::String)::Union{Nothing,String}`.
  - `IO.init_data(params::Dict, input::PeriLabInput, path::String, comm::MPI.Comm)` (`params` still needed by `contact_basis`, slice 2b-5).
  - `IO.load_and_evaluate_mesh(params::Dict, input::PeriLabInput, path::String, ranksize::Int64)` (`params` still needed by `apply_bond_filters` until Task 4).
  - `IO.read_mesh(input::PeriLabInput, path::String)`; every importer: `read_mesh(input::PeriLabInput, filename::String)`.
  - `IO.create_neighborhoodlist(mesh::DataFrame, blocks::Dict{String,BlockParams}, scaling::Vector{Float64}, dof::Int64)`, `IO.neighbors(mesh::DataFrame, blocks::Dict{String,BlockParams}, scaling::Vector{Float64}, coor)`.
  - `IO.extrude_surface_mesh(mesh::DataFrame, extrusion::Union{Nothing,SurfaceExtrusionParams})`.
  - `IO.create_consistent_neighborhoodlist(external_topology::DataFrame, add_neighbor_search::Bool, nlist, dof)`.

- [ ] **Step 1: Write the failing tests**

`test/unit_tests/Support/Parameters/Input/ut_node_sets_typed.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
using DataFrames
const PS = PeriLab.ParameterSpec
const ID = PeriLab.InputDeck
const PH = PeriLab.Parameter_Handling

@testset "read_node_sets equals get_node_sets" begin
    dir = mktempdir()
    write(joinpath(dir, "ns.txt"), "header: global_id\n2\n3\n")
    sets = Dict{String,Any}("id" => 2, "list" => "1 3", "range" => "1:3",
                            "expr" => "x > 0.5", "all" => "All", "file" => "ns.txt")
    raw = Dict{String,Any}("Type" => "Text File", "Input Mesh File" => "m.txt",
                           "Node Sets" => sets)
    ctx = PS.ParseContext()
    d = PS.convert_value(ID.DiscretizationParams, raw, "D", ctx)
    @test isempty(ctx.errors)
    mesh = DataFrame(x = [0.0, 1.0, 2.0], y = [0.0, 0.0, 0.0])
    typed = PH.read_node_sets(d, dir, mesh)
    reference = PH.get_node_sets(Dict("Discretization" => raw), dir, mesh)
    @test keys(typed) == keys(reference)
    for key in keys(reference)
        @test typed[key] == reference[key]
    end
end

@testset "external topology file" begin
    dir = mktempdir()
    write(joinpath(dir, "topo.txt"), "1 2 3\n")
    ctx = PS.ParseContext()
    d = PS.convert_value(ID.DiscretizationParams,
                         Dict{String,Any}("Type" => "Text File", "Input Mesh File" => "m.txt",
                                          "Input External Topology" => Dict{String,Any}("File" => "topo.txt")),
                         "D", ctx)
    @test PH.external_topology_file(d, dir) == "topo.txt"
    @test_throws PeriLab.PeriLabError PH.external_topology_file(d, mktempdir())
    d2 = PS.convert_value(ID.DiscretizationParams,
                          Dict{String,Any}("Type" => "Text File", "Input Mesh File" => "m.txt"),
                          "D", ctx)
    @test PH.external_topology_file(d2, dir) === nothing
end
```

Add `"ut_node_sets_typed.jl"` to the list in `input_tests.jl` (after `"ut_blocks_mesh_typed.jl"`).

In `test/unit_tests/IO/ut_mesh_data.jl`, convert the mesh-import, extrusion and topology tests:
- every `params = Dict("Discretization" => Dict("Type" => <T>, "Input Mesh File" => <F>))` used with `PeriLab.IO.read_mesh` becomes `params = typed_input(Dict("Discretization" => Dict("Type" => <T>, "Input Mesh File" => <F>)))` (same `<T>`, `<F>`; the calls `PeriLab.IO.read_mesh(params, …)` stay);
- `params = Dict()` before the first `create_consistent_neighborhoodlist` call becomes `params = false`; `params = Dict("Add Neighbor Search" => false)` becomes `params = false`; `params = Dict("Add Neighbor Search" => true)` becomes `params = true`;
- `params = Dict("Discretization" => Dict("Bla" => "Bla"))` becomes `params = nothing`;
- every `params = Dict("Discretization" => Dict("Surface Extrusion" => Dict(<entries>)))` becomes `params = typed_section(PeriLab.InputDeck.SurfaceExtrusionParams, Dict(<entries>))`.

After editing, `grep -n 'params = Dict' test/unit_tests/IO/ut_mesh_data.jl` lists only dicts not passed to the functions above (check each remaining hit by its next call).

- [ ] **Step 2: Run tests to verify they fail**

Run: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl`
Expected: FAIL in `ut_node_sets_typed.jl` with `UndefVarError: read_node_sets not defined in PeriLab.Parameter_Handling`.
Run the single-file command for `unit_tests/IO/ut_mesh_data.jl`.
Expected: FAIL with `MethodError`s for `read_mesh(::PeriLab.InputDeck.PeriLabInput, ::String)`, `extrude_surface_mesh(::DataFrame, ::Nothing)` and `create_consistent_neighborhoodlist(::DataFrame, ::Bool, …)`.

- [ ] **Step 3: Node sets and external topology (Parameter_Handling)**

In `src/Support/Parameters/parameter_handling.jl`, extend `using ..InputDeck: read_input` with `, DiscretizationParams`. In `parameter_handling_mesh.jl`, add at the top (after its `using` lines) `export read_node_sets, external_topology_file`, and append:

```julia
"""
    external_topology_file(d, path)

File name of the external element topology, or `nothing`. Aborts if the file
does not exist in `path`.
"""
function external_topology_file(d::DiscretizationParams, path::String)
    topology = d.input_external_topology
    topology === nothing && return nothing
    filename = joinpath(path, topology.file)
    if !isfile(filename)
        @abort "External topology file: ''$filename'' does not exist"
        return
    end
    return topology.file
end

"""
    read_node_sets(d, path, mesh_df)

Node sets of the discretization: from an Exodus mesh, or from `Node Sets`
entries (node id, id list, range `a:b`, coordinate expression in `x`/`y`/`z`,
`All`, or a node set file).
"""
function read_node_sets(d::DiscretizationParams, path::String, mesh_df::DataFrame)
    nsets = Dict{String,Vector{Int64}}()
    if d.type == "Exodus"
        exo = ExodusDatabase(joinpath(path, d.input_mesh_file), "r")
        nset_names = read_names(exo, NodeSet)
        conn = collect_element_connectivities(exo)
        for (id, entry) in enumerate(nset_names)
            nset_nodes = Vector{Int64}(read_set(exo, NodeSet, id).nodes)
            name = length(entry) == 0 ? "Set-" * string(id) : entry
            nsets[name] = findall(row -> all(val -> any(val .== nset_nodes), row), conn)
        end
        @info "Found $(length(nsets)) node sets"
        close(exo)
        return nsets
    end
    for (entry, value) in d.node_sets
        if value isa Int64
            nsets[entry] = [value]
        elseif occursin(".txt", value)
            if isnothing(get_header(joinpath(path, value)))
                @warn "Node set file " * value *
                      " is not correctly specified. Please check the examples. The node set is excluded."
                continue
            end
            header_line, header = get_header(joinpath(path, value))
            nodes = CSV.read(joinpath(path, value), DataFrame; delim = " ", header = false,
                             skipto = header_line + 1,)
            if size(nodes) == (0, 0)
                @abort "Node set file is empty " * value * ". The node set is excluded."
                return
            end
            nsets[entry] = nodes.Column1
        elseif occursin(":", value)
            nsets[entry] = collect(eval(Meta.parse(value)))
        elseif occursin("x", value) || occursin("y", value) || occursin("z", value)
            nodes = []
            for id in 1:size(mesh_df, 1)
                global x = mesh_df[!, "x"][id]
                global y = mesh_df[!, "y"][id]
                if occursin("z", value)
                    global z = mesh_df[!, "z"][id]
                end
                try
                    if eval(Meta.parse(value))
                        push!(nodes, id)
                    end
                catch UndefVarError
                    @abort "Failed to eval nodeset value: '$(value)', $UndefVarError"
                    return
                end
            end
            nsets[entry] = nodes
        elseif value == "All"
            nsets[entry] = collect(1:size(mesh_df, 1))
        else
            nsets[entry] = parse.(Int, split(value))
        end
    end
    return nsets
end
```

(This is `get_node_sets` with the dict lookups replaced by fields; the evaluation logic is unchanged.)

- [ ] **Step 4: Mesh import reads `PeriLabInput`**

`src/IO/Mesh_Import/Mesh_Import.jl`: replace `using ...Parameter_Handling: get_mesh_name` with `using ...InputDeck: PeriLabInput`; in `read_mesh`, change `function read_mesh(params::Dict, path::String)` to `function read_mesh(input::PeriLabInput, path::String)`, `mesh_name = get_mesh_name(params)` to `mesh_name = input.sections.discretization.input_mesh_file`, `type = params["Discretization"]["Type"]` to `type = input.sections.discretization.type`, `return importer.read_mesh(params, mesh_path)` to `return importer.read_mesh(input, mesh_path)`; update the docstring (`input::PeriLabInput`: "The typed input deck"; dispatch on `Discretization.Type`).

In each importer (`Text_Mesh.jl`, `Exodus_Mesh.jl`, `Gmsh_Mesh.jl`, `Abaqus_Mesh.jl`, `Gcode_Mesh.jl`): add `using ....InputDeck: PeriLabInput` next to `using ....PeriLabExceptions` (Gcode: also `gcode_block_ids`), change `function read_mesh(params::Dict, filename::String)` to `function read_mesh(input::PeriLabInput, filename::String)` and update its docstring argument. Then:

`Abaqus_Mesh.jl` — replace

```julia
    for boundary_condtion in keys(params["Boundary Conditions"])
        if haskey(params["Boundary Conditions"][boundary_condtion], "Node Set")
            push!(nset_names,
                  params["Boundary Conditions"][boundary_condtion]["Node Set"])
        end
    end
```

with

```julia
    for bc in values(input.sections.boundary_conditions)
        push!(nset_names, bc.node_set)
    end
```

`Gcode_Mesh.jl` — replace the block from `sampling = params["Discretization"]["Gcode"]["Sampling"]` through `commands_dict["End"] = get(params["Discretization"]["Gcode"], "End Command", nothing)` with

```julia
    gcode = input.sections.discretization.gcode
    gcode === nothing && @abort "Discretization type \"Gcode\" needs a \"Gcode\" section"
    sampling = gcode.sampling
    scale = gcode.scale
    width = gcode.width
    height = gcode.height
    blocks = gcode_block_ids(gcode)

    commands_dict = Dict{String,Any}()
    commands_dict["Start"] = gcode.start_command
    commands_dict["Stop"] = gcode.stop_command
    commands_dict["End"] = gcode.end_command
```

(`gcode_block_ids` restores integer block IDs for the `block_id = block[1]` assignment.)

After the edits `grep -n 'params' src/IO/Mesh_Import/{Text_Mesh,Exodus_Mesh,Gmsh_Mesh,Abaqus_Mesh,Gcode_Mesh}.jl` must show no `params[` reads (docstring mentions are fine only if they describe the new argument).

- [ ] **Step 5: `mesh_data.jl` pipeline**

- `using ..Parameter_Handling:` list: remove `get_mesh_name`, `get_node_sets`, `get_external_topology_name`, `get_horizon`, `get_mesh_scaling`; add `read_node_sets, external_topology_file`; extend `using ..InputDeck: BlockParams, block_by_id, block_angles` with `, PeriLabInput, SurfaceExtrusionParams, mesh_scaling`.
- `init_data`: signature `function init_data(params::Dict, input::PeriLabInput, path::String, comm::MPI.Comm)`; in it, `load_and_evaluate_mesh(params,` + `path,` + `size)` becomes `load_and_evaluate_mesh(params, input, path, size)` (keep the line layout), and `Data_Manager.set_horizon_mesh_scaling(get_mesh_scaling(params))` becomes `Data_Manager.set_horizon_mesh_scaling(mesh_scaling(input.sections.discretization))`. Update the docstring.
- Replace the head of `load_and_evaluate_mesh` from its signature through `mesh, surface_ns = extrude_surface_mesh(mesh, params)` with:

```julia
function load_and_evaluate_mesh(params::Dict,
                                input::PeriLabInput,
                                path::String,
                                ranksize::Int64)
    discretization = input.sections.discretization
    if discretization.type == "Abaqus"
        @timeit "read_mesh" mesh, nsets=read_mesh(input, path)
    else
        # Gmsh, Gcode, Text File, Exodus. (For Gcode the old code compared
        # Type == "Text File" instead of assigning it, so an existing .txt mesh
        # is never reused; that behaviour is kept.)
        @debug "Read node sets"
        @timeit "read_mesh" mesh=read_mesh(input, path)
        nsets = read_node_sets(discretization, path, mesh)
    end
    nnodes = size(mesh, 1) + 1
    mesh, surface_ns = extrude_surface_mesh(mesh, discretization.surface_extrusion)
```

- In the same function replace

```julia
    if !isnothing(get_external_topology_name(params, path))
        external_topology = read_external_topology(joinpath(path,
                                                            get_external_topology_name(params,
                                                                                       path)))
    end
```

with

```julia
    topology_file = external_topology_file(discretization, path)
    if !isnothing(topology_file)
        external_topology = read_external_topology(joinpath(path, topology_file))
    end
```

replace `create_neighborhoodlist(mesh, params, dof)` with `create_neighborhoodlist(mesh, input.sections.blocks, mesh_scaling(discretization), dof)`, replace

```julia
        topology = create_consistent_neighborhoodlist(external_topology,
                                                      params["Discretization"]["Input External Topology"],
                                                      nlist,
                                                      dof)
```

with

```julia
        topology = create_consistent_neighborhoodlist(external_topology,
                                                      something(discretization.input_external_topology.add_neighbor_search,
                                                                false),
                                                      nlist,
                                                      dof)
```

replace `if haskey(params["Discretization"], "Distribution Type")` with `if discretization.distribution_type !== nothing`, the argument `params["Discretization"]["Distribution Type"])` with `discretization.distribution_type)`, and `if haskey(params, "FEM") && !isnothing(external_topology)` with `if input.sections.fem !== nothing && !isnothing(external_topology)`. (`apply_bond_filters(nlist, mesh, params, dof)` stays until Task 4.) Update the docstring.

- `create_neighborhoodlist`:

```julia
function create_neighborhoodlist(mesh::DataFrame, blocks::Dict{String,BlockParams},
                                 scaling::Vector{Float64}, dof::Int64)
    coor = names(mesh)
    return neighbors(mesh, blocks, scaling, coor[1:dof])
end
```

- `neighbors`: signature `function neighbors(mesh::DataFrame, blocks::Dict{String,BlockParams}, scaling::Vector{Float64}, coor::Union{Vector{Int64},Vector{String}})`; replace

```julia
    radius = zeros(Float64, nnodes)
    for iID in 1:nnodes
        radius[iID] = get_horizon(params, mesh[!, "block_id"][iID])
    end
    mesh_scaling = get_mesh_scaling(params)

    return get_nearest_neighbors(1:nnodes, dof, data, data, radius, neighborList;
                                 mesh_scaling = mesh_scaling)
```

with

```julia
    horizon = Dict(id => block_by_id(blocks, id)[2].horizon for id in block_ids)
    radius = [horizon[id] for id in mesh[!, "block_id"]]

    return get_nearest_neighbors(1:nnodes, dof, data, data, radius, neighborList;
                                 mesh_scaling = scaling)
```

- `extrude_surface_mesh`: signature `function extrude_surface_mesh(mesh::DataFrame, extrusion::Union{Nothing,SurfaceExtrusionParams})`; replace

```julia
    if !("Surface Extrusion" in keys(params["Discretization"]))
        return mesh, nothing
    end
    direction = params["Discretization"]["Surface Extrusion"]["Direction"]
    step_x = params["Discretization"]["Surface Extrusion"]["Step_X"]
    step_y = params["Discretization"]["Surface Extrusion"]["Step_Y"]
    step_z = params["Discretization"]["Surface Extrusion"]["Step_Z"]
    number = params["Discretization"]["Surface Extrusion"]["Number"]
```

with

```julia
    extrusion === nothing && return mesh, nothing
    direction = extrusion.direction
    step_x = extrusion.step_x
    step_y = extrusion.step_y
    step_z = extrusion.step_z
    number = extrusion.number
```

- `create_consistent_neighborhoodlist`: change the second argument `params::Dict,` to `add_neighbor_search::Bool,` and replace

```julia
    pd_neighbors::Bool = false
    if haskey(params, "Add Neighbor Search")
        pd_neighbors = params["Add Neighbor Search"]
    end
```

with

```julia
    pd_neighbors::Bool = add_neighbor_search
```

Update the docstrings of all changed functions.

`src/IO/IO.jl`, `initialize_data`: `params=init_data(deck, filedirectory, comm)` → `params=init_data(deck, input, filedirectory, comm)`.

- [ ] **Step 6: Run tests to verify they pass**

Run: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl` — Expected: PASS.
Run the single-file command for `unit_tests/IO/ut_mesh_data.jl` — Expected: PASS.
Run (background, ≈35 min): `julia --project=. -e 'using Pkg; Pkg.test()' > /tmp/full_2b2_t3.log 2>&1; tail -5 /tmp/full_2b2_t3.log` — Expected: `Testing PeriLab tests passed`.

- [ ] **Step 7: Commit**

```bash
git add src test
git -c user.name="Jan-Timo Hesse" -c user.email="jthesse@gmx.de" commit -m "Mesh import, node sets, neighbours, extrusion and topology read typed input"
```

---

### Task 4: Bond filters take `BondFilterParams`

**Files:**
- Modify: `src/IO/bond_filter.jl`, `src/IO/Bond_Filter/disk_filter.jl`, `rectangular_filter.jl`, `filter_template.jl`
- Modify: `src/IO/mesh_data.jl` (`load_and_evaluate_mesh`, `init_data`)
- Modify: `test/unit_tests/IO/ut_bond_filter.jl`

**Interfaces:**
- Consumes: Task 1 `BondFilterParams` (+ its `check!`); Task 3 pipeline.
- Produces:
  - `IO.apply_bond_filters(nlist::BondScalarState{Int64}, mesh::DataFrame, filters::Dict{String,BondFilterParams}, dof::Int64)`.
  - `run_bond_filter(nnodes::Int64, data::Matrix{Float64}, filter::BondFilterParams, nlist::BondScalarState{Int64}, dof::Int64)` in every filter module (an extension point: new filters read typed fields).
  - `IO.load_and_evaluate_mesh(input::PeriLabInput, path::String, ranksize::Int64)` (no `params`).

- [ ] **Step 1: Update the bond filter tests (they become the failing tests)**

In `test/unit_tests/IO/ut_bond_filter.jl`, wrap the two filter dicts: replace `filter = Dict("Center X" => 0.0,` with `filter = typed_section(PeriLab.InputDeck.BondFilterParams, Dict("Type" => "Disk", "Center X" => 0.0,` and `filter = Dict("Lower Left Corner X" => -0.5,` with `filter = typed_section(PeriLab.InputDeck.BondFilterParams, Dict("Type" => "Rectangular_Plane", "Lower Left Corner X" => -0.5,`, and add one closing `)` at the end of each of the two dict expressions. Any other `Dict(` passed as `filter` or as `params` to `apply_bond_filters` in this file is converted the same way (bond filter dicts → `typed_section(PeriLab.InputDeck.BondFilterParams, …)` with a `"Type"`; a `params` dict with `"Discretization" => Dict("Bond Filters" => F)` → `typed_section(Dict{String,PeriLab.InputDeck.BondFilterParams}, F)`).

- [ ] **Step 2: Run the test to verify it fails**

Run the single-file command for `unit_tests/IO/ut_bond_filter.jl`.
Expected: FAIL with `MethodError: no method matching run_bond_filter(::Int64, ::Matrix{Float64}, ::PeriLab.InputDeck.BondFilterParams, …)`.

- [ ] **Step 3: Write the implementation**

`src/IO/bond_filter.jl`: replace `using ...Parameter_Handling: get_bond_filters` with `using ....InputDeck: BondFilterParams` (dots of its `PeriLabExceptions` import); change the signature's `params::Dict,` to `filters::Dict{String,BondFilterParams},`; replace

```julia
    bond_filters = get_bond_filters(params)
```

with nothing (delete), `if bond_filters[1]` with `if !isempty(filters)`, both `bond_filters[2]` with `filters`, both `get(filter, "Allow Contact", false)` with `filter.allow_contact`, and `filter["Type"]` (both occurrences) with `filter.type`. Update the docstring.

Filter modules — add `using .....InputDeck: BondFilterParams` next to their `Data_Manager` import, change `filter::Dict,` to `filter::BondFilterParams,` in `run_bond_filter` (and in its docstring), and replace field reads:
- `disk_filter.jl`: `filter["Center X"]` → `filter.center_x`, `filter["Center Y"]` → `filter.center_y`, `filter["Center Z"]` → `filter.center_z`, `filter["Normal X"]` → `filter.normal_x`, `filter["Normal Y"]` → `filter.normal_y`, `filter["Normal Z"]` → `filter.normal_z`, `filter["Radius"]` → `filter.radius`.
- `rectangular_filter.jl`: `filter["Normal X"]` → `filter.normal_x`, `filter["Normal Y"]` → `filter.normal_y`, `filter["Normal Z"]` → `filter.normal_z`, `filter["Lower Left Corner X"]` → `filter.lower_left_corner_x` (Y, Z likewise: `lower_left_corner_y`, `lower_left_corner_z`), `filter["Bottom Unit Vector X"]` → `filter.bottom_unit_vector_x` (Y, Z likewise), `filter["Bottom Length"]` → `filter.bottom_length`, `filter["Side Length"]` → `filter.side_length`.
- `filter_template.jl`: signature and docstring only.

After the edits `grep -n 'filter\["' src/IO/Bond_Filter/*.jl src/IO/bond_filter.jl` must print nothing.

`src/IO/mesh_data.jl`: in `load_and_evaluate_mesh` remove the `params::Dict,` parameter (signature becomes `(input::PeriLabInput, path::String, ranksize::Int64)`) and replace the `apply_bond_filters(nlist, mesh, params, dof)` arguments with `apply_bond_filters(nlist, mesh, discretization.bond_filters, dof)`; in `init_data` call it as `load_and_evaluate_mesh(input, path, size)`. Confirm with `grep -n 'params' src/IO/mesh_data.jl` that the only remaining `params` uses are `init_data`'s signature, `contact_basis(params, …)`, `return params` and docstrings.

- [ ] **Step 4: Run tests to verify they pass**

Run the single-file command for `unit_tests/IO/ut_bond_filter.jl` and `unit_tests/IO/ut_mesh_data.jl` — Expected: PASS.
Run: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl` — Expected: PASS.
Run (background, ≈35 min): `julia --project=. -e 'using Pkg; Pkg.test()' > /tmp/full_2b2_t4.log 2>&1; tail -5 /tmp/full_2b2_t4.log` — Expected: `Testing PeriLab tests passed`.

- [ ] **Step 5: Commit**

```bash
git add src test
git -c user.name="Jan-Timo Hesse" -c user.email="jthesse@gmx.de" commit -m "Bond filters read typed filter parameters"
```
