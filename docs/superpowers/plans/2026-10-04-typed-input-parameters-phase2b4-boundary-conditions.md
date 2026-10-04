# Typed Input Parameters — Phase 2b-4 (Boundary Conditions Consume Typed Input) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Boundary conditions are built from the typed `BoundaryConditionParams` into a small runtime struct `BoundaryCondition` (resolved field, local nodes, cached value) instead of copying and mutating the input deck dicts; `Step ID` and `+`-joined node set lists are parsed once.

**Architecture:** Input logic (`Step ID` list, node set names) lives next to the type in `InputDeck`, with an input check for malformed `Step ID` lists. `BC_manager.jl` gets `mutable struct BoundaryCondition` — runtime state, not input: the value starts as a number or expression string and is replaced by the compiled expression after the first evaluation, as today. `init_BCs` takes `Dict{String,BoundaryConditionParams}`; the solvers and the solver manager pass `Dict{String,BoundaryCondition}` instead of `Dict{Any,Any}`. Same rule as before: read fields directly; functions only for real logic.

**Tech Stack:** Julia 1.12, existing dependencies only.

**Spec:** `docs/superpowers/specs/2026-10-03-typed-input-parameters-design.md` — §3 ("Boundary conditions: `Step ID` is parsed into a `Vector{Int64}` once …; `Node Set` `+`-lists are parsed into `Vector{String}` once"), phase 2 of §5. Previous slices: 2b-1 solvers, 2b-2 mesh / blocks, 2b-3 outputs / computes. Remaining after this one: 2b-5 contact, FEM / coupling, surface correction.

## Global Constraints

- Julia `1.12`; no new packages.
- Branch `feature/typed-input-parameters`; one commit per task, message ending with the session's `Co-Authored-By` line; commit identity `-c user.name="Jan-Timo Hesse" -c user.email="jthesse@gmx.de"`. Never merge.
- Simulation results must not change: every fullscale test must still pass.
- `get_bc_definitions` (Dict getter) stays unchanged (phase 4); it is the reference in the equivalence test.
- No thin wrappers around single fields.
- Read the full-suite result before recording a task as complete or committing it.
- Relative imports: write the `InputDeck` import with the same number of dots as the file's existing `PeriLabExceptions` import.
- Test commands (from the repository root):
  - Input tests: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl`
  - Single existing test file as `runtests.jl` runs it:
    `cd test && julia --project=.. -e 'using Test, Logging, MPI; import PeriLab; Logging.disable_logging(Logging.Warn); include("helper.jl"); MPI.Init(); @testset "t" begin include("unit_tests/<path>.jl") end'`
  - Full suite (≈35 min, background, output to a file): `julia --project=. -e 'using Pkg; Pkg.test()'`

## Review Focus

1. `Step ID` given as an integer (`2`) or a list (`"1,3"`, also with spaces `"1, 3"`) must restrict the BC to those steps; a single-solver run (step `-1`) applies every BC — tests in Task 1 ("step ids") and Task 2 ("step filter").
2. A `Node Set` joining several sets (`"top + left"`) must apply to the nodes of all of them, in order — tests in Task 1 ("node set names") and Task 2 (existing `ut_boundary_condition` expectations).
3. A BC without `Type` must still default to Dirichlet with the warning — test in Task 2 (existing `ut_apply_bc` decks have no `Type`; new testset "type defaults to Dirichlet").
4. A BC value evaluated once as an expression must be reused as the compiled function in later steps (no re-parse) — test in Task 2 ("value is compiled once").
5. A malformed `Step ID` (`"1,a"`) must be an input error, not a parse crash at the step change — test in Task 1 ("step ids").

Behaviour change (intended bug fix, decided in this plan): a BC on coordinate `z` in a 2D run used to `break` out of the BC loop, silently dropping every BC that came after it in dictionary order. It now skips only that BC (with the same warning). Dictionary order was not defined, so the old outcome was arbitrary.

---

### Task 1: Typed BC input helpers and `Step ID` check, proven equivalent on every shipped deck

**Files:**
- Modify: `src/Support/Parameters/Input/conditions.jl`, `InputDeck.jl`
- Create: `test/unit_tests/Support/Parameters/Input/ut_bc_typed.jl`
- Modify: `test/unit_tests/Support/Parameters/Input/input_tests.jl`

**Interfaces:**
- Consumes: phase 2a `BoundaryConditionParams`, `check!`.
- Produces (exported from `InputDeck`):
  - `bc_node_set_names(bc::BoundaryConditionParams)::Vector{String}` — the `+`-separated names, stripped.
  - `bc_step_ids(bc::BoundaryConditionParams)::Union{Nothing,Vector{Int64}}` — `nothing` if not given; `[id]` for an integer; parsed comma-separated list for a string.
  - `check!(::BoundaryConditionParams)`: a string `Step ID` must be a comma-separated list of integers.
  - `BoundaryConditionParams` added to the `export` statement.

- [ ] **Step 1: Write the failing tests**

`test/unit_tests/Support/Parameters/Input/ut_bc_typed.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const PS = PeriLab.ParameterSpec
const ID = PeriLab.InputDeck

function ut_bc(raw)
    ctx = PS.ParseContext()
    value = PS.convert_value(ID.BoundaryConditionParams, raw, "BC", ctx)
    return value, ctx
end

ut_bc_raw(extra...) = Dict{String,Any}("Variable" => "Displacements", "Node Set" => "Set-1",
                                       "Value" => 0.0, extra...)

@testset "node set names" begin
    bc, _ = ut_bc(ut_bc_raw("Node Set" => "top + left+right"))
    @test ID.bc_node_set_names(bc) == ["top", "left", "right"]
    bc, _ = ut_bc(ut_bc_raw())
    @test ID.bc_node_set_names(bc) == ["Set-1"]
end

@testset "step ids" begin
    @test ID.bc_step_ids(first(ut_bc(ut_bc_raw()))) === nothing
    @test ID.bc_step_ids(first(ut_bc(ut_bc_raw("Step ID" => 2)))) == [2]
    @test ID.bc_step_ids(first(ut_bc(ut_bc_raw("Step ID" => "1,3")))) == [1, 3]
    @test ID.bc_step_ids(first(ut_bc(ut_bc_raw("Step ID" => "1, 3")))) == [1, 3]
    _, ctx = ut_bc(ut_bc_raw("Step ID" => "1,a"))
    @test ctx.errors[1].path == "BC.\"Step ID\""
    @test ctx.errors[1].message ==
          "expected an integer or a comma-separated list of integers, got \"1,a\""
end

const UT_BC_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", "..", ".."))
const UT_BC_SKIP = ("test/fullscale_tests/test_PD_Solid_Elastic/strain_xx_external.yaml",)

function ut_bc_decks()
    decks = String[]
    for root in ("examples", "test"), (dir, _, files) in walkdir(joinpath(UT_BC_ROOT, root))
        append!(decks, [joinpath(dir, f) for f in files if endswith(f, ".yaml")])
    end
    return sort!(decks)
end

function ut_bc_compare_all_decks()
    for file in ut_bc_decks()
        relpath(file, UT_BC_ROOT) in UT_BC_SKIP && continue
        raw = PeriLab.IO.read_input(file)
        (raw isa AbstractDict && haskey(raw, "PeriLab")) || continue
        deck = raw["PeriLab"]
        input, _ = ID.read_input(deck, dirname(file))
        reference = PeriLab.Parameter_Handling.get_bc_definitions(deck)
        @test sort(collect(keys(input.sections.boundary_conditions))) ==
              sort(string.(collect(keys(reference))))
        for (name, bc) in input.sections.boundary_conditions
            old = reference[name]
            @test ID.bc_node_set_names(bc) == String.(strip.(split(old["Node Set"], "+")))
            if haskey(old, "Step ID")
                @test ID.bc_step_ids(bc) ==
                      parse.(Int64, strip.(split(string(old["Step ID"]), ",")))
            else
                @test ID.bc_step_ids(bc) === nothing
            end
        end
    end
end

@testset "BC helpers equal the Dict definitions on every shipped deck" begin
    ut_bc_compare_all_decks()
end
```

Update `input_tests.jl`: add `"ut_bc_typed.jl"` after `"ut_outputs_typed.jl"` in the list.

- [ ] **Step 2: Run tests to verify they fail**

Run: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl`
Expected: FAIL in `ut_bc_typed.jl` with `UndefVarError: bc_node_set_names not defined in PeriLab.InputDeck`.

- [ ] **Step 3: Write the implementation**

Append to `src/Support/Parameters/Input/conditions.jl`:

```julia
function check!(bc::BoundaryConditionParams, path::String, ctx::ParseContext)
    if bc.step_id isa String &&
       any(part -> tryparse(Int64, strip(part)) === nothing, split(bc.step_id, ","))
        add_error!(ctx, join_path(path, "Step ID"),
                   "expected an integer or a comma-separated list of integers, got \"$(bc.step_id)\"")
    end
    return nothing
end

"Names of the node sets of a boundary condition (`Node Set` joined with `+`)."
bc_node_set_names(bc::BoundaryConditionParams) = String.(strip.(split(bc.node_set, "+")))

"Solver steps a boundary condition applies to, or `nothing` (all steps)."
function bc_step_ids(bc::BoundaryConditionParams)
    bc.step_id === nothing && return nothing
    bc.step_id isa Int64 && return [bc.step_id]
    return parse.(Int64, strip.(split(bc.step_id, ",")))
end
```

In `src/Support/Parameters/Input/InputDeck.jl`, append `, BoundaryConditionParams, bc_node_set_names, bc_step_ids` to the `export` statement (after `ComputeClassParams`).

- [ ] **Step 4: Run tests to verify they pass**

Run: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl`
Expected: PASS — including the golden deck test, which now applies the `Step ID` check to every shipped deck.

- [ ] **Step 5: Commit**

```bash
git add src/Support/Parameters/Input test/unit_tests/Support/Parameters/Input
git -c user.name="Jan-Timo Hesse" -c user.email="jthesse@gmx.de" commit -m "Typed BC input helpers; Step ID check"
```

---

### Task 2: `BoundaryCondition` runtime struct; BC manager, solvers and solver manager use it

**Files:**
- Modify: `src/Core/BC_manager.jl`
- Modify: `src/Core/Solver/Solver_manager.jl`, `Verlet_solver.jl`, `Static_solver.jl`, `Matrix_linear_static.jl`, `Matrix_Verlet.jl`, `Newmark.jl`
- Modify: `test/unit_tests/Core/ut_BC_manager.jl`

**Interfaces:**
- Consumes: Task 1 `bc_node_set_names`, `bc_step_ids`, `BoundaryConditionParams`; `typed_section` (tests).
- Produces (in `Boundary_Conditions`, exported: `BoundaryCondition`):
  - `mutable struct BoundaryCondition` with fields `variable::String`, `time::String` (`"NP1"` or `"Constant"`), `type::String`, `initial::Bool`, `coordinate::Union{Nothing,String}`, `node_set::Vector{Int64}` (local node ids), `value::Any` (number, expression string, or compiled function).
  - `boundary_condition(bcs_in::Dict{String,BoundaryConditionParams})::Dict{String,Vector{Int64}}` — local node ids per BC; aborts `"Node Set '<name>' is missing"`.
  - `check_valid_bcs(bcs_in::Dict{String,BoundaryConditionParams}, node_sets::Dict{String,Vector{Int64}})::Dict{String,BoundaryCondition}` — step filter, 2D `z` check, `Type` default, field resolution.
  - `init_BCs(bcs_in::Dict{String,BoundaryConditionParams})::Dict{String,BoundaryCondition}`.
  - `find_bc_free_dof(bcs::Dict{String,BoundaryCondition})`, `apply_bc_dirichlet(allowed_variables::Vector{String}, bcs::Dict{String,BoundaryCondition}, time, step_time)`, `apply_bc_neumann(bcs::Dict{String,BoundaryCondition}, time, step_time)`.
  - Solvers: `bcs::Dict{String,BoundaryCondition}` instead of `bcs::Dict{Any,Any}` in `init_solver` / `run_solver`; `Solver_Manager.solver(…, bcs::Dict{String,BoundaryCondition}, …)`.

- [ ] **Step 1: Update `ut_BC_manager.jl` to the new API (it becomes the failing test)**

Mechanical conversions in `test/unit_tests/Core/ut_BC_manager.jl`:

```bash
F=test/unit_tests/Core/ut_BC_manager.jl
perl -0pi -e 's/params = Dict\("Boundary Conditions" => Dict\(((?:(?!\n    \S).)*?)\)\)\n/params = typed_section(Dict{String,PeriLab.InputDeck.BoundaryConditionParams},\n                           Dict($1))\n/gs; s/^    params = Dict\(\)$/    params = typed_section(Dict{String,PeriLab.InputDeck.BoundaryConditionParams}, Dict())/m; s/bcs\["(BC_\d)"\]\["Variable"\]/bcs["$1"].variable/g; s/bcs\["(BC_\d)"\]\["Coordinate"\]/bcs["$1"].coordinate/g; s/bcs\["(BC_\d)"\]\["Value"\]/bcs["$1"].value/g; s/bcs\["(BC_\d)"\]\["Time"\]/bcs["$1"].time/g; s/bcs\["(BC_\d)"\]\["Node Set"\]/bcs["$1"].node_set/g; s/check_valid_bcs\(bcs\)/check_valid_bcs(params, bcs)/g' $F
grep -n 'Dict("Boundary Conditions"\|\]\["' $F
```

The final `grep` must print nothing. Then, in `@testset "ut_boundary_condition"`, `boundary_condition` now returns only the local node ids per BC (the other entries stay in the typed input), so replace

```julia
    # params representation
    @test bcs["BC_1"].variable == "Forces"
    @test bcs["BC_1"].coordinate == "x"
    @test bcs["BC_1"].value == "20*t"
    @test bcs["BC_1"].node_set == [1, 3, 4]
    @test bcs["BC_2"].variable == "Displacements"
    @test bcs["BC_2"].coordinate == "z"
    @test bcs["BC_2"].value == "0"
    @test bcs["BC_2"].node_set == [4, 2, 7, 10]
```

with

```julia
    # local node ids per boundary condition, in node set order
    @test bcs["BC_1"] == [1, 3, 4]
    @test bcs["BC_2"] == [4, 2, 7, 10]
```

and append to the file:

```julia
@testset "type defaults to Dirichlet" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(3)
    PeriLab.Data_Manager.set_dof(2)
    PeriLab.Data_Manager.set_nset("Nset_1", [1, 2])
    PeriLab.Data_Manager.set_glob_to_loc(Dict(1 => 1, 2 => 2, 3 => 3))
    PeriLab.Data_Manager.create_node_vector_field("Displacements", Float64, 2)
    params = typed_section(Dict{String,PeriLab.InputDeck.BoundaryConditionParams},
                           Dict("BC_1" => Dict("Variable" => "Displacements",
                                               "Node Set" => "Nset_1",
                                               "Coordinate" => "x", "Value" => 1.0)))
    bcs = PeriLab.Solver_Manager.Boundary_Conditions.init_BCs(params)
    @test bcs["BC_1"].type == "Dirichlet"
    @test !bcs["BC_1"].initial
    @test bcs["BC_1"].time == "NP1"
end

@testset "step filter" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(3)
    PeriLab.Data_Manager.set_dof(2)
    PeriLab.Data_Manager.set_nset("Nset_1", [1, 2])
    PeriLab.Data_Manager.set_glob_to_loc(Dict(1 => 1, 2 => 2, 3 => 3))
    PeriLab.Data_Manager.create_node_vector_field("Displacements", Float64, 2)
    params = typed_section(Dict{String,PeriLab.InputDeck.BoundaryConditionParams},
                           Dict("always" => Dict("Variable" => "Displacements",
                                                 "Node Set" => "Nset_1",
                                                 "Coordinate" => "x", "Value" => 0.0),
                                "steps_1_3" => Dict("Variable" => "Displacements",
                                                    "Node Set" => "Nset_1",
                                                    "Coordinate" => "y", "Value" => 0.0,
                                                    "Step ID" => "1, 3")))
    for (step, expected) in ((-1, ["always", "steps_1_3"]), (1, ["always", "steps_1_3"]),
                             (2, ["always"]), (3, ["always", "steps_1_3"]))
        PeriLab.Data_Manager.set_step(step)
        bcs = PeriLab.Solver_Manager.Boundary_Conditions.init_BCs(params)
        @test sort(collect(keys(bcs))) == expected
    end
end

@testset "value is compiled once" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(3)
    PeriLab.Data_Manager.set_dof(2)
    PeriLab.Data_Manager.set_nset("Nset_1", [1, 2])
    PeriLab.Data_Manager.set_glob_to_loc(Dict(1 => 1, 2 => 2, 3 => 3))
    PeriLab.Data_Manager.create_constant_node_vector_field("Coordinates", Float64, 2)
    displacements_N, displacements_NP1 = PeriLab.Data_Manager.create_node_vector_field("Displacements",
                                                                                      Float64,
                                                                                      2)
    params = typed_section(Dict{String,PeriLab.InputDeck.BoundaryConditionParams},
                           Dict("BC_1" => Dict("Variable" => "Displacements",
                                               "Node Set" => "Nset_1",
                                               "Coordinate" => "x", "Value" => "10*t")))
    bcs = PeriLab.Solver_Manager.Boundary_Conditions.init_BCs(params)
    @test bcs["BC_1"].value == "10*t"
    PeriLab.Solver_Manager.Boundary_Conditions.apply_bc_dirichlet(["Displacements"], bcs, 1.0,
                                                                  1.0)
    @test bcs["BC_1"].value isa Function
    compiled = bcs["BC_1"].value
    PeriLab.Solver_Manager.Boundary_Conditions.apply_bc_dirichlet(["Displacements"], bcs, 2.0,
                                                                  2.0)
    @test bcs["BC_1"].value === compiled
    @test displacements_NP1[1:2, 1] == [20.0, 20.0]
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run the single-file command for `unit_tests/Core/ut_BC_manager.jl`.
Expected: FAIL with `MethodError`s for `boundary_condition(::Dict{String, PeriLab.InputDeck.BoundaryConditionParams})`, `init_BCs(::Dict{String, PeriLab.InputDeck.BoundaryConditionParams})` and `check_valid_bcs(::Dict{String, …}, …)`.

- [ ] **Step 3: Rewrite the BC manager**

In `src/Core/BC_manager.jl`:
- replace `using ...Parameter_Handling: get_bc_definitions` with `using ...InputDeck: BoundaryConditionParams, bc_node_set_names, bc_step_ids` (same dots as its `PeriLabExceptions` import), and add `export BoundaryCondition` next to the existing exports;
- replace the functions `find_bc_free_dof`, `check_valid_bcs`, `init_BCs` and `boundary_condition` (each with its docstring) with:

```julia
"""
    BoundaryCondition

A boundary condition prepared for one solver step: the resolved field
(`variable`, `time` = `"NP1"` or `"Constant"`), its `type`, the local nodes it
acts on, and its `value` — a number or an expression string from the input
deck, replaced by the compiled expression after the first evaluation.
"""
mutable struct BoundaryCondition
    variable::String
    time::String
    type::String
    initial::Bool
    coordinate::Union{Nothing,String}
    node_set::Vector{Int64}
    value::Any
end

"""
    find_bc_free_dof(bcs)

Finds all dof without a displacement boundary condition and stores them in the
Data_Manager.
"""
function find_bc_free_dof(bcs::Dict{String,BoundaryCondition})
    nnodes = Data_Manager.get_nnodes()
    dof = Data_Manager.get_dof()
    bc_free_dof = vec([(i, j) for i in 1:nnodes, j in 1:dof])
    dof_mapping = Dict{String,Int8}("x" => 1, "y" => 2, "z" => 3)
    for bc in values(bcs)
        if bc.variable == "Displacements" && bc.type == "Dirichlet"
            act = Vector{Tuple{Int64,Int64}}([(node, dof_mapping[bc.coordinate])
                                              for node in bc.node_set])
            bc_free_dof = setdiff(bc_free_dof, act)
        end
    end
    Data_Manager.set_bc_free_dof([(t[1] + (t[2] - 1) * nnodes) for t in bc_free_dof])
end

"""
    boundary_condition(bcs_in)

Local node ids of every boundary condition (all node sets of its `Node Set`
list, in order). Aborts if a node set does not exist.
"""
function boundary_condition(bcs_in::Dict{String,BoundaryConditionParams})
    nsets = Data_Manager.get_nsets()
    node_sets = Dict{String,Vector{Int64}}()
    for (name, bc) in bcs_in
        nodes = Int64[]
        for node_set_name in bc_node_set_names(bc)
            if !haskey(nsets, node_set_name)
                @abort "Node Set '$node_set_name' is missing"
                return
            end
            append!(nodes, Data_Manager.get_local_nodes(nsets[node_set_name]))
        end
        node_sets[name] = nodes
    end
    return node_sets
end

"""
    check_valid_bcs(bcs_in, node_sets)

The boundary conditions active in the current solver step, with their field
resolved. A `z` condition in a 2D run is skipped with a warning; a missing
`Type` means Dirichlet. Aborts if the field does not exist.
"""
function check_valid_bcs(bcs_in::Dict{String,BoundaryConditionParams},
                         node_sets::Dict{String,Vector{Int64}})
    working_bcs = Dict{String,BoundaryCondition}()
    step = Data_Manager.get_step()
    dof = Data_Manager.get_dof()
    for (name, bc) in bcs_in
        steps = bc_step_ids(bc)
        if steps !== nothing && step != -1 && !(step in steps)
            continue
        end
        if bc.coordinate == "z" && dof < 3
            @warn "Boundary condition $name is not possible with $dof DOF"
            continue
        end
        type = bc.type
        if type === nothing
            type = "Dirichlet"
            @warn "Missing boundary condition type for $name. Assuming Dirichlet."
        end
        time = nothing
        for data_entry in Data_Manager.get_all_field_keys()
            if bc.variable * "NP1" == data_entry
                time = "NP1"
                break
            elseif bc.variable == data_entry
                time = "Constant"
                break
            end
        end
        if time === nothing
            @abort "Boundary condition $name is not valid: Variable $(bc.variable) not found. Please check if the physical model is activated."
            return
        end
        working_bcs[name] = BoundaryCondition(bc.variable, time, type, type == "Initial",
                                              bc.coordinate, node_sets[name], bc.value)
    end
    return working_bcs
end

"""
    init_BCs(bcs_in)

The boundary conditions of the current solver step.
"""
function init_BCs(bcs_in::Dict{String,BoundaryConditionParams})
    return check_valid_bcs(bcs_in, boundary_condition(bcs_in))
end
```

- in `apply_bc_dirichlet` and `apply_bc_neumann`: change the `bcs::Dict` parameter to `bcs::Dict{String,BoundaryCondition}`, replace `for name in keys(bcs)` + `bc = bcs[name]` with `for (name, bc) in bcs`, and replace every `bc["Type"]` → `bc.type`, `bc["Variable"]` → `bc.variable`, `bc["Time"]` → `bc.time`, `bc["Coordinate"]` → `bc.coordinate`, `bc["Node Set"]` → `bc.node_set`, `bc["Initial"]` → `bc.initial`, and `bc["Value"] = eval_bc!(… bc["Value"], …)` → `bc.value = eval_bc!(… bc.value, …)`. (`haskey(dof_mapping, bc.coordinate)` with `nothing` is `false` and reaches the existing "Coordinate … must be x,y or z" abort.) Update both docstrings (`bcs::Dict{String,BoundaryCondition}`).

After the edits `grep -n 'bc\["\|bcs\[bc\]' src/Core/BC_manager.jl` must print nothing.

- [ ] **Step 4: Solvers and solver manager pass `Dict{String,BoundaryCondition}`**

In each of `Verlet_solver.jl`, `Static_solver.jl`, `Matrix_linear_static.jl`, `Matrix_Verlet.jl`, `Newmark.jl`: add `BoundaryCondition` to the `using ..Boundary_Conditions:` import, and replace every `bcs::Dict{Any,Any}` (signatures and docstrings) with `bcs::Dict{String,BoundaryCondition}`.

In `Solver_manager.jl`: change `using .Boundary_Conditions: init_BCs` to `using .Boundary_Conditions: init_BCs, BoundaryCondition`; replace `@timeit "init_BCs" bcs=init_BCs(params)` with `@timeit "init_BCs" bcs=init_BCs(input.sections.boundary_conditions)`; replace every `bcs::Dict{Any,Any}` (the `solver` signature and the docstrings) with `bcs::Dict{String,BoundaryCondition}`.

Check: `grep -rn 'bcs::Dict{Any,Any}' src` prints nothing.

- [ ] **Step 5: Run tests to verify they pass**

Run the single-file command for `unit_tests/Core/ut_BC_manager.jl` — Expected: PASS.
Run: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl` — Expected: PASS.
Run (background, ≈35 min): `julia --project=. -e 'using Pkg; Pkg.test()' > /tmp/full_2b4_t2.log 2>&1; tail -5 /tmp/full_2b4_t2.log` — read the result; Expected: `Testing PeriLab tests passed` (the fullscale tests run every BC kind — Dirichlet / Neumann / Initial, expressions, multistep `Step ID`, `Forces` / `Force Densities` — against reference results).

- [ ] **Step 6: Commit**

```bash
git add src test
git -c user.name="Jan-Timo Hesse" -c user.email="jthesse@gmx.de" commit -m "Boundary conditions built from typed input into a runtime struct"
```
