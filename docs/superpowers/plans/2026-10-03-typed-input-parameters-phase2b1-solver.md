# Typed Input Parameters — Phase 2b-1 (Solvers Consume Typed Solver Options) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** The five solvers, the solver manager and model reduction read their options directly from the typed `SolverParams` built in phase 2a instead of the params `Dict`, and the typed `PeriLabInput` is passed from `PeriLab.run` down to the solver manager.

**Architecture:** Typed structs are read by field access; functions exist only where there is real logic. The five solver option structs share an abstract supertype `SolverOptions` whose documented common fields are `safety_factor`, `fixed_dt` and `numerical_damping`; `active_options(s)` returns the selected one and `solver_name` is derived from its type. Start and end time, which depend on the current simulation time, become pure functions `start_time(s, t)` / `end_time(s, t)`. Step selection (`solver_steps`, `solver_step`) and `model_options` live next to the types in `InputDeck`. The old `Dict` getters stay untouched (deleted in phase 4) and serve as the reference in an equivalence test over every shipped deck.

**Tech Stack:** Julia 1.12, existing dependencies only.

**Spec:** `docs/superpowers/specs/2026-10-03-typed-input-parameters-design.md` — §3 step 4 ("Use") for the solver part; phase 2 of §5. Phase 2a (structs, `read_input`) is in `src/Support/Parameters/Input/`.

**Phase 2b slices** (each its own plan, in this order): **2b-1 solvers (this plan)** → 2b-2 mesh, discretization and blocks → 2b-3 outputs and compute classes → 2b-4 boundary conditions → 2b-5 contact, FEM / coupling, surface correction, influence function. Rule for all slices: read fields directly; add a function only for real logic (derived values, runtime-dependent values, selection); BCs, outputs and computes get runtime structs built from the immutable parameters.

## Global Constraints

- Julia `1.12`; no new packages.
- Branch `feature/typed-input-parameters`; one commit per task, message ending with the session's `Co-Authored-By` line; commit identity `-c user.name="Jan-Timo Hesse" -c user.email="jthesse@gmx.de"` (the container has no git identity). Never merge.
- Simulation results must not change: every fullscale test must still pass.
- The `Dict` getters in `parameter_handling_solver.jl` stay unchanged (deleted in phase 4); no typed methods are added to them.
- No thin wrapper functions around single fields: solver code reads `params.maximum_damage`, `params.verlet.safety_factor`, etc. directly.
- Relative imports: in each file, use the same number of dots for `InputDeck` as that file already uses for `Data_Manager` (e.g. `using ...Data_Manager` → `using ...InputDeck: SolverParams`).
- Every new source file starts with the SPDX header (`# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>` / `#` / `# SPDX-License-Identifier: BSD-3-Clause`).
- Test commands (from the repository root):
  - Spec tests: `julia --project=. test/unit_tests/Support/Parameters/Spec/run_spec_tests.jl`
  - Input tests: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl`
  - Full suite (≈35 min, run in the background, output to a file): `julia --project=. -e 'using Pkg; Pkg.test()'`

## Review Focus

1. A Static solver deck that gives `Fixed dt` and `Number of Steps` must still warn and use `Fixed dt`; one with neither must still run exactly one step — test in Task 1 ("number of steps given or not").
2. A multistep deck whose later steps have no `Initial Time` must start each step at the current time; `Additional Time` must count from the current time — test in Task 1 ("start and end time depend on the current time").
3. `Reduction Blocks` given as an integer (`2`) or a list (`"2, 3"`) must both work — test in Task 2 ("reduction blocks forms").
4. `Model Reduction: true` under Newmark (a path that crashed before, calling an unimported `init_reduce_model` with the wrong argument) must be rejected during input validation with a clear message — test in Task 2 ("Newmark model reduction is rejected clearly").
5. Multistep `Step ID`s given out of order (3, 1, 2) must run 1, 2, 3 — test in Task 1 ("solver steps are sorted").

---

## File Structure

| File | Change |
|---|---|
| `src/Support/Parameters/Spec/params_macro.jl` | `@params` accepts a supertype |
| `src/Support/Parameters/Input/solver.jl` | `SolverOptions`; `number_of_steps::Union{Nothing,Int64}`; `reduction_blocks::Union{Nothing,Int64,String}`; `active_options`, `solver_name`, `start_time`, `end_time`, `model_options`; Newmark check |
| `src/Support/Parameters/Input/input.jl` | `solver_steps`, `solver_step` |
| `src/Support/Parameters/Input/InputDeck.jl` | imports, exports |
| `src/Core/Model_reduction/Model_reduction.jl` | `parse_reduction_blocks` / `init_reduce_model` take `ModelReductionParams` |
| `src/Core/Solver/{Verlet_solver,Static_solver,Matrix_linear_static,Matrix_Verlet,Newmark}.jl` | `init_solver(solver_options, params::SolverParams, bcs, block_nodes)`, field access |
| `src/Core/Solver/Solver_manager.jl` | typed solver selection; `init(params, input, step_id)` |
| `src/Support/Parameters/parameter_handling.jl`, `src/IO/read_inputdeck.jl`, `src/IO/IO.jl`, `src/PeriLab.jl` | typed input plumbing |
| `test/unit_tests/Support/Parameters/Spec/ut_extensions.jl` | supertype test |
| `test/unit_tests/Support/Parameters/Input/ut_solver_typed.jl`, `ut_model_reduction_params.jl` | new tests |
| `test/unit_tests/Support/Parameters/Input/ut_solver.jl`, `ut_validate_yaml.jl`, `input_tests.jl` | updated tests |

---

### Task 1: Typed solver helpers next to the types, proven equivalent on every shipped deck

**Files:**
- Modify: `src/Support/Parameters/Spec/params_macro.jl`
- Modify: `src/Support/Parameters/Input/solver.jl`
- Modify: `src/Support/Parameters/Input/input.jl`
- Modify: `src/Support/Parameters/Input/InputDeck.jl`
- Modify: `test/unit_tests/Support/Parameters/Spec/ut_extensions.jl`
- Modify: `test/unit_tests/Support/Parameters/Input/ut_solver.jl`
- Create: `test/unit_tests/Support/Parameters/Input/ut_solver_typed.jl`
- Modify: `test/unit_tests/Support/Parameters/Input/input_tests.jl`

**Interfaces:**
- Consumes: phase 2a `SolverParams`, the five option structs, `PeriLabInput`, `read_input`; existing `Dict` getters (reference only).
- Produces:
  - `@params struct Name <: Super … end` is supported.
  - `abstract type SolverOptions end`; `VerletParams`, `StaticParams`, `LinearStaticMatrixParams`, `VerletMatrixParams`, `NewmarkParams` `<: SolverOptions`, each with fields `safety_factor::Float64`, `fixed_dt::Float64`, `numerical_damping::Float64`.
  - `SolverParams.number_of_steps::Union{Nothing,Int64}` (default `nothing`; 1 step if not given).
  - `ModelReductionParams.reduction_blocks::Union{Nothing,Int64,String}`.
  - In `InputDeck` (exported): `active_options(s::SolverParams)::SolverOptions`; `solver_name(o::SolverOptions)::String` and `solver_name(s::SolverParams)::String`; `start_time(s::SolverParams, current_time::Float64)::Float64`; `end_time(s::SolverParams, current_time::Float64)::Float64`; `model_options(s::SolverParams)::Vector{String}`; `solver_steps(input::PeriLabInput)::Vector{Int64}` (`[-1]` for a single `Solver`); `solver_step(input::PeriLabInput, step_id::Int64)::SolverParams` (`-1` → `Solver`). `start_time` / `end_time` / `solver_step` abort (`@abort`) with the same messages as the `Dict` getters.

- [ ] **Step 1: Write the failing tests**

Append to `test/unit_tests/Support/Parameters/Spec/ut_extensions.jl`:

```julia
abstract type UTAbstractOptions end

PS.@params struct UTConcreteOptions <: UTAbstractOptions
    factor::Float64 = opt("Factor"; default = 1.0)
end

PS.@params struct UTDependentOptions <: UTAbstractOptions
    value::Dependent = req("Value")
end

@testset "@params accepts a supertype" begin
    @test UTConcreteOptions <: UTAbstractOptions
    @test UTDependentOptions{PS.Constant} <: UTAbstractOptions
    ctx = PS.ParseContext()
    p = PS.parse_section(UTConcreteOptions, Dict{String,Any}("Factor" => 2), "O", ctx)
    @test isempty(ctx.errors) && p.factor === 2.0
end
```

In `test/unit_tests/Support/Parameters/Input/ut_solver.jl`, in `@testset "Verlet solver with defaults"`, replace

```julia
    @test s.final_time === 1.0 && s.number_of_steps === 1 && s.maximum_damage === Inf
```

with

```julia
    @test s.final_time === 1.0 && s.number_of_steps === nothing && s.maximum_damage === Inf
```

`test/unit_tests/Support/Parameters/Input/ut_solver_typed.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const PS = PeriLab.ParameterSpec
const ID = PeriLab.InputDeck
const PH = PeriLab.Parameter_Handling

function ut_st_solver(dict)
    ctx = PS.ParseContext()
    s = PS.parse_section(ID.SolverParams, dict, "Solver", ctx)
    @test isempty(ctx.errors)
    return s
end

@testset "active options and solver name" begin
    s = ut_st_solver(Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                      "Static" => Dict{String,Any}("m" => 3)))
    @test ID.active_options(s) === s.static
    @test ID.active_options(s) isa ID.SolverOptions
    @test ID.solver_name(s) == "Static"
    @test ID.active_options(s).safety_factor === 1.0
    names = [ID.solver_name(ut_st_solver(Dict{String,Any}("Initial Time" => 0.0,
                                                          "Final Time" => 1.0,
                                                          key => Dict{String,Any}())))
             for key in ("Verlet", "Static", "Linear Static Matrix Based",
                         "Verlet Matrix Based", "Newmark")]
    @test names ==
          ["Verlet", "Static", "Linear Static Matrix Based", "Verlet Matrix Based", "Newmark"]
end

@testset "number of steps given or not" begin
    s = ut_st_solver(Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                      "Static" => Dict{String,Any}()))
    @test s.number_of_steps === nothing
    s = ut_st_solver(Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                      "Number of Steps" => 7,
                                      "Static" => Dict{String,Any}("Fixed dt" => 0.1)))
    @test s.number_of_steps === 7 && s.static.fixed_dt === 0.1
end

@testset "start and end time depend on the current time" begin
    s = ut_st_solver(Dict{String,Any}("Initial Time" => 0.5, "Final Time" => 2.0,
                                      "Verlet" => Dict{String,Any}()))
    @test ID.start_time(s, 0.0) === 0.5
    @test ID.start_time(s, 1.0) === 1.0          # already past the initial time
    @test ID.end_time(s, 1.0) === 2.0
    s = ut_st_solver(Dict{String,Any}("Additional Time" => 3.0, "Step ID" => 2,
                                      "Verlet" => Dict{String,Any}()))
    @test ID.start_time(s, 1.5) === 1.5          # later multistep step
    @test ID.end_time(s, 1.5) === 4.5
    @test_throws PeriLab.PeriLabError ID.start_time(s, 0.0)
end

@testset "model options keep the old order" begin
    s = ut_st_solver(Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                      "Damage Models" => true, "Thermal Models" => true,
                                      "Additive Models" => true,
                                      "Verlet" => Dict{String,Any}()))
    @test ID.model_options(s) ==
          ["Additive", "Damage", "Pre_Calculation", "Thermal", "Material"]
end

function ut_st_deck(; multistep = nothing)
    deck = Dict{String,Any}("Discretization" => Dict{String,Any}("Type" => "Text File",
                                                                 "Input Mesh File" => "m.txt"),
                            "Blocks" => Dict{String,Any}("b" => Dict{String,Any}("Block ID" => 1,
                                                                                 "Density" => 1.0,
                                                                                 "Horizon" => 1.0)),
                            "Models" => Dict{String,Any}())
    if multistep === nothing
        deck["Solver"] = Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                          "Verlet" => Dict{String,Any}())
    else
        deck["Multistep Solver"] = multistep
    end
    input, ctx = ID.read_input(deck)
    @test isempty(ctx.errors)
    return input
end

@testset "solver steps are sorted" begin
    step(id) = Dict{String,Any}("Step ID" => id, "Initial Time" => 0.0, "Final Time" => 1.0,
                                "Verlet" => Dict{String,Any}())
    input = ut_st_deck(multistep = Dict{String,Any}("c" => step(3), "a" => step(1),
                                                    "b" => step(2)))
    @test ID.solver_steps(input) == [1, 2, 3]
    @test ID.solver_step(input, 2).step_id === 2
    @test_throws PeriLab.PeriLabError ID.solver_step(input, 9)
    input = ut_st_deck()
    @test ID.solver_steps(input) == [-1]
    @test ID.solver_step(input, -1) === input.sections.solver
end

const UT_ST_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", "..", ".."))
const UT_ST_SKIP = ("test/fullscale_tests/test_PD_Solid_Elastic/strain_xx_external.yaml",)

function ut_st_decks()
    decks = String[]
    for root in ("examples", "test"), (dir, _, files) in walkdir(joinpath(UT_ST_ROOT, root))
        append!(decks, [joinpath(dir, f) for f in files if endswith(f, ".yaml")])
    end
    return sort!(decks)
end

@testset "typed solver options equal the Dict getters on every shipped deck" begin
    t = PeriLab.Data_Manager.get_current_time()
    for file in ut_st_decks()
        relpath(file, UT_ST_ROOT) in UT_ST_SKIP && continue
        raw = PeriLab.IO.read_input(file)
        (raw isa AbstractDict && haskey(raw, "PeriLab")) || continue
        deck = raw["PeriLab"]
        input, _ = ID.read_input(deck, dirname(file))
        steps = PH.get_solver_steps(deck)
        @test ID.solver_steps(input) == steps
        for step in steps
            d = step == -1 ? deck["Solver"] : PH.get_solver_params(deck, step)
            s = ID.solver_step(input, step)
            options = ID.active_options(s)
            calc = PH.get_calculation_options(d)
            @test ID.solver_name(s) == PH.get_solver_name(d)
            @test something(s.number_of_steps, 1) == PH.get_nsteps(d)
            @test options.safety_factor == PH.get_safety_factor(d)
            @test options.fixed_dt == PH.get_fixed_dt(d)
            @test options.numerical_damping == PH.get_numerical_damping(d)
            @test s.maximum_damage == PH.get_max_damage(d)
            @test ID.model_options(s) == PH.get_model_options(d)
            @test s.calculate_cauchy == calc["Calculate Cauchy"]
            @test s.calculate_von_mises_stress == calc["Calculate von Mises stress"]
            @test s.calculate_strain == calc["Calculate Strain"]
            @test ID.end_time(s, t) == PH.get_final_time(d)
            if haskey(d, "Initial Time") || t != 0.0     # otherwise both abort
                @test ID.start_time(s, t) == PH.get_initial_time(d)
            end
        end
    end
end
```

Update `input_tests.jl`:

```julia
for file in ["ut_mesh_blocks.jl", "ut_solver.jl", "ut_outputs_conditions_contact.jl",
             "ut_read_input.jl", "ut_golden_decks.jl", "ut_validate_yaml.jl",
             "ut_solver_typed.jl"]
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `julia --project=. test/unit_tests/Support/Parameters/Spec/run_spec_tests.jl`
Expected: FAIL in `ut_extensions.jl` with `ParamsDefinitionError: @params struct UTConcreteOptions <: UTAbstractOptions: write a plain name without type parameters or supertype`.
Run: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl`
Expected: FAIL — `ut_solver.jl` "Verlet solver with defaults" (`number_of_steps` is `1`) and `ut_solver_typed.jl` with `UndefVarError: active_options not defined in PeriLab.InputDeck`.

- [ ] **Step 3: `@params` accepts a supertype**

In `src/Support/Parameters/Spec/params_macro.jl`, inside `macro params`, replace

```julia
    name = structdef.args[2]
    if !(name isa Symbol)
        return _definition_error("@params struct $(name): write a plain name without type parameters or supertype")
    end
```

with

```julia
    name = structdef.args[2]
    supertype = nothing
    if name isa Expr && name.head === :<: && name.args[1] isa Symbol
        name, supertype = name.args[1], name.args[2]
    end
    if !(name isa Symbol)
        return _definition_error("@params struct $(name): write a plain name without type parameters")
    end
```

and replace

```julia
    head = isempty(typeparams) ? name : Expr(:curly, name, typeparams...)
```

with

```julia
    head = isempty(typeparams) ? name : Expr(:curly, name, typeparams...)
    supertype === nothing || (head = Expr(:<:, head, supertype))
```

Run: `julia --project=. test/unit_tests/Support/Parameters/Spec/run_spec_tests.jl`
Expected: PASS.

- [ ] **Step 4: Solver types and helpers**

In `src/Support/Parameters/Input/solver.jl`, replace the comment at the top

```julia
# Every solver option section carries "Safety Factor", "Fixed dt" and
# "Numerical Damping": the solver getters read them from the active section.
```

with

```julia
"""
    SolverOptions

Options section of one solver. Every subtype has the fields `safety_factor`,
`fixed_dt` (-1 if not given) and `numerical_damping`.
"""
abstract type SolverOptions end
```

and change the five option struct headers to:

```julia
@params struct VerletParams <: SolverOptions
@params struct StaticParams <: SolverOptions
@params struct LinearStaticMatrixParams <: SolverOptions
@params struct VerletMatrixParams <: SolverOptions
@params struct NewmarkParams <: SolverOptions
```

(`ModelReductionParams` and `SolverParams` get no supertype.) Replace

```julia
    number_of_steps::Int64 = opt("Number of Steps"; default = 1, min = 1)
```

with

```julia
    number_of_steps::Union{Nothing,Int64} = opt("Number of Steps"; default = nothing, min = 1,
                                                description = "1 if not given")
```

replace

```julia
    reduction_blocks::Union{Nothing,String} = opt("Reduction Blocks"; default = nothing)
```

with

```julia
    reduction_blocks::Union{Nothing,Int64,String} = opt("Reduction Blocks"; default = nothing,
                                                        description = "Block id or list, e.g. \"2, 3\"")
```

and append at the end of the file:

```julia
"The options section of the solver `s` selects (exactly one is given, see `check!`)."
function active_options(s::SolverParams)
    for options in (s.verlet, s.static, s.linear_static_matrix_based, s.verlet_matrix_based,
                    s.newmark)
        options === nothing || return options
    end
    throw(ArgumentError("no solver section given"))
end

solver_name(::VerletParams) = "Verlet"
solver_name(::StaticParams) = "Static"
solver_name(::LinearStaticMatrixParams) = "Linear Static Matrix Based"
solver_name(::VerletMatrixParams) = "Verlet Matrix Based"
solver_name(::NewmarkParams) = "Newmark"
solver_name(s::SolverParams) = solver_name(active_options(s))

"""
    start_time(s, current_time)

Start time of the solver step: `Initial Time`, or the current time if that is
later or `Initial Time` is not given (later steps of a multistep run).
"""
function start_time(s::SolverParams, current_time::Float64)
    if s.initial_time !== nothing
        return max(s.initial_time, current_time)
    end
    current_time != 0.0 && return current_time
    @abort "No initial time defined"
end

"""
    end_time(s, current_time)

End time of the solver step: `Final Time`, or `Additional Time` after the
current time.
"""
function end_time(s::SolverParams, current_time::Float64)
    s.final_time !== nothing && return s.final_time
    s.additional_time !== nothing && return current_time + s.additional_time
    @abort "No final time defined"
end

"Active model categories, in the order the models are evaluated."
function model_options(s::SolverParams)
    return [name
            for (name, used) in (("Additive", s.additive_models), ("Damage", s.damage_models),
                                 ("Pre_Calculation", s.pre_calculation_models),
                                 ("Thermal", s.thermal_models),
                                 ("Degradation", s.degradation_models),
                                 ("Material", s.material_models)) if used]
end
```

In `src/Support/Parameters/Input/input.jl`, append:

```julia
"Step IDs to run, sorted; `[-1]` for a single `Solver`."
function solver_steps(input::PeriLabInput)
    isempty(input.sections.multistep_solver) && return [-1]
    return sort!([step.step_id for step in values(input.sections.multistep_solver)])
end

"Solver parameters of step `step_id` (`-1`: the single `Solver`)."
function solver_step(input::PeriLabInput, step_id::Int64)
    step_id == -1 && return input.sections.solver
    for step in values(input.sections.multistep_solver)
        step.step_id == step_id && return step
    end
    @abort "Step ID $step_id not found"
end
```

In `src/Support/Parameters/Input/InputDeck.jl`, add after the `import ..ParameterSpec: check!` line:

```julia
using ..PeriLabExceptions: @abort
```

and replace

```julia
export read_input, PeriLabInput
```

with

```julia
export read_input, PeriLabInput, SolverParams, SolverOptions, ModelReductionParams,
       active_options, solver_name, start_time, end_time, model_options, solver_steps,
       solver_step
```

`PeriLab.jl` includes `IO/exceptions.jl` (module `PeriLabExceptions`) before `InputDeck.jl`, so the import resolves.

- [ ] **Step 5: Run tests to verify they pass**

Run: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl`
Expected: PASS. If the equivalence testset fails for a deck, the typed helper (not the test) is wrong unless the `Dict` getter throws for that deck — record any such case as a ruling.

- [ ] **Step 6: Commit**

```bash
git add src/Support/Parameters test/unit_tests/Support/Parameters
git -c user.name="Jan-Timo Hesse" -c user.email="jthesse@gmx.de" commit -m "Typed solver options: SolverOptions, time and step helpers"
```

---

### Task 2: Solvers and model reduction read `SolverParams` fields (temporary adapter in the solver manager)

**Files:**
- Modify: `src/Support/Parameters/Input/solver.jl` (Newmark check)
- Modify: `src/Core/Model_reduction/Model_reduction.jl`
- Modify: `src/Core/Solver/Verlet_solver.jl`, `Static_solver.jl`, `Matrix_linear_static.jl`, `Matrix_Verlet.jl`, `Newmark.jl`
- Modify: `src/Core/Solver/Solver_manager.jl`
- Create: `test/unit_tests/Support/Parameters/Input/ut_model_reduction_params.jl`
- Modify: `test/unit_tests/Support/Parameters/Input/input_tests.jl`

**Interfaces:**
- Consumes: Task 1 helpers and types; phase 1 `ParseContext`, `parse_section`, `report!`.
- Produces:
  - `init_solver(solver_options::Dict{Any,Any}, params::SolverParams, bcs::Dict{Any,Any}, block_nodes::Dict{Int64,Vector{Int64}})` in all five solver modules (argument keeps its name `params`).
  - `Model_reduction.parse_reduction_blocks(model_param::ModelReductionParams)`, `Model_reduction.init_reduce_model(model_param::ModelReductionParams, block_nodes, density)`.
  - `Solver_Manager._typed_solver(solver_params::Dict)::SolverParams` — temporary, removed in Task 3.
  - `Solver_Manager._calculation_options(s::SolverParams)::Dict{String,Any}` — builds the runtime `solver_options["Calculation"]` dict the models read.
  - `solver_options["Model Reduction"]` keeps its values (`false`, or the reduction parameters — now a `ModelReductionParams`), so `setup_reduced_state` is unchanged.

- [ ] **Step 1: Write the failing test**

`test/unit_tests/Support/Parameters/Input/ut_model_reduction_params.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const PS = PeriLab.ParameterSpec
const ID = PeriLab.InputDeck

# Model_reduction.jl is included by Matrix_Verlet.jl
const UT_MR = PeriLab.Solver_Manager.Matrix_Verlet.Model_reduction

function ut_mr(dict)
    ctx = PS.ParseContext()
    p = PS.parse_section(ID.ModelReductionParams, dict, "MR", ctx)
    @test isempty(ctx.errors)
    return p
end

@testset "reduction blocks forms" begin
    @test UT_MR.parse_reduction_blocks(ut_mr(Dict{String,Any}("Type" => "Guyan",
                                                              "Reduction Blocks" => 2))) == [2]
    @test UT_MR.parse_reduction_blocks(ut_mr(Dict{String,Any}("Type" => "Guyan",
                                                              "Reduction Blocks" => "2, 3"))) ==
          [2, 3]
    @test UT_MR.parse_reduction_blocks(ut_mr(Dict{String,Any}("Type" => "Guyan"))) === nothing
end

@testset "Newmark model reduction is rejected clearly" begin
    ctx = PS.ParseContext()
    PS.parse_section(ID.SolverParams,
                     Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                      "Newmark" => Dict{String,Any}("Model Reduction" => true)),
                     "Solver", ctx)
    @test ctx.errors[1].path == "Solver.Newmark"
    @test ctx.errors[1].message ==
          "\"Model Reduction\" is not supported by the Newmark solver; use \"Verlet Matrix Based\" with a \"Model Reduction\" section"
    ctx = PS.ParseContext()
    PS.parse_section(ID.SolverParams,
                     Dict{String,Any}("Initial Time" => 0.0, "Final Time" => 1.0,
                                      "Newmark" => Dict{String,Any}("Model Reduction" => false)),
                     "Solver", ctx)
    @test isempty(ctx.errors)
end
```

Update `input_tests.jl`:

```julia
for file in ["ut_mesh_blocks.jl", "ut_solver.jl", "ut_outputs_conditions_contact.jl",
             "ut_read_input.jl", "ut_golden_decks.jl", "ut_validate_yaml.jl",
             "ut_solver_typed.jl", "ut_model_reduction_params.jl"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl`
Expected: FAIL in `ut_model_reduction_params.jl` — `MethodError: no method matching parse_reduction_blocks(::PeriLab.InputDeck.ModelReductionParams)`, and "Newmark model reduction is rejected clearly" (no error reported).

- [ ] **Step 3: Model reduction takes `ModelReductionParams`; Newmark rejects it in validation**

In `src/Support/Parameters/Input/solver.jl`, at the end of `check!(p::SolverParams, …)` before `return nothing`, add:

```julia
    if p.newmark !== nothing && p.newmark.model_reduction
        add_error!(ctx, join_path(path, "Newmark"),
                   "\"Model Reduction\" is not supported by the Newmark solver; use \"Verlet Matrix Based\" with a \"Model Reduction\" section")
    end
```

In `src/Core/Model_reduction/Model_reduction.jl`, add after `using ...Data_Manager`:

```julia
using ...InputDeck: ModelReductionParams
```

Replace

```julia
function parse_reduction_blocks(model_param::Dict)
    reduction_blocks = get(model_param, "Reduction Blocks", nothing)
```

with

```julia
function parse_reduction_blocks(model_param::ModelReductionParams)
    reduction_blocks = model_param.reduction_blocks
```

and in its docstring replace ``- `model_param::Dict`: The `"Model Reduction"` solver parameters`` with ``- `model_param::ModelReductionParams`: The `"Model Reduction"` solver parameters``. In `init_reduce_model`, replace

```julia
function init_reduce_model(model_param::Dict, block_nodes::Dict{Int64,Vector{Int64}},
```

with

```julia
function init_reduce_model(model_param::ModelReductionParams,
                           block_nodes::Dict{Int64,Vector{Int64}},
```

and in its body replace `model_param["Type"]` (every occurrence) with `model_param.type`, `get(model_param, "Number of Modes", 1)` with `model_param.number_of_modes`, and `get(model_param, "Material Point Region", true)` with `model_param.material_point_region`; update its docstring argument type the same way.

- [ ] **Step 4: Run the test**

Run: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl`
Expected: PASS.

- [ ] **Step 5: Solvers read fields**

In each of `Verlet_solver.jl`, `Static_solver.jl`, `Matrix_linear_static.jl`, `Matrix_Verlet.jl`, `Newmark.jl`:
- add `using ...InputDeck: SolverParams, start_time, end_time` next to the file's `using ...Data_Manager`;
- in `function init_solver(solver_options::Dict{Any,Any},` change the next line `params::Dict,` to `params::SolverParams,` and update the docstring's `params` type to `SolverParams`;
- remove `get_initial_time`, `get_final_time`, `get_safety_factor`, `get_fixed_dt`, `get_nsteps`, `get_numerical_damping`, `get_max_damage` from the file's `using ...Parameter_Handling:` list (delete the whole `using` statement if nothing remains).

Line replacements (current line → new line):

`Verlet_solver.jl`:

| current | new |
|---|---|
| `initial_time = get_initial_time(params)` | `initial_time = start_time(params, Data_Manager.get_current_time())` |
| `final_time = get_final_time(params)` | `final_time = end_time(params, Data_Manager.get_current_time())` |
| `safety_factor = get_safety_factor(params)` | `safety_factor = params.verlet.safety_factor` |
| `fixed_dt = get_fixed_dt(params)` | `fixed_dt = params.verlet.fixed_dt` |
| `numerical_damping = get_numerical_damping(params)` | `numerical_damping = params.verlet.numerical_damping` |
| `max_damage = get_max_damage(params)` | `max_damage = params.maximum_damage` |

`Static_solver.jl`:

| current | new |
|---|---|
| `initial_time = get_initial_time(params)` | `initial_time = start_time(params, Data_Manager.get_current_time())` |
| `final_time = get_final_time(params)` | `final_time = end_time(params, Data_Manager.get_current_time())` |
| `fixed_dt = get_fixed_dt(params)` | `fixed_dt = params.static.fixed_dt` |
| `numerical_damping = get_numerical_damping(params)` | `numerical_damping = params.static.numerical_damping` |
| `max_damage = get_max_damage(params)` | `max_damage = params.maximum_damage` |

and in `Static_solver.jl` replace

```julia
    if fixed_dt == -1.0
        if haskey(params, "Number of Steps")
            nsteps = params["Number of Steps"]
            dt = (final_time - initial_time) / nsteps
        else
            nsteps = Int64(1)
            dt = final_time - initial_time
        end
    else
        if haskey(params, "Number of Steps")
```

with

```julia
    if fixed_dt == -1.0
        if params.number_of_steps !== nothing
            nsteps = params.number_of_steps
            dt = (final_time - initial_time) / nsteps
        else
            nsteps = Int64(1)
            dt = final_time - initial_time
        end
    else
        if params.number_of_steps !== nothing
```

and replace the block from `solver_specifics = Dict("Solution tolerance" => 1e-7,` through the `end` of `if haskey(params["Static"], "Linear Start Value")` with

```julia
    static = params.static
    solver_specifics = Dict("Solution tolerance" => static.solution_tolerance,
                            "Residual tolerance" => static.residual_tolerance,
                            "Maximum number of iterations" => static.maximum_number_of_iterations,
                            "Show solver iteration" => static.show_solver_iteration,
                            "Residual scaling" => static.residual_scaling,
                            "m" => static.m,
                            "Linear Start Value" => static.linear_start_value === nothing ?
                                                    zeros(2 * dof) :
                                                    parse.(Float64,
                                                           split(static.linear_start_value)))
```

`Matrix_linear_static.jl`:

| current | new |
|---|---|
| `solver_options["Initial Time"] = get_initial_time(params)` | `solver_options["Initial Time"] = start_time(params, Data_Manager.get_current_time())` |
| `solver_options["Final Time"] = get_final_time(params)` | `solver_options["Final Time"] = end_time(params, Data_Manager.get_current_time())` |
| `solver_options["Number of Steps"] = get_nsteps(params)` | `solver_options["Number of Steps"] = something(params.number_of_steps, 1)` |

and replace the statement starting `solver_options["Matrix Update"] = get(params["Linear Static Matrix Based"],` (through its closing parenthesis) with

```julia
    solver_options["Matrix Update"] = params.linear_static_matrix_based.matrix_update
```

`Matrix_Verlet.jl`:

| current | new |
|---|---|
| `initial_time = get_initial_time(params)` | `initial_time = start_time(params, Data_Manager.get_current_time())` |
| `final_time = get_final_time(params)` | `final_time = end_time(params, Data_Manager.get_current_time())` |
| `nsteps = get_nsteps(params)` | `nsteps = something(params.number_of_steps, 1)` |
| `safety_factor = get_safety_factor(params)` | `safety_factor = params.verlet_matrix_based.safety_factor` |
| `fixed_dt = get_fixed_dt(params)` | `fixed_dt = params.verlet_matrix_based.fixed_dt` |
| `solver_options["Numerical Damping"] = get_numerical_damping(params)` | `solver_options["Numerical Damping"] = params.verlet_matrix_based.numerical_damping` |

and replace

```julia
    model_reduction = get(params["Verlet Matrix Based"], "Model Reduction", false)
    reduce = true
    if model_reduction == false
        reduce = false
    end
    solver_options["Model Reduction"] = model_reduction
```

with

```julia
    model_reduction = something(params.verlet_matrix_based.model_reduction, false)
    reduce = model_reduction !== false
    solver_options["Model Reduction"] = model_reduction
```

`Newmark.jl`:

| current | new |
|---|---|
| `solver_options["Initial Time"] = get_initial_time(params)` | `solver_options["Initial Time"] = start_time(params, Data_Manager.get_current_time())` |
| `solver_options["Final Time"] = get_final_time(params)` | `solver_options["Final Time"] = end_time(params, Data_Manager.get_current_time())` |
| `solver_options["Number of Steps"] = get_nsteps(params)` | `solver_options["Number of Steps"] = something(params.number_of_steps, 1)` |

and replace

```julia
    solver_options["Newmark Beta"] = get(params["Newmark"], "Newmark Delta", 0.5) # name is delta because of beta. Internally it is beta.
    solver_options["Newmark Alpha"] = get(params["Newmark"], "Newmark Alpha",
                                          0.25 * (0.5 + solver_options["Newmark Beta"])^2)
    solver_options["Matrix Update"] = get(params["Newmark"], "Matrix Update", false)
```

with

```julia
    solver_options["Newmark Beta"] = params.newmark.newmark_delta # name is delta because of beta. Internally it is beta.
    solver_options["Newmark Alpha"] = something(params.newmark.newmark_alpha,
                                                0.25 * (0.5 + solver_options["Newmark Beta"])^2)
    solver_options["Matrix Update"] = params.newmark.matrix_update
```

and delete

```julia
    model_reduction = get(params["Newmark"], "Model Reduction", false)
    if model_reduction
        init_reduce_model(solver_options, block_nodes, density)
        return
    end
```

(input validation now rejects `Model Reduction: true` for Newmark, so the branch is unreachable).

After the edits, `grep -n 'params\["' src/Core/Solver/{Verlet_solver,Static_solver,Matrix_linear_static,Matrix_Verlet,Newmark}.jl` must print nothing.

`Solver_manager.jl` — add the imports

```julia
using ..InputDeck: SolverParams, solver_name, model_options
using ..ParameterSpec: ParseContext, parse_section, report!
```

remove `get_solver_name`, `get_model_options`, `get_calculation_options` from its `using ..Parameter_Handling:` list, and add before `function init(`:

```julia
# Temporary until the typed input is passed in (phase 2b-1, Task 3).
function _typed_solver(solver_params::Dict)
    ctx = ParseContext()
    solver = parse_section(SolverParams, solver_params, "Solver", ctx)
    report!(ctx)
    return solver
end

# Runtime format read by the models
function _calculation_options(s::SolverParams)
    return Dict{String,Any}("Calculate Cauchy" => s.calculate_cauchy,
                            "Calculate von Mises stress" => s.calculate_von_mises_stress,
                            "Calculate Strain" => s.calculate_strain)
end
```

In `init`, make these replacements:

| current | new |
|---|---|
| `solver_params = step_id == -1 ? params["Solver"] :` and the next line `get_solver_params(params, step_id)` | `solver_params = _typed_solver(step_id == -1 ? params["Solver"] : get_solver_params(params, step_id))` |
| `solver_options["Models"] = get_model_options(solver_params)` | `solver_options["Models"] = model_options(solver_params)` |
| `solver_options["All Models"] = get_model_options(solver_params)` | `solver_options["All Models"] = model_options(solver_params)` |
| `solver_options["Calculation"] = get_calculation_options(solver_params)` | `solver_options["Calculation"] = _calculation_options(solver_params)` |
| `step_solver_params = get_solver_params(params, step)` | `step_solver_params = _typed_solver(get_solver_params(params, step))` |
| `get_model_options(step_solver_params))` | `model_options(step_solver_params))` |
| `calc_options = get_calculation_options(step_solver_params)` | `calc_options = _calculation_options(step_solver_params)` |

and replace every remaining `get_solver_name(solver_params)` in `init` (five occurrences: `solver_options["Solver"] = …`, the `@info`, `create_module_specifics(…`, the `@abort`, `Data_Manager.set_model_module(…`) with `solver_name(solver_params)`.

- [ ] **Step 6: Run everything**

Run: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl`
Expected: PASS.
Run (background, ≈35 min): `julia --project=. -e 'using Pkg; Pkg.test()' > /tmp/full_2b1_t2.log 2>&1; tail -5 /tmp/full_2b1_t2.log`
Expected: `Testing PeriLab tests passed` — the fullscale tests run Verlet, Static, Linear Static Matrix Based, Verlet Matrix Based with Craig-Bampton reduction, Newmark and multistep decks against reference results.

- [ ] **Step 7: Commit**

```bash
git add src/Core src/Support test/unit_tests/Support/Parameters/Input
git -c user.name="Jan-Timo Hesse" -c user.email="jthesse@gmx.de" commit -m "Solvers and model reduction read typed solver options"
```

---

### Task 3: Pass the typed input from `run` to the solver manager; remove the adapter

**Files:**
- Modify: `src/Support/Parameters/parameter_handling.jl` (`validate_input`, `validate_yaml`)
- Modify: `src/IO/read_inputdeck.jl` (`read_input_deck`)
- Modify: `src/IO/IO.jl` (`initialize_data`)
- Modify: `src/Core/Solver/Solver_manager.jl` (`init`)
- Modify: `src/PeriLab.jl` (`run`)
- Modify: `test/unit_tests/Support/Parameters/Input/ut_validate_yaml.jl`

**Interfaces:**
- Consumes: Task 1 `solver_steps`, `solver_step`; Task 2 solvers; phase 2a `read_input`.
- Produces:
  - `Parameter_Handling.validate_input(params::Dict; directory = "", no_strict = false)::Tuple{Dict,PeriLabInput}`; `validate_yaml` returns `first(validate_input(...))` (unchanged API).
  - `IO.read_input_deck(filename::String; directory = dirname(filename), no_strict = false)::Tuple{Dict,PeriLabInput}`.
  - `IO.initialize_data(filename, filedirectory, comm; no_strict = false)` returns `(params, input, steps)` with `steps = solver_steps(input)`.
  - `Solver_Manager.init(params::Dict, input::PeriLabInput, step_id::Int64)`.

- [ ] **Step 1: Write the failing test**

Append to `test/unit_tests/Support/Parameters/Input/ut_validate_yaml.jl`:

```julia
@testset "validate_input returns the typed input" begin
    params = ut_valid_params()
    deck, input = PeriLab.Parameter_Handling.validate_input(params)
    @test deck === params["PeriLab"]
    @test input isa PeriLab.InputDeck.PeriLabInput
    @test input.sections.solver.verlet !== nothing
end

@testset "read_input_deck" begin
    dir = mktempdir()
    file = joinpath(dir, "deck.yaml")
    write(file,
          """
          PeriLab:
            Discretization:
              Type: Text File
              Input Mesh File: m.txt
            Blocks:
              block_1:
                Block ID: 1
                Density: 1.0
                Horizon: 1.0
            Models:
              Material Models:
                mat_1:
                  Material Model: a
            Solver:
              Initial Time: 0.0
              Final Time: 1.0
              Number of Steps: 4
              Verlet:
                Safety Factor: 0.5
          """)
    deck, input = PeriLab.IO.read_input_deck(file)
    @test deck["Solver"]["Number of Steps"] == 4
    @test input.sections.solver.number_of_steps === 4
    @test PeriLab.InputDeck.solver_steps(input) == [-1]
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl`
Expected: FAIL with `UndefVarError: validate_input not defined in PeriLab.Parameter_Handling` and `UndefVarError: read_input_deck not defined in PeriLab.IO`.

- [ ] **Step 3: Write the implementation**

In `parameter_handling.jl`, replace the whole `validate_yaml` function (docstring included) with:

```julia
"""
    validate_input(params; directory = "", no_strict = false) -> (deck, input)

Validates a loaded input deck against the typed input declarations
(`InputDeck.read_input`) and, for `Models`, the legacy structure. Reports
every problem at once and aborts if there is an error. Returns the unchanged
`params["PeriLab"]` dict and the typed `PeriLabInput`.
"""
function validate_input(params::Dict; directory::AbstractString = "", no_strict::Bool = false)
    if !haskey(params, "PeriLab") || !(params["PeriLab"] isa AbstractDict) ||
       length(params["PeriLab"]) < 2
        @abort "Yaml file is not valid."
        return
    end
    deck = params["PeriLab"]
    input, ctx = read_input(deck, directory;
                            strict = strict_mode(deck; no_strict_flag = no_strict))
    validate_models(deck) ||
        add_error!(ctx, "Models", "invalid model parameters (see the warnings above)")
    report!(ctx)
    return deck, input
end

"""
    validate_yaml(params; directory = "", no_strict = false)

Like `validate_input`, returning only the deck dict.
"""
function validate_yaml(params::Dict; directory::AbstractString = "", no_strict::Bool = false)
    return first(validate_input(params; directory = directory, no_strict = no_strict))
end
```

and change the module's `export validate_yaml` line to `export validate_yaml, validate_input`.

In `src/IO/read_inputdeck.jl`, replace

```julia
using ...Parameter_Handling: validate_yaml

export read_input_file
```

with

```julia
using ...Parameter_Handling: validate_yaml, validate_input

export read_input_file, read_input_deck
```

and append:

```julia
"""
    read_input_deck(filename; directory = dirname(filename), no_strict = false)

Reads and validates the input deck. Returns the deck dict (for consumers not
yet switched to typed input) and the typed `PeriLabInput`.
"""
function read_input_deck(filename::String; directory::AbstractString = dirname(filename),
                         no_strict::Bool = false)
    if !isfile(filename)
        @abort "$(filename) can not be found. Make sure the file exist and is readable."
        return
    end
    if !occursin("yaml", filename)
        @abort "Not a supported filetype $filename"
        return
    end
    @info "Read input file $filename"
    return validate_input(read_input(filename); directory = directory, no_strict = no_strict)
end
```

In `src/IO/IO.jl`, add `using ..InputDeck: solver_steps` next to its `Parameter_Handling` import, remove `get_solver_steps` from that import list, and in `initialize_data` replace

```julia
    @timeit "init_data" params=init_data(read_input_file(filename;
                                                         directory = filedirectory,
                                                         no_strict = no_strict),
                                         filedirectory, comm)
    steps = get_solver_steps(params)
    return params, steps
```

with

```julia
    deck, input = read_input_deck(filename; directory = filedirectory, no_strict = no_strict)
    @timeit "init_data" params=init_data(deck, filedirectory, comm)
    steps = solver_steps(input)
    return params, input, steps
```

and update its docstring's return list to `params::Dict, input::PeriLabInput, steps::Vector{Int64}`.

In `src/Core/Solver/Solver_manager.jl`, delete the `_typed_solver` function and its comment; replace

```julia
using ..InputDeck: SolverParams, solver_name, model_options
using ..ParameterSpec: ParseContext, parse_section, report!
```

with

```julia
using ..InputDeck: PeriLabInput, SolverParams, solver_name, model_options, solver_step
```

remove `get_solver_params` from the `using ..Parameter_Handling:` list; replace

```julia
function init(params::Dict,
              step_id::Int64)
```

with

```julia
function init(params::Dict,
              input::PeriLabInput,
              step_id::Int64)
```

replace

```julia
    solver_params = _typed_solver(step_id == -1 ? params["Solver"] : get_solver_params(params, step_id))
```

with

```julia
    solver_params = solver_step(input, step_id)
```

and

```julia
            step_solver_params = _typed_solver(get_solver_params(params, step))
```

with

```julia
            step_solver_params = solver_step(input, step)
```

Update the docstring of `init` to list `input::PeriLabInput`.

In `src/PeriLab.jl`, in `run`, replace

```julia
            @timeit "IO.initialize_data" params,
                                         steps=IO.initialize_data(filename,
```

with

```julia
            @timeit "IO.initialize_data" params, input,
                                         steps=IO.initialize_data(filename,
```

and

```julia
                                              solver_options=Solver_Manager.init(params,
                                                                                 step_id)
```

with

```julia
                                              solver_options=Solver_Manager.init(params,
                                                                                 input,
                                                                                 step_id)
```

Finally run `grep -rn "Solver_Manager.init(\|initialize_data(filename" src test` and confirm no caller with the old arity remains.

- [ ] **Step 4: Run tests to verify they pass**

Run: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl`
Expected: PASS.
Run (background, ≈35 min): `julia --project=. -e 'using Pkg; Pkg.test()' > /tmp/full_2b1_t3.log 2>&1; tail -5 /tmp/full_2b1_t3.log`
Expected: `Testing PeriLab tests passed`.

- [ ] **Step 5: Commit**

```bash
git add src test
git -c user.name="Jan-Timo Hesse" -c user.email="jthesse@gmx.de" commit -m "Pass typed input from run to the solver manager"
```
