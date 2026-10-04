# Typed Input Parameters — Phase 2b-3 (Outputs and Compute Classes Consume Typed Input) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Output setup (file names, output frequencies, field selection, result mapping) and compute classes read the typed `OutputParams` / `ComputeClassParams` instead of the params `Dict`; the input deck dict is no longer written to at runtime (`"fieldnames"`, rewritten `"Variable"`).

**Architecture:** Same rule as 2b-1/2b-2: read fields directly; functions only for real logic. The output/compute logic (file names with duplicate check, per-step frequencies, field selection, filtering computes to existing fields) becomes typed functions in `Parameter_Handling` next to the old `Dict` getters, which stay as the reference in an equivalence test over every shipped deck (deleted in phase 4). The per-output runtime mapping built in `IO.get_results_mapping` (`"Fields"`, `"flush_file"`, `"Bond Export"`, …) is runtime state read by about 15 places in IO / exodus / CSV export; it stays a `Dict`, now filled from the typed input, and the compute parameters stored in it become `ComputeClassParams`. Turning the runtime mapping into a struct is a separate refactor (no input handling, no hot path) and is out of scope.

**Tech Stack:** Julia 1.12, existing dependencies only.

**Spec:** `docs/superpowers/specs/2026-10-03-typed-input-parameters-design.md` (§3 step 4, phase 2 of §5). Previous slices: 2b-1 (solvers), 2b-2 (mesh, blocks). Remaining after this one: 2b-4 boundary conditions, 2b-5 contact, FEM / coupling, surface correction.

## Global Constraints

- Julia `1.12`; no new packages.
- Branch `feature/typed-input-parameters`; one commit per task, message ending with the session's `Co-Authored-By` line; commit identity `-c user.name="Jan-Timo Hesse" -c user.email="jthesse@gmx.de"`. Never merge.
- Simulation results and written output files must not change: every fullscale test must still pass.
- Output order invariant: file names, result mapping entries and output frequencies are matched by position. All three must iterate the same collection (`values(input.sections.outputs)` / `enumerate(input.sections.outputs)`), never a mix of the typed dict and the deck dict.
- The `Dict` getters in `parameter_handling_output.jl` / `parameter_handling_computes.jl` stay unchanged; no typed methods are added to them; new typed functions get new names.
- No thin wrappers around single fields.
- Read the full-suite result before recording a task as complete or committing it (2b-2 lesson).
- Test files that run standalone under `mpiexec` (e.g. `test/unit_tests/MPI_communication/ut_MPI.jl`) cannot use `test/helper.jl` functions.
- Test commands (from the repository root):
  - Input tests: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl`
  - Single existing test file as `runtests.jl` runs it:
    `cd test && julia --project=.. -e 'using Test, Logging, MPI; import PeriLab; Logging.disable_logging(Logging.Warn); include("helper.jl"); MPI.Init(); @testset "t" begin include("unit_tests/<path>.jl") end'`
  - Full suite (≈35 min, background, output to a file): `julia --project=. -e 'using Pkg; Pkg.test()'`

## Review Focus

1. A per-step frequency string (`Output Frequency: "10 100"` in a multistep run) must pick the entry of the current step — test in Task 1 ("output frequencies").
2. With both `Output Frequency` and `Number of Output Steps` given, the old code warned and used `Number of Output Steps`; that must stay — test in Task 1 ("output frequencies").
3. `Number of Output Steps` must give `ceil(nsteps / value)` and every frequency must be clamped to `nsteps` — test in Task 1 ("output frequencies").
4. File names, result mapping and frequencies must line up output by output — test in Task 1 ("filenames and frequencies follow the same output order").
5. A compute class whose `Variable` is not a field must be skipped with a warning, not crash — test in Task 1 ("active computes").

---

### Task 1: Typed output and compute functions, proven equivalent on every shipped deck

**Files:**
- Modify: `src/Support/Parameters/parameter_handling.jl` (import), `parameter_handling_output.jl`, `parameter_handling_computes.jl`
- Create: `test/unit_tests/Support/Parameters/Input/ut_outputs_typed.jl`
- Modify: `test/unit_tests/Support/Parameters/Input/input_tests.jl`

**Interfaces:**
- Consumes: phase 2a `OutputParams`, `ComputeClassParams`; `check_for_duplicates` (existing); `Dict` getters (reference only).
- Produces (exported from `Parameter_Handling`):
  - `output_filenames(outputs::Dict{String,OutputParams}, output_dir::String)::Vector{String}` — `<filename>.csv` for CSV, `<filename>.e` otherwise, in `values(outputs)` order; aborts on duplicates.
  - `output_frequencies(outputs::Dict{String,OutputParams}, nsteps::Int64, step_id::Int64)::Vector{Int64}` — in `values(outputs)` order.
  - `output_fieldnames(variables::Dict{String,Bool}, field_keys::Vector{String}, compute_names::Vector{String}, output_type::String)::Vector{Vector{String}}` — `[name, "Constant"|"NP1"]` entries.
  - `compute_names(computes::Dict{String,ComputeClassParams})::Vector{String}` — sorted.
  - `active_computes(computes::Dict{String,ComputeClassParams}, field_keys::Vector{String})::Dict{String,ComputeClassParams}` — computes whose `Variable` is a field (or `Variable * "NP1"` is); others warned and skipped.
- `InputDeck` must export `OutputParams` and `ComputeClassParams` (add them to the `export` statement).

- [ ] **Step 1: Write the failing tests**

`test/unit_tests/Support/Parameters/Input/ut_outputs_typed.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using Test
const PS = PeriLab.ParameterSpec
const ID = PeriLab.InputDeck
const PH = PeriLab.Parameter_Handling

function ut_ot(T, raw)
    ctx = PS.ParseContext()
    value = PS.convert_value(T, raw, "X", ctx)
    @test isempty(ctx.errors)
    return value
end

ut_ot_outputs(raw) = ut_ot(Dict{String,ID.OutputParams}, raw)
ut_ot_vars() = Dict{String,Any}("Displacements" => true)

@testset "output frequencies" begin
    outputs = ut_ot_outputs(Dict{String,Any}("freq" => Dict{String,Any}("Output Filename" => "a",
                                                                         "Output Frequency" => 5,
                                                                         "Output Variables" => ut_ot_vars()),
                                             "steps" => Dict{String,Any}("Output Filename" => "b",
                                                                          "Number of Output Steps" => 3,
                                                                          "Output Variables" => ut_ot_vars()),
                                             "both" => Dict{String,Any}("Output Filename" => "c",
                                                                         "Output Frequency" => 1,
                                                                         "Number of Output Steps" => 2,
                                                                         "Output Variables" => ut_ot_vars()),
                                             "per_step" => Dict{String,Any}("Output Filename" => "d",
                                                                             "Output Frequency" => "10 100",
                                                                             "Output Variables" => ut_ot_vars()),
                                             "clamped" => Dict{String,Any}("Output Filename" => "e",
                                                                            "Output Frequency" => 50,
                                                                            "Output Variables" => ut_ot_vars())))
    by_name(nsteps, step) = Dict(zip(keys(outputs),
                                     PH.output_frequencies(outputs, nsteps, step)))
    f = by_name(20, 1)
    @test f["freq"] == 5
    @test f["steps"] == 7          # ceil(20 / 3)
    @test f["both"] == 10          # Number of Output Steps wins: ceil(20 / 2)
    @test f["per_step"] == 10      # first entry for step 1
    @test f["clamped"] == 20       # clamped to nsteps
    @test by_name(200, 2)["per_step"] == 100
end

@testset "filenames and frequencies follow the same output order" begin
    outputs = ut_ot_outputs(Dict{String,Any}("o$i" => Dict{String,Any}("Output Filename" => "file_$i",
                                                                       "Output Frequency" => i,
                                                                       "Output File Type" => isodd(i) ?
                                                                                            "Exodus" :
                                                                                            "CSV",
                                                                       "Output Variables" => ut_ot_vars())
                                             for i in 1:6))
    filenames = PH.output_filenames(outputs, "out")
    frequencies = PH.output_frequencies(outputs, 100, 1)
    for (k, output) in enumerate(values(outputs))
        expected = joinpath("out",
                            output.output_filename *
                            (output.output_file_type == "CSV" ? ".csv" : ".e"))
        @test filenames[k] == expected
        @test frequencies[k] == output.output_frequency
    end
    duplicate = ut_ot_outputs(Dict{String,Any}("a" => Dict{String,Any}("Output Filename" => "same",
                                                                       "Output Frequency" => 1,
                                                                       "Output Variables" => ut_ot_vars()),
                                               "b" => Dict{String,Any}("Output Filename" => "same",
                                                                       "Output Frequency" => 1,
                                                                       "Output Variables" => ut_ot_vars())))
    @test_throws PeriLab.PeriLabError PH.output_filenames(duplicate, "out")
end

@testset "output fieldnames" begin
    variables = Dict("Displacements" => true, "Forces" => true, "Temperature" => false,
                     "Missing" => true, "Reaction" => true)
    field_keys = ["DisplacementsNP1", "Forces"]
    fields = PH.output_fieldnames(variables, field_keys, ["Reaction"], "Exodus")
    @test sort(fields) == sort([["Displacements", "NP1"], ["Forces", "Constant"],
                                ["Reaction", "Constant"]])
    @test PH.output_fieldnames(variables, field_keys, ["Reaction"], "CSV") ==
          [["Reaction", "Constant"]]
    @test PH.output_fieldnames(variables, field_keys, ["Reaction"], "Exodus") ==
          PH.get_output_fieldnames(variables, field_keys, ["Reaction"], "Exodus")
end

@testset "active computes" begin
    computes = ut_ot(Dict{String,ID.ComputeClassParams},
                     Dict{String,Any}("b_max" => Dict{String,Any}("Compute Class" => "Block_Data",
                                                                  "Variable" => "Displacements",
                                                                  "Block" => "block_1"),
                                      "a_sum" => Dict{String,Any}("Compute Class" => "Block_Data",
                                                                  "Variable" => "Forces",
                                                                  "Block" => "block_1"),
                                      "bad" => Dict{String,Any}("Compute Class" => "Block_Data",
                                                                "Variable" => "Nope",
                                                                "Block" => "block_1")))
    @test PH.compute_names(computes) == ["a_sum", "b_max", "bad"]
    active = PH.active_computes(computes, ["DisplacementsNP1", "Forces"])
    @test sort(collect(keys(active))) == ["a_sum", "b_max"]
    @test active["b_max"] === computes["b_max"]
end

const UT_OT_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", "..", ".."))
const UT_OT_SKIP = ("test/fullscale_tests/test_PD_Solid_Elastic/strain_xx_external.yaml",)

function ut_ot_decks()
    decks = String[]
    for root in ("examples", "test"), (dir, _, files) in walkdir(joinpath(UT_OT_ROOT, root))
        append!(decks, [joinpath(dir, f) for f in files if endswith(f, ".yaml")])
    end
    return sort!(decks)
end

function ut_ot_compare_all_decks()
    for file in ut_ot_decks()
        relpath(file, UT_OT_ROOT) in UT_OT_SKIP && continue
        raw = PeriLab.IO.read_input(file)
        (raw isa AbstractDict && haskey(raw, "PeriLab")) || continue
        deck = raw["PeriLab"]
        input, _ = ID.read_input(deck, dirname(file))
        outputs = input.sections.outputs
        computes = input.sections.compute_class_parameters
        @test sort(PH.output_filenames(outputs, "out")) ==
              sort(PH.get_output_filenames(deck, "out"))
        if haskey(deck, "Outputs")
            for nsteps in (1, 7, 1000)
                typed = Dict(zip(keys(outputs), PH.output_frequencies(outputs, nsteps, 1)))
                reference = Dict(zip(string.(keys(deck["Outputs"])),
                                     PH.get_output_frequencies(deck, nsteps, 1)))
                @test typed == reference
            end
        end
        @test PH.compute_names(computes) == PH.get_computes_names(deck)
        variables = String[c.variable for c in values(computes)]
        @test sort(collect(keys(PH.active_computes(computes, variables)))) ==
              sort(collect(keys(PH.get_computes(deck, variables))))
    end
end

@testset "output and compute functions equal the Dict getters on every shipped deck" begin
    ut_ot_compare_all_decks()
end
```

Update `input_tests.jl`: add `"ut_outputs_typed.jl"` after `"ut_node_sets_typed.jl"` in the list.

- [ ] **Step 2: Run tests to verify they fail**

Run: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl`
Expected: FAIL in `ut_outputs_typed.jl` with `UndefVarError: OutputParams not defined in PeriLab.InputDeck` (or `output_frequencies not defined in PeriLab.Parameter_Handling` once exported).

- [ ] **Step 3: Write the implementation**

In `src/Support/Parameters/Input/InputDeck.jl`, add `OutputParams, ComputeClassParams` to the `export` statement (append after `gcode_block_ids`).

In `src/Support/Parameters/parameter_handling.jl`, extend `using ..InputDeck: read_input, DiscretizationParams` with `, OutputParams, ComputeClassParams`.

In `parameter_handling_output.jl`, add after the existing `export` lines:

```julia
export output_filenames, output_frequencies, output_fieldnames
```

and append:

```julia
# Typed functions (phase 2b). The Dict getters above stay until phase 4.

"""
    output_filenames(outputs, output_dir)

Result file names in `values(outputs)` order (`.csv` for CSV outputs, `.e`
otherwise). Aborts if a name is used twice.
"""
function output_filenames(outputs::Dict{String,OutputParams}, output_dir::String)
    filenames = String[]
    for output in values(outputs)
        extension = output.output_file_type == "CSV" ? ".csv" : ".e"
        push!(filenames, joinpath(output_dir, output.output_filename * extension))
    end
    check_for_duplicates(filenames)
    return filenames
end

"""
    output_frequencies(outputs, nsteps, step_id)

Output frequency (write every n-th step) per output, in `values(outputs)`
order. `Number of Output Steps` wins over `Output Frequency`; a value given as
a string holds one entry per solver step; the result is clamped to `nsteps`.
"""
function output_frequencies(outputs::Dict{String,OutputParams}, nsteps::Int64,
                            step_id::Int64)
    frequencies = zeros(Int64, length(outputs))
    for (id, output) in enumerate(values(outputs))
        use_frequency = output.number_of_output_steps === nothing
        if !use_frequency && output.output_frequency !== nothing
            @warn "Double output step / frequency definition. First option is used. ''Output Frequency'' is ignored."
        end
        value = use_frequency ? output.output_frequency : output.number_of_output_steps
        if value isa String
            value = parse(Int, split(value)[step_id])
        end
        frequencies[id] = use_frequency ? value : Int64(ceil(nsteps / value))
        frequencies[id] = min(frequencies[id], nsteps)
    end
    return frequencies
end

"""
    output_fieldnames(variables, field_keys, compute_names, output_type)

The selected output variables as `[name, "Constant"]` or `[name, "NP1"]`.
CSV outputs only take compute classes.
"""
function output_fieldnames(variables::Dict{String,Bool}, field_keys::Vector{String},
                           compute_names::Vector{String}, output_type::String)
    fieldnames = Vector{Vector{String}}()
    for (name, selected) in variables
        selected || continue
        if output_type == "CSV"
            if name in compute_names
                push!(fieldnames, [name, "Constant"])
            else
                @warn '"' * name * '"' * " is not defined as global variable"
            end
        elseif name in field_keys || name in compute_names
            push!(fieldnames, [name, "Constant"])
        elseif name * "NP1" in field_keys
            push!(fieldnames, [name, "NP1"])
        else
            @warn '"' * name * '"' * " is not defined as variable"
        end
    end
    return fieldnames
end
```

In `parameter_handling_computes.jl`, add after the existing `export` lines:

```julia
export compute_names, active_computes
```

and append:

```julia
# Typed functions (phase 2b). The Dict getters above stay until phase 4.

"Names of the compute classes, sorted."
compute_names(computes::Dict{String,ComputeClassParams}) = sort!(collect(keys(computes)))

"""
    active_computes(computes, field_keys)

The compute classes whose `Variable` is an existing field (directly or as
`<Variable>NP1`); the others are skipped with a warning.
"""
function active_computes(computes::Dict{String,ComputeClassParams},
                         field_keys::Vector{String})
    active = Dict{String,ComputeClassParams}()
    for (name, compute) in computes
        if compute.variable in field_keys || compute.variable * "NP1" in field_keys
            active[name] = compute
        else
            @warn '"' * compute.variable * '"' * " is not defined as variable"
        end
    end
    return active
end
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl`
Expected: PASS. If the equivalence testset fails for a deck, the typed function is wrong unless the `Dict` getter throws for that deck (record a ruling).

- [ ] **Step 5: Commit**

```bash
git add src/Support/Parameters test/unit_tests/Support/Parameters/Input
git -c user.name="Jan-Timo Hesse" -c user.email="jthesse@gmx.de" commit -m "Typed output and compute functions"
```

---

### Task 2: Output setup and compute evaluation read typed input

**Files:**
- Modify: `src/IO/IO.jl` (`get_results_mapping`, `init_write_results`, `set_output_frequency`, `get_global_values`, imports)
- Modify: `src/PeriLab.jl` (`run`)
- Modify: `test/unit_tests/IO/ut_IO.jl`

**Interfaces:**
- Consumes: Task 1 functions; `PeriLabInput`; `typed_input` (tests).
- Produces:
  - `IO.get_results_mapping(input::PeriLabInput, path::String)::Dict{Int64,Dict}` — same runtime structure as before; `"compute_params"` holds a `ComputeClassParams` (or `nothing` for non-compute fields).
  - `IO.init_write_results(input::PeriLabInput, output_dir::String, path::String, PERILAB_VERSION::String, qa_vector::Vector{String}, reuse::Bool)`.
  - `IO.set_output_frequency(input::PeriLabInput, nsteps::Int64, step_id::Int64, reuse::Bool)`.
  - `IO.get_global_values(output::Dict)` reads `ComputeClassParams` fields.

- [ ] **Step 1: Update the existing test (it becomes the failing test)**

In `test/unit_tests/IO/ut_IO.jl`, change the top-level `params = Dict("Outputs" => …, "Compute Class Parameters" => …)` assignment: wrap the whole `Dict(...)` in `typed_input(...)`, and add `"Output Frequency" => 1,` as the first entry of each of the three output dicts (`Output1`, `Output2`, `Output3`) — the typed input requires a frequency. The result:

```julia
params = typed_input(Dict("Outputs" => Dict("Output1" => Dict("Output Frequency" => 1,
                                                              "Output Filename" => filename1,
                                                              "Flush File" => false,
                                                              "Output Variables" => Dict("Forces" => true)),
                                            "Output2" => Dict("Output Frequency" => 1,
                                                              "Output Filename" => filename2,
                                                              "Flush File" => false,
                                                              "Output Variables" => Dict("Displacements" => true,
                                                                                         "Forces" => true)),
                                            "Output3" => Dict("Output Frequency" => 1,
                                                              "Output Filename" => filename3,
                                                              "Output File Type" => "CSV",
                                                              "Output Variables" => Dict("External_Displacement" => true))),
                          "Compute Class Parameters" => Dict("External_Displacement" => Dict("Block" => "block_1",
                                                                                             "Calculation Type" => "Maximum",
                                                                                             "Compute Class" => "Block_Data",
                                                                                             "Variable" => "Displacements"))))
```

The calls `PeriLab.IO.get_results_mapping(params, "")` and `PeriLab.IO.init_write_results(params, "", "", …)` stay. Leave `ut_exodus_export.jl` unchanged: its `compute_params` dicts are never read by the export functions it calls.

- [ ] **Step 2: Run the test to verify it fails**

Run the single-file command for `unit_tests/IO/ut_IO.jl`.
Expected: FAIL with `MethodError: no method matching get_results_mapping(::PeriLab.InputDeck.PeriLabInput, ::String)` (and the same for `init_write_results`).

- [ ] **Step 3: Write the implementation**

`src/IO/IO.jl` imports: delete the `using ..Parameter_Handling: get_flush_file, get_write_after_damage, get_start_time, get_end_time, get_outputs, get_output_frequencies, get_output_filenames, get_computes_names, get_computes` statement (all names replaced) and add

```julia
using ..Parameter_Handling: output_filenames, output_frequencies, output_fieldnames,
                            active_computes, compute_names
```

and extend `using ..InputDeck: solver_steps, BlockParams` with `, PeriLabInput`.

`get_results_mapping` — replace the head of the function from its signature through `output_mapping[id]["end_time"] = get_end_time(outputs, output)` and the following bond-export lines (through the `end` of `if haskey(outputs[output], "Bond Blocks")`) with:

```julia
function get_results_mapping(input::PeriLabInput, path::String)
    field_keys = Data_Manager.get_all_field_keys()
    all_compute_names = compute_names(input.sections.compute_class_parameters)
    computes = active_computes(input.sections.compute_class_parameters, field_keys)
    output_mapping = Dict{Int64,Dict{}}()
    nsets = Data_Manager.get_nsets()

    for (id, (output_name, output)) in enumerate(input.sections.outputs)
        output_mapping[id] = Dict{}()
        output_mapping[id]["Fields"] = Dict{}()

        isempty(output.output_variables) &&
            @warn "No output variables are defined for " * output_name * "."
        fieldnames = output_fieldnames(output.output_variables, field_keys, all_compute_names,
                                       output.output_file_type)
        output_mapping[id]["flush_file"] = output.flush_file
        output_mapping[id]["write_after_damage"] = output.write_after_damage
        output_mapping[id]["start_time"] = output.start_time
        output_mapping[id]["end_time"] = output.end_time

        # Bond export settings have to be carried over, otherwise
        # init_bond_information_export never sees them and silently exports every block.
        output_mapping[id]["Bond Export"] = output.bond_export
        if output.bond_blocks !== nothing
            output_mapping[id]["Bond Blocks"] = output.bond_blocks
        end
```

In the field loop of the same function, replace `compute_params = Dict{}` with `compute_params = nothing`, and replace the compute lookup block

```julia
            for key in keys(computes)
                if fieldname[1] == key
                    fieldname[1] = computes[key]["Variable"]
```

through the end of its `elseif computes[key]["Compute Class"] == "Nearest_Point_Data"` branch with the same logic on the typed compute:

```julia
            for (key, compute) in computes
                if fieldname[1] == key
                    fieldname[1] = compute.variable
                    fieldname[2] = "Constant"
                    if Data_Manager.has_key(fieldname[1] * "NP1")
                        fieldname[2] = "NP1"
                    end
                    compute_name = string(key)
                    compute_params = compute
                    global_var = true
                    if compute.compute_class == "Node_Set_Data"
                        nodeset = compute.node_set
                        num_nodes = length(nsets[nodeset])
                        if num_nodes > 1 && compute.calculation_type === nothing
                            if num_nodes > 10
                                @warn "Compute $key references $num_nodes nodes and will create output entries for each of them in $output_name, make sure this is intendend!"
                            end
                            multi_ids = true
                            node_ids = nsets[nodeset]
                        end
                    elseif compute.compute_class == "Nearest_Point_Data"
                        #find neares_point_id and reduce over cores

                        coor = Data_Manager.get_field("Coordinates")
                        point = [compute.x, compute.y, compute.z]
```

(the rest of the `Nearest_Point_Data` branch — `tree = KDTree(...)` onwards — is unchanged). Update the docstring of `get_results_mapping` (`input::PeriLabInput`: "The typed input deck").

`get_global_values` — replace

```julia
        compute_class = output[varname]["compute_params"]["Compute Class"]
        calculation_type = get(output[varname]["compute_params"], "Calculation Type",
                               "Single_Point")
        fieldname = output[varname]["compute_params"]["Variable"]
        extra_equation = get(output[varname]["compute_params"], "Equation", nothing)
```

with

```julia
        compute = output[varname]["compute_params"]
        compute_class = compute.compute_class
        calculation_type = something(compute.calculation_type, "Single_Point")
        fieldname = compute.variable
        extra_equation = compute.equation
```

and

```julia
            if !haskey(output[varname]["compute_params"], "Block")
                @abort "Missing Block for compute class $varname"
            end
            block = output[varname]["compute_params"]["Block"]
```

with

```julia
            if compute.block === nothing
                @abort "Missing Block for compute class $varname"
            end
            block = compute.block
```

`init_write_results` — change the first argument `params::Dict,` to `input::PeriLabInput,`, `filenames = get_output_filenames(params, output_dir)` to `filenames = output_filenames(input.sections.outputs, output_dir)`, and `outputs = get_results_mapping(params, path)` to `outputs = get_results_mapping(input, path)`; update the docstring.

`set_output_frequency` — change `function set_output_frequency(params::Dict,` to `function set_output_frequency(input::PeriLabInput,` and `output_frequencies = get_output_frequencies(params, nsteps, step_id)` to `output_frequencies = output_frequencies(input.sections.outputs, nsteps, step_id)` — rename the local to avoid shadowing the function: write `frequencies = output_frequencies(input.sections.outputs, nsteps, step_id)` and replace the two later uses of the local `output_frequencies` (`eachindex(output_frequencies)`, `output_frequencies[id]`, in both branches) with `frequencies`. Update the docstring.

After the edits `grep -n 'params' src/IO/IO.jl` must show no reads of the input deck in these four functions (other functions are other slices).

`src/PeriLab.jl`, in `run`: in `IO.init_write_results(params,` replace `params` with `input`, and in `IO.set_output_frequency(params,` replace `params` with `input`.

- [ ] **Step 4: Run tests to verify they pass**

Run the single-file command for `unit_tests/IO/ut_IO.jl` — Expected: `ut_get_results_mapping` passes; `ut_init_write_result_and_write_results` errors only with the standalone-only `Data_Manager` error it already had before this slice ("Field ''Neighborhoodlist'' does not exist"), which the full suite does not show.
Run: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl` — Expected: PASS.
Run (background, ≈35 min): `julia --project=. -e 'using Pkg; Pkg.test()' > /tmp/full_2b3_t2.log 2>&1; tail -5 /tmp/full_2b3_t2.log` — read the result; Expected: `Testing PeriLab tests passed` (the fullscale tests compare written Exodus/CSV results, including compute classes, against references).

- [ ] **Step 5: Commit**

```bash
git add src test
git -c user.name="Jan-Timo Hesse" -c user.email="jthesse@gmx.de" commit -m "Output setup and compute evaluation read typed input"
```
