# Parameters

The input deck is read into typed structs. Every section and every model declares
its YAML keys as an `@params` struct (`src/Support/Parameters/Spec`): type,
required or default, `min` / `max`, allowed values, a quantity (documentation only;
PeriLab has no fixed units) and a description.

`PeriLab.InputDeck.read_input` validates a deck against these declarations and
reports every problem at once (unknown keys with a suggestion, wrong types, values
out of range). With `Strict Validation: false` (or `--no_strict`) unknown keys are
warnings instead of errors.

A model module declares its own keys and registers them under its model name, e.g.

    @params struct MyMaterialParams
        "Stress at which the material starts to yield"
        yield_stress::Dependent = req("Yield Stress"; min = 0, quantity = :stress)
    end
    __init__() = register_material("My Material", MyMaterialParams)

A docstring above a field is its description (the same as
`description = "..."`; giving both is an error).

See the templates under `src/Models/*/…_template` for every category.

A solver is registered the same way. Its options struct is a subtype of
`SolverOptions`, and the input deck names the solver by its key inside
`Solver` (or a `Multistep Solver` step):

    @params struct MySolverParams <: SolverOptions
        safety_factor::Float64 = opt("Safety Factor"; default = 1.0, min = 0)
    end
    __init__() = register_solver("My Solver", MySolverParams)

The solver module defines `init_solver(options, params, models, bcs, block_nodes)`,
which returns a `TimeWindow` and the solver's state, and two functions for the
time steps:

- `prepare_solver(run::SolverRun{MySolverParams}, ctx::SolverContext)` is called
  once before the first step and returns the solver's workspace (fields, buffers,
  start values);
- `time_step!(ws, run, ctx, clock::StepClock)` computes one step from N to NP1 at
  `clock.time` with increment `clock.dt`, and returns `StepFailed(reason)` if the
  step cannot be computed.

The time loop is shared by all solvers (`time_loop`): it writes the results,
applies `Maximum Damage`, cancels, switches `NP1 → N` and advances the time. A
solver with its own time control may still define
`run_solver(run::SolverRun{MySolverParams}, ctx::SolverContext)`. See
`src/Core/Solver/Solver_template/solver_template.jl`.

A failed step is discarded. If the solver section has a `Fallback`, that solver
recomputes the step from the last converged state and runs to the final time of
the solver step; otherwise the run is cancelled:

    Solver:
      Final Time: 1.0e-3
      Static:
        Maximum number of iterations: 50
      Fallback:
        Verlet:
          Safety Factor: 0.9

When the fallback takes over, each output keeps the time between two outputs, so
a fallback with many short steps does not write every step. `Fallback Output
Frequency` and `Fallback Number of Output Steps` in an output set the frequency
for the fallback steps instead, like `Output Frequency` and `Number of Output
Steps` do for the solver steps.

`SolverParams` selects the solver with a `one_of` field:
`options::SolverOptions = one_of(:solver)` holds exactly one sub-section whose
key is a registered solver name. Errors for a missing, a second or a misspelled
solver come from `one_of`; schema, `describe` and the generated pages list the
registered solvers.

The [input reference](../../generated/input_sections.md) lists every section and model with its
keys; it is generated from the declarations (`PeriLab.generate_parameter_docs`). In a
Julia session, `PeriLab.describe("Correspondence Plastic")` prints the keys of one
model and `PeriLab.describe("Correspondence Plastic"; template = true)` a YAML block
to start from. `PeriLab.to_json_schema(PeriLab.InputDeck.PeriLabInput)` returns a
JSON Schema of the whole deck for editors.

!!! note "Good start"
    Please check some of the full scale tests. There are several yaml files with parameter definitions.
