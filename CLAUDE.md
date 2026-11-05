<!--
SPDX-FileCopyrightText: 2026 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>

SPDX-License-Identifier: BSD-3-Clause
-->

# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**PeriLab** is a Julia package (v3.0.0) for peridynamic simulations — a non-local continuum mechanics framework. It supports bond-based, bond-associated, ordinary and non-ordinary state-based formulations with multiple material models, damage models, solvers, and multi-physics capabilities. Runs on CPU via MPI for distributed computing.

**Requirements:** Julia 1.12+

## Commands

### Installation / Environment Setup

```julia
julia
julia> ]
pkg> app dev .
```

This activates the dev environment and installs dependencies (see `Project.toml`). Also run `Pkg.instantiate()` to ensure `Manifest.toml` is up to date.

### Running the Application

```sh
# Run an example YAML input deck
PeriLab examples/DCB/DCBmodel.yaml

# Dry run (estimates runtime/file size, runs only 10 steps)
# Always do a dry run before starting a large simulation
PeriLab --dry_run examples/DCB/DCBmodel.yaml

# Unknown keys in the input deck are errors by default; --no_strict (or `Strict Validation: false` in the deck) only warns
PeriLab --no_strict examples/DCB/DCBmodel.yaml

# Interactive example picker
PeriLab -e

# Parallel MPI run (requires julia MPI.install_mpiexecjl()); -bind-to core pins each rank to a core
mpiexecjl -n 2 -bind-to core julia -e 'using PeriLab' examples/DCB/DCBmodel.yaml

# Hybrid MPI + threads: one core per thread for every rank
mpiexecjl -n 4 -bind-to core:4 julia -t 4 -e 'using PeriLab' examples/DCB/DCBmodel.yaml
```

### Running Tests

```julia
julia> ]
pkg> test PeriLab
```

Or directly:
```sh
julia --project -e 'using Pkg; Pkg.test()'
```

Tests include ~90 unit tests (`test/unit_tests/`) and ~45 fullscale integration tests (`test/fullscale_tests/`). The full suite initializes and finalizes MPI.

To run one test file on its own, mirror the setup of `test/runtests.jl`: `using PeriLab, MPI, Test, JSON3, Logging`, `MPI.Init()`, `Logging.disable_logging(Logging.Warn)` (otherwise `@test_logs` checks capture extra warnings), and `include("test/helper.jl")` (defines e.g. `typed_section`, `run_perilab`). The test-only dependencies (JSON3, TestSetExtensions, ...) are in `[extras]` of `Project.toml`, not in the main environment.

### Formatting

Code style follows SciML style (`.JuliaFormatter.toml`). Pre-commit hooks include JuliaFormatter.jl, YAML checks, trailing-whitespace, and REUSE license compliance:

```sh
pre-commit install
pre-commit run --all-files
```

## Architecture

### Entry Points

- **`src/main.jl`** — CLI argument parsing via ArgParse; copies examples when `-e` is passed.
- **`src/PeriLab.jl`** — Main module (`PeriLab`). Exports `main` function. Contains `run()` which orchestrates a full simulation: logging init → data init → IO initialization → solver loop (per step) → result writing.

### High-Level Flow (one `run()` call)

1. Initialize logging and MPI
2. `Data_Manager.initialize_data()` / `set_silent()` / `set_verbose()`
3. `IO.read_input(filename)` parses the YAML once; `Additive.load_licensed_models(raw_deck)` loads licensed models (only for decks with additive models)
4. `IO.initialize_data(...)` validates the deck into the typed `PeriLabInput` and returns it with the step list
5. For each step:
   - `IO.init_orientations(input.sections.blocks)`
   - `Solver_Manager.init(input, step_id)` → block_nodes, bcs, `SolverRun`
   - `IO.init_write_results()` (first step only) → result_files, outputs
   - `Solver_Manager.solver(run, ctx::SolverContext)` → shared `time_loop` (see Solver_Common), writes results
6. Close/merge/delete result files (Exodus output)

### Source Directory Structure (`src/`)

| Directory | Purpose |
|---|---|
| `Core/` | Simulation core: `Data_Manager` (central state holding everything), `Solver/` (Verlet, Static, Newmark, Matrix-based solvers; `Solver_Common.jl` holds the shared time loop), `Model_reduction/` (Craig-Bampton, Guyan), `BC_manager.jl` (boundary conditions), `Influence_function.jl` |
| `Models/` | Pluggable physics models via Factory pattern. Subdirectories: `Material/` (bond-based, correspondence, ordinary, UMAT/VUMAT), `Damage/` (critical stretch, energy release), `Contact/` (penalty, short-range + search), `Thermal/` (heat transfer, expansion, HETVAL), `Additive/` (additive manufacturing), `Degradation/` (corrosion, thermal decomposition), `Pre_calculation/`, `Surface_correction/` |
| `IO/` | Input/output: YAML reading (`read_inputdeck.jl`), node sets, Exodus export, Abaqus import, Gmsh import, CSV export, G-code, bond filters (disk, rectangular), volume calculation, logging |
| `FEM/` | FEM integration: Lagrange elements, Arlequin FEM/PD coupling, FEM factory/basis/template |
| `Compute/` | Computation helpers: field values, global values |
| `MPI_communication/` | MPI communication utilities |
| `Support/` | Utilities: `Helpers.jl`, `Geometry.jl`, `Parameters/Spec/` (`ParameterSpec`: the `@params` DSL, conversion, validation errors, model registry, JSON Schema/docs generators), `Parameters/Input/` (`InputDeck`: typed sections of the input deck, `PeriLabInput`), `ModuleLoader.jl` (finds model files by their register call, loads licensed modules) |

### Key Patterns

- **Typed input**: The input deck is validated into immutable structs declared with `@params` (`field::T = req("YAML key"; ...)` / `opt("YAML key"; default = ...)`), with optional `check!` (cross-field rules) and `derive` hooks. No raw parameter dicts travel through the code; models read their typed parameters (e.g. `Data_Manager.get_block_models(block)`).
- **Factory pattern + registration**: Every model type (Material, Damage, Contact, Thermal, Additive, Degradation) and every solver has a `*_Factory.jl`/manager that dispatches on the typed parameters. A model file registers itself in `__init__()` (e.g. `register_material("PD Solid Elastic", PDSolidElasticParams)`); the factory includes every file that contains such a call.
- **Template pattern**: Each model category has a `*_template.jl` defining the interface (e.g., `material_template.jl`, `damage_template.jl`).
- **Data_Manager**: Central singleton-like module that holds all simulation state (particles, bonds, fields, blocks, boundary conditions, status variables). Sub-modules in `Data_manager/` split responsibilities (FEM, MPI, contact, fields, solving process, status vars, utilities).
- **Multi-step**: Input decks can define multiple solver steps with different conditions; the main loop iterates over steps.
- **Solver time loop**: Solvers implement `init_solver`, `prepare_solver(run, ctx)` and `time_step!(ws, run, ctx, clock)`; `Solver_Common.time_loop` does steps, damage limit, output and `NP1 → N`. A step returning `StepFailed` is discarded (`Data_Manager.reset_NP1()`) and, if a `Fallback` solver is given, recomputed by it up to the end of the solver step.
- **Update flags**: Damage models must set `Update[iID] = true` for every node whose bonds they change. `compute_models` recomputes only flagged nodes in its update pass, and the shape tensor pre-calculation skips unflagged nodes (`ONLY_FLAGGED_BONDS_CHANGED`).
- **Bond fields**: Bond states are CSR `ArraysOfArrays.VectorOfVectors` (vector/tensor states with `SVector`/`SMatrix` elements), see `docs/src/man/dev/bond_storage.md`.
- **Threading**: Hot node loops use `Polyester.@batch minbatch=MINBATCH` (`Helpers.MINBATCH`). Only thread loops whose writes are node-local; with more than one thread `@batch` passes arrays as `PtrArray`, so kernel signatures take `AbstractArray`, not `Matrix`/`Vector`.

### Test Organization

- **`test/unit_tests/`** — Mirror the `src/` directory structure; small targeted tests for individual functions/modules.
- **`test/fullscale_tests/`** — Each subdirectory has a YAML input deck and a `.jl` test runner; these run full simulations and compare against reference results.
- **`test/helper.jl`** — Shared test utilities.
- **`test/multifunctional_test/`** — Cross-cutting integration tests (e.g., MPI data manager).
- **`test/tryout/`** — Experimental/development scratchpad files (not loaded by `runtests.jl`).

### Documentation

Built via Documenter.jl from `docs/make.jl`. The docs site is at https://perihub.github.io/PeriLab.jl/. The docs project has its own `Project.toml`/`Manifest.toml`.

### Build

The `build/` directory contains build scripts. See `docs/build/` for the doc build pipeline.
