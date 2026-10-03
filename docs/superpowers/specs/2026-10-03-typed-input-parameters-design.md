# Typed, Validated Input Parameters — Design

Date: 2026-10-03
Status: Draft for review

## 1. Goal

Replace PeriLab's untyped `Dict{String,Any}` input handling with typed,
immutable parameter structs that are declared once per module and validated
against the YAML input deck before a run starts.

Success criteria:

1. Every input value has one declaration: Julia type, YAML key, required or
   optional, default, and optional `min` / `max` / `allowed` / `quantity` /
   `description`.
2. The YAML deck is validated against these declarations before any
   computation; all errors are reported together, with paths and near-match
   suggestions.
3. Modules receive concrete structs. Hot loops (`compute_stresses`, damage,
   thermal, ...) are type-stable; no string-keyed lookups in model code.
4. New modules (local, or licensed via `ModuleLoader`) declare their own
   parameters; no central schema has to be edited.
5. Every existing input deck in `test/` and `examples/` runs unchanged.

### Decisions taken during brainstorming

| Topic | Decision |
|---|---|
| Foundation | Struct-first: the struct declaration is the single source of truth. `src/Support/Parameters/Schema.jl` and `PeriLabInputSchema.jl` are an abandoned tryout and are deleted. |
| Scope | The whole input tree (Discretization, Blocks, Models, Solver, Outputs, Boundary Conditions, Compute Classes, Contact, FEM), implemented in phases. |
| Backward compatibility | Full. Same YAML keys (including spaces/apostrophes), same flat layout, same `"A + B"` model composition. |
| Generators | JSON Schema export, documentation generation and `describe()` / template output are in scope, in the last phase. |
| Unknown keys | Error with near-match suggestions. Opt-out downgrades them to warnings. `Globals` passes through unvalidated. |
| Audience | Domain scientists writing modules by copying a template. Declarations must read like a plain annotated struct and produce friendly errors. Licensed modules are maintained by the same team; no frozen versioned API is required. |
| Mechanism | In-house thin `@params` macro (no external dependency). |
| Units | PeriLab has no fixed unit system; values only need to be consistent. Declarations carry an optional, documentation-only `quantity` (e.g. `:stress`), never a unit, and nothing is converted. |
| Dependent values | A value may depend on exactly one node field (any field available in `Data_Manager`, named by the data file header). Multi-field tables and time/global dependence are out of scope. |
| Licensing | All basic modules are always registered. With a valid license, all licensed modules it grants are registered at startup. Without one, PeriLab runs with basic modules only. |

## 2. Architecture

Five units, each with a single responsibility.

### 2.1 `ParameterSpec` — declaring parameters

Location: `src/Support/Parameters/Spec/`.

Provides `@params`, `req`, `opt`, `parameter_spec(::Type)`, and the
`Dependent` types. It knows nothing about YAML files, modules or the registry.

`@params` rewrites a struct with annotated fields into:

- a plain struct (parametric only over `Dependent` fields, see 2.6),
- `parameter_spec(::Type{T})` returning an ordered collection of `FieldSpec`
  (field name, Julia type, YAML alias, required flag, default, `min`, `max`,
  `allowed`, `quantity`, `description`),
- `from_dict(::Type{T}, dict, path, ctx)` building the struct from a validated
  dict.

The macro output is equivalent to a hand-written plain struct plus a
hand-written `parameter_spec` method; authors can write that form directly if
they ever need to. `@macroexpand` shows exactly that.

Supported field types: `Float64`, `Int64`, `Bool`, `String`, any `@enum`
(YAML string matched against the enum instance names), `Vector{Float64}`,
`Vector{Int64}`, `Vector{String}`, `Dependent`, a nested `@params` struct, and
`Dict{String,T}` of a nested `@params` struct (for user-named entries such as
blocks, outputs, bond filters).

Definition-time checks (errors raised from the macro, pointing at the field):

- field without `req(...)` / `opt(...)`,
- unsupported field type,
- default value not convertible to the field type,
- default value outside its own `min` / `max` / `allowed`,
- duplicate YAML alias within one struct,
- `opt` without `default` (use `opt(...; default = nothing)` with a
  `Union{Nothing,T}` field only where absence is meaningful).

### 2.2 `ParameterRegistry` — available modules

Maps `(category, model name) → parameter struct type`. Categories: Material,
Damage, Thermal, Additive, Degradation, Pre Calculation, Surface Correction,
Contact, FEM element / coupling.

```julia
register_material("Bondbased Elastic", BondbasedElastic)
```

Fixed sections (Discretization, Solver, Outputs, ...) are ordinary `@params`
structs and are not registered.

Registration-time checks:

- duplicate model name within a category,
- composite compatibility: when two structs in the same category declare the
  same YAML alias with different types, this is reported when a composite
  using both is parsed (see 2.4).

Optional: a licensed-module manifest may register model names as
`unavailable` stubs, so a deck referencing them without a license gets a
specific error message (see 2.5).

### 2.3 `InputReader` — YAML dict to `PeriLabInput`

Replaces `validate_yaml` and most of `parameter_handling*.jl`.

`InputReader.parse(dict::Dict; strict::Bool)::PeriLabInput`

- Walks the tree top-down, building structs via `from_dict`.
- For model sections, splits the model name on `+`, looks up each part in the
  registry, builds each part's struct from the same flat YAML block.
- Collects all errors in a `Vector{InputError}` (path, message) instead of
  aborting at the first one; after the walk, reports all and calls `@abort`
  if any are errors.
- Unknown keys: a key is unknown if no struct consuming that YAML object
  declares it. In strict mode an error, otherwise a warning. Suggestions use
  edit distance against the declared aliases of that object (case- and
  punctuation-insensitive).
- `Globals` is copied unvalidated into `PeriLabInput.globals::Dict{String,Any}`.
- Strict mode: top-level YAML key `Strict Validation: false` or command line
  flag `--no-strict`. Default is strict.

Error format (path segments that are not plain identifiers are quoted):

```
Input errors (3):
  Models."Material Models".Steel."Poissons Ratio": unknown key — did you mean "Poisson's Ratio"?
  Models."Material Models".Steel."Young's Modulus": missing (required by Bondbased Elastic)
  Blocks.block_1.Horizon: -0.1 is below minimum 0
```

### 2.4 Model composition

`"Correspondence Elastic + Correspondence Plastic"` yields
`Composite{Tuple{CorrespondenceElastic, CorrespondencePlastic}}`.

- Each part reads the keys it declares from the shared flat YAML block. A part
  needing a value another part also uses (e.g. `Young's Modulus`) declares it
  itself; both read the same key.
- A key in the block is unknown only if no part declares it.
- If two parts declare the same alias with different Julia types, parsing
  fails with an error naming both modules.
- Model functions receive only their own part's struct. The composite
  dispatch (iterating over the tuple) lives in the factory / correspondence
  base, replacing the current `split(material_parameter["Material Model"], "+")`
  calls in `Material_Factory.jl`, `Thermal_Factory.jl`, `Correspondence.jl`
  and `Bond_Associated_Correspondence.jl`.
- Blocks without a model of a category get `NoModel()`.

### 2.5 License-gated module availability

At startup, before the input file is read:

1. All local (basic) modules register (when the factories include them, as
   today).
2. If a license (or `ModuleLoader` local dev mode) is configured,
   `check_license_on_startup` runs; on success all licensed modules the
   license grants are loaded and register.
3. Without a license, or with an invalid one: info message, PeriLab continues
   with basic modules only.
4. A deck referencing an unregistered model fails during parsing:
   - with a near-match suggestion if the name resembles a registered model,
   - with "requires a license that is not available" if the name is a
     registered `unavailable` stub,
   - otherwise "not found; it may require a licensed module".

`ModuleLoader`'s public API is unchanged.

### 2.6 Dependent values

A parameter declared as `Dependent` accepts in YAML either a number or a path
to a `.txt` data file (current file format: first column is the field the
value depends on, further columns named after parameters with spaces replaced
by underscores). The string exists only while reading; the struct never holds
a path.

```julia
@params struct MyMaterial
    youngs_modulus::Dependent = req("Young's Modulus"; min = 0, quantity = :stress)
    poissons_ratio::Float64   = req("Poisson's Ratio")
end

E = p.youngs_modulus
for iID in nodes
    E_i = value(E, iID)
end
```

- `abstract type Dependent end` with concrete `Constant <: Dependent`
  (`value::Float64`) and `Table1D{S} <: Dependent` (field name, spline,
  min/max of the input range, bound field array, warn-once flag).
- `@params` makes the struct parametric over each `Dependent` field
  (`MyMaterial{E<:Dependent}`), so every block's struct is concrete:
  `MyMaterial{Constant}` or `MyMaterial{Table1D{...}}`. Only combinations that
  occur in a deck are compiled.
- `value(::Constant, iID)` returns the constant; `value(::Table1D, iID)`
  evaluates the spline at the bound field's value for `iID`, warning once if
  outside the data range (current behaviour).
- Validation: file exists (relative to the input deck directory), header
  contains a column for this parameter, data is numeric, input column is
  sorted; `min` / `max` apply to every table value. Errors name the file and
  line.
- Binding: after `Data_Manager` fields exist, `bind!` resolves each table's
  field name to the `NP1` field array. A missing field aborts with parameter
  path, file and field name. Evaluating an unbound table throws.
  A table holds a reference to the array it was bound to, and
  `switch_NP1_to_N` swaps which array is `NP1` every step, so the binding must
  be renewed after each switch (phase 3 calls `bind_dependents!` from the
  switch, or binds to a stable holder instead of the raw array).
- A `.txt` path given for a non-`Dependent` field is a type error: "expected
  number, got file path — this parameter does not support dependent values".
- Values derived from other parameters (e.g. shear modulus from Young's
  modulus and Poisson's ratio) are produced in `derive` (2.7). If an input is
  a `Table1D`, the derived value is a `Table1D` computed point by point on the
  same input grid.

Replaces `find_data_files`, `csv_reader_temporary` usage in
`get_model_parameter`, `is_dependent`, `get_dependent_value`,
`ConstantValue` / `InterpolatedValue` in `Helpers.jl`.

### 2.7 Module author interface

What a template author writes:

```julia
module My_Material

using .......ParameterSpec: @params, req, opt
import .......ParameterSpec: derive   # required to extend the hook
using .......ParameterSpec: register_material

@enum Hardening Linear Exponential

@params struct MyMaterial
    youngs_modulus::Float64 = req("Young's Modulus"; min = 0, quantity = :stress,
                                  description = "Elastic stiffness")
    poissons_ratio::Float64 = req("Poisson's Ratio"; min = -1, max = 0.5)
    yield_stress::Float64   = opt("Yield Stress"; default = Inf, min = 0,
                                  quantity = :stress)
    hardening::Hardening    = opt("Hardening"; default = Linear)
end

derive(p::MyMaterial) = p          # optional

# registration runs at load time, never during precompilation
__init__() = register_material("My Material", MyMaterial)

init_model(nodes, p::MyMaterial, block) = ...
compute_stresses(iID, dof, p::MyMaterial, time, dt,
                 strain_inc, stress_N, stress_NP1) = ...
end
```

- Existing model function names stay; `material_parameter::Dict` becomes the
  module's struct type.
- `derive(p)` runs once after parsing and returns a struct of the same
  struct type (type parameters may differ, e.g. a derived field becoming a
  `Table1D`). Derived quantities are declared as optional fields (e.g.
  `bulk_modulus`, `shear_modulus`) that `derive` fills in. It replaces
  init-time dict mutation: `get_all_elastic_moduli` (becomes a shared helper
  returning the completed moduli), UMAT/VUMAT/HETVAL file path normalisation,
  code-level defaults such as `Penalty_model.jl`'s contact stiffness (which
  become declared defaults instead), `Correspondence.jl`'s injected
  `"Bond Associated"`.
- String-valued switches compared in compute code (`thermal_flow.jl`
  `"Type"`, `Penalty_model.jl` `"Symmetry"`, `Material_Basis` `"Symmetry"`)
  become `@enum` fields or dispatch types. A YAML string matches an enum
  instance ignoring case, spaces and punctuation ("plane stress" →
  `PlaneStress`); spellings that are not valid Julia names ("3D") are added via
  `ParameterSpec.enum_aliases(::Type{Symmetry}) = Dict("3D" => Full3D)`.
- Per-node state stays in `Data_Manager` fields; structs are immutable.
- All templates (`Material_template`, `FEM_template`, and the templates of the
  other model categories) are rewritten to this form.

## 3. Lifecycle and integration

1. **Register** (startup): basic modules; licensed modules if licensed (2.5).
2. **Load**: `read_input_file` → YAML `Dict` (unchanged) →
   `InputReader.parse` → `PeriLabInput`. Each MPI rank parses on its own, as
   each rank reads the YAML today (`IO.initialize_data`); nothing is
   broadcast.
3. **Derive**: `derive` on every model struct.
4. **Bind**: after `Data_Manager` fields are created and `init_model` has run,
   `bind!` all `Dependent` tables.
5. **Use**: `Data_Manager` stores a `BlockModels` per block, replacing
   `data["properties"][block][model_name]::Dict{String,Any}`:

   ```julia
   struct BlockModels{M,D,T,A,G,P,S}
       material::M
       damage::D
       thermal::T
       additive::A
       degradation::G
       pre_calculation::P
       surface_correction::S
   end
   ```

   `Model_Factory` retrieves a block's `BlockModels` once and calls model
   functions with the concrete parts — the function barrier that makes hot
   loops type-stable. `get_properties` / `get_property` / `check_property` /
   `set_property(ies)` are removed; "block has no damage model" is dispatch on
   `NoModel`.

Other sections:

- **Multistep solver**: each step is a `@params` struct.
- **Boundary conditions**: `"Step ID"` is parsed into a `Vector{Int64}` once
  (replacing `split(string(...), ",")` per evaluation in `BC_manager.jl`);
  `"Node Set"` `+`-lists are parsed into `Vector{String}` once.
- **Contact**: `Global Search Frequency` moves from `Globals` into the Contact
  struct (as an alias also accepted under `Globals` for compatibility).
- **Globals**: plain `Dict{String,Any}`, unvalidated escape hatch.

## 4. Generators (phase 5)

All are pure functions over `parameter_spec` and the registry, so they
reflect exactly the registered (installed and licensed) modules:

- `to_json_schema(PeriLabInput)` — JSON Schema document for PeriHub and
  editors: types, required, defaults, `minimum` / `maximum`, `enum`,
  `description`; `Dependent` becomes `number | string (data file path)`;
  composites and user-named entries map to `additionalProperties`.
- `generate_parameter_docs(dir)` — Documenter pages under `docs/src/` listing,
  per section and per model: YAML key, type, required, default, range,
  quantity ("in your consistent unit system"), description.
- `describe(name)` — prints a model's or section's parameters;
  `describe(name; template = true)` prints a commented YAML block.

## 5. Phasing

Each phase leaves the code building and the full test suite passing.

1. **Core machinery** — `@params`, `req` / `opt`, `Dependent` (`Constant`,
   `Table1D`), registry, `InputReader` with strict mode and suggestions,
   `NoModel`, `Composite`. Unit tests only; no solver changes.
2. **Fixed sections** — Discretization, Blocks, Solver / Multistep, Outputs,
   Boundary Conditions, Compute Classes, Contact, FEM as `@params` structs;
   their consumers (mesh import, `BC_manager`, IO, solvers) switch to fields.
3. **Models by category** — Material first (elastic, plastic, correspondence,
   UMAT / VUMAT), then Damage, Thermal (incl. HETVAL), Additive, Degradation,
   Pre Calculation, Surface Correction. Each category's factory moves to
   `BlockModels` dispatch; templates and licensed modules of that category are
   updated in the same phase.
4. **Remove the old path** — delete `parameter_handling*.jl`, `validate_yaml`,
   `get_properties` & co., `find_data_files`, old dependent-value helpers,
   `Schema.jl`, `PeriLabInputSchema.jl`, `example_usage.jl`,
   `apply_defaults.jl`.
5. **Generators** — JSON Schema, docs, `describe`.

Bridge during phases 2–3: `InputReader` also produces the validated `Dict` so
unmigrated consumers keep working. Removed in phase 4.

## 6. Testing

- **Unit tests** per feature: aliases, defaults, `min` / `max` / `allowed`,
  enums, nested structs, `Dict` of named entries, `Dependent` constant and
  table (incl. range warnings, binding errors), composites incl. alias type
  conflicts, unknown keys with suggestions, strict on/off, `Globals`
  pass-through, and the macro's definition-time error messages.
- **Golden compatibility test**: every YAML under `test/` and `examples/`
  parses in strict mode without errors, or is on an explicit allowlist with a
  reason. Allowlisted decks are fixed or documented before phase 4 ends.
- **Regression**: existing test suite after every phase; numerical results
  unchanged.
- **Type stability**: `@inferred` / JET checks on `compute_stresses`, damage,
  thermal and contact entry points for representative concrete parameter
  types; a benchmark of one material loop before and after.
- **License gating**: basic-only run; local-dev-mode licensed modules; a deck
  referencing a licensed model without a license (asserting the message).

## 7. Out of scope

- Dependent values over more than one field, or over time / global scalars
  (the `Dependent` abstraction allows adding `Table2D` etc. later without
  changing module code).
- Unit conversion or unit checking.
- Changes to the YAML format.
- A frozen, versioned public API for external module authors.
