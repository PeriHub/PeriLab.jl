# Typed Input Parameters — Phase 3c-1: Correspondence Family on Typed Block Materials

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Correspondence materials compute from the typed `BlockMaterial` instead of the material dict. This covers the dispatcher, Elastic, Plastic, UMAT, VUMAT, bond-associated and zero energy control.

**Architecture:**
- **Hooke matrix.**
  - `get_Hooke_matrix`'s formulas move into `_hooke_matrix(c, symmetry, dof, ID)`, which gets its constants from a function `c(key, ID)`.
  - Both the legacy dict path and the new `hooke_matrix(material, dof, ID)` use it, so the formulas exist once.
- **Typed methods are added next to the dict methods.** Every correspondence module and the zero energy control get them. In the last task the Material factory switches correspondence blocks to the typed methods.
- **Dict methods stay as dead code until phase 3c-2.** 3c-2 removes them together with the material dict, `get_all_elastic_moduli` and the legacy tests. So each task stays green and the old tests keep pinning the old behaviour, which the typed path is compared against.
- **Indexed keys (`Property_N`)** are kept by the parser in `WithBase.extras` and passed on in `BlockMaterial.extras`.

**Tech Stack:** Julia 1.12, `PeriLab.ParameterSpec`, `PeriLab.Data_Manager`.

**Spec:** `docs/superpowers/specs/2026-10-03-typed-input-parameters-design.md` (§2.4 composites, §2.7, §3)

## Scope split (decided while planning)

Phase 3c is split.
- **3c-1 (this plan):** the correspondence family.
- **3c-2:**
  - the matrix-based correspondence path and the matrix solvers;
  - `compute_field_values` (strain via the Hooke matrix), `Pre_Calculation.check_dependencies`, the critical time step, `Accuracy Order` in `Model_Factory`, local damping;
  - then removal of the material dict and of every dict method (including those left dead here), `write_moduli!`, `get_all_elastic_moduli`, the dict `get_Hooke_matrix` and the legacy tests.

## Global Constraints

- Numerical results unchanged: the full suite must stay green. That includes the fullscale correspondence decks:
  - `test_Correspondence_Elastic`, `test_Correspondence_Elastic_Plastic`, `test_correspondence_elastic_3D`;
  - `test_correspondence_elastic_with_zero_E_control`, `test_3D_aniso_material`, `test_symmetry`;
  - `test_Umat`, `test_DCB`, `test_Dogbone`.
- Read fields directly; functions only for real logic.
- Structs are immutable; per-node state stays in `Data_Manager` fields.
- Dependent tables are bound before use (`Material.bind_material!`). Correspondence blocks call it before init and before every compute.
- Branch `feature/typed-input-parameters`; one commit per task; never merge. Read the full-suite result before committing a task that runs it.

## Review Focus

1. **An orthotropic material whose `Young's Modulus X` comes from a temperature table** (`test_symmetry/symmetry_eng_dep.yaml`) gets a per-node Hooke matrix following the current temperature. Pinned in Task 2 (`Hooke matrix from a table`) and by the symmetry fullscale test.
2. **Every symmetry branch of the Hooke matrix** (isotropic 3D / plane strain / plane stress / unknown symmetry, orthotropic, transverse isotropic 3D / plane strain / plane stress, anisotropic) gives the same matrix typed and legacy. Pinned in Task 2 (`typed Hooke matrix equals legacy`).
3. **An active `Flaw Function` with missing size, magnitude or location** aborts with a clear message instead of a `KeyError`. Pinned in Task 2 (`typed flaw function`).
4. **A UMAT block with `Property_1 … Property_N`** passes the same properties to the UMAT as before, and a missing `Property_k` is 0.0. Pinned in Task 4 (`UMAT properties from extras`).
5. **Zero energy control is skipped for UMAT materials and applied for VUMAT and the others**, as before. Pinned in Task 3 (`zero energy control skips UMAT`).

---

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
- Spec tests: `julia --project=. test/unit_tests/Support/Parameters/Spec/run_spec_tests.jl 2>&1 | tail -15`
- Input tests: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl 2>&1 | tail -15`
- Full suite (~30 min, background, read the tail): `julia --project=. -e 'using Pkg; Pkg.test()' > <workspace>/suite.log 2>&1`
- `ut_MPI.jl` runs standalone under mpiexec without `helper.jl`. Never use `typed_block_material` there.

## Module paths and imports

| Module | Path | ParameterSpec import |
|---|---|---|
| Material_Basis | `PeriLab.Solver_Manager.Material_Basis` | `using ......ParameterSpec: value` |
| Material factory | `PeriLab.Solver_Manager.Model_Factory.Material` | (has it) |
| Correspondence | `….Material.Correspondence` | — |
| Elastic/Plastic/UMAT/VUMAT | `….Material.Correspondence.Correspondence_*` | `.......ParameterSpec` |
| Bond_Associated_Correspondence | `….Material.Correspondence.Bond_Associated_Correspondence` | — |
| Zero_Energy_Control | `PeriLab.Solver_Manager.Zero_Energy_Control` (+ `.Global_Zero_Energy_Control`) | — |

`BlockMaterial` is defined in the Material factory *after* the model modules are included. Modules therefore take `material` untyped and dispatch on their own parameter struct.

---

### Task 1: Parser keeps indexed-key values (`WithBase.extras`)

**Files:**
- Modify: `src/Support/Parameters/Spec/model.jl` (`WithBase`, `_check_key_patterns!`, `_parse_model`)
- Test: `test/unit_tests/Support/Parameters/Spec/ut_model.jl` (append)

**Interfaces:**
- Produces:
  - `struct WithBase{B,M}` with fields `base::B`, `model::M`, `extras::Dict{String,Any}`.
  - A two-argument constructor `WithBase(base, model)` gives empty extras.
  - `_check_key_patterns!(dict, known, types, path, ctx)` returns `Dict{String,Any}` of the converted pattern values (failed conversions are left out).

- [ ] **Step 1: Write the failing test** (append to `ut_model.jl`)

```julia
@testset "indexed key values are kept in extras" begin
    m, ctx = ut_parse_based(Dict{String,Any}(UT_KEY => "UT Patterned", "File" => "a.so",
                                             "Property_1" => 1, "Property_27" => 2.5))
    @test isempty(ctx.errors)
    @test m.extras == Dict{String,Any}("Property_1" => 1.0, "Property_27" => 2.5)
    @test m.extras["Property_1"] isa Float64
    m, ctx = ut_parse_based(Dict{String,Any}(UT_KEY => "UT Empty"))
    @test isempty(m.extras)
    @test PS.WithBase(1, 2).extras == Dict{String,Any}()
end
```

- [ ] **Step 2: Run to verify it fails**

Run: Spec tests
Expected: `type WithBase has no field extras`.

- [ ] **Step 3: Implement** in `src/Support/Parameters/Spec/model.jl`

Replace

```julia
struct WithBase{B,M}
    base::B
    model::M
end
```

with

```julia
struct WithBase{B,M}
    base::B
    model::M
    extras::Dict{String,Any}   # values of indexed keys (see `key_patterns`), e.g. Property_1
end
WithBase(base, model) = WithBase(base, model, Dict{String,Any}())
```

In `_check_key_patterns!`:
- Add `extras = Dict{String,Any}()` before the loop.
- Replace the line `convert_value(last(key_patterns(T)[pattern]), dict[k], join_path(path, key), ctx)` with

```julia
            v = convert_value(last(key_patterns(T)[pattern]), dict[k], join_path(path, key),
                              ctx)
            v === FAILED || (extras[key] = v)
```

- Replace the final `return nothing` with `return extras`.

In `_parse_model`:
- `_check_key_patterns!(dict, known, all_types, path, ctx)` becomes `extras = _check_key_patterns!(dict, known, all_types, path, ctx)`.
- `WithBase(base_part, model)` becomes `WithBase(base_part, model, extras)`.

- [ ] **Step 4: Run to verify it passes**

Run: Spec tests, then input tests
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add src/Support/Parameters/Spec/model.jl test/unit_tests/Support/Parameters/Spec/ut_model.jl
git commit -m "ParameterSpec keeps the values of indexed keys (WithBase.extras)"
```

---

### Task 2: One Hooke-matrix implementation for dict and typed materials; typed flaw function

**Files:**
- Modify: `src/Models/Material/Material_Basis.jl` (`_hooke_matrix`, `get_Hooke_matrix`, `hooke_matrix`, typed `flaw_function`)
- Modify: `src/Models/Material/Material_Factory.jl` (`BlockMaterial` gains `hooke_symmetry` and `extras`; `hooke_symmetry`; `block_material`)
- Test: `test/unit_tests/Models/Material/ut_block_material.jl` (append)

**Interfaces:**
- Consumes (Task 1): `WithBase.extras`.
- Produces:
  - Material factory: `BlockMaterial{B,M,E}` fields in this order: `base, model, symmetry::String, hooke_symmetry::String, moduli::E, tables::Vector{Table1D}, extras::Dict{String,Any}`.
  - `hooke_symmetry(symmetry::Union{Nothing,String}, dof::Int64)::String`: the symmetry string the Hooke matrix uses. In 3D, a trailing `plane strain` / `plane stress` is stripped (as `check_symmetry` does to the dict); a missing symmetry gives `"isotropic"` (as `write_moduli!` sets in the dict).
  - Material_Basis:
    - `hooke_matrix(material, dof::Int64, ID::Int64 = 1)`: the Hooke matrix of a `BlockMaterial` at node `ID`.
    - `get_Hooke_matrix(parameter::Dict, symmetry, dof, ID = 1)`: unchanged behaviour.
    - `flaw_function(flaw, coor, stress)` for `flaw::Nothing` or a `FlawFunctionParams`.

- [ ] **Step 1: Write the failing tests** (append to `ut_block_material.jl`)

```julia
function ut_legacy_dict(raw; dof)
    ut_reset(dof)
    legacy = Dict{String,Any}(raw)
    BBASIS.get_all_elastic_moduli(legacy)
    return legacy
end

@testset "typed Hooke matrix equals legacy" begin
    ortho = Dict("Young's Modulus X" => 2.0e3, "Young's Modulus Y" => 1.5e3,
                 "Young's Modulus Z" => 1.0e3, "Poisson's Ratio XY" => 0.3,
                 "Poisson's Ratio YZ" => 0.25, "Poisson's Ratio XZ" => 0.2,
                 "Shear Modulus XY" => 700.0, "Shear Modulus YZ" => 600.0,
                 "Shear Modulus XZ" => 500.0)
    aniso = Dict("C$i$j" => (i == j ? 100.0 * i : 1.0 * i + j) for i in 1:6 for j in i:6)
    cases = [(3, Dict("Symmetry" => "isotropic", "Bulk Modulus" => 10.0, "Shear Modulus" => 4.0)),
             (2, Dict("Symmetry" => "isotropic plane strain", "Bulk Modulus" => 10.0,
                      "Shear Modulus" => 4.0)),
             (2, Dict("Symmetry" => "isotropic plane stress", "Young's Modulus" => 10.0,
                      "Poisson's Ratio" => 0.3)),
             (2, Dict("Symmetry" => "something else", "Young's Modulus" => 10.0,
                      "Poisson's Ratio" => 0.3)),
             (3, merge(Dict("Symmetry" => "orthotropic"), ortho)),
             (2, merge(Dict("Symmetry" => "orthotropic plane stress"), ortho)),
             (3, merge(Dict("Symmetry" => "transverse isotropic"), ortho)),
             (2, merge(Dict("Symmetry" => "transverse isotropic plane strain"), ortho)),
             (2, merge(Dict("Symmetry" => "transverse isotropic plane stress"), ortho)),
             (3, merge(Dict{String,Any}("Symmetry" => "anisotropic"), aniso)),
             (2, merge(Dict{String,Any}("Symmetry" => "anisotropic plane strain"), aniso))]
    for (dof, raw) in cases
        raw = merge(Dict{String,Any}("Material Model" => "Correspondence Elastic"), raw)
        legacy = ut_legacy_dict(raw; dof = dof)
        ut_reset(dof)
        m = typed_block_material(raw; dof = dof)
        @test Matrix(BBASIS.hooke_matrix(m, dof, 2)) ≈
              Matrix(BBASIS.get_Hooke_matrix(legacy, legacy["Symmetry"], dof, 2))
    end
end

@testset "hooke_symmetry follows the dict the legacy path used" begin
    @test BMAT.hooke_symmetry(nothing, 3) == "isotropic"
    @test BMAT.hooke_symmetry("isotropic plane strain", 3) == "isotropic "
    @test BMAT.hooke_symmetry("isotropic plane strain", 2) == "isotropic plane strain"
    @test BMAT.hooke_symmetry("orthotropic", 3) == "orthotropic"
end

@testset "Hooke matrix from a table" begin
    ut_reset(3; nnodes = 2)
    dir = mktempdir()
    write(joinpath(dir, "ex.txt"), "header: Temperature Young's_Modulus_X\n0 1000\n100 3000\n")
    N, NP1 = PeriLab.Data_Manager.create_node_scalar_field("Temperature", Float64)
    NP1 .= [0.0, 100.0]
    raw = Dict{String,Any}("Material Model" => "Correspondence Elastic",
                           "Symmetry" => "orthotropic", "Young's Modulus X" => "ex.txt",
                           "Young's Modulus Y" => 1.5e3, "Young's Modulus Z" => 1.0e3,
                           "Poisson's Ratio XY" => 0.3, "Poisson's Ratio YZ" => 0.25,
                           "Poisson's Ratio XZ" => 0.2, "Shear Modulus XY" => 700.0,
                           "Shear Modulus YZ" => 600.0, "Shear Modulus XZ" => 500.0)
    ctx = PeriLab.ParameterSpec.ParseContext(directory = dir)
    wb = PeriLab.ParameterSpec.parse_model(:material, raw, "M", ctx;
                                           name_key = "Material Model")
    @test isempty(ctx.errors)
    m = BMAT.block_material(wb, "Correspondence Elastic", 3)
    BMAT.bind_material!(m)
    constant_x(E) = begin
        r = copy(raw)
        r["Young's Modulus X"] = E
        ut_reset(3; nnodes = 2)
        BBASIS.hooke_matrix(typed_block_material(r; dof = 3), 3, 1)
    end
    @test Matrix(BBASIS.hooke_matrix(m, 3, 1)) ≈ Matrix(constant_x(1000.0))
    @test Matrix(BBASIS.hooke_matrix(m, 3, 2)) ≈ Matrix(constant_x(3000.0))
end

@testset "typed flaw function" begin
    ut_reset(3)
    flaw = Dict("Active" => true, "Function" => "Pre-defined", "Flaw Size" => 0.2,
                "Flaw Magnitude" => 0.5, "Flaw Location X" => 1.0, "Flaw Location Y" => 0.5)
    m = typed_block_material(Dict("Material Model" => "PD Solid Elastic",
                                  "Bulk Modulus" => 1.0, "Shear Modulus" => 1.0,
                                  "Flaw Function" => flaw))
    for coor in ([1.0, 0.5], [1.1, 0.4, 0.2], [3.0, 3.0])
        @test BBASIS.flaw_function(m.base.flaw_function, coor, 10.0) ≈
              BBASIS.flaw_function(Dict("Flaw Function" => flaw), coor, 10.0)
    end
    @test BBASIS.flaw_function(nothing, [0.0, 0.0], 10.0) == 10.0
    inactive = typed_block_material(Dict("Material Model" => "PD Solid Elastic",
                                         "Bulk Modulus" => 1.0, "Shear Modulus" => 1.0,
                                         "Flaw Function" => Dict("Active" => false,
                                                                 "Function" => "Pre-defined")))
    @test BBASIS.flaw_function(inactive.base.flaw_function, [0.0, 0.0], 10.0) == 10.0
    incomplete = typed_block_material(Dict("Material Model" => "PD Solid Elastic",
                                           "Bulk Modulus" => 1.0, "Shear Modulus" => 1.0,
                                           "Flaw Function" => Dict("Active" => true,
                                                                   "Function" => "Pre-defined")))
    @test_logs (:error,
                "An active Flaw Function needs Flaw Size, Flaw Magnitude, Flaw Location X and Flaw Location Y.") @test_throws PeriLab.PeriLabError begin
        BBASIS.flaw_function(incomplete.base.flaw_function, [0.0, 0.0], 10.0)
    end
end

@testset "extras reach the block material" begin
    ut_reset(3)
    m = typed_block_material(Dict("Material Model" => "Correspondence UMAT", "File" => "x.so",
                                  "Number of Properties" => 2, "Property_1" => 3.0,
                                  "Young's Modulus" => 1.0, "Poisson's Ratio" => 0.3))
    @test m.extras == Dict{String,Any}("Property_1" => 3.0)
    @test m.hooke_symmetry == "isotropic"
end
```

- [ ] **Step 2: Run to verify they fail**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Models/Material/ut_block_material`
Expected: `UndefVarError: hooke_matrix` / `hooke_symmetry`; the flaw tests fail with a `MethodError` for `flaw_function(::FlawFunctionParams…)`.

- [ ] **Step 3: Implement**

`src/Models/Material/Material_Basis.jl`:
- Add `using ......ParameterSpec: value` after the Data_Manager import.
- Refactor `get_Hooke_matrix` mechanically with this script. Run it once from the repo root; it fails loudly if a pattern is missing.

```julia
f = "src/Models/Material/Material_Basis.jl"
s = read(f, String)
function sub!(a, b; all = false)
    occursin(a, s) || error("missing: " * first(string(a), 80))
    global s = all ? replace(s, a => b) : replace(s, a => b; count = 1)
end
sub!("function get_Hooke_matrix(parameter::Dict,\n                          symmetry::String,\n                          dof::Int64,\n                          ID::Int64 = 1)",
     "function _hooke_matrix(c, symmetry::String, dof::Int64, ID::Int64)")
sub!(r"get_dependent_value_with_ID\(\"C\" \* string\(iID\) \* string\(jID\),\s*parameter\)",
     "c(\"C\" * string(iID) * string(jID), 1)")
sub!(r"get_dependent_value_with_ID\((\"[^\"]+\"), parameter, ID\)", s"c(\1, ID)"; all = true)
sub!("    iID = ID\n    if parameter[\"Poisson's Ratio\"] isa Float64\n        iID = 1\n    end\n", "")
sub!("parameter[\"Poisson's Ratio\"][iID]", "c(\"Poisson's Ratio\", ID)"; all = true)
sub!("parameter[\"Young's Modulus\"][iID]", "c(\"Young's Modulus\", ID)"; all = true)
sub!("parameter[\"Shear Modulus\"][iID]", "c(\"Shear Modulus\", ID)"; all = true)
occursin(r"_hooke_matrix[\s\S]*parameter\[", s[findfirst("function _hooke_matrix", s)[1]:end][1:min(end, 9000)]) &&
    error("a parameter[...] read is left in _hooke_matrix")
write(f, s)
```

  Then add directly above the (renamed) `function _hooke_matrix`:

```julia
# constant of a material dict at node `id` (the legacy reads of get_Hooke_matrix)
function _dict_constant(parameter::Dict, key::String, id::Int64)
    if key in ("Poisson's Ratio", "Young's Modulus", "Shear Modulus")
        iID = parameter["Poisson's Ratio"] isa Float64 ? 1 : id
        return parameter[key][iID]
    end
    return get_dependent_value_with_ID(key, parameter, id)
end

"""
	get_Hooke_matrix(parameter::Dict, symmetry::String, dof::Int64, ID::Int64=1)

Returns the Hooke matrix of the material (material dict; see `hooke_matrix`
for typed block materials).
"""
function get_Hooke_matrix(parameter::Dict, symmetry::String, dof::Int64, ID::Int64 = 1)
    return _hooke_matrix((key, id) -> _dict_constant(parameter, key, id), symmetry, dof, ID)
end

const _HOOKE_BASE_FIELDS = Dict("Young's Modulus X" => :youngs_modulus_x,
                                "Young's Modulus Y" => :youngs_modulus_y,
                                "Young's Modulus Z" => :youngs_modulus_z,
                                "Poisson's Ratio XY" => :poissons_ratio_xy,
                                "Poisson's Ratio YZ" => :poissons_ratio_yz,
                                "Poisson's Ratio XZ" => :poissons_ratio_xz,
                                "Shear Modulus XY" => :shear_modulus_xy,
                                "Shear Modulus YZ" => :shear_modulus_yz,
                                "Shear Modulus XZ" => :shear_modulus_xz)

@inline _at(x::Real, ::Int64) = x
@inline _at(x::AbstractVector, id::Int64) = x[id]

# constant of a typed block material at node `id`
function _typed_constant(material, key::String, id::Int64)
    key == "Poisson's Ratio" && return _at(material.moduli.poissons_ratio, id)
    key == "Young's Modulus" && return _at(material.moduli.youngs_modulus, id)
    key == "Shear Modulus" && return _at(material.moduli.shear_modulus, id)
    startswith(key, "C") && return getfield(material.base, Symbol(lowercase(key)))
    return value(getfield(material.base, _HOOKE_BASE_FIELDS[key]), id)
end

"""
	hooke_matrix(material, dof, ID = 1)

Hooke matrix of a typed block material (`BlockMaterial`) at node `ID`. Dependent
constants must be bound (`Material.bind_material!`).
"""
function hooke_matrix(material, dof::Int64, ID::Int64 = 1)
    return _hooke_matrix((key, id) -> _typed_constant(material, key, id),
                         material.hooke_symmetry, dof, ID)
end
```

  Delete the old docstring directly above `_hooke_matrix` if it still names `get_Hooke_matrix`. Replace it with a one-line comment: `# formulas of the Hooke matrix; c(key, id) returns a constant at node id`.
- Add `export hooke_matrix` next to the existing exports.
- After the legacy `flaw_function(params::Dict, …)` add:

```julia
flaw_function(::Nothing, coor::AbstractVector{<:Real}, stress::Union{Int64,Float64}) = Float64(stress)

# typed flaw function (`FlawFunctionParams`); same formula as the dict version
function flaw_function(flaw, coor::AbstractVector{<:Real},
                       stress::T)::Float64 where {T<:Union{Int64,Float64}}
    flaw.active || return stress
    if flaw.flaw_size === nothing || flaw.flaw_magnitude === nothing ||
       flaw.flaw_location_x === nothing || flaw.flaw_location_y === nothing
        @abort "An active Flaw Function needs Flaw Size, Flaw Magnitude, Flaw Location X and Flaw Location Y."
    end
    flaw_size::Float64 = flaw.flaw_size
    flaw_magnitude::Float64 = flaw.flaw_magnitude
    if !(0 < flaw_magnitude <= 1)
        @abort "Flaw Magnitude should be between 0 and 1"
    end
    if flaw_size <= 0
        @abort "Flaw Size must be positive."
    end
    dx = Float64(coor[1]) - flaw.flaw_location_x
    dy = Float64(coor[2]) - flaw.flaw_location_y
    distance_squared = dx * dx + dy * dy
    if length(coor) == 3
        dz = Float64(coor[3]) - something(flaw.flaw_location_z, 0.0)
        distance_squared += dz * dz
    end
    return stress *
           (1 - flaw_magnitude * exp(-distance_squared / (flaw_size * flaw_size)))
end
```

`src/Models/Material/Material_Factory.jl`:
- Replace the `BlockMaterial` struct with

```julia
struct BlockMaterial{B,M,E}
    base::B
    model::M
    symmetry::String          # used by the force models: "plane strain", "plane stress" or "3D"
    hooke_symmetry::String    # used by the Hooke matrix (see hooke_symmetry)
    moduli::E
    tables::Vector{Table1D}   # dependent tables of base and model, re-bound every step
    extras::Dict{String,Any}  # values of indexed keys, e.g. Property_1
end
```

- Add after `material_symmetry`:

```julia
"""
    hooke_symmetry(symmetry, dof)

The symmetry string the Hooke matrix uses: as given, with a trailing plane strain /
plane stress removed in 3D (as `check_symmetry` does), `"isotropic"` if missing.
"""
function hooke_symmetry(symmetry::Union{Nothing,String}, dof::Int64)
    symmetry === nothing && return "isotropic"
    dof == 3 || return symmetry
    return replace(replace(symmetry, r"plane strain$" => ""), r"plane stress$" => "")
end
```

- In `block_material`, the `return BlockMaterial(…)` becomes

```julia
    return BlockMaterial(wb.base, wb.model, material_symmetry(wb.base.symmetry, dof),
                         hooke_symmetry(wb.base.symmetry, dof),
                         elastic_moduli(wb.base, bond_based, dof),
                         vcat(_tables(wb.base), _tables(wb.model)), wb.extras)
```

- [ ] **Step 4: Run to verify they pass**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Models/Material/ut_block_material unit_tests/Models/Material/ut_material_basis`
Expected: all pass. `ut_material_basis`'s existing `get_Hooke_matrix` and `flaw_function` tests prove that the refactor kept the legacy results.

- [ ] **Step 5: Full suite, then commit**

Full suite in background; read the tail. Expected: 0 failures, 0 errors.

```bash
git add src/Models/Material/Material_Basis.jl src/Models/Material/Material_Factory.jl \
        test/unit_tests/Models/Material/ut_block_material.jl
git commit -m "One Hooke-matrix implementation for dict and typed materials; typed flaw function"
```

---

### Task 3: Typed methods for Correspondence Elastic, Correspondence Plastic and zero energy control

**Files:**
- Modify: `src/Models/Material/Material_Models/Correspondence/Correspondence_Elastic.jl`, `Correspondence_Plastic.jl`
- Modify: `src/Models/Material/Material_Models/Zero_Energy_Control/Zero_Energy_Control.jl`, `Global_Zero_Energy_Control.jl`
- Test: `test/unit_tests/Models/Material/ut_block_material.jl` (append)

**Interfaces:**
- Consumes (Task 2): `Material_Basis.hooke_matrix(material, dof, ID)`, `flaw_function(flaw, coor, stress)`, `BlockMaterial` fields.
- Produces the typed correspondence model interface (the dict methods stay):
  - `init_model(nodes::AbstractVector{Int64}, p::<Params>, material)`
  - `compute_stresses(nodes::AbstractVector{Int64}, dof::Int64, p::<Params>, material, time::Float64, dt::Float64, strain_increment, stress_N, stress_NP1)`
  - `compute_stresses_ba(nodes, nlist, dof::Int64, p::<Params>, material, time::Float64, dt::Float64, strain_increment, stress_N, stress_NP1)`
  - `Zero_Energy_Control.init_model(nodes, material, block::Int64)` and `Zero_Energy_Control.compute_zero_energy_control(nodes, material, block::Int64, time::Float64, dt::Float64)`.
  - `Global_Zero_Energy_Control.init_model(nodes, material)` and `compute_control(nodes, material, time::Float64, dt::Float64)`.
  - In each case the dict method has a `Dict` argument type and so wins dispatch for dicts.

- [ ] **Step 1: Write the failing tests** (append to `ut_block_material.jl`)

```julia
const BCORR = BMAT.Correspondence

@testset "Correspondence Elastic typed init equals legacy" begin
    for (dof, raw) in ((3, Dict("Symmetry" => "isotropic", "Bulk Modulus" => 10.0,
                                "Shear Modulus" => 4.0)),
                       (2, Dict("Symmetry" => "isotropic plane strain",
                                "Young's Modulus" => 10.0, "Poisson's Ratio" => 0.3)))
        raw = merge(Dict{String,Any}("Material Model" => "Correspondence Elastic"), raw)
        legacy = ut_legacy_dict(raw; dof = dof)
        ut_reset(dof; nnodes = 2)
        m = typed_block_material(raw; dof = dof)
        BCORR.Correspondence_Elastic.init_model([1, 2], m.model, m)
        C = PeriLab.Data_Manager.get_field("Material Gradient")
        for iID in 1:2
            @test C[iID, :, :] ≈ Matrix(BBASIS.get_Hooke_matrix(legacy, legacy["Symmetry"],
                                                                dof, iID))
        end
        @test hasmethod(BCORR.Correspondence_Elastic.compute_stresses,
                        Tuple{Vector{Int64},Int64,typeof(m.model),typeof(m),Float64,Float64,
                              Array{Float64,3},Array{Float64,3},Array{Float64,3}})
    end
end

@testset "Correspondence Plastic typed compute equals legacy" begin
    dof = 3
    nnodes = 2
    raw = Dict{String,Any}("Material Model" => "Correspondence Elastic + Correspondence Plastic",
                           "Symmetry" => "isotropic", "Bulk Modulus" => 10.0,
                           "Shear Modulus" => 4.0, "Yield Stress" => 0.01)
    results = []
    for typed in (false, true)
        ut_reset(dof; nnodes = nnodes)
        coor = PeriLab.Data_Manager.create_constant_node_vector_field("Coordinates", Float64,
                                                                      dof)
        coor .= [0.0 0.0 0.0; 1.0 0.0 0.0]
        if typed
            m = typed_block_material(raw; dof = dof)
            p = m.model.parts[2]
            BCORR.Correspondence_Plastic.init_model(collect(1:nnodes), p, m)
        else
            legacy = Dict{String,Any}(raw)
            BBASIS.get_all_elastic_moduli(legacy)
            BCORR.Correspondence_Plastic.init_model(collect(1:nnodes), legacy)
        end
        strain_inc = zeros(nnodes, dof, dof)
        strain_inc[:, 1, 1] .= 0.02
        strain_inc[:, 1, 2] .= 0.01
        strain_inc[:, 2, 1] .= 0.01
        stress_N = zeros(nnodes, dof, dof)
        stress_NP1 = zeros(nnodes, dof, dof)
        stress_NP1[:, 1, 1] .= 0.5
        if typed
            BCORR.Correspondence_Plastic.compute_stresses(collect(1:nnodes), dof, p, m, 0.0,
                                                          1.0, strain_inc, stress_N,
                                                          stress_NP1)
        else
            BCORR.Correspondence_Plastic.compute_stresses(collect(1:nnodes), dof, legacy,
                                                          0.0, 1.0, strain_inc, stress_N,
                                                          stress_NP1)
        end
        push!(results, copy(stress_NP1))
    end
    @test results[1] ≈ results[2]
end

@testset "zero energy control skips UMAT" begin
    ut_reset(3; nnodes = 2)
    elastic = typed_block_material(Dict("Material Model" => "Correspondence Elastic",
                                        "Bulk Modulus" => 1.0, "Shear Modulus" => 1.0))
    umat = typed_block_material(Dict("Material Model" => "Correspondence UMAT",
                                     "File" => "x.so", "Number of Properties" => 1,
                                     "Young's Modulus" => 1.0, "Poisson's Ratio" => 0.3))
    vumat = typed_block_material(Dict("Material Model" => "Correspondence VUMAT",
                                      "File" => "x.so", "Number of Properties" => 1,
                                      "Young's Modulus" => 1.0, "Poisson's Ratio" => 0.3))
    GZEC = PeriLab.Solver_Manager.Zero_Energy_Control.Global_Zero_Energy_Control
    @test !GZEC.is_umat(elastic)
    @test GZEC.is_umat(umat)
    @test !GZEC.is_umat(vumat)
end

@testset "typed zero energy control init" begin
    ut_reset(3; nnodes = 2)
    m = typed_block_material(Dict("Material Model" => "Correspondence Elastic",
                                  "Bulk Modulus" => 10.0, "Shear Modulus" => 4.0,
                                  "Zero Energy Control" => "Global"))
    ZEC = PeriLab.Solver_Manager.Zero_Energy_Control
    ZEC.init_model([1, 2], m, 1)
    @test PeriLab.Data_Manager.get_analysis_model("Zero Energy Control Model", 1) == ["Global"]
    @test PeriLab.Data_Manager.get_field("Material Gradient")[2, :, :] ≈
          Matrix(BBASIS.hooke_matrix(m, 3, 2))
    ut_reset(3; nnodes = 2)
    plain = typed_block_material(Dict("Material Model" => "Correspondence Elastic",
                                      "Bulk Modulus" => 10.0, "Shear Modulus" => 4.0))
    ZEC.init_model([1, 2], plain, 1)
    @test PeriLab.Data_Manager.get_analysis_model("Zero Energy Control Model", 1) == [""]
end
```

- [ ] **Step 2: Run to verify they fail**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Models/Material/ut_block_material`
Expected: `MethodError`s for the typed `init_model` / `compute_stresses`, and `UndefVarError: is_umat`.

If `get_analysis_model` returns a value of another shape than a one-element vector, print it and adapt the two `==` assertions to that shape. Ledger it.

- [ ] **Step 3: Implement**

`Correspondence_Elastic.jl`:
- Change `using .....Material_Basis: get_Hooke_matrix` to `using .....Material_Basis: get_Hooke_matrix, hooke_matrix`.
- Rename the body of the dict `compute_stresses(nodes::AbstractVector{Int64}, dof::Int64, material_parameter::Dict, …)` into a shared kernel, without changing its body. Its signature becomes

```julia
function _elastic_stresses!(nodes::AbstractVector{Int64},
                            dof::Int64,
                            strain_increment::NodeTensorField{Float64},
                            stress_N::NodeTensorField{Float64},
                            stress_NP1::NodeTensorField{Float64})
```

- Add after it:

```julia
function compute_stresses(nodes::AbstractVector{Int64}, dof::Int64, material_parameter::Dict,
                          time::Float64, dt::Float64,
                          strain_increment::NodeTensorField{Float64},
                          stress_N::NodeTensorField{Float64},
                          stress_NP1::NodeTensorField{Float64})
    return _elastic_stresses!(nodes, dof, strain_increment, stress_N, stress_NP1)
end

function compute_stresses(nodes::AbstractVector{Int64}, dof::Int64,
                          p::CorrespondenceElasticParams, material,
                          time::Float64, dt::Float64,
                          strain_increment::NodeTensorField{Float64},
                          stress_N::NodeTensorField{Float64},
                          stress_NP1::NodeTensorField{Float64})
    return _elastic_stresses!(nodes, dof, strain_increment, stress_N, stress_NP1)
end

function init_model(nodes::AbstractVector{Int64}, p::CorrespondenceElasticParams, material)
    dof::Int64 = Data_Manager.get_dof()
    hooke::NodeTensorField{Float64} = Data_Manager.create_constant_node_tensor_field("Material Gradient",
                                                                                     Float64,
                                                                                     Int64((dof *
                                                                                            (dof +
                                                                                             1)) /
                                                                                           2))
    for iID in nodes
        @views hooke[iID, :, :] = hooke_matrix(material, dof, iID)
    end
end

function compute_stresses_ba(nodes, nlist, dof::Int64, p::CorrespondenceElasticParams,
                             material, time::Float64, dt::Float64, strain_increment,
                             stress_N, stress_NP1)
    @views mapping = get_mapping(dof)
    for iID in nodes
        @views hookeMatrix = hooke_matrix(material, dof, iID)
        @fastmath @inbounds @simd for jID in eachindex(nlist[iID])
            @views sNP1 = stress_NP1[iID][jID, :, :]
            @views sInc = strain_increment[iID][jID, :, :]
            @views sN = stress_N[iID][jID, :, :]
            fast_mul!(sNP1, hookeMatrix, sInc, sN, mapping)
        end
    end
end
```

  `CorrespondenceElasticParams` is declared after `correspondence_name()`, so move its `@params` block and `__init__` line to just after the module's `export` lines, as was done for the non-correspondence modules in 3b. The same move applies to Plastic, UMAT and VUMAT in this task and in Task 4.

`Correspondence_Plastic.jl`:
- Import `flaw_function` is already there. Add `using .......ParameterSpec: value`, which may be merged into the existing ParameterSpec import.
- Move the `@params` block up as described above.
- Add:

```julia
function init_model(nodes::AbstractVector{Int64}, p::CorrespondencePlasticParams, material)
    if material.moduli === nothing
        @abort "Shear Modulus must be defined to be able to run this plastic material"
        return
    end
    Data_Manager.create_node_scalar_field("von Mises Yield Stress", Float64)
    Data_Manager.create_node_scalar_field("Plastic Strain", Float64)
    if material.base.bond_associated
        Data_Manager.create_bond_scalar_state("von Mises Bond Yield Stress", Float64)
        Data_Manager.create_bond_scalar_state("Plastic Bond Strain", Float64)
    end
end
```

- For `compute_stresses` and `compute_stresses_ba`: copy each dict method and change only these lines in the copy:
  - The signature's `material_parameter::Dict` becomes `p::CorrespondencePlasticParams, material`.
  - `yield_stress_fn = get_dependent_value("Yield Stress", material_parameter)` is deleted.
  - `yield_stress = yield_stress_fn(iID)` becomes `yield_stress = value(p.yield_stress, iID)`.
  - `flaw_function(material_parameter, ` becomes `flaw_function(material.base.flaw_function, `.
  - `material_parameter["Shear Modulus"]` becomes `material.moduli.shear_modulus`.

  The two copies are the typed methods. Run `grep -n 'material_parameter' ` on the typed methods only: none may remain.

`Zero_Energy_Control.jl`: add below the dict methods:

```julia
function init_model(nodes::AbstractVector{Int64}, material, block::Int64)
    zero_energy_model = material.base.zero_energy_control
    if zero_energy_model !== nothing
        @debug "Init zero energy control model ''$zero_energy_model'' at block $block."
        Data_Manager.set_analysis_model("Zero Energy Control Model", block,
                                        zero_energy_model)
        mod = create_module_specifics(zero_energy_model,
                                      module_list,
                                      @__MODULE__,
                                      "control_name")
        Data_Manager.set_model_module(zero_energy_model, mod)
        mod.init_model(nodes, material)
    else
        Data_Manager.set_analysis_model("Zero Energy Control Model", block, "")
        @warn "No zero energy control activated for corresponcence in block $block. This might cause errors."
    end
end

function compute_zero_energy_control(nodes::AbstractVector{Int64}, material, block::Int64,
                                     time::Float64, dt::Float64)
    for zero_energy_model in Data_Manager.get_analysis_model("Zero Energy Control Model",
                                                             block)
        zero_energy_model == "" && continue
        mod = Data_Manager.get_model_module(zero_energy_model)
        mod.compute_control(nodes, material, time, dt)
    end
end
```

`Global_Zero_Energy_Control.jl`:
- Change the Material_Basis import to `using ...Material_Basis: get_Hooke_matrix, hooke_matrix`.
- Add:

```julia
# the legacy dict marked UMAT materials by the key "UMAT Material Name"
_model_parts(model) = hasfield(typeof(model), :parts) ? model.parts : (model,)
is_umat(material) = any(part -> nameof(typeof(part)) === :CorrespondenceUMATParams,
                        _model_parts(material.model))

function init_model(nodes::AbstractVector{Int64}, material)
    dof::Int64 = Data_Manager.get_dof()
    Data_Manager.create_constant_node_tensor_field("Zero Energy Stiffness", Float64, dof)
    "Material Gradient" in Data_Manager.get_all_field_keys() && return
    hooke::NodeTensorField{Float64} = Data_Manager.create_constant_node_tensor_field("Material Gradient",
                                                                                     Float64,
                                                                                     Int64((dof *
                                                                                            (dof +
                                                                                             1)) /
                                                                                           2))
    for iID in nodes
        @views hooke[iID, :, :] = hooke_matrix(material, dof, iID)
    end
end
```

- Rename the body of the dict `compute_control` into `_compute_control!(nodes, apply_zero_energy::Bool, time, dt)`. Inside it, replace `if !haskey(material_parameter, "UMAT Material Name")` with `if apply_zero_energy` and leave everything else as is. Then add:

```julia
compute_control(nodes::AbstractVector{Int64}, material_parameter::Dict{String,Any},
                time::Float64, dt::Float64) = _compute_control!(nodes,
                                                                !haskey(material_parameter,
                                                                        "UMAT Material Name"),
                                                                time, dt)
compute_control(nodes::AbstractVector{Int64}, material, time::Float64, dt::Float64) = _compute_control!(nodes,
                                                                                                        !is_umat(material),
                                                                                                        time,
                                                                                                        dt)
```

- [ ] **Step 4: Run to verify they pass**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Models/Material/ut_block_material unit_tests/Models/Material/Material_Models/Correspondence/ut_Correspondence_Plastic unit_tests/Models/Material/Zero_Energy_Control/ut_Global_Zero_Energy_Control`
Expected: all pass. The legacy tests still pass because the dict methods are unchanged in behaviour.

- [ ] **Step 5: Commit** (the typed methods are not used by a run yet, so the unit runs above are the gate)

```bash
git add src/Models/Material/Material_Models/Correspondence/Correspondence_Elastic.jl \
        src/Models/Material/Material_Models/Correspondence/Correspondence_Plastic.jl \
        src/Models/Material/Material_Models/Zero_Energy_Control \
        test/unit_tests/Models/Material/ut_block_material.jl
git commit -m "Typed methods for correspondence elastic/plastic and zero energy control"
```

---

### Task 4: Typed methods for Correspondence UMAT and VUMAT

**Files:**
- Modify: `src/Models/Material/Material_Models/Correspondence/Correspondence_UMAT.jl`, `Correspondence_VUMAT.jl`
- Test: `test/unit_tests/Models/Material/ut_block_material.jl` (append)

**Interfaces:**
- Consumes (Tasks 1–3): `material.extras`, `hooke_matrix`, the typed correspondence model interface.
- Produces:
  - `init_model(nodes, p::CorrespondenceUMATParams, material)` and `compute_stresses(nodes, dof, p::CorrespondenceUMATParams, material, time, dt, strain_increment, stress_N, stress_NP1)`, plus `compute_stresses_ba(…)`, which aborts as the dict version does. The same three methods exist for `CorrespondenceVUMATParams`.

- [ ] **Step 1: Write the failing tests** (append to `ut_block_material.jl`)

```julia
function ut_umat_file()
    file = "./src/Models/Material/UMATs/libperuser.so"
    isfile(file) || (file = "../src/Models/Material/UMATs/libperuser.so")
    return file
end

@testset "UMAT properties from extras" begin
    ut_reset(3; nnodes = 2)
    UMAT = BCORR.Correspondence_UMAT
    file = ut_umat_file()
    m = typed_block_material(Dict("Material Model" => "Correspondence UMAT", "File" => file,
                                  "Number of Properties" => 3, "Property_1" => 2,
                                  "Property_3" => 2.4, "Young's Modulus" => 2.0,
                                  "Poisson's Ratio" => 0.1))
    UMAT.init_model([1, 2], m.model, m)
    props = PeriLab.Data_Manager.get_field("Properties")
    @test props[1] == 2.0 && props[2] == 0.0 && props[3] == 2.4
    @test PeriLab.Data_Manager.get_field("Material Gradient")[1, :, :] ≈
          Matrix(BBASIS.hooke_matrix(m, 3, 1))
    @test UMAT.umat_file_path == joinpath(pwd(), PeriLab.Data_Manager.get_directory(), file)
end

@testset "UMAT and VUMAT typed init errors" begin
    for (UM, name_key, model) in ((BCORR.Correspondence_UMAT, "UMAT Material Name",
                                   "Correspondence UMAT"),
                                  (BCORR.Correspondence_VUMAT, "VUMAT Material Name",
                                   "Correspondence VUMAT"))
        ut_reset(3; nnodes = 2)
        file = ut_umat_file()
        missing_file = typed_block_material(Dict("Material Model" => model,
                                                 "File" => file * "_not_there",
                                                 "Number of Properties" => 1,
                                                 "Young's Modulus" => 2.0,
                                                 "Poisson's Ratio" => 0.1))
        @test_logs (:error,
                    "File $(joinpath(pwd(), PeriLab.Data_Manager.get_directory(), file * "_not_there")) does not exist, please check name and directory.") @test_throws PeriLab.PeriLabError begin
            UM.init_model([1, 2], missing_file.model, missing_file)
        end
        long_name = typed_block_material(Dict("Material Model" => model, "File" => file,
                                              "Number of Properties" => 1,
                                              name_key => "a"^81,
                                              "Young's Modulus" => 2.0,
                                              "Poisson's Ratio" => 0.1))
        @test_logs (:error,
                    "Due to old Fortran standards only a name length of 80 is supported") @test_throws PeriLab.PeriLabError begin
            UM.init_model([1, 2], long_name.model, long_name)
        end
    end
end

@testset "UMAT predefined fields (typed)" begin
    ut_reset(3; nnodes = 2)
    t2 = PeriLab.Data_Manager.create_constant_node_scalar_field("test_field_2", Float64)
    t2[1] = 7.3
    t3 = PeriLab.Data_Manager.create_constant_node_scalar_field("test_field_3", Float64)
    t3 .= 3
    m = typed_block_material(Dict("Material Model" => "Correspondence UMAT",
                                  "File" => ut_umat_file(), "Number of Properties" => 1,
                                  "Predefined Field Names" => "test_field_2 test_field_3",
                                  "Young's Modulus" => 2.0, "Poisson's Ratio" => 0.1))
    BCORR.Correspondence_UMAT.init_model([1, 2], m.model, m)
    fields = PeriLab.Data_Manager.get_field("Predefined Fields")
    @test fields[1, 1] == 7.3 && fields[2, 2] == 3.0
end
```

- [ ] **Step 2: Run to verify they fail**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Models/Material/ut_block_material`
Expected: `MethodError` for `init_model(::Vector{Int64}, ::CorrespondenceUMATParams, ::BlockMaterial…)` (and the VUMAT equivalent).

- [ ] **Step 3: Implement**

For both modules:
- Move the `@params` block, `key_patterns` and `__init__` to just after the `export` lines.
- Import `hooke_matrix` with `get_Hooke_matrix` (UMAT only).

`Correspondence_UMAT.jl`, typed init (add after the dict `init_model`):

```julia
function init_model(nodes::AbstractVector{Int64}, p::CorrespondenceUMATParams, material)
    num_state_vars::Int64 = something(p.number_of_state_variables, 1)
    file = joinpath(pwd(), Data_Manager.get_directory(), p.file)
    global umat_file_path = file
    if !isfile(file)
        @abort "File $file does not exist, please check name and directory."
        return 1
    end
    if num_state_vars == 1
        Data_Manager.create_constant_node_scalar_field("State Variables", Float64)
    else
        Data_Manager.create_constant_node_vector_field("State Variables", Float64,
                                                       num_state_vars)
    end
    properties = Data_Manager.create_constant_free_size_field("Properties", Float64,
                                                              (p.number_of_properties, 1))
    for iID in 1:p.number_of_properties
        if !haskey(material.extras, "Property_$iID")
            @warn "Property_$iID is missing. Make sure that all properties are defined."
            properties[iID] = 0.0
        else
            properties[iID] = material.extras["Property_$iID"]
        end
    end
    if p.umat_material_name === nothing
        @warn "No UMAT Material Name is defined. Please check if you use it as method to check different material in your UMAT."
    elseif length(p.umat_material_name) > 80
        @abort "Due to old Fortran standards only a name length of 80 is supported"
    end
    _init_umat_fields!(nodes, p.predefined_field_names)
    for iID in nodes
        @views Data_Manager.get_field("Material Gradient")[iID, :, :] = hooke_matrix(_state_scaled(material),
                                                                                    Data_Manager.get_dof(),
                                                                                    iID)
    end
end

# moduli scaled by the state variable named by State Factor ID (the legacy init re-ran
# get_all_elastic_moduli after creating the State Variables field)
function _state_scaled(material)
    id = material.base.state_factor_id
    (id === nothing || material.moduli === nothing) && return material
    factor = Data_Manager.get_field("State Variables")[:, id]
    m = material.moduli
    moduli = (bulk_modulus = m.bulk_modulus .* factor,
              youngs_modulus = m.youngs_modulus .* factor,
              shear_modulus = m.shear_modulus .* factor,
              poissons_ratio = m.poissons_ratio)
    Data_Manager.get_field("Bulk_Modulus") .= moduli.bulk_modulus
    Data_Manager.get_field("Young's_Modulus") .= moduli.youngs_modulus
    Data_Manager.get_field("Shear_Modulus") .= moduli.shear_modulus
    return (base = material.base, moduli = moduli, hooke_symmetry = material.hooke_symmetry)
end
```

- Refactor the dict `init_model`:
  - Everything from `dof = Data_Manager.get_dof()` (after the `UMAT name` default) down to and including the `zStiff = Data_Manager.create_constant_node_tensor_field("Zero Energy Stiffness", …)` statement moves, unchanged, into `function _init_umat_fields!(nodes::AbstractVector{Int64}, predefined_field_names::Union{Nothing,String})`.
  - Inside it, the condition `if haskey(material_parameter, "Predefined Field Names")` becomes `if predefined_field_names !== nothing`, and `split(material_parameter["Predefined Field Names"], " ")` becomes `split(predefined_field_names, " ")`.
  - The dict `init_model` then calls `_init_umat_fields!(nodes, get(material_parameter, "Predefined Field Names", nothing))` at that place.
  - Keep its `get_all_elastic_moduli` and `get_Hooke_matrix` lines after the call.
- Refactor the dict `compute_stresses` to a kernel `_umat_stresses!(nodes, dof, nstatev::Int64, nprops::Int64, cmname::String, time, dt, strain_increment, stress_N, stress_NP1)`:
  - `nstatev = material_parameter["Number of State Variables"]` and `nprops = material_parameter["Number of Properties"]` are deleted, because they are now arguments.
  - `malloc_cstring(material_parameter["UMAT Material Name"])` becomes `malloc_cstring(cmname)`.
- Then add the two thin methods:

```julia
compute_stresses(nodes::AbstractVector{Int64}, dof::Int64, material_parameter::Dict,
                 time::Float64, dt::Float64, strain_increment::AbstractArray{Float64,3},
                 stress_N::AbstractArray{Float64,3}, stress_NP1::AbstractArray{Float64,3}) = _umat_stresses!(nodes,
                                                                                                             dof,
                                                                                                             material_parameter["Number of State Variables"],
                                                                                                             material_parameter["Number of Properties"],
                                                                                                             material_parameter["UMAT Material Name"],
                                                                                                             time,
                                                                                                             dt,
                                                                                                             strain_increment,
                                                                                                             stress_N,
                                                                                                             stress_NP1)
compute_stresses(nodes::AbstractVector{Int64}, dof::Int64, p::CorrespondenceUMATParams,
                 material, time::Float64, dt::Float64,
                 strain_increment::AbstractArray{Float64,3},
                 stress_N::AbstractArray{Float64,3}, stress_NP1::AbstractArray{Float64,3}) = _umat_stresses!(nodes,
                                                                                                             dof,
                                                                                                             something(p.number_of_state_variables,
                                                                                                                       1),
                                                                                                             p.number_of_properties,
                                                                                                             something(p.umat_material_name,
                                                                                                                       ""),
                                                                                                             time,
                                                                                                             dt,
                                                                                                             strain_increment,
                                                                                                             stress_N,
                                                                                                             stress_NP1)
compute_stresses_ba(nodes, nlist, dof::Int64, p::CorrespondenceUMATParams, material,
                    time::Float64, dt::Float64, strain_increment, stress_N, stress_NP1) = @abort "$(correspondence_name()) not yet implemented for bond associated."
```

  `Number of State Variables` defaults to 1 in the typed compute, as in init. The dict compute read the key directly and failed with a `KeyError` without it. Ledger this.

`Correspondence_VUMAT.jl`: the same pattern.
- **Typed init:**
  - Use the same `File`, state variables and properties lines as above, with `global vumat_file_path = file`.
  - The name check uses `p.vumat_material_name` and the warning text "No VUMAT Material Name is defined. Please check if you use it as method to check different material in your VUMAT.".
  - The tail matches the dict init: everything after the `VUMAT name` default (`dof = …` through `create_node_scalar_field("Temperature", Float64)`). Move it into `_init_vumat_fields!()`, which both inits call.
- **Compute kernel:** `_vumat_stresses!(nodes, dof, nstatev, nprops, cmname::String, time, dt, strain_increment, stress_N, stress_NP1)`. The three dict reads (`"Number of State Variables"`, `"Number of Properties"`, `"VUMAT Material Name"`) become arguments. Add thin dict and typed `compute_stresses` methods as for UMAT, plus the typed `compute_stresses_ba` that aborts.

- [ ] **Step 4: Run to verify they pass**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Models/Material/ut_block_material unit_tests/Models/Material/Material_Models/Correspondence/ut_Correspondence_UMAT unit_tests/Models/Material/Material_Models/Correspondence/ut_Correspondence_VUMAT`
Expected: all pass. The legacy UMAT/VUMAT tests still pass on the refactored dict methods.

- [ ] **Step 5: Commit**

```bash
git add src/Models/Material/Material_Models/Correspondence/Correspondence_UMAT.jl \
        src/Models/Material/Material_Models/Correspondence/Correspondence_VUMAT.jl \
        test/unit_tests/Models/Material/ut_block_material.jl
git commit -m "Typed methods for correspondence UMAT and VUMAT"
```

---

### Task 5: Correspondence blocks run on the typed path

**Files:**
- Modify: `src/Models/Material/Material_Models/Correspondence/Correspondence.jl` (typed `init_model`, `compute_model`, `compute_correspondence_model`, `fields_for_local_synchronization`)
- Modify: `src/Models/Material/Material_Models/Correspondence/Bond_Associated_Correspondence.jl` (typed `init_model`, `compute_model`)
- Modify: `src/Models/Material/Material_Factory.jl` (correspondence branches of `init_model`, `fields_for_local_synchronization`, `compute_model`; `compute_correspondence_bond_forces`)
- Modify: `src/Models/Model_Factory.jl` (`compute_matrix_based_bond_forces`)
- Modify: `src/Models/Material/Material_Models/Correspondence/Correspondence_matrix_based.jl:622` (zero energy control init)
- Test: `test/unit_tests/Models/Material/ut_block_material.jl` (append)

**Interfaces:**
- Consumes (Tasks 2–4): the typed model interface, typed ZEC, `bind_material!`, `hooke_matrix`.
- Produces:
  - `Correspondence.init_model(nodes, block::Int64, material)`
  - `Correspondence.compute_model(nodes, material, block::Int64, time::Float64, dt::Float64)`
  - `Correspondence.compute_correspondence_model(nodes, material, block::Int64, time::Float64, dt::Float64)`
  - `Correspondence.fields_for_local_synchronization(model::String, block::Int64, material)`
  - `Bond_Associated_Correspondence.init_model(nodes, material)` and `compute_model(nodes, material, block::Int64, time::Float64, dt::Float64)`
  - `Material.compute_correspondence_bond_forces(nodes, material, block::Int64, time::Float64, dt::Float64)`

- [ ] **Step 1: Write the failing tests** (append to `ut_block_material.jl`)

```julia
@testset "typed correspondence dispatcher" begin
    ut_reset(3; nnodes = 2)
    PeriLab.Data_Manager.set_rotation(false)
    m = typed_block_material(Dict("Material Model" => "Correspondence Elastic",
                                  "Symmetry" => "isotropic", "Bulk Modulus" => 10.0,
                                  "Shear Modulus" => 4.0))
    BCORR.init_model([1, 2], 1, m)
    @test PeriLab.Data_Manager.get_analysis_model("Correspondence Model", 1) ==
          ["Correspondence Elastic"]
    @test PeriLab.Data_Manager.get_field("Material Gradient")[1, :, :] ≈
          Matrix(BBASIS.hooke_matrix(m, 3, 1))
    ut_reset(3; nnodes = 2)
    nosym = typed_block_material(Dict("Material Model" => "Correspondence Elastic",
                                      "Bulk Modulus" => 10.0, "Shear Modulus" => 4.0))
    @test_throws PeriLab.PeriLabError BCORR.init_model([1, 2], 1, nosym)
    @test hasmethod(BCORR.compute_model, Tuple{Vector{Int64},typeof(m),Int64,Float64,Float64})
    @test hasmethod(BCORR.Bond_Associated_Correspondence.compute_model,
                    Tuple{Vector{Int64},typeof(m),Int64,Float64,Float64})
end
```

If `Data_Manager.set_rotation` does not exist, remove that line: the init does not need it. Ledger it if you remove it.

- [ ] **Step 2: Run to verify it fails**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Models/Material/ut_block_material`
Expected: `MethodError: no method matching init_model(::Vector{Int64}, ::Int64, ::BlockMaterial…)`.

- [ ] **Step 3: Implement**

`Correspondence.jl`: add below the dict methods:

```julia
_model_parts(model) = hasfield(typeof(model), :parts) ? model.parts : (model,)

function init_model(nodes::AbstractVector{Int64}, block::Int64, material)
    if material.base.symmetry === nothing
        @abort "Symmetry for correspondence material is missing; options are 'isotropic plane strain', 'isotropic plane stress', 'anisotropic plane stress', 'anisotropic plane stress','isotropic' and 'anisotropic'. For 3D the plane stress or plane strain option is ignored."
        return
    end
    dof = Data_Manager.get_dof()
    Data_Manager.create_node_tensor_field("Strain", Float64, dof)
    Data_Manager.create_constant_node_tensor_field("Strain Increment", Float64, dof)
    Data_Manager.create_node_tensor_field("Cauchy Stress", Float64, dof)
    Data_Manager.create_node_scalar_field("von Mises Stress", Float64)
    for part in _model_parts(material.model)
        mod = parentmodule(typeof(part))
        Data_Manager.set_analysis_model("Correspondence Model", block,
                                        mod.correspondence_name())
        Data_Manager.set_model_module(mod.correspondence_name(), mod)
        mod.init_model(nodes, part, material)
    end
    if material.base.bond_associated
        return Bond_Associated_Correspondence.init_model(nodes, material)
    end
    Zero_Energy_Control.init_model(nodes, material, block)
end

function fields_for_local_synchronization(model::String, block::Int64, material)
    for material_model in Data_Manager.get_analysis_model("Correspondence Model", block)
        mod = Data_Manager.get_model_module(material_model)
        mod.fields_for_local_synchronization(model)
        if material.base.bond_associated
            Bond_Associated_Correspondence.fields_for_local_synchronization(model)
        end
    end
end

function compute_model(nodes::AbstractVector{Int64}, material, block::Int64, time::Float64,
                       dt::Float64)
    if material.base.bond_associated
        return Bond_Associated_Correspondence.compute_model(nodes, material, block, time, dt)
    end
    return compute_correspondence_model(nodes, material, block, time, dt)
end
```

- Copy the dict `compute_correspondence_model` into a typed method and change only:
  - The signature's `material_parameter::Dict{String,Any}` becomes `material`.
  - `if haskey(material_parameter, "Linear Strain") && material_parameter["Linear Strain"]` becomes `if material.base.linear_strain`.
  - The `@timeit "compute material"` block becomes

```julia
    @timeit "compute material" begin
        for part in _model_parts(material.model)
            parentmodule(typeof(part)).compute_stresses(nodes, dof, part, material, time, dt,
                                                        strain_increment, stress_N,
                                                        stress_NP1)
        end
    end
```

  - The zero energy call passes `material` instead of `material_parameter`.

`Bond_Associated_Correspondence.jl`: add typed methods.
- `init_model(nodes::AbstractVector{Int64}, material)` is a copy of the dict init with these changes:
  - The symmetry check becomes `material.base.symmetry === nothing && @abort "<same message>"`.
  - The accuracy-order lines become `material.base.accuracy_order === nothing || Data_Manager.set_accuracy_order(material.base.accuracy_order)`.
- `compute_model(nodes::AbstractVector{Int64}, material, block::Int64, time::Float64, dt::Float64)` is a copy of the dict compute, where the `material_models = split(…)` block through its loop becomes

```julia
    for part in (hasfield(typeof(material.model), :parts) ? material.model.parts :
                 (material.model,))
        parentmodule(typeof(part)).compute_stresses_ba(nodes, nlist, dof, part, material,
                                                       time, dt, strain_increment, stress_N,
                                                       stress_NP1)
    end
```

`Material_Factory.jl`:
- In `init_model`, replace

```julia
    if occursin("Correspondence", model_param["Material Model"])
        Data_Manager.set_model_module("Correspondence", Correspondence)
        return Correspondence.init_model(nodes, block, model_param)
    end
```

  with

```julia
    if occursin("Correspondence", model_param["Material Model"])
        Data_Manager.set_model_module("Correspondence", Correspondence)
        material = Data_Manager.get_block_material(block)
        if material === nothing
            @abort "Block $block has no typed material parameters."
            return
        end
        bind_material!(material)
        return Correspondence.init_model(nodes, block, material)
    end
```

- In `fields_for_local_synchronization`, the correspondence call becomes `Correspondence.fields_for_local_synchronization(model, block, Data_Manager.get_block_material(block))`.
- In `compute_model`, the correspondence call `Correspondence.compute_model(nodes, model_param, block, time, dt)` becomes `Correspondence.compute_model(nodes, bind_material!(Data_Manager.get_block_material(block)), block, time, dt)`.
- `compute_correspondence_bond_forces` becomes

```julia
function compute_correspondence_bond_forces(nodes::AbstractVector{Int64}, material,
                                            block::Int64, time::Float64, dt::Float64)
    Correspondence.compute_correspondence_model(nodes, bind_material!(material), block,
                                                time, dt)
end
```

`Model_Factory.jl`, `compute_matrix_based_bond_forces`: in the call `Material.compute_correspondence_bond_forces(active_nodes, material_parameter, block, time, dt)`, the argument `material_parameter` becomes `Data_Manager.get_block_material(block)`. The `occursin("Correspondence", material_parameter["Material Model"])` test stays until 3c-2.

`Correspondence_matrix_based.jl:622`: `Zero_Energy_Control.init_model(nodes, material_parameter, block_id)` becomes `Zero_Energy_Control.init_model(nodes, Data_Manager.get_block_material(block_id), block_id)`.

Do not delete any dict method. They are unused after this task and are removed in 3c-2 with the material dict. Ledger this.

- [ ] **Step 4: Run to verify it passes**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Models/Material/ut_block_material unit_tests/Models/Material/ut_Material_Factory unit_tests/Models/ut_Model_Factory`
Expected: all pass.

Then the correspondence fullscale tests: find their runners with `ls test/fullscale_tests/{test_Correspondence_Elastic,test_Correspondence_Elastic_Plastic,test_correspondence_elastic_3D,test_correspondence_elastic_with_zero_E_control,test_3D_aniso_material,test_symmetry,test_Umat,test_DCB}/*.jl`. Run them with the unit runner (`fullscale_tests/<dir>/<file>`) and `JULIA_PROJECT` set.
Expected: pass.

- [ ] **Step 5: Full suite, then commit**

Full suite in background; read the tail. Expected: 0 failures, 0 errors. `ut_MPI.jl` does not use correspondence; if it fails, read the error before touching it.

```bash
git add src/Models/Material src/Models/Model_Factory.jl test/unit_tests/Models/Material
git commit -m "Correspondence blocks run on the typed block material"
```
