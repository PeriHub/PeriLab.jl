# Typed Input Parameters — Phase 3b: Material Runtime on Typed Structs (core, non-correspondence)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Every block gets a typed `BlockMaterial` (base parameters, model struct, normalised symmetry, runtime elastic moduli). The bond-based, PD Solid and Rigid models compute from it instead of the material dict.

**Architecture:**
- `elastic_moduli` (in the Material factory) is a faithful port of `Material_Basis.get_all_elastic_moduli` to `MaterialBaseParams`.
  - It returns `ElasticModuli` whose entries are `Float64`, or node-field vectors for per-node and state-scaled moduli.
  - It is now the only place K/E/G/ν are completed.
  - `read_properties` builds one `BlockMaterial` per block, stores it in `Data_Manager`, and writes the completed moduli back into the block's material dict.
  - The write-back keeps correspondence, the strain compute class, the critical time step and local damping working unchanged until phase 3c.
- Material factory dispatch:
  - Non-correspondence blocks call each model part's module (found with `parentmodule(typeof(part))`) with the typed part and the `BlockMaterial`.
  - Dependent tables are re-bound to the current `NP1` field before every compute call.
- Parse-time completeness: `check!(::MaterialBaseParams)` rejects anisotropic, orthotropic and transverse isotropic materials that lack constants.

**Tech Stack:** Julia 1.12, `PeriLab.ParameterSpec`, `PeriLab.InputDeck`, `PeriLab.Data_Manager`.

**Spec:** `docs/superpowers/specs/2026-10-03-typed-input-parameters-design.md` (§2.6 binding, §2.7, §3 Use/BlockModels)

## Decisions taken with the user for this phase

- **Runtime moduli struct.** Moduli can come from mesh columns (`Bulk_Modulus` …), from state-variable scaling, or be constants. They are therefore computed at initialisation into `ElasticModuli` (`Float64` or node-field vector per entry), not at parse time. `check!` performs only the parse-time completeness checks that do not depend on the mesh.
- **Split:**
  - 3b: core infrastructure plus bond-based, PD Solid and Rigid.
  - 3c: Correspondence Elastic/Plastic/UMAT/VUMAT, bond-associated, zero energy control, and the solvers, compute classes and pre-calculation that read the dict. Then the material dict is removed.
  - Until then the dict stays, filled with the typed moduli.

## Global Constraints

- Numerical results unchanged: the full suite (`Pkg.test()`) must stay green after every task that runs it, including the fullscale decks (`test_Bond_Based_Elastic`, `test_PD_Solid_Elastic`, `test_PD_solid_plastic`, `test_point_wise_material`, `test_symmetry`, …).
- Full YAML backward compatibility: every shipped deck parses in strict mode (golden test).
- Read fields directly; functions only for real logic.
- Structs are immutable; per-node state stays in `Data_Manager` fields.
- A `Table1D` keeps a reference to the array it was bound to, and `switch_NP1_to_N` swaps arrays every step. So bind right before use, in every compute call (spec §2.6).
- Branch `feature/typed-input-parameters`; one commit per task; never merge. Read the full-suite result before committing a task that runs it.

## Review Focus

1. **Per-node moduli from a mesh column (`Bulk_Modulus` in `test_point_wise_material/test_mesh.txt`)** give node vectors identical to the legacy completion, and the fields are written back. Pinned in Task 1 (`moduli from a mesh field`), plus the point-wise fullscale tests.
2. **A bond-based material with an explicit Poisson's ratio** gets the fixed ν (1/3 in 2D, 1/4 in 3D) and its E unchanged. Pinned in Task 1 (`bond-based fixed Poisson's ratio`).
3. **A `Yield Stress` table bound to a field whose N/NP1 arrays swap every step** follows the current NP1 array. Pinned in Task 5 (`table follows the NP1 switch`).
4. **`PD Solid Elastic + PD Solid Plastic`** computes elastic, then plastic, as before. Pinned in Task 5 (`composite parts in order`) and the fullscale test.
5. **An orthotropic material missing `Shear Modulus XZ`** is an input error naming the missing constant, not a runtime abort. Pinned in Task 2 (`orthotropic completeness`).

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
- Fullscale file: the same runner with `fullscale_tests/<dir>/<file>`. Set `JULIA_PROJECT` so MPI subprocesses find PeriLab.
- Input tests: `julia --project=. test/unit_tests/Support/Parameters/Input/run_input_tests.jl 2>&1 | tail -15`
- Full suite (~30 min, background): `julia --project=. -e 'using Pkg; Pkg.test()' > <workspace>/suite.log 2>&1`

`Logging.disable_logging(Logging.Warn)` hides `:warn` records; assert only `:error` logs.

## Module paths

| Module | Path |
|---|---|
| Material factory | `PeriLab.Solver_Manager.Model_Factory.Material` (`src/Models/Material/Material_Factory.jl`) |
| Material_Basis | `PeriLab.Solver_Manager.Material_Basis` |
| model modules | `….Material.Bondbased_Elastic`, `.OneD_Bond_Based_Elastic`, `.Unified_Bondbased_Elastic`, `.PD_Solid_Elastic`, `.PD_Solid_Plastic`, `.Rigid`, `.Material_template` |

Imports of `ParameterSpec` from a model module use `......ParameterSpec` (extra dots stay at the package root).

---

### Task 1: `ElasticModuli`, `BlockMaterial` and the typed moduli completion

**Files:**
- Modify: `src/Models/Material/Material_Factory.jl` (types and functions after the `module_list` include loop, before `export init_model`)
- Modify: `test/helper.jl` (`typed_block_material`)
- Create: `test/unit_tests/Models/Material/ut_block_material.jl`
- Modify: `test/runtests.jl` (include the new file next to `ut_material_basis.jl`)

**Interfaces:**
- Consumes: `MaterialBaseParams` (phase 3a), `ParameterSpec.WithBase`, `Data_Manager.has_key/get_field/create_constant_node_scalar_field/get_dof`.
- Produces (in `Material`):
  - `struct ElasticModuli{K,E,G,N}` with fields `bulk_modulus::K, youngs_modulus::E, shear_modulus::G, poissons_ratio::N`. Each entry is `Float64` or a node-field `Vector{Float64}`.
  - `modulus(x::Real, iID) = x`, `modulus(x::AbstractVector, iID) = x[iID]`.
  - `struct BlockMaterial{B,M,E}` with fields:
    - `base::B` (`MaterialBaseParams`);
    - `model::M` (model struct or `Composite`);
    - `symmetry::String` (`"plane strain"`, `"plane stress"` or `"3D"`);
    - `moduli::E` (`ElasticModuli` or `nothing` for anisotropic, orthotropic and transverse isotropic).
  - `material_symmetry(symmetry::Union{Nothing,String}, dof::Int64)::String`.
  - `elastic_moduli(base::MaterialBaseParams, bond_based::Bool, dof::Int64)`.
  - `block_material(wb::WithBase, model_name::String, dof::Int64)::BlockMaterial`.
  - `write_moduli!(dict::Dict{String,Any}, material::BlockMaterial)`.
  - `model_parts(model)`: the tuple of model structs.
  - Test helper `typed_block_material(raw::AbstractDict; dof::Int64 = PeriLab.Data_Manager.get_dof())::BlockMaterial`.

- [ ] **Step 1: Add the test helper** (append to `test/helper.jl`)

```julia
"""
    typed_block_material(raw; dof)

A `BlockMaterial` from a raw material block (with `"Material Model"`), parsed the
way the input reader does; aborts on input errors.
"""
function typed_block_material(raw::AbstractDict;
                              dof::Int64 = PeriLab.Data_Manager.get_dof())
    ctx = PeriLab.ParameterSpec.ParseContext()
    dict = Dict{String,Any}(string(k) => v for (k, v) in raw)
    wb = PeriLab.ParameterSpec.parse_model(:material, dict, "test", ctx;
                                           name_key = "Material Model")
    PeriLab.ParameterSpec.report!(ctx)
    return PeriLab.Solver_Manager.Model_Factory.Material.block_material(wb,
                                                                        String(dict["Material Model"]),
                                                                        dof)
end
```

- [ ] **Step 2: Write the failing tests**

Create `test/unit_tests/Models/Material/ut_block_material.jl`:

```julia
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

const BMAT = PeriLab.Solver_Manager.Model_Factory.Material
const BBASIS = PeriLab.Solver_Manager.Material_Basis

function ut_reset(dof; nnodes = 3)
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(nnodes)
    PeriLab.Data_Manager.set_dof(dof)
end

# legacy completion and typed completion on the same raw block, each on a fresh Data_Manager
function ut_both(raw; dof = 3, setup = () -> nothing)
    ut_reset(dof)
    setup()
    legacy = Dict{String,Any}(raw)
    BBASIS.get_all_elastic_moduli(legacy)
    legacy_values = Dict(k => copy(legacy[k])
                         for k in ("Bulk Modulus", "Young's Modulus", "Shear Modulus",
                                   "Poisson's Ratio"))
    ut_reset(dof)
    setup()
    typed = typed_block_material(raw; dof = dof)
    return legacy_values, typed
end

function ut_same(legacy, m::BMAT.ElasticModuli)
    return isapprox(legacy["Bulk Modulus"], m.bulk_modulus) &&
           isapprox(legacy["Young's Modulus"], m.youngs_modulus) &&
           isapprox(legacy["Shear Modulus"], m.shear_modulus) &&
           isapprox(legacy["Poisson's Ratio"], m.poissons_ratio)
end

@testset "isotropic completion matches the legacy completion" begin
    for raw in [Dict("Material Model" => "PD Solid Elastic", "Bulk Modulus" => 10.0,
                     "Shear Modulus" => 10.0),
                Dict("Material Model" => "PD Solid Elastic", "Bulk Modulus" => 5.0,
                     "Young's Modulus" => 1.25),
                Dict("Material Model" => "PD Solid Elastic", "Poisson's Ratio" => 0.45,
                     "Shear Modulus" => 1.25),
                Dict("Material Model" => "PD Solid Elastic", "Young's Modulus" => 5.0,
                     "Poisson's Ratio" => 0.125, "Symmetry" => "isotropic"),
                Dict("Material Model" => "PD Solid Elastic", "Bulk Modulus" => 1.0,
                     "Shear Modulus" => 10.0, "Poisson's Ratio" => 0.2),
                Dict("Material Model" => "Unified Bond-based Elastic",
                     "Young's Modulus" => 5.0, "Poisson's Ratio" => 0.125,
                     "Symmetry" => "isotropic plane strain")]
        legacy, typed = ut_both(raw)
        @test typed isa BMAT.BlockMaterial
        @test ut_same(legacy, typed.moduli)
    end
end

@testset "bond-based fixed Poisson's ratio" begin
    for (dof, nu) in ((2, 1 / 3), (3, 1 / 4))
        legacy, typed = ut_both(Dict("Material Model" => "Bond-based Elastic",
                                     "Young's Modulus" => 5.0, "Poisson's Ratio" => 0.125,
                                     "Symmetry" => dof == 2 ? "isotropic plane stress" :
                                                   "isotropic"); dof = dof)
        @test typed.moduli.poissons_ratio == nu
        @test typed.moduli.youngs_modulus == 5.0
        @test ut_same(legacy, typed.moduli)
    end
end

@testset "moduli from a mesh field" begin
    setup = () -> PeriLab.Data_Manager.create_constant_node_scalar_field("Bulk_Modulus",
                                                                        Float64;
                                                                        default_value = 10)
    legacy, typed = ut_both(Dict("Material Model" => "PD Solid Elastic",
                                 "Shear Modulus" => 10.0); setup = setup)
    @test typed.moduli.youngs_modulus == [22.5, 22.5, 22.5]
    @test ut_same(legacy, typed.moduli)
    @test PeriLab.Data_Manager.get_field("Young's_Modulus") == [22.5, 22.5, 22.5]
    @test BMAT.modulus(typed.moduli.youngs_modulus, 2) == 22.5
    @test BMAT.modulus(5.0, 2) == 5.0
end

@testset "Hooke-matrix symmetries have no isotropic moduli" begin
    ut_reset(3)
    m = typed_block_material(Dict("Material Model" => "PD Solid Elastic",
                                  "Symmetry" => "orthotropic",
                                  "Young's Modulus X" => 1.0, "Young's Modulus Y" => 1.0,
                                  "Young's Modulus Z" => 1.0, "Poisson's Ratio XY" => 0.3,
                                  "Poisson's Ratio YZ" => 0.3, "Poisson's Ratio XZ" => 0.3,
                                  "Shear Modulus XY" => 1.0, "Shear Modulus YZ" => 1.0,
                                  "Shear Modulus XZ" => 1.0))
    @test m.moduli === nothing
end

@testset "too few isotropic constants" begin
    ut_reset(3)
    @test_logs (:error,
                "Minimum of two parameters are needed for isotropic material") @test_throws PeriLab.PeriLabError begin
        typed_block_material(Dict("Material Model" => "PD Solid Elastic",
                                  "Bulk Modulus" => 10.0))
    end
end

@testset "material_symmetry follows check_symmetry and get_symmetry" begin
    for (sym, dof) in (("isotropic plane strain", 2), ("isotropic plane stress", 2),
                       ("isotropic plane strain", 3), ("Isotropic Plane Stress", 2),
                       ("isotropic", 3), (nothing, 2), (nothing, 3),
                       ("iso Plane Strain", 3))
        legacy = sym === nothing ? Dict{String,Any}() : Dict{String,Any}("Symmetry" => sym)
        if dof == 3 && haskey(legacy, "Symmetry")
            legacy["Symmetry"] = replace(replace(legacy["Symmetry"], r"plane strain$" => ""),
                                         r"plane stress$" => "")
        end
        @test BMAT.material_symmetry(sym, dof) == BBASIS.get_symmetry(legacy)
    end
end

@testset "write_moduli! fills the material dict" begin
    ut_reset(3)
    m = typed_block_material(Dict("Material Model" => "PD Solid Elastic",
                                  "Bulk Modulus" => 10.0, "Shear Modulus" => 10.0))
    dict = Dict{String,Any}("Material Model" => "PD Solid Elastic", "Bulk Modulus" => 10.0,
                            "Shear Modulus" => 10.0)
    BMAT.write_moduli!(dict, m)
    @test dict["Young's Modulus"] == 22.5 && dict["Poisson's Ratio"] == 0.125
    @test dict["Computed"] === true
    @test dict["Symmetry"] == "isotropic"
    @test BMAT.model_parts(m.model) == (m.model,)
end
```

In `test/runtests.jl`, next to the line that includes `unit_tests/Models/Material/ut_material_basis.jl`, add a testset including `unit_tests/Models/Material/ut_block_material.jl`, in the same style as its neighbours.

- [ ] **Step 3: Run to verify they fail**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Models/Material/ut_block_material`
Expected: errors such as `UndefVarError: block_material` / `material_symmetry` / `ElasticModuli` not defined in `Material`.

- [ ] **Step 4: Implement** in `src/Models/Material/Material_Factory.jl`

Extend the ParameterSpec import to `using .....ParameterSpec: @params, Dependent, register_base!, WithBase, Composite`. After the `using ...Material_Basis: …` block, add:

```julia
"""
    ElasticModuli

Completed isotropic elastic constants of a block. Each entry is a `Float64`, or
a node field (`Vector{Float64}`) when moduli are given per node (mesh columns
`Bulk_Modulus`, …) or scaled by a state variable. Read one with `modulus(x, iID)`.
"""
struct ElasticModuli{K,E,G,N}
    bulk_modulus::K
    youngs_modulus::E
    shear_modulus::G
    poissons_ratio::N
end

@inline modulus(x::Real, ::Int64) = x
@inline modulus(x::AbstractVector, iID::Int64) = x[iID]

"""
    BlockMaterial

Typed material of one block: the shared base parameters, the model struct (or
`Composite`), the symmetry used by the force models (`"plane strain"`,
`"plane stress"` or `"3D"`), and the completed elastic moduli (`nothing` for
materials defined by a stiffness matrix).
"""
struct BlockMaterial{B,M,E}
    base::B
    model::M
    symmetry::String
    moduli::E
end

model_parts(model::Composite) = model.parts
model_parts(model) = (model,)

"""
    material_symmetry(symmetry, dof)

The symmetry the force models use. Plane strain / plane stress are ignored in 3D
(as `check_symmetry` does), and the rest follows `get_symmetry`.
"""
function material_symmetry(symmetry::Union{Nothing,String}, dof::Int64)
    symmetry === nothing && return "3D"
    s = symmetry
    if dof == 3
        s = replace(replace(s, r"plane strain$" => ""), r"plane stress$" => "")
    end
    s = lowercase(s)
    occursin("plane strain", s) && return "plane strain"
    occursin("plane stress", s) && return "plane stress"
    return "3D"
end

const _MODULI = (("Bulk Modulus", :bulk_modulus), ("Young's Modulus", :youngs_modulus),
                 ("Shear Modulus", :shear_modulus), ("Poisson's Ratio", :poissons_ratio))

_modulus_field(key::String) = replace(key, " " => "_")

# value of one modulus: the node field if it exists, a new node field if any
# modulus is per node, otherwise the given constant (0.0 if not given)
function _modulus_value(given, field_allocated::Bool, any_field_allocated::Bool,
                        key::String)
    field_allocated && return Data_Manager.get_field(_modulus_field(key))
    if any_field_allocated
        return given === nothing ?
               Data_Manager.create_constant_node_scalar_field(_modulus_field(key), Float64) :
               Data_Manager.create_constant_node_scalar_field(_modulus_field(key), Float64;
                                                              default_value = given)
    end
    return given === nothing ? 0.0 : given
end

"""
    elastic_moduli(base, bond_based, dof)

Completes the isotropic elastic constants from any two of bulk modulus, Young's
modulus, shear modulus and Poisson's ratio (bond-based models: Poisson's ratio
is fixed). Moduli given per node (fields `Bulk_Modulus`, …) or a `State Factor
ID` make the result per node, and the node fields are updated. Returns `nothing`
for anisotropic, orthotropic and transverse isotropic materials (stiffness
matrix; completeness is checked when the input is read).
"""
function elastic_moduli(base::MaterialBaseParams, bond_based::Bool, dof::Int64)
    state_factor_defined = base.state_factor_id !== nothing
    allocated = Dict(key => Data_Manager.has_key(_modulus_field(key)) for (key, _) in _MODULI)
    any_field_allocated = any(values(allocated)) || state_factor_defined
    given = Dict(key => getfield(base, name) for (key, name) in _MODULI)
    has = Dict(key => given[key] !== nothing || allocated[key] for (key, _) in _MODULI)

    K = _modulus_value(given["Bulk Modulus"], allocated["Bulk Modulus"],
                       any_field_allocated, "Bulk Modulus")
    E = _modulus_value(given["Young's Modulus"], allocated["Young's Modulus"],
                       any_field_allocated, "Young's Modulus")
    G = _modulus_value(given["Shear Modulus"], allocated["Shear Modulus"],
                       any_field_allocated, "Shear Modulus")
    nu = _modulus_value(given["Poisson's Ratio"], allocated["Poisson's Ratio"],
                        any_field_allocated, "Poisson's Ratio")
    bulk, youngs, shear, poissons = has["Bulk Modulus"], has["Young's Modulus"],
                                    has["Shear Modulus"], has["Poisson's Ratio"]

    if bond_based
        nu_fixed = dof == 2 ? 1 / 3 : 1 / 4
        if nu != 0.0 && nu != nu_fixed
            @warn "Chosen Bond-based model only supports a fixed Poisson's ratio of " *
                  string(nu_fixed)
        end
        nu = nu_fixed
        poissons = true
    end
    if base.symmetry !== nothing
        symmetry = lowercase(base.symmetry)
        if occursin("anisotropic", symmetry) || occursin("transverse isotropic", symmetry) ||
           occursin("orthotropic", symmetry)
            return nothing
        end
    else
        @warn "Material symmetry is not defined, assuming isotropic material"
    end

    if bulk + youngs + shear + poissons < 2
        @abort "Minimum of two parameters are needed for isotropic material"
    elseif bulk + youngs + shear + poissons > 2
        @warn "Only two parameters are needed for isotropic material, ignoring additional parameters"
    end

    if bulk && poissons
        E = 3 .* K .* (1 .- 2 .* nu)
        G = 3 .* K .* (1 .- 2 .* nu) ./ (2 .+ 2 .* nu)
    end
    if shear && poissons
        E = 2 .* G .* (1 .+ nu)
        K = 2 .* G .* (1 .+ nu) ./ (3 .- 6 .* nu)
    end
    if bulk && shear
        E = 9 .* K .* G ./ (3 .* K .+ G)
        nu = (3 .* K .- 2 .* G) ./ (6 .* K .+ 2 .* G)
    end
    if youngs && shear
        K = E .* G ./ (9 .* G .- 3 .* E)
        nu = E ./ (2 .* G) .- 1
    end
    if youngs && bulk
        G = 3 .* K .* E ./ (9 .* K .- E)
        nu = (3 .* K .- E) ./ (6 .* K)
    end
    if youngs && poissons
        K = E ./ (3 .- 6 .* nu)
        G = E ./ (2 .+ 2 .* nu)
    end

    if state_factor_defined && Data_Manager.has_key("State Variables")
        state_factor = Data_Manager.get_field("State Variables")[:, base.state_factor_id]
        K = K .* state_factor
        E = E .* state_factor
        G = G .* state_factor
    end
    if any_field_allocated
        Data_Manager.get_field("Bulk_Modulus") .= K
        Data_Manager.get_field("Young's_Modulus") .= E
        Data_Manager.get_field("Shear_Modulus") .= G
        Data_Manager.get_field("Poisson's_Ratio") .= nu
    end
    return ElasticModuli(K, E, G, nu)
end
```

```julia
"""
    block_material(wb, model_name, dof)

The `BlockMaterial` of a parsed material block. `model_name` is the block's
`Material Model` string (bond-based models fix Poisson's ratio).
"""
function block_material(wb::WithBase, model_name::String, dof::Int64)
    bond_based = occursin("Bond-based", model_name) &&
                 !occursin("Unified Bond-based", model_name)
    return BlockMaterial(wb.base, wb.model, material_symmetry(wb.base.symmetry, dof),
                         elastic_moduli(wb.base, bond_based, dof))
end

"""
    write_moduli!(dict, material)

Writes the completed moduli into a block's material dict, for the code that
still reads the dict (correspondence, compute classes; phase 3c).
"""
function write_moduli!(dict::Dict{String,Any}, material::BlockMaterial)
    material.moduli === nothing && return dict
    dict["Bulk Modulus"] = material.moduli.bulk_modulus
    dict["Young's Modulus"] = material.moduli.youngs_modulus
    dict["Shear Modulus"] = material.moduli.shear_modulus
    dict["Poisson's Ratio"] = material.moduli.poissons_ratio
    dict["Computed"] = true
    haskey(dict, "Symmetry") || (dict["Symmetry"] = "isotropic")
    return dict
end
```

In the legacy function, the symmetry check (with its early `return` and the "assuming isotropic" warning) sits *between* the bond-based fix and the two-parameter check. The code above keeps that order: bond-based fix, then symmetry, then counts.

- [ ] **Step 5: Run to verify they pass**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Models/Material/ut_block_material unit_tests/Models/Material/ut_material_basis`
Expected: all pass.

- [ ] **Step 6: Commit**

```bash
git add src/Models/Material/Material_Factory.jl test/helper.jl \
        test/unit_tests/Models/Material/ut_block_material.jl test/runtests.jl
git commit -m "Typed block material: elastic moduli, symmetry, dict write-back"
```

---

### Task 2: Parse-time completeness of stiffness-matrix materials

**Files:**
- Modify: `src/Models/Material/Material_Factory.jl` (`check!(::MaterialBaseParams, …)`)
- Test: `test/unit_tests/Support/Parameters/Input/ut_material_params.jl` (append)

**Interfaces:**
- Consumes: `ParameterSpec.check!`, `add_error!`, `join_path`, `ParseContext`.
- Produces: input errors at `<material path>.Symmetry`:
  - anisotropic: `"\"<symmetry>\" requires C12, …"`;
  - orthotropic and transverse isotropic: `"\"<symmetry>\" requires Shear Modulus XZ, …"`.
  - Requirements are copied from the legacy runtime checks in `get_all_elastic_moduli`.

- [ ] **Step 1: Write the failing tests** (append to `ut_material_params.jl`)

```julia
@testset "orthotropic completeness" begin
    full = Dict{String,Any}("Material Model" => "PD Solid Elastic", "Symmetry" => "Orthotropic",
                            "Young's Modulus X" => 1.0, "Young's Modulus Y" => 1.0,
                            "Young's Modulus Z" => 1.0, "Poisson's Ratio XY" => 0.3,
                            "Poisson's Ratio YZ" => 0.3, "Poisson's Ratio XZ" => 0.3,
                            "Shear Modulus XY" => 1.0, "Shear Modulus YZ" => 1.0,
                            "Shear Modulus XZ" => 1.0)
    m, ctx = ut_material(full)
    @test isempty(ctx.errors)
    delete!(full, "Shear Modulus XZ")
    m, ctx = ut_material(full)
    @test only(ctx.errors).path == "Models.\"Material Models\".M.Symmetry"
    @test only(ctx.errors).message == "\"Orthotropic\" requires Shear Modulus XZ"
end

@testset "anisotropic and transverse isotropic completeness" begin
    m, ctx = ut_material(Dict("Material Model" => "PD Solid Elastic",
                              "Symmetry" => "anisotropic", "C11" => 1.0))
    @test startswith(only(ctx.errors).message, "\"anisotropic\" requires C12, C13")
    m, ctx = ut_material(Dict("Material Model" => "PD Solid Elastic",
                              "Symmetry" => "transverse isotropic plane stress",
                              "Young's Modulus X" => 1.0, "Young's Modulus Y" => 1.0,
                              "Poisson's Ratio XY" => 0.3))
    @test only(ctx.errors).message ==
          "\"transverse isotropic plane stress\" requires Shear Modulus XY"
    m, ctx = ut_material(Dict("Material Model" => "PD Solid Elastic",
                              "Symmetry" => "isotropic plane strain",
                              "Bulk Modulus" => 1.0, "Shear Modulus" => 1.0))
    @test isempty(ctx.errors)
end
```

- [ ] **Step 2: Run to verify they fail**

Run: input tests
Expected: the incomplete cases report no errors, so the `only(ctx.errors)` assertions fail.

- [ ] **Step 3: Implement** in `src/Models/Material/Material_Factory.jl`

Change the ParameterSpec import line to

```julia
using .....ParameterSpec: @params, Dependent, register_base!, WithBase, Composite,
                          ParseContext, add_error!, join_path
import .....ParameterSpec: check!
```

After the `MaterialBaseParams` struct, add:

```julia
const _ORTHOTROPIC_KEYS = ((:youngs_modulus_x, "Young's Modulus X"),
                           (:youngs_modulus_y, "Young's Modulus Y"),
                           (:youngs_modulus_z, "Young's Modulus Z"),
                           (:poissons_ratio_xy, "Poisson's Ratio XY"),
                           (:poissons_ratio_yz, "Poisson's Ratio YZ"),
                           (:poissons_ratio_xz, "Poisson's Ratio XZ"),
                           (:shear_modulus_xy, "Shear Modulus XY"),
                           (:shear_modulus_yz, "Shear Modulus YZ"),
                           (:shear_modulus_xz, "Shear Modulus XZ"))

function _transverse_keys(symmetry::String)
    keys = [(:youngs_modulus_x, "Young's Modulus X"), (:youngs_modulus_y, "Young's Modulus Y"),
            (:poissons_ratio_xy, "Poisson's Ratio XY")]
    occursin("plane stress", symmetry) || push!(keys, (:poissons_ratio_yz, "Poisson's Ratio YZ"))
    push!(keys, (:shear_modulus_xy, "Shear Modulus XY"))
    if !occursin("plane strain", symmetry) && !occursin("plane stress", symmetry)
        push!(keys, (:shear_modulus_yz, "Shear Modulus YZ"))
    end
    return keys
end

# stiffness-matrix materials must define all their constants
function check!(p::MaterialBaseParams, path::String, ctx::ParseContext)
    p.symmetry === nothing && return nothing
    symmetry = lowercase(p.symmetry)
    required = if occursin("anisotropic", symmetry)
        [(Symbol("c$i$j"), "C$i$j") for i in 1:6 for j in i:6]
    elseif occursin("transverse isotropic", symmetry)
        _transverse_keys(symmetry)
    elseif occursin("orthotropic", symmetry)
        collect(_ORTHOTROPIC_KEYS)
    else
        Tuple{Symbol,String}[]
    end
    missing_keys = [alias for (field, alias) in required if getfield(p, field) === nothing]
    isempty(missing_keys) ||
        add_error!(ctx, join_path(path, "Symmetry"),
                   "\"$(p.symmetry)\" requires $(join(missing_keys, ", "))")
    return nothing
end
```

- [ ] **Step 4: Run to verify they pass**

Run: input tests
Expected: all pass, including the golden decks. If a shipped deck fails the new check, the legacy runtime would have aborted on it too. Check whether a fullscale test runs it:
- If one does, the check is wrong: fix it and ledger the finding.
- If none does, leave the check as is. Then fix the deck's material (if the missing constant is obvious from the deck) or add the deck to `UT_DECK_ALLOWLIST` with the reason "incomplete orthotropic constants; not run by any test". Ledger which you chose.

- [ ] **Step 5: Commit**

```bash
git add src/Models/Material/Material_Factory.jl test/unit_tests/Support/Parameters/Input/ut_material_params.jl
git commit -m "Stiffness-matrix materials must define all constants (input error)"
```

---

### Task 3: `read_properties` builds the block materials

**Files:**
- Modify: `src/Core/Data_manager.jl` (`Block Materials` slot plus `set_block_material`/`get_block_material`, next to `set_fem_params`)
- Modify: `src/Models/Model_Factory.jl` (`read_properties(params, input, material_model)`; import `block_by_id` is not needed, the block name list is used)
- Modify: `src/Core/Solver/Solver_manager.jl:124`
- Test: `test/unit_tests/Models/ut_Model_Factory.jl` (the `read_properties` testset), plus a new testset

**Interfaces:**
- Consumes (Task 1): `Material.block_material`, `Material.write_moduli!`, `Material.check_material_symmetry`.
- Produces:
  - `Data_Manager.set_block_material(block::Int64, material)` / `Data_Manager.get_block_material(block::Int64)`: a `BlockMaterial` or `nothing`.
  - `Model_Factory.read_properties(params::Dict, input::PeriLabInput, material_model::Bool)`.
  - For each block whose `Material Model` names an entry of `input.materials`, it builds the `BlockMaterial`, stores it, and writes the moduli into the block's material dict.
  - A block without a typed material falls back to the legacy completion. The input reader guarantees typed materials, so this happens only in tests that pass raw dicts.

- [ ] **Step 1: Write the failing test** (append to `test/unit_tests/Models/ut_Model_Factory.jl`)

```julia
@testset "read_properties builds typed block materials" begin
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(3)
    PeriLab.Data_Manager.set_dof(2)
    input = typed_input(Dict("Blocks" => Dict("block_1" => Dict("Block ID" => 1,
                                                                 "Density" => 1.0,
                                                                 "Horizon" => 1.0,
                                                                 "Material Model" => "Mat")),
                             "Models" => Dict("Material Models" => Dict("Mat" => Dict("Material Model" => "PD Solid Elastic",
                                                                                       "Symmetry" => "isotropic plane strain",
                                                                                       "Bulk Modulus" => 10.0,
                                                                                       "Shear Modulus" => 10.0)))))
    params = Dict{String,Any}("Blocks" => Dict{String,Any}("block_1" => Dict{String,Any}("Block ID" => 1,
                                                                                         "Material Model" => "Mat")),
                              "Models" => input.models)
    PeriLab.Data_Manager.set_block_name_list(["block_1"])
    PeriLab.Data_Manager.set_block_id_list([1])
    PeriLab.Solver_Manager.Model_Factory.read_properties(params, input, true)
    m = PeriLab.Data_Manager.get_block_material(1)
    @test m isa PeriLab.Solver_Manager.Model_Factory.Material.BlockMaterial
    @test m.symmetry == "plane strain"
    @test m.moduli.youngs_modulus == 22.5
    @test PeriLab.Data_Manager.get_property(1, "Material Model", "Young's Modulus") == 22.5
    @test PeriLab.Data_Manager.get_property(1, "Material Model", "Computed") === true
end
```

In the existing `read_properties` testset of `ut_Model_Factory.jl`, change both calls `read_properties(params, false)` / `read_properties(params, true)` to `read_properties(params, typed_input(Dict()), false)` / `read_properties(params, typed_input(Dict()), true)`. `typed_input(Dict())` has no material named "a" or "c", so those blocks use the legacy fallback as before.

- [ ] **Step 2: Run to verify it fails**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Models/ut_Model_Factory`
Expected: `MethodError: no method matching read_properties(::Dict…, ::PeriLabInput, ::Bool)`.

- [ ] **Step 3: Implement**

`src/Core/Data_manager.jl`:
- In `initialize_data`, next to `data["FEM Parameters"] = nothing`, add `data["Block Materials"] = Dict{Int64,Any}()`.
- Next to `set_fem_params`, add:

```julia
"""
	set_block_material(block, material)

Stores the typed material (`BlockMaterial`) of a block.
"""
function set_block_material(block::Int64, material)
    data["Block Materials"][block] = material
end

"""
	get_block_material(block)

The typed material of a block, or `nothing`.
"""
function get_block_material(block::Int64)
    return get(data["Block Materials"], block, nothing)
end
```

- Export both next to `export get_fem_params`.

`src/Models/Model_Factory.jl`, `read_properties`:
- In the docstring, change the signature to `read_properties(params::Dict, input::PeriLabInput, material_model::Bool)`; it already documents `input`.
- Replace the function with:

```julia
function read_properties(params::Dict, input::PeriLabInput, material_model::Bool)
    Data_Manager.init_properties()
    block_name_list = Data_Manager.get_block_name_list()
    block_id_list = Data_Manager.get_block_id_list()
    prop_keys = Data_Manager.init_properties()
    directory = Data_Manager.get_directory()
    get_block_model_definition(params,
                               block_name_list,
                               block_id_list,
                               prop_keys,
                               Data_Manager.set_properties,
                               directory,
                               material_model)
    if material_model
        dof = Data_Manager.get_dof()
        for (block_name, block) in zip(block_name_list, block_id_list)
            Material.check_material_symmetry(block)
            properties = Data_Manager.get_properties(block, "Material Model")
            block_params = get(input.sections.blocks, block_name, nothing)
            material_name = block_params === nothing ? nothing : block_params.material_model
            if material_name !== nothing && haskey(input.materials, material_name)
                material = Material.block_material(input.materials[material_name],
                                                   String(properties["Material Model"]), dof)
                Data_Manager.set_block_material(block, material)
                Material.write_moduli!(properties, material)
            else
                # raw dicts without typed materials (unit tests)
                Material.determine_isotropic_parameter(properties)
            end
        end
    end
end
```

`Data_Manager.get_properties` returns the stored dict itself (a `convert` to `Dict{String,Any}` of a `Dict{String,Any}` is a no-op), so `write_moduli!` updates the stored properties, as the legacy `determine_isotropic_parameter` did.

`src/Core/Solver/Solver_manager.jl:124`: `read_properties(params, "Material" in solver_options["Models"])` becomes `read_properties(params, input, "Material" in solver_options["Models"])`.

- [ ] **Step 4: Run to verify it passes**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Models/ut_Model_Factory unit_tests/Models/Material/ut_block_material`
Expected: all pass.

- [ ] **Step 5: Full suite, then commit**

Full suite in background; read the tail. Expected: 0 failures, 0 errors. This proves the typed completion plus write-back reproduces every fullscale result, including point-wise and symmetry decks.

```bash
git add src/Core/Data_manager.jl src/Models/Model_Factory.jl src/Core/Solver/Solver_manager.jl \
        test/unit_tests/Models/ut_Model_Factory.jl
git commit -m "Block materials are built from the typed input; moduli written back to the dict"
```

---

### Task 4: Material factory dispatch; bond-based and Rigid models on `BlockMaterial`

**Files:**
- Modify: `src/Models/Material/Material_Factory.jl` (`init_model`, `compute_model`, new `compute_block_material`, `bind_material!`)
- Modify: `src/Models/Material/Material_Models/BondBased/Bondbased_Elastic.jl`, `1D_Bondbased_Elastic.jl`, `Unified_Bondbased_Elastic.jl`
- Modify: `src/Models/Material/Material_Models/Rigid/Rigid.jl`, `src/Models/Material/Material_Models/Material_template/material_template.jl`
- Test:
  - `test/unit_tests/Models/Material/Material_Models/BondBased/ut_Bondbased_Elastic.jl`, `ut_1D_Bondbased_Elastic.jl`, `ut_Unified_Bondbased_Elastic.jl`
  - `test/unit_tests/Models/Material/ut_block_material.jl` (append)

**Interfaces:**
- Consumes (Tasks 1, 3): `BlockMaterial`, `model_parts`, `Data_Manager.get_block_material`.
- Produces: the model module interface for non-correspondence materials:
  - `init_model(nodes::AbstractVector{Int64}, p::<ModelParams>, material)`.
  - `compute_model(nodes::AbstractVector{Int64}, p::<ModelParams>, material, block::Int64, time::Float64, dt::Float64)`.
  - `material` is the block's `BlockMaterial` (left untyped in the modules, because the type is defined in the parent module after they are included).
  - The factory resolves a part's module with `parentmodule(typeof(part))` and registers it under `mod.material_name()`, so `fields_for_local_synchronization` keeps working by name.
  - `Material.bind_material!(material)` binds every `Table1D` of the material to the current node field (NP1 if it has states). Errors abort with the parameter path.

- [ ] **Step 1: Write the failing tests**

Append to `test/unit_tests/Models/Material/ut_block_material.jl`:

```julia
@testset "factory dispatches typed parts to their modules" begin
    ut_reset(3; nnodes = 2)
    m = typed_block_material(Dict("Material Model" => "Bond-based Elastic",
                                  "Young's Modulus" => 1.0))
    @test parentmodule(typeof(m.model)) === BMAT.Bondbased_Elastic
    @test hasmethod(BMAT.Bondbased_Elastic.compute_model,
                    Tuple{Vector{Int64},typeof(m.model),typeof(m),Int64,Float64,Float64})
    @test hasmethod(BMAT.Rigid.init_model,
                    Tuple{Vector{Int64},BMAT.Rigid.RigidParams,Any})
end
```

In `ut_Bondbased_Elastic.jl`, replace the two init/compute pairs:

```julia
    material = typed_block_material(Dict("Material Model" => "Bond-based Elastic",
                                         "Bulk Modulus" => 1.0, "Young's Modulus" => 1.0))
    PeriLab.Solver_Manager.Model_Factory.Material.Bondbased_Elastic.init_model(Vector{Int64}(1:nodes),
                                                                               material.model,
                                                                               material)
    PeriLab.Solver_Manager.Model_Factory.Material.Bondbased_Elastic.compute_model(Vector{Int64}(1:nodes),
                                                                                  material.model,
                                                                                  material,
                                                                                  1,
                                                                                  0.0,
                                                                                  0.0)
```

and for the second pair the same with `"Symmetry" => "here is something"` added to the dict. Keep all the `isapprox` assertions unchanged.

In `ut_1D_Bondbased_Elastic.jl`, the `init_model` call becomes

```julia
    material = typed_block_material(Dict("Material Model" => "1D Bond-based Elastic",
                                         "Bulk Modulus" => 1.0, "Young's Modulus" => 1.0,
                                         "Id1" => 1, "Id2" => 2))
    PeriLab.Solver_Manager.Model_Factory.Material.OneD_Bond_Based_Elastic.init_model(Vector{Int64}(1:nodes),
                                                                                     material.model,
                                                                                     material)
```

Leave its commented-out compute call commented.

In `ut_Unified_Bondbased_Elastic.jl`, run `grep -n 'init_model\|compute_model' ` on the file. Every call taking a `Dict` gets the same conversion. Build `typed_block_material(Dict("Material Model" => "Unified Bond-based Elastic", <the dict's moduli and Symmetry>))` and pass `material.model, material`, plus `1, time, dt` for compute. If the file calls neither function (it may only test helpers), leave it unchanged.

- [ ] **Step 2: Run to verify they fail**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Models/Material/ut_block_material unit_tests/Models/Material/Material_Models/BondBased/ut_Bondbased_Elastic unit_tests/Models/Material/Material_Models/BondBased/ut_1D_Bondbased_Elastic`
Expected:
- `hasmethod` assertions fail.
- The model tests error with `MethodError: no method matching init_model(::Vector{Int64}, ::BondbasedElasticParams, ::BlockMaterial…)`.

- [ ] **Step 3: Implement**

`src/Models/Material/Material_Factory.jl`:
- Extend the ParameterSpec import with `bind_dependents!, report!`.
- Replace `init_model(nodes, block)` from `material_models = split(...)` through the end of its `for material_model in material_models … end` loop with:

```julia
    material = Data_Manager.get_block_material(block)
    if material === nothing
        @abort "Block $block has no typed material parameters."
        return
    end
    for part in model_parts(material.model)
        mod = parentmodule(typeof(part))
        Data_Manager.set_analysis_model("Material Model", block, mod.material_name())
        Data_Manager.set_model_module(mod.material_name(), mod)
        mod.init_model(nodes, part, material)
    end
```

  Keep everything before it (the `haskey` abort and the `Correspondence` branch) and after it (the filtered-nlist bond-norm block) unchanged.
- In `compute_model(nodes, model_param, block, time, dt)`, replace the `for material_model in Data_Manager.get_analysis_model("Material Model", block) … end` loop with

```julia
        @timeit "material" compute_block_material(nodes,
                                                  Data_Manager.get_block_material(block),
                                                  block, time, dt)
```

- Add after `compute_model`:

```julia
# node field a dependent table reads: the NP1 state if the field has states
function _dependent_field(name::String)
    Data_Manager.has_key(name * "NP1") && return Data_Manager.get_field(name, "NP1")
    Data_Manager.has_key(name) && return Data_Manager.get_field(name)
    return nothing
end

"""
    bind_material!(material)

Binds the dependent tables of a block material to the current node fields. Call
it before every evaluation, because the N/NP1 field arrays are swapped every step.
"""
function bind_material!(material::BlockMaterial)
    ctx = ParseContext()
    bind_dependents!(material.base, _dependent_field, "Material", ctx)
    bind_dependents!(material.model, _dependent_field, "Material", ctx)
    report!(ctx)
    return material
end

# function barrier: `material` has a concrete type here
function compute_block_material(nodes::AbstractVector{Int64}, material::BlockMaterial,
                                block::Int64, time::Float64, dt::Float64)
    bind_material!(material)
    for part in model_parts(material.model)
        parentmodule(typeof(part)).compute_model(nodes, part, material, block, time, dt)
    end
    return nothing
end
```

Model modules (only signatures and parameter reads change; the numerics stay):

`Bondbased_Elastic.jl`:
- `init_model(nodes::AbstractVector{Int64}, material_parameter::Dict)` becomes `init_model(nodes::AbstractVector{Int64}, p::BondbasedElasticParams, material)`.
- `symmetry::String = get_symmetry(material_parameter)` becomes `symmetry::String = material.symmetry`.
- In `compute_model`:
  - The signature becomes `compute_model(nodes::AbstractVector{Int64}, p::BondbasedElasticParams, material, block::Int64, time::Float64, dt::Float64)`.
  - `E = material_parameter["Young's Modulus"]` becomes `E = material.moduli.youngs_modulus`.
  - Delete the `dependend_value, dependent_field = is_dependent(…)` lines. Replace the final `if dependend_value … else … end` with `@timeit "apply_pointwise_E" apply_pointwise_E(nodes, E, bond_force)`. The dependent branch could not be reached, because a `.txt` Young's Modulus was never accepted. Ledger this as a ruling.
- Drop `is_dependent` and `get_symmetry` from the imports if they become unused.
- Update both docstrings' argument lines: `p::BondbasedElasticParams` and `material::BlockMaterial` (typed block material).

`1D_Bondbased_Elastic.jl`:
- `init_model(nodes, material_parameter::Dict)` becomes `init_model(nodes::AbstractVector{Int64}, p::OneDBondbasedElasticParams, material)`.
- The `compute_model` signature becomes `(nodes::AbstractVector{Int64}, p::OneDBondbasedElasticParams, material, block::Int64, time::Float64, dt::Float64)`.
- `E = material_parameter["Young's Modulus"]` becomes `E = material.moduli.youngs_modulus`.
- `id1 = material_parameter["Id1"]` becomes `id1 = p.id1`, and `id2 = material_parameter["Id2"]` becomes `id2 = p.id2`.

`Unified_Bondbased_Elastic.jl`:
- `init_model` becomes `(nodes::AbstractVector{Int64}, p::UnifiedBondbasedElasticParams, material)`.
- `symmetry::String = get_symmetry(material_parameter)` becomes `material.symmetry`.
- `nu = material_parameter["Poisson's Ratio"]` becomes `nu = material.moduli.poissons_ratio`.
- `E = material_parameter["Young's Modulus"]` becomes `E = material.moduli.youngs_modulus`.
- `compute_model` becomes `(nodes::AbstractVector{Int64}, p::UnifiedBondbasedElasticParams, material, block::Int64, time::Float64, dt::Float64)`, with the same two replacements for `E` and `symmetry`.

`Rigid.jl`:
- `init_model(nodes, material_parameter::Dict{String,Any})` becomes `init_model(nodes::AbstractVector{Int64}, p::RigidParams, material)`.
- `compute_model(nodes, material_parameter::Dict{String,Any}, block, time, dt)` becomes `compute_model(nodes::AbstractVector{Int64}, p::RigidParams, material, block::Int64, time::Float64, dt::Float64)`.

`material_template.jl`:
- `init_model` becomes `(nodes::AbstractVector{Int64}, p::MaterialTemplateParams, material)`.
- `compute_model` becomes `(nodes::AbstractVector{Int64}, p::MaterialTemplateParams, material, block::Int64, time::Float64, dt::Float64)`.
- Update the docstrings: `p` holds the model's own parameters; `material.base`, `material.moduli` and `material.symmetry` hold the shared ones.

Then run `grep -rn 'material_parameter' src/Models/Material/Material_Models/BondBased src/Models/Material/Material_Models/Rigid src/Models/Material/Material_Models/Material_template/material_template.jl`. Expected: no code reads remain (docstring mentions may be updated or removed).

- [ ] **Step 4: Run to verify they pass**

Run: the Step 2 command, plus `unit_tests/Models/Material/Material_Models/BondBased/ut_Unified_Bondbased_Elastic` and `unit_tests/Models/Material/ut_Material_Factory`
Expected: all pass, with the same force values as before.

Then the bond-based fullscale tests: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl fullscale_tests/test_Bond_Based_Elastic/test_Bond_Based_Elastic` (check the runner name with `ls test/fullscale_tests/test_Bond_Based_Elastic/*.jl`).
Expected: pass.

- [ ] **Step 5: Full suite, then commit**

Full suite in background; read the tail. Expected: 0 failures, 0 errors.

```bash
git add src/Models/Material test/unit_tests/Models/Material
git commit -m "Material factory dispatches typed parts; bond-based and rigid models read the block material"
```

---

### Task 5: PD Solid Elastic and PD Solid Plastic on `BlockMaterial`

**Files:**
- Modify: `src/Models/Material/Material_Models/Ordinary/PD_Solid_Elastic.jl`, `PD_Solid_Plastic.jl`
- Test:
  - `test/unit_tests/Models/Material/Material_Models/Ordinary/ut_PD_Solid_Plastic.jl` (`ut_init_model`)
  - `ut_PD_Solid_Elastic.jl` (if it calls init/compute with a Dict)
  - `test/unit_tests/Models/Material/ut_block_material.jl` (append)

**Interfaces:**
- Consumes (Task 4): the module interface `init_model(nodes, p, material)` and `compute_model(nodes, p, material, block, time, dt)`; `Material.bind_material!`.
- Produces: `PDSolidElasticParams` / `PDSolidPlasticParams` models computing from `material.moduli`, `material.symmetry` and `p.yield_stress` (`value(p.yield_stress, iID)`).

- [ ] **Step 1: Write the failing tests**

Append to `ut_block_material.jl`:

```julia
@testset "table follows the NP1 switch" begin
    ut_reset(3; nnodes = 2)
    dir = mktempdir()
    write(joinpath(dir, "ys.txt"), "header: Temperature Yield_Stress\n0 10\n100 20\n")
    N, NP1 = PeriLab.Data_Manager.create_node_scalar_field("Temperature", Float64)
    NP1 .= [0.0, 100.0]
    ctx = PeriLab.ParameterSpec.ParseContext(directory = dir)
    wb = PeriLab.ParameterSpec.parse_model(:material,
                                           Dict{String,Any}("Material Model" => "PD Solid Plastic",
                                                            "Bulk Modulus" => 1.0,
                                                            "Shear Modulus" => 1.0,
                                                            "Yield Stress" => "ys.txt"),
                                           "M", ctx; name_key = "Material Model")
    @test isempty(ctx.errors)
    m = BMAT.block_material(wb, "PD Solid Plastic", 3)
    BMAT.bind_material!(m)
    @test PeriLab.ParameterSpec.value(m.model.yield_stress, 2) ≈ 20.0
    PeriLab.Data_Manager.switch_NP1_to_N()
    PeriLab.Data_Manager.get_field("Temperature", "NP1") .= [100.0, 0.0]
    BMAT.bind_material!(m)
    @test PeriLab.ParameterSpec.value(m.model.yield_stress, 1) ≈ 20.0
    @test PeriLab.ParameterSpec.value(m.model.yield_stress, 2) ≈ 10.0
end

@testset "composite parts in order" begin
    ut_reset(3; nnodes = 2)
    m = typed_block_material(Dict("Material Model" => "PD Solid Elastic + PD Solid Plastic",
                                  "Bulk Modulus" => 1.0, "Shear Modulus" => 1.0,
                                  "Yield Stress" => 2.0))
    @test parentmodule.(typeof.(BMAT.model_parts(m.model))) ==
          (BMAT.PD_Solid_Elastic, BMAT.PD_Solid_Plastic)
end
```

In `ut_PD_Solid_Plastic.jl`, `@testset "ut_init_model"`:
- Replace the `@test_logs (:error, "Yield Stress is not defined in input deck") … Dict()` block with

```julia
    @test_throws PeriLab.PeriLabError typed_block_material(Dict("Material Model" => "PD Solid Plastic",
                                                                "Bulk Modulus" => 1.0,
                                                                "Shear Modulus" => 1.0))
```

- Replace the two `init_model(…, Dict("Yield Stress" => 5.3))` and `init_model(…, Dict("Yield Stress" => 2.2, "Symmetry" => "plane stress"))` calls with

```julia
    material = typed_block_material(Dict("Material Model" => "PD Solid Plastic",
                                         "Bulk Modulus" => 1.0, "Shear Modulus" => 1.0,
                                         "Yield Stress" => 5.3))
    PeriLab.Solver_Manager.Model_Factory.Material.PD_Solid_Plastic.init_model(Vector{Int64}(1:nodes),
                                                                              material.model,
                                                                              material)
```

  and, for the second (2D) case,

```julia
    PeriLab.Data_Manager.set_dof(2)
    material = typed_block_material(Dict("Material Model" => "PD Solid Plastic",
                                         "Bulk Modulus" => 1.0, "Shear Modulus" => 1.0,
                                         "Yield Stress" => 2.2,
                                         "Symmetry" => "plane stress"))
    PeriLab.Solver_Manager.Model_Factory.Material.PD_Solid_Plastic.init_model(Vector{Int64}(1:nodes),
                                                                              material.model,
                                                                              material)
```

  The old test kept dof 3 with `"Symmetry" => "plane stress"`. The typed symmetry ignores plane stress in 3D, as `check_symmetry` does in a real run. The 2D formula is therefore tested with dof 2, and the expected values stay unchanged. Ledger this.
- Keep the expected `isapprox` values unchanged.

In `ut_PD_Solid_Elastic.jl`: `grep -n 'init_model\|compute_model'`. Convert any Dict call the same way (`"Material Model" => "PD Solid Elastic"` with the dict's moduli). If there are none, leave the file unchanged.

- [ ] **Step 2: Run to verify they fail**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Models/Material/ut_block_material unit_tests/Models/Material/Material_Models/Ordinary/ut_PD_Solid_Plastic unit_tests/Models/Material/Material_Models/Ordinary/ut_PD_Solid_Elastic`
Expected:
- "table follows the NP1 switch" and "composite parts in order" pass already: they test Task 1/4 infrastructure. They are pins for this task's code paths; ledger that.
- The `ut_PD_Solid_Plastic` `init_model` calls fail with a `MethodError`.

- [ ] **Step 3: Implement**

`PD_Solid_Elastic.jl`:
- `init_model(nodes, material_parameter::Dict)` becomes `init_model(nodes::AbstractVector{Int64}, p::PDSolidElasticParams, material)`.
- The `compute_model` signature becomes `(nodes::AbstractVector{Int64}, p::PDSolidElasticParams, material, block::Int64, time::Float64, dt::Float64)`.
- Replace

```julia
    shear_modulus::Union{Float64,Vector{Float64}} = material_parameter["Shear Modulus"]
    bulk_modulus::Union{Float64,Vector{Float64}} = material_parameter["Bulk Modulus"]

    symmetry::String = get_symmetry(material_parameter)
```

with

```julia
    shear_modulus = material.moduli.shear_modulus
    bulk_modulus = material.moduli.bulk_modulus
    symmetry::String = material.symmetry
```

- Drop `get_symmetry` from the imports if unused.

`PD_Solid_Plastic.jl`:
- Imports: replace `is_dependent, interpol_data, get_dependent_value` in the Helpers import (keep the others), and add `using ......ParameterSpec: value, Table1D`. Drop `get_symmetry` if unused.
- Replace `init_model` with:

```julia
function init_model(nodes::AbstractVector{Int64}, p::PDSolidPlasticParams, material)
    horizon = Data_Manager.get_field("Horizon")
    yield = Data_Manager.create_constant_node_scalar_field("Yield Value", Float64)
    set_yield_value!(yield, nodes, p.yield_stress, material.symmetry, horizon)

    Data_Manager.create_constant_bond_scalar_state("Deviatoric Plastic Extension State",
                                                   Float64)
    Data_Manager.create_node_scalar_field("Lambda Plastic", Float64)
    Data_Manager.create_constant_node_scalar_field("TD Norm", Float64)
    Data_Manager.create_constant_bond_scalar_state("Bond Forces Deviatoric", Float64)
    Data_Manager.create_constant_bond_scalar_state("Bond Forces Isotropic", Float64)
end

"Yield value per node from the (possibly field dependent) yield stress."
function set_yield_value!(yield, nodes::AbstractVector{Int64}, yield_stress, symmetry::String,
                          horizon)
    if symmetry == "3D"
        for iID in nodes
            ys = value(yield_stress, iID)
            yield[iID] = 25 * ys * ys / (8 * pi * horizon[iID]^5)
        end
    else
        thickness::Float64 = 1 # is a placeholder
        for iID in nodes
            ys = value(yield_stress, iID)
            yield[iID] = 225 * ys * ys / (24 * thickness * pi * horizon[iID]^4)
        end
    end
    return yield
end
```

  In `init_model`, the yield stress table is not yet bound when it comes from a file, because the factory binds before compute only. So the factory's `init_model` path must bind first: in `Material_Factory.jl` `init_model`, call `bind_material!(material)` right before the `for part in model_parts(material.model)` loop. If the dependent field does not exist yet at init, `bind_material!` aborts with a message naming the field. The legacy `get_dependent_value` aborted in that case too.
- In `compute_model`:
  - The signature becomes `(nodes::AbstractVector{Int64}, p::PDSolidPlasticParams, material, block::Int64, time::Float64, dt::Float64)`.
  - `symmetry::String = get_symmetry(material_parameter)` becomes `material.symmetry`.
  - `shear_modulus = material_parameter["Shear Modulus"]` becomes `material.moduli.shear_modulus`, and `bulk_modulus = material_parameter["Bulk Modulus"]` becomes `material.moduli.bulk_modulus`.
  - Replace the whole `dependend_value, dependent_field = is_dependent(…)` block (through its `end`) with

```julia
    if p.yield_stress isa Table1D
        set_yield_value!(yield_value, nodes, p.yield_stress, symmetry,
                         Data_Manager.get_field("Horizon"))
    end
```

  The legacy dependent branch referenced an undefined `horizon`, so it could not run. This version recomputes the yield value from the table each step, as that branch intended. Ledger it.

Run `grep -rn 'material_parameter' src/Models/Material/Material_Models/Ordinary`. Expected: no code reads remain.

- [ ] **Step 4: Run to verify they pass**

Run: the Step 2 command
Expected: all pass.

Then the PD Solid fullscale tests and the point-wise material test. Find their runners with `ls test/fullscale_tests/test_PD_Solid_Elastic/*.jl test/fullscale_tests/test_PD_solid_plastic/*.jl test/fullscale_tests/test_PD_solid_elastic_3D/*.jl test/fullscale_tests/test_point_wise_material/*.jl`, and run them with the unit runner and `JULIA_PROJECT` set.
Expected: pass.

- [ ] **Step 5: Full suite, then commit**

Full suite in background; read the tail. Expected: 0 failures, 0 errors.

```bash
git add src/Models/Material test/unit_tests/Models/Material
git commit -m "PD Solid Elastic and Plastic read the block material; yield stress tables bound per step"
```
