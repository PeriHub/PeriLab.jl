# Typed Input Parameters — Phase 3c-2: Last Material-Dict Readers, Then Remove the Material Dict

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Nothing reads `Data_Manager.get_properties(block, "Material Model")` any more. The material properties dict is no longer stored. All dict-based material code is deleted: `get_all_elastic_moduli`, the dict Hooke matrix, `write_moduli!` and every dict method of the material modules.

**Architecture:**
- **Correspondence flag.** `BlockMaterial` gains `correspondence::Bool`, set from the block's model name exactly like the old `occursin("Correspondence", …)` tests. Every remaining dict reader switches to the block material:
  - the matrix-based correspondence path and the three matrix solvers;
  - the strain compute class (via `hooke_matrix`), `Pre_Calculation.check_dependencies`;
  - the critical time step, `Accuracy Order`, local damping;
  - the FEM material lookup.
- **Generic dispatch.** `Model_Factory`'s generic per-block dispatch asks `Data_Manager.get_block_material` for the material category instead of `check_property` / `get_properties`.
- **No material dict.** `read_properties` stops storing it.
- **Legacy removal.** The legacy functions are deleted from `src`. The equivalence tests keep comparing against a test-only copy of the old functions (`legacy_material_oracle.jl`), so the typed code stays pinned to the old numbers.

**Tech Stack:** Julia 1.12, `PeriLab.Data_Manager`, `PeriLab.ParameterSpec`.

**Spec:** `docs/superpowers/specs/2026-10-03-typed-input-parameters-design.md` (§3 Use: "`get_properties` / … are removed"; §5 phase 3 material category)

## Global Constraints

- Numerical results unchanged: the full suite must stay green. That includes the matrix-based decks:
  - `test_linear_static_matrix_based_solver`, `test_large_model_matrix_based`, `test_model_reduction_Verlet_matrix`, `test_Newmark`;
  - `test_calculation` (Calculate Strain), the bond-based, PD solid and correspondence decks, and `test_DCB` (local damping, bond-associated).
- `Data_Manager.get_properties` / `check_property` stay for the other categories (damage, thermal, additive, degradation, pre-calculation; later phases). Only the material category stops using them.
- Read fields directly; functions only for real logic.
- `ut_MPI.jl` runs standalone under mpiexec without `helper.jl`. Never use `typed_block_material` there.
- `/dev/shm` is 64 MB. If MPI tests die with "Bus error (signal 7)", check that no MPI process runs, then delete the stale `/dev/shm/mpich_shm_*` files.
- Branch `feature/typed-input-parameters`; one commit per task; never merge. Read the full-suite result before committing a task that runs it.

## Review Focus

1. **The critical time step of every symmetry** equals the legacy value:
   - isotropic: the bulk modulus, also per node;
   - orthotropic: the compliance formula;
   - anisotropic: C44/C55/C66;
   - transverse isotropic: G_xy / 2;
   - no constants: abort.

   Pinned in Task 1 (`critical bulk modulus`).
2. **A bond-associated correspondence block** gets the pre-calculations `Bond Associated Correspondence`; any other correspondence block gets `Shape Tensor` + `Deformation Gradient`; a non-correspondence block gets only `Deformed Bond Geometry`. Pinned in Task 1 (`pre-calculation dependencies`).
3. **A 2D material without plane strain / plane stress in `Symmetry`** still aborts with "Model definition is missing; plane stress or plane strain has to be defined for 2D". Pinned in Task 3 (`2D symmetry check`).
4. **`Calculate Strain` for a non-correspondence block** uses the inverse of the same Hooke matrix as before. Pinned by `test_calculation` (fullscale) and Task 1 (`strain Hooke matrix`).
5. **A block without a material** is skipped by the material dispatch, and a block whose material is missing aborts with "Block N has no material model defined." Pinned in Task 3 (`dispatch without the material dict`).

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
- Full suite (~30 min, background, read the tail): `julia --project=. -e 'using Pkg; Pkg.test()' > <workspace>/suite.log 2>&1`

---

### Task 1: Correspondence flag; typed critical time step, accuracy order, local damping, pre-calculation dependencies, strain compute

**Files:**
- Modify: `src/Models/Material/Material_Factory.jl` (`BlockMaterial.correspondence`, `block_material`, `critical_bulk_modulus`)
- Modify: `src/Models/Model_Factory.jl` (accuracy order in `init_models`; local damping in `init_models`; `compute_crititical_time_step`)
- Modify: `src/Models/Material/Material_Basis.jl` (`init_local_damping_due_to_damage(nodes, symmetry::String, damage_parameter)`)
- Modify: `src/Models/Material/Material_Factory.jl` (`init_local_damping(nodes, symmetry, damage_parameter)`)
- Modify: `src/Models/Pre_calculation/Pre_Calculation_Factory.jl` (`check_dependencies`)
- Modify: `src/Compute/compute_field_values.jl` (`calculate_stresses`)
- Test: `test/unit_tests/Models/Material/ut_block_material.jl` (append), `test/unit_tests/Models/Material/ut_material_basis.jl` (`ut_init_local_damping_due_to_damage`)

**Interfaces:**
- Produces:
  - `BlockMaterial` field `correspondence::Bool`, placed after `extras`, so it is the last field.
  - `Material.critical_bulk_modulus(material)`: the bulk modulus the critical time step uses (`Float64`, a node vector, or `nothing`).
  - `Material_Basis.init_local_damping_due_to_damage(nodes, symmetry::String, damage_parameter)` and `Material.init_local_damping(nodes, symmetry::String, damage_parameter)`.

- [ ] **Step 1: Write the failing tests** (append to `ut_block_material.jl`)

```julia
@testset "correspondence flag" begin
    ut_reset(3)
    @test typed_block_material(Dict("Material Model" => "Correspondence Elastic",
                                    "Bulk Modulus" => 1.0, "Shear Modulus" => 1.0)).correspondence
    @test !typed_block_material(Dict("Material Model" => "PD Solid Elastic",
                                     "Bulk Modulus" => 1.0, "Shear Modulus" => 1.0)).correspondence
end

@testset "critical bulk modulus" begin
    ut_reset(3)
    iso = typed_block_material(Dict("Material Model" => "PD Solid Elastic",
                                    "Bulk Modulus" => 7.0, "Shear Modulus" => 2.0))
    @test BMAT.critical_bulk_modulus(iso) == 7.0
    ortho = typed_block_material(Dict("Material Model" => "Correspondence Elastic",
                                      "Symmetry" => "orthotropic",
                                      "Young's Modulus X" => 2.0, "Young's Modulus Y" => 1.5,
                                      "Young's Modulus Z" => 1.0, "Poisson's Ratio XY" => 0.3,
                                      "Poisson's Ratio YZ" => 0.25, "Poisson's Ratio XZ" => 0.2,
                                      "Shear Modulus XY" => 0.7, "Shear Modulus YZ" => 0.6,
                                      "Shear Modulus XZ" => 0.5))
    s11, s22, s33 = 1 / 2.0, 1 / 1.5, 1 / 1.0
    s12, s23, s13 = -0.3 / 2.0, -0.25 / 1.0, -0.2 / 1.0
    @test BMAT.critical_bulk_modulus(ortho) ≈
          1 / (s11 + s22 + s33 + 2 * (s12 + s23 + s13))
    aniso = typed_block_material(merge(Dict{String,Any}("Material Model" => "Correspondence Elastic",
                                                        "Symmetry" => "anisotropic"),
                                       Dict("C$i$j" => (i == j ? 10.0 * i : 1.0)
                                            for i in 1:6 for j in i:6)))
    @test BMAT.critical_bulk_modulus(aniso) == maximum([40.0 / 2, 50.0 / 2, 60.0 / 2])
    transverse = typed_block_material(Dict("Material Model" => "Correspondence Elastic",
                                           "Symmetry" => "transverse isotropic plane stress",
                                           "Young's Modulus X" => 2.0,
                                           "Young's Modulus Y" => 1.5,
                                           "Poisson's Ratio XY" => 0.3,
                                           "Shear Modulus XY" => 0.8))
    @test BMAT.critical_bulk_modulus(transverse) == 0.4
end

@testset "pre-calculation dependencies" begin
    function ut_dependencies(raw)
        ut_reset(3; nnodes = 2)
        PeriLab.Data_Manager.set_block_id_list([1])
        PeriLab.Data_Manager.init_properties()
        PeriLab.Data_Manager.set_block_material(1, typed_block_material(raw))
        PeriLab.Solver_Manager.Model_Factory.Pre_Calculation.check_dependencies(Dict(1 => [1,
                                                                                            2]))
        return PeriLab.Data_Manager.get_properties(1, "Pre Calculation Model")
    end
    corr = Dict{String,Any}("Material Model" => "Correspondence Elastic",
                            "Symmetry" => "isotropic", "Bulk Modulus" => 1.0,
                            "Shear Modulus" => 1.0)
    p = ut_dependencies(merge(corr, Dict{String,Any}("Bond Associated" => true)))
    @test p["Bond Associated Correspondence"] && p["Deformed Bond Geometry"]
    @test !haskey(p, "Shape Tensor")
    p = ut_dependencies(corr)
    @test p["Shape Tensor"] && p["Deformation Gradient"]
    p = ut_dependencies(Dict{String,Any}("Material Model" => "PD Solid Elastic",
                                         "Bulk Modulus" => 1.0, "Shear Modulus" => 1.0))
    @test p["Deformed Bond Geometry"] && !haskey(p, "Shape Tensor")
end

@testset "strain Hooke matrix" begin
    ut_reset(2)
    m = typed_block_material(Dict("Material Model" => "PD Solid Elastic",
                                  "Symmetry" => "isotropic plane strain",
                                  "Bulk Modulus" => 10.0, "Shear Modulus" => 4.0))
    legacy = ut_legacy_dict(Dict{String,Any}("Material Model" => "PD Solid Elastic",
                                             "Symmetry" => "isotropic plane strain",
                                             "Bulk Modulus" => 10.0,
                                             "Shear Modulus" => 4.0); dof = 2)
    @test Matrix(BBASIS.hooke_matrix(m, 2)) ≈
          Matrix(BBASIS.get_Hooke_matrix(legacy, legacy["Symmetry"], 2))
end
```

In `ut_material_basis.jl`, `@testset "ut_init_local_damping_due_to_damage"`: in both calls, replace the second argument `Dict()` with `"3D"`.

If the "pre-calculation dependencies" testset reveals that `check_dependencies` (lines after the shown part) also adds `Shape Tensor` for bond-associated blocks through its own dependency resolution, adapt only the `!haskey(p, "Shape Tensor")` line to what the legacy function produced on the dict path. Find that out by running the legacy function once on the same dict before you change it, and ledger it.

- [ ] **Step 2: Run to verify they fail**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Models/Material/ut_block_material unit_tests/Models/Material/ut_material_basis`
Expected:
- `type BlockMaterial has no field correspondence`;
- `UndefVarError: critical_bulk_modulus`;
- pre-calculation dependencies: a `KeyError` or a missing `Bond Associated Correspondence` (the dict is empty);
- local damping: a `MethodError`.

"strain Hooke matrix" passes already; it pins the replacement used in Step 3.

- [ ] **Step 3: Implement**

`Material_Factory.jl`:
- Add `correspondence::Bool   # model name contains "Correspondence" (as the legacy dict tests)` as the last field of `BlockMaterial`.
- In `block_material`, append the argument `occursin("Correspondence", model_name)` to the `BlockMaterial(…)` call.
- Add after `block_material`:

```julia
_constant(x) = x === nothing ? nothing : value(x, 1)

"""
    critical_bulk_modulus(material)

Bulk modulus used for the critical time step (legacy rules): the isotropic bulk
modulus; for orthotropic constants the compliance-based estimate; else C44/C55/C66
or Shear Modulus XY as estimates; `nothing` if none is defined.
"""
function critical_bulk_modulus(material::BlockMaterial)
    material.moduli === nothing || return material.moduli.bulk_modulus
    base = material.base
    nu_xy, nu_yz, nu_xz = _constant(base.poissons_ratio_xy), _constant(base.poissons_ratio_yz),
                          _constant(base.poissons_ratio_xz)
    if nu_xy !== nothing && nu_yz !== nothing && nu_xz !== nothing
        E_x, E_y, E_z = _constant(base.youngs_modulus_x), _constant(base.youngs_modulus_y),
                        _constant(base.youngs_modulus_z)
        s11 = 1 / E_x
        s22 = 1 / E_y
        s33 = 1 / E_z
        s12 = -nu_xy / E_x
        s23 = -nu_yz / E_z
        s13 = -nu_xz / E_z
        return 1 / (s11 + s22 + s33 + 2 * (s12 + s23 + s13))
    elseif base.c44 !== nothing && base.c55 !== nothing && base.c66 !== nothing
        return maximum([base.c44 / 2, base.c55 / 2, base.c66 / 2])
    elseif base.shear_modulus_xy !== nothing
        return _constant(base.shear_modulus_xy) / 2
    end
    return nothing
end
```

- Extend the factory's ParameterSpec import with `value`.
- Change `init_local_damping(nodes, material_parameter, damage_parameter)` to `init_local_damping(nodes, symmetry::String, damage_parameter)`; it passes `symmetry` on.

`Material_Basis.jl`, `init_local_damping_due_to_damage`:
- The signature becomes `(nodes::AbstractVector{Int64}, symmetry::String, damage_parameter)`.
- Delete the line `symmetry::String = get_symmetry(material_parameter)`; the argument is used directly.

`Model_Factory.jl`:
- **Accuracy order:** in `init_models`, replace the `if haskey(params["Models"], "Material Models") … end` block (inside `"Pre_Calculation" in solver_options["Models"]`) with

```julia
        for material in values(input.materials)
            material.base.accuracy_order === nothing ||
                Data_Manager.set_accuracy_order(material.base.accuracy_order)
        end
```

- **Local damping:** in `init_models`, the call `Material.init_local_damping(block_nodes[block], Data_Manager.get_properties(block, "Material Model"), Data_Manager.get_properties(block, "Damage Model"))` becomes `Material.init_local_damping(block_nodes[block], Data_Manager.get_block_material(block).symmetry, Data_Manager.get_properties(block, "Damage Model"))`.
- **Critical time step:** in `compute_crititical_time_step`, replace everything inside `if mechanical` from `bulk_modulus = Data_Manager.get_property(iblock, "Material Model", "Bulk Modulus")` down to (and including) the `else @abort "No time step for material is determined because of missing properties." return nothing end` chain with

```julia
            material = Data_Manager.get_block_material(iblock)
            bulk_modulus = material === nothing ? nothing :
                           Material.critical_bulk_modulus(material)
            if isnothing(bulk_modulus)
                @abort "No time step for material is determined because of missing properties."
                return nothing
            end
```

  Keep the following `t = compute_mechanical_critical_time_step(…)` lines.

`Pre_Calculation_Factory.jl`, `check_dependencies`:
- `if !Data_Manager.check_property(block_id, "Material Model") continue end` and `model_param = …` become

```julia
        material = Data_Manager.get_block_material(block_id)
        material === nothing && continue
```

- `if occursin("Correspondence", model_param["Material Model"])` becomes `if material.correspondence`.
- `if haskey(model_param, "Bond Associated") && model_param["Bond Associated"]` becomes `if material.base.bond_associated`.

`compute_field_values.jl`, `calculate_stresses`:
- Add `hooke_matrix` to the `using ..Material_Basis:` import list.
- `correspondence = occursin("Correspondence", Data_Manager.get_properties(block, "Material Model")["Material Model"])` becomes

```julia
        material = Data_Manager.get_block_material(block)
        correspondence = material.correspondence
```

- The `material_parameter = …` and `hookeMatrix = get_Hooke_matrix(material_parameter, material_parameter["Symmetry"], Data_Manager.get_dof())` lines become `hookeMatrix = hooke_matrix(material, Data_Manager.get_dof())`.

Then run `grep -n '"Material Model"' src/Models/Model_Factory.jl src/Compute/compute_field_values.jl src/Models/Pre_calculation/Pre_Calculation_Factory.jl`. Expected: only the generic-dispatch name checks (`active_model_name == "Material Model"`) and the `read_properties` / `get_block_model_definition` lines, which change in Task 3.

- [ ] **Step 4: Run to verify they pass**

Run: the Step 2 command, plus `unit_tests/Models/ut_Model_Factory`
Expected: all pass.

- [ ] **Step 5: Full suite, then commit**

Full suite in background; read the tail. Expected: 0 failures, 0 errors (`test_calculation`, `test_DCB`, critical-time-step decks).

```bash
git add src/Models/Material src/Models/Model_Factory.jl src/Models/Pre_calculation/Pre_Calculation_Factory.jl \
        src/Compute/compute_field_values.jl test/unit_tests/Models/Material
git commit -m "Critical time step, accuracy order, local damping, pre-calculation and strain read the block material"
```

---

### Task 2: Matrix-based correspondence, matrix solvers and FEM read the block material

**Files:**
- Modify: `src/Models/Material/Material_Models/Correspondence/Correspondence_matrix_based.jl` (`init_model`, `_setup_zero_energy`)
- Modify: `src/Core/Solver/Matrix_linear_static.jl`, `Matrix_Verlet.jl`, `Newmark.jl` (block loop calling `init_model`)
- Modify: `src/Core/Solver/Solver_manager.jl` (`init_FEM` material names from `input.materials`)
- Modify: `src/Models/Model_Factory.jl` (`compute_matrix_based_bond_forces`)
- Modify: `src/FEM/FEM_basis.jl` (delete the unused `get_FE_material_model`)
- Test: `test/unit_tests/Models/Material/ut_block_material.jl` (append); `test/unit_tests/FEM/ut_FEM_routines.jl` (delete `@testset "ut_get_FE_material_model"`)

**Interfaces:**
- Consumes (Task 1): `BlockMaterial.correspondence`.
- Produces: `Correspondence_matrix_based.init_model(nodes::AbstractVector{Int64}, material, block_id::Int64)`.

- [ ] **Step 1: Write the failing test** (append to `ut_block_material.jl`)

```julia
@testset "matrix-based correspondence init takes the block material" begin
    MB = PeriLab.Solver_Manager.Correspondence_matrix_based
    ut_reset(3; nnodes = 2)
    m = typed_block_material(Dict("Material Model" => "Correspondence Elastic",
                                  "Symmetry" => "isotropic", "Bulk Modulus" => 1.0,
                                  "Shear Modulus" => 1.0, "Zero Energy Control" => "Global"))
    @test hasmethod(MB.init_model, Tuple{Vector{Int64},typeof(m),Int64})
    @test !hasmethod(MB.init_model, Tuple{Vector{Int64},Dict{String,Any},Int64})
end
```

If `PeriLab.Solver_Manager.Correspondence_matrix_based` is not the module path, find it with `grep -rn 'module Correspondence_matrix_based' src` and how `Solver_manager.jl` includes it, and use that path. Ledger it.

- [ ] **Step 2: Run to verify it fails**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Models/Material/ut_block_material`
Expected: both assertions fail — only the Dict method exists (no method for the block material, and the Dict method is still there).

- [ ] **Step 3: Implement**

`Correspondence_matrix_based.jl`:
- `function init_model(nodes::AbstractVector{Int64}, material_parameter::Dict, block_id::Int64)` becomes `function init_model(nodes::AbstractVector{Int64}, material, block_id::Int64)`.
- Its first line becomes `Zero_Energy_Control.init_model(nodes, material, block_id)`.
- In `_setup_zero_energy`, replace

```julia
    if haskey(Data_Manager.get_properties(1, "Material Model"), "Zero Energy Control")
        if Data_Manager.get_properties(1, "Material Model")["Zero Energy Control"] ==
           "Global"
```

  with

```julia
    material = Data_Manager.get_block_material(1)
    if material !== nothing && material.base.zero_energy_control !== nothing
        if material.base.zero_energy_control == "Global"
```

  Keep the rest. Block 1 as before: a pre-existing limitation, ledgered in 3c-1.

`Matrix_linear_static.jl`: the block loop becomes

```julia
    for (block, nodes) in pairs(block_nodes)
        material = Data_Manager.get_block_material(block)
        if !material.correspondence
            @abort "Only Correspondence Models are supported with the Linear Static Matrix based solver"
        end
        init_model(nodes, material, block)
    end
```

`Matrix_Verlet.jl` and `Newmark.jl`: the block loop becomes

```julia
    for (block, nodes) in pairs(block_nodes)
        init_model(nodes, Data_Manager.get_block_material(block), block)
    end
```

`Solver_manager.jl`: `get(input.models, "Material Models", Dict{String,Any}())` in the `init_FEM` call becomes `input.materials`.

`Model_Factory.jl`, `compute_matrix_based_bond_forces`:
- `material_parameter = Data_Manager.get_properties(block, "Material Model")` becomes `material = Data_Manager.get_block_material(block)`.
- `if occursin("Correspondence", material_parameter["Material Model"])` becomes `if material.correspondence`.
- The `compute_correspondence_bond_forces` call keeps passing the block material (now `material`).

`FEM_basis.jl`: delete `get_FE_material_model` (with its docstring); `ut_FEM_routines.jl`: delete `@testset "ut_get_FE_material_model"`.

- [ ] **Step 4: Run to verify it passes**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Models/Material/ut_block_material unit_tests/FEM/ut_FEM_routines unit_tests/FEM/ut_FEM_Factory`
Expected: all pass.

Then the matrix-based fullscale tests. Find their runners with `ls test/fullscale_tests/{test_linear_static_matrix_based_solver,test_large_model_matrix_based,test_model_reduction_Verlet_matrix,test_Newmark,test_FEM,test_FEM_Coupling}/*.jl`, and run them with `JULIA_PROJECT` set.
Expected: pass.

- [ ] **Step 5: Full suite, then commit**

```bash
git add src test/unit_tests/Models/Material/ut_block_material.jl test/unit_tests/FEM/ut_FEM_routines.jl
git commit -m "Matrix-based correspondence, matrix solvers and FEM read the block material"
```

---

### Task 3: Material dispatch without the material dict; the dict is no longer stored

**Files:**
- Modify: `src/Models/Model_Factory.jl` (`has_block_model`, `block_model_parameters`, the four generic dispatch sites, `read_properties`, `get_block_model_definition`)
- Modify: `src/Models/Material/Material_Factory.jl` (`init_model`, `fields_for_local_synchronization`, `compute_model`, `check_material_symmetry(material, dof)`)
- Test: `test/unit_tests/Models/Material/ut_block_material.jl` (append), `test/unit_tests/Models/ut_Model_Factory.jl` (`ut_read_properties`, `read_properties builds typed block materials`), `test/unit_tests/Models/Material/ut_Material_Factory.jl`

**Interfaces:**
- Consumes (Tasks 1–2): every reader uses `Data_Manager.get_block_material`.
- Produces:
  - `Model_Factory.has_block_model(block::Int64, name::String)::Bool` and `Model_Factory.block_model_parameters(block::Int64, name::String)`, which for `"Material Model"` returns the `BlockMaterial`.
  - `Material.compute_model(nodes, material, block::Int64, time::Float64, dt::Float64)`.
  - `Material.check_material_symmetry(material::BlockMaterial, dof::Int64)`: aborts in 2D without plane strain / plane stress, warns in 3D (legacy `check_symmetry` messages).
  - After this task, `Data_Manager.get_properties(block, "Material Model")` is an empty dict in a run.

- [ ] **Step 1: Write the failing tests**

Append to `ut_block_material.jl`:

```julia
@testset "dispatch without the material dict" begin
    MF = PeriLab.Solver_Manager.Model_Factory
    ut_reset(3; nnodes = 2)
    PeriLab.Data_Manager.set_block_id_list([1, 2])
    PeriLab.Data_Manager.init_properties()
    m = typed_block_material(Dict("Material Model" => "PD Solid Elastic",
                                  "Bulk Modulus" => 1.0, "Shear Modulus" => 1.0))
    PeriLab.Data_Manager.set_block_material(1, m)
    @test MF.has_block_model(1, "Material Model")
    @test !MF.has_block_model(2, "Material Model")
    @test MF.block_model_parameters(1, "Material Model") === m
    @test !MF.has_block_model(1, "Damage Model")
    @test_logs (:error, "Block 2 has no material model defined.") @test_throws PeriLab.PeriLabError begin
        BMAT.init_model([1, 2], 2)
    end
    @test hasmethod(BMAT.compute_model, Tuple{Vector{Int64},typeof(m),Int64,Float64,Float64})
end

@testset "2D symmetry check" begin
    ut_reset(2)
    m = typed_block_material(Dict("Material Model" => "PD Solid Elastic",
                                  "Symmetry" => "isotropic", "Bulk Modulus" => 1.0,
                                  "Shear Modulus" => 1.0); dof = 2)
    @test_logs (:error,
                "Model definition is missing; plane stress or plane strain has to be defined for 2D") @test_throws PeriLab.PeriLabError begin
        BMAT.check_material_symmetry(m, 2)
    end
    ok = typed_block_material(Dict("Material Model" => "PD Solid Elastic",
                                   "Symmetry" => "isotropic plane strain",
                                   "Bulk Modulus" => 1.0, "Shear Modulus" => 1.0); dof = 2)
    @test isnothing(BMAT.check_material_symmetry(ok, 2))
end
```

In `ut_Model_Factory.jl`:
- In `@testset "read_properties builds typed block materials"`, replace the two assertions reading `get_property(1, "Material Model", …)` with

```julia
    @test isempty(PeriLab.Data_Manager.get_properties(1, "Material Model"))
```

- In `@testset "ut_read_properties"`, the second call `read_properties(params, typed_input(Dict()), true)` and the assertions after it that read `"Material Model"` properties (the four `get_property(…, "Material Model", …)` tests) are replaced by

```julia
    PeriLab.Solver_Manager.Model_Factory.read_properties(params, typed_input(Dict()), true)
    @test isempty(PeriLab.Data_Manager.get_properties(1, "Material Model"))
    @test PeriLab.Data_Manager.get_property(3, "Damage Model", "value") ==
          params["Models"]["Damage Models"]["a"]["value"]
```

In `ut_Material_Factory.jl`, `@testset "ut_init_model"`: `set_property(1, "Material Model", "E", 1)` no longer makes a block "have" a material. Keep the two abort assertions; their message is unchanged.

- [ ] **Step 2: Run to verify they fail**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Models/Material/ut_block_material unit_tests/Models/ut_Model_Factory unit_tests/Models/Material/ut_Material_Factory`
Expected:
- `UndefVarError: has_block_model`;
- `MethodError` for `check_material_symmetry(::BlockMaterial, ::Int64)`;
- the read_properties assertions fail because the dict is still stored.

- [ ] **Step 3: Implement**

`Material_Factory.jl`:
- Replace `check_material_symmetry(block::Int64)` with

```julia
"""
    check_material_symmetry(material, dof)

2D needs plane strain or plane stress in `Symmetry`; in 3D these are ignored
(with a warning). A missing Symmetry is not checked here (isotropic).
"""
function check_material_symmetry(material::BlockMaterial, dof::Int64)
    symmetry = material.base.symmetry
    symmetry === nothing && return nothing
    if dof == 2 && !occursin("plane strain", symmetry) && !occursin("plane stress", symmetry)
        @abort "Model definition is missing; plane stress or plane strain has to be defined for 2D"
        return
    end
    if dof == 3 && occursin("plane strain", symmetry)
        @warn "Plane strain symmetry is not supported for 3D, going to ignore it"
    end
    if dof == 3 && occursin("plane stress", symmetry)
        @warn "Plane stress symmetry is not supported for 3D, going to ignore it"
    end
    return nothing
end
```

- In `init_model(nodes, block)`, replace the first lines

```julia
    model_param = Data_Manager.get_properties(block, "Material Model")::Dict{String,Any}
    if !haskey(model_param, "Material Model")
        @abort "Block " * string(block) * " has no material model defined."
        return
    end

    if occursin("Correspondence", model_param["Material Model"])
        Data_Manager.set_model_module("Correspondence", Correspondence)
        material = Data_Manager.get_block_material(block)
        if material === nothing
            @abort "Block $block has no typed material parameters."
            return
        end
```

  with

```julia
    material = Data_Manager.get_block_material(block)
    if material === nothing
        @abort "Block " * string(block) * " has no material model defined."
        return
    end

    if material.correspondence
        Data_Manager.set_model_module("Correspondence", Correspondence)
```

  Keep the rest of the correspondence branch (`bind_material!`, `return Correspondence.init_model(…)`). In the non-correspondence part, delete the second `material = Data_Manager.get_block_material(block)` and its `=== nothing` abort, which are now redundant.
- `fields_for_local_synchronization(model, block)`:
  - `model_param = …; if occursin("Correspondence", model_param["Material Model"])` becomes `material = Data_Manager.get_block_material(block); if material.correspondence`.
  - The call passes `material`.
- `compute_model`: the signature becomes `function compute_model(nodes::AbstractVector{Int64}, material, block::Int64, time::Float64, dt::Float64)`, and the body becomes

```julia
    @timeit "all" begin
        if material.correspondence
            @timeit "corresponcence" begin
                Correspondence.compute_model(nodes, bind_material!(material), block, time, dt)
                return
            end
        end
        @timeit "material" compute_block_material(nodes, material, block, time, dt)
    end
```

  Update its docstring argument line to `material::BlockMaterial`.

`Model_Factory.jl`:
- Add near the top (after the `using` lines):

```julia
# the material category is typed (Data_Manager.get_block_material); the other
# categories still use the property dicts
has_block_model(block::Int64, name::String) = name == "Material Model" ?
                                              Data_Manager.get_block_material(block) !==
                                              nothing :
                                              Data_Manager.check_property(block, name)
block_model_parameters(block::Int64, name::String) = name == "Material Model" ?
                                                     Data_Manager.get_block_material(block) :
                                                     Data_Manager.get_properties(block, name)
```

- At the four generic dispatch sites (`init_models` model loop, both `compute_models` loops, `compute_stiff_matrix_compatible_models`): `Data_Manager.check_property(block, active_model_name)` becomes `has_block_model(block, active_model_name)`. In the two `compute_models` loops and in `compute_stiff_matrix_compatible_models`, `Data_Manager.get_properties(block, active_model_name)` becomes `block_model_parameters(block, active_model_name)`.
- `get_block_model_definition`: in `for model in prop_keys`, replace `if model == "Material Model" && !material_model continue end` with `model == "Material Model" && continue   # typed: input.materials`. Drop the now unused `material_model` argument from its signature and from the call in `read_properties`.
- `read_properties`: the material block becomes

```julia
    if material_model
        dof = Data_Manager.get_dof()
        for (block_name, block) in zip(block_name_list, block_id_list)
            block_params = get(input.sections.blocks, block_name, nothing)
            material_name = block_params === nothing ? nothing : block_params.material_model
            (material_name === nothing || !haskey(input.materials, material_name)) && continue
            model_name = String(input.models["Material Models"][material_name]["Material Model"])
            material = Material.block_material(input.materials[material_name], model_name, dof)
            Material.check_material_symmetry(material, dof)
            Data_Manager.set_block_material(block, material)
        end
    end
```

  Here `write_moduli!`, `determine_isotropic_parameter` and the dict `check_material_symmetry(block)` are no longer called.

Then run `grep -rn 'get_properties([^)]*"Material Model"\|get_property([^)]*"Material Model"\|check_property([^)]*"Material Model"' src`. Expected: no matches.

- [ ] **Step 4: Run to verify they pass**

Run: the Step 2 command.
Expected: all pass.

- [ ] **Step 5: Full suite, then commit**

Full suite in background; read the tail. Expected: 0 failures, 0 errors.

```bash
git add src/Models test/unit_tests/Models
git commit -m "Material dispatch reads the block material; the material dict is no longer stored"
```

---

### Task 4: Remove the legacy material code; correspondence template on the typed interface

**Files:**
- Create: `test/unit_tests/Models/Material/legacy_material_oracle.jl`
- Modify: `test/runtests.jl` (include the oracle before `ut_block_material.jl`)
- Modify: `src/Models/Material/Material_Basis.jl`, `src/Models/Material/Material_Factory.jl`
- Modify: the correspondence modules, `Bond_Associated_Correspondence.jl`, `Zero_Energy_Control.jl`, `Global_Zero_Energy_Control.jl`
- Modify: `src/Models/Material/Material_Models/Material_template/correspondence_template.jl`
- Modify tests:
  - `test/unit_tests/Models/Material/ut_block_material.jl`, `ut_material_basis.jl`;
  - `Material_Models/Correspondence/ut_Correspondence.jl`, `ut_Correspondence_Plastic.jl`, `ut_Correspondence_UMAT.jl`, `ut_Correspondence_VUMAT.jl`, `ut_Bond_Associated_Correspondence.jl`.

**Interfaces:**
- Consumes (Task 3): nothing in `src` calls a material dict method any more.
- Produces: test module `LegacyMaterialOracle` with `get_all_elastic_moduli(dict)`, `get_Hooke_matrix(dict, symmetry, dof, ID = 1)`, `get_symmetry(dict)` and `flaw_function(dict, coor, stress)`. These are verbatim copies of the deleted legacy functions, so the equivalence tests keep their reference.

- [ ] **Step 1: Freeze the legacy reference before deleting anything**

Create the oracle by extracting the legacy functions verbatim (run from the repo root while they still exist):

```julia
src = read("src/Models/Material/Material_Basis.jl", String)
function grab(header)
    i = findfirst(header, src); i === nothing && error("missing " * header)
    j = findnext("\nend\n", src, last(i))
    return src[first(i):last(j)]
end
parts = [grab("function get_value(parameter::Union{Dict{Any,Any},Dict{String,Any}},"),
         grab("function get_all_elastic_moduli(parameter::Union{Dict{Any,Any},Dict{String,Any}})"),
         grab("function get_symmetry(material::Dict)"),
         grab("function _dict_constant(parameter::Dict, key::String, id::Int64)"),
         grab("function get_Hooke_matrix(parameter::Dict, symmetry::String, dof::Int64, ID::Int64 = 1)"),
         grab("function flaw_function(params::Dict,")]
header = """
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# Test-only copies of the legacy (material dict) functions removed in phase 3c-2.
# The typed implementation is compared against them, so its numbers stay pinned.
module LegacyMaterialOracle
using PeriLab.Data_Manager
using PeriLab.PeriLabExceptions: @abort
using PeriLab.Helpers: get_dependent_value_with_ID
using PeriLab.Solver_Manager.Material_Basis: _hooke_matrix

"""
write("test/unit_tests/Models/Material/legacy_material_oracle.jl",
      header * join(parts, "\n") * "end\n")
```

  Then:
  - In `test/runtests.jl`, include `unit_tests/Models/Material/legacy_material_oracle.jl` once, directly before the `ut_block_material` testset (outside any testset, so the module is defined once).
  - In `ut_block_material.jl`, add at the top `const LEGACY = isdefined(@__MODULE__, :LegacyMaterialOracle) ? LegacyMaterialOracle : (include(joinpath(@__DIR__, "legacy_material_oracle.jl")); LegacyMaterialOracle)`.
  - Replace every `BBASIS.get_all_elastic_moduli` / `BBASIS.get_Hooke_matrix` / `BBASIS.get_symmetry` / `BBASIS.flaw_function(Dict(` with `LEGACY.…`.
  - Leave typed calls (`BBASIS.hooke_matrix`, `BBASIS.flaw_function(m.base…` / `(nothing…`) unchanged.

- [ ] **Step 2: Snapshot the two equivalence tests that call dict module methods**

`"Correspondence Plastic typed compute equals legacy"` compares against the dict methods of `Correspondence_Plastic`. Run it once with this added at the end of its loop body: `typed || println(repr(stress_NP1))`. Paste the printed array as `const UT_PLASTIC_REFERENCE = <array>` above the testset, change the loop to `for typed in (true,)` and compare `results[1] ≈ UT_PLASTIC_REFERENCE`. Remove the `println`.

`"Correspondence Elastic typed init equals legacy"` compares against `get_Hooke_matrix`; it now uses `LEGACY.get_Hooke_matrix` and stays as is.

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Models/Material/ut_block_material`
Expected: all pass (still against the legacy code, now via the oracle).

- [ ] **Step 3: Delete the legacy code**

`Material_Basis.jl`: delete `get_value`, `get_all_elastic_moduli`, `check_symmetry(block)`, `get_symmetry(material::Dict)`, `_dict_constant`, the dict `get_Hooke_matrix`, and the dict `flaw_function(params::Dict, …)`, with their docstrings and `export` lines. Keep `_hooke_matrix`, `hooke_matrix`, `_typed_constant`, the typed `flaw_function` methods, `get_2D_Hooke_matrix` and the rest.

`Material_Factory.jl`:
- Delete `determine_isotropic_parameter`, `write_moduli!` and the removed names from the `using ...Material_Basis:` list.
- Delete `get_all_elastic_moduli`, which is listed twice there, and `check_symmetry`.

Correspondence family. Delete the dict methods (the ones whose material argument is a `Dict`):
- `Correspondence.jl`: dict `init_model(nodes, block, material_parameter::Dict{String,Any})`, dict `fields_for_local_synchronization(model, block, model_param::Dict)`, dict `compute_model`, dict `compute_correspondence_model`.
- `Bond_Associated_Correspondence.jl`: dict `init_model`, dict `compute_model`.
- `Correspondence_Elastic.jl`: dict `init_model`, dict `compute_stresses` wrapper, dict `compute_stresses_ba`, single-node `compute_stresses(dof::Int64, material_parameter::Dict, …)`.
- `Correspondence_Plastic.jl`: dict `init_model`, dict `compute_stresses`, dict `compute_stresses_ba`.
- `Correspondence_UMAT.jl` / `Correspondence_VUMAT.jl`: dict `init_model`, dict `compute_stresses` one-liner, dict `compute_stresses_ba`.
- `Zero_Energy_Control.jl`: dict `init_model`, dict `compute_zero_energy_control`. `Global_Zero_Energy_Control.jl`: dict `init_model`, dict `compute_control` one-liner.

Afterwards run `grep -rn 'material_parameter::Dict\|model_param::Dict\|get_dependent_value("Yield Stress"' src/Models/Material src/Models/Material/Material_Models`. Expected: no matches.

Correspondence template: rewrite `correspondence_template.jl` to the typed interface.
- Add after its `using` lines:

```julia
using .......ParameterSpec: @params, register_material

"""
    CorrespondenceTemplateParams

Declare the YAML keys your model needs beyond the shared material keys (Symmetry,
moduli, … are in `material.base` / `material.moduli`). Register it under your model
name by uncommenting `__init__` (the template stays unregistered so that a copy never
collides with it).
"""
@params struct CorrespondenceTemplateParams
end
# __init__() = register_material("Correspondence Template", CorrespondenceTemplateParams)
```

- Delete the comment block that showed this as a comment (the lines starting with `# Declare the YAML keys of your model …` through `#     __init__() = …`).
- `init_model(nodes, material_parameter::Dict, block::Int64)` becomes `init_model(nodes::AbstractVector{Int64}, p::CorrespondenceTemplateParams, material)`, with the docstring argument lines `p` (model parameters) and `material::BlockMaterial` (base, moduli, symmetry).
- The node `compute_stresses(iID::Int64, dof, material_parameter::Dict, …)` becomes `compute_stresses(nodes::AbstractVector{Int64}, dof::Int64, p::CorrespondenceTemplateParams, material, time::Float64, dt::Float64, strain_increment, stress_N, stress_NP1)`, with docstring lines updated the same way. Its `@info` line about `material_parameter` becomes "The Data_Manager, p and material hold all you need to solve your problem on material level."
- Delete the single-node `compute_stresses(dof::Int64, material_parameter::Dict, …)`.
- `compute_stresses_ba(nodes, nlist, dof, material_parameter::Dict, …)` becomes `compute_stresses_ba(nodes, nlist, dof::Int64, p::CorrespondenceTemplateParams, material, time::Float64, dt::Float64, strain_increment, stress_N, stress_NP1)`.

Tests:
- `ut_material_basis.jl`: delete the testsets `"ut_flaw_function"`, `"check_symmetry"`, `"get_symmetry"`, `"get_all_elastic_moduli"` and `"get_Hooke_matrix"`. Their behaviour is pinned through the oracle comparisons in `ut_block_material.jl`.
- `ut_Correspondence.jl`: delete `@testset "ut_init_model"`. The missing-symmetry abort is gone since 3c-1, and an unknown correspondence model name is a parse error since 3a. Delete the file and its `runtests.jl` include if no testset remains.
- `ut_Correspondence_UMAT.jl` / `ut_Correspondence_VUMAT.jl`: delete `@testset "init exceptions"`; the typed equivalents are in `ut_block_material.jl`.
- `ut_Correspondence_Plastic.jl`: replace `@testset "ut_init_model"` with

```julia
@testset "ut_init_model" begin
    nodes = 2
    PeriLab.Data_Manager.initialize_data()
    PeriLab.Data_Manager.set_num_controller(nodes)
    PeriLab.Data_Manager.set_dof(3)
    PLASTIC = PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.Correspondence_Plastic
    m = typed_block_material(Dict("Material Model" => "Correspondence Elastic + Correspondence Plastic",
                                  "Symmetry" => "isotropic", "Bulk Modulus" => 10.0,
                                  "Shear Modulus" => 10.5, "Yield Stress" => 3.4))
    PLASTIC.init_model(Vector{Int64}(1:nodes), m.model.parts[2], m)
    @test PeriLab.Data_Manager.has_key("von Mises Yield StressN")
    @test PeriLab.Data_Manager.has_key("Plastic StrainN")
    ba = typed_block_material(Dict("Material Model" => "Correspondence Elastic + Correspondence Plastic",
                                   "Symmetry" => "isotropic", "Bulk Modulus" => 10.0,
                                   "Shear Modulus" => 10.5, "Yield Stress" => 3.4,
                                   "Bond Associated" => true))
    PLASTIC.init_model(Vector{Int64}(1:nodes), ba.model.parts[2], ba)
    @test PeriLab.Data_Manager.has_key("von Mises Bond Yield StressN")
    @test PeriLab.Data_Manager.has_key("Plastic Bond StrainN")
end
```

- `ut_Bond_Associated_Correspondence.jl`: in `@testset "ut_init_Bond-Associated"`, delete the missing-symmetry `@test_logs … Dict()` block. Replace the two `init_model(nodes, Dict(…))` calls with

```julia
    m = typed_block_material(Dict("Material Model" => "Correspondence Elastic",
                                  "Symmetry" => "isotropic", "Bulk Modulus" => 1.0,
                                  "Shear Modulus" => 1.0))
    PeriLab.Solver_Manager.Model_Factory.Material.Correspondence.Bond_Associated_Correspondence.init_model(nodes,
                                                                                                           m)
```

  and, for the second call, the same with `"Accuracy Order" => 2` added. Keep the two `get_accuracy_order()` assertions.

Ledger one ruling listing the deleted testsets and why each is covered.

- [ ] **Step 4: Run to verify**

Run: `JULIA_PROJECT=/home/PeriLab.jl julia <workspace>/unit.jl unit_tests/Models/Material/ut_block_material unit_tests/Models/Material/ut_material_basis unit_tests/Models/Material/Material_Models/Correspondence/ut_Correspondence_Plastic unit_tests/Models/Material/Material_Models/Correspondence/ut_Correspondence_UMAT unit_tests/Models/Material/Material_Models/Correspondence/ut_Correspondence_VUMAT unit_tests/Models/Material/Material_Models/Correspondence/ut_Bond_Associated_Correspondence unit_tests/Models/Material/Zero_Energy_Control/ut_Global_Zero_Energy_Control`
Expected: all pass.

- [ ] **Step 5: Full suite, then commit**

Full suite in background; read the tail. Expected: 0 failures, 0 errors.

```bash
git add src/Models/Material test/unit_tests/Models/Material test/runtests.jl
git commit -m "Remove the legacy material dict code; correspondence template on the typed interface"
```
