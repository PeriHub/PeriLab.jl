# Input File

A PeriLab input deck is a YAML file with a single top-level key `PeriLab`. Its sections are
read and checked against the parameter declarations in the source; unknown keys are errors
(run with `--no_strict` or set `Strict Validation: false` to get warnings instead). Every key,
its type, default and allowed range is listed in the generated
[input sections](../generated/input_sections.md) and [input models](../generated/input_models.md)
reference.

```yaml
PeriLab:
  Discretization:            # mesh, node sets, bond filters
    Type: "Text File"
    Input Mesh File: "mesh.txt"
    Node Sets:
      Node Set 1: "ns_bottom.txt"
      Node Set 2: "ns_top.txt"
  Models:                    # named models; blocks refer to them by name
    Material Models:
      Aluminium:
        Material Model: "PD Solid Elastic"
        Symmetry: "isotropic plane strain"
        Young's Modulus: 70000.0
        Poisson's Ratio: 0.3
    Damage Models:
      Damage:
        Damage Model: "Critical Stretch"
        Critical Value: 0.01
  Blocks:                    # one entry per block of the mesh
    block_1:
      Block ID: 1
      Material Model: "Aluminium"
      Damage Model: "Damage"
      Density: 2.7e-9
      Horizon: 0.3
  Boundary Conditions:
    BC_1:
      Type: "Dirichlet"
      Variable: "Displacements"
      Node Set: "Node Set 1"
      Coordinate: "y"
      Value: "0"
    BC_2:
      Type: "Dirichlet"
      Variable: "Displacements"
      Node Set: "Node Set 2"
      Coordinate: "y"
      Value: "0.1*t"
  Outputs:
    Output1:
      Output Filename: "result"
      Number of Output Steps: 50
      Output Variables:
        Displacements: true
        Damage: true
  Solver:
    Initial Time: 0.0
    Final Time: 1.0e-3
    Damage Models: true
    Verlet:
      Safety Factor: 0.9
```

| Section | Content | Reference |
|---|---|---|
| `Discretization` | mesh file, node sets, bond filters, surface extrusion, influence function | [Discretization](../generated/input_sections.md#Discretization), [mesh input](mesh_input.md), [bond filters](bond_filter.md) |
| `Models` | material, damage, thermal, additive and degradation models, pre-calculations | [input models](../generated/input_models.md), [models](models/overview.md) |
| `Blocks` | per mesh block: density, horizon and the names of its models | [Blocks](../generated/input_sections.md#Blocks) |
| `Boundary Conditions` | Dirichlet / Neumann conditions on node sets, values as expressions in `t`, `x`, `y`, `z` | [boundary conditions](boundary_conditions.md), [Boundary Conditions](../generated/input_sections.md#Boundary-Conditions) |
| `Solver` or `Multistep Solver` | time range, active model types, one solver | [Solver](../generated/input_sections.md#Solver), [solver options](solver/solver.md) |
| `Outputs` | Exodus or CSV output sets and their variables | [Outputs](../generated/input_sections.md#Outputs), [output](output.md) |
| `Compute Class Parameters` | global values (sum, maximum, ...) of blocks or node sets for CSV output | [Compute Class Parameters](../generated/input_sections.md#Compute-Class-Parameters) |
| `Contact` | contact models between blocks | [contact](models/contact.md) |
| `FEM` | finite element blocks and FE-PD coupling | [FEM](../generated/input_sections.md#FEM) |
| `Surface Correction` | correction of the surface effect | [Surface Correction](../generated/input_sections.md#Surface-Correction) |

In a Julia session, `PeriLab.describe("Blocks")` prints the keys of a section or a model and
`PeriLab.describe("Blocks"; template = true)` a YAML block to start from. The typed structs the
deck is read into are described in the [developer guide](@ref "Parameters").
