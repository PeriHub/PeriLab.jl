# Boundary Conditions

Every entry of `Boundary Conditions` prescribes one field on one or more node sets:

```yaml
Boundary Conditions:
  BC_1:
    Type: "Dirichlet"
    Variable: "Displacements"
    Node Set: "Node Set 1"
    Coordinate: "y"
    Value: "0"
  BC_2:
    Type: "Neumann"
    Variable: "Force Densities"
    Node Set: "Node Set 2 + Node Set 3"
    Coordinate: "x"
    Value: "1e3*sin(t)"
  Temperature:
    Type: "Initial"
    Variable: "Temperature"
    Node Set: "Node Set 4"
    Value: "20 + 2*x"
```

- **Type**: `Dirichlet` (default) sets the value in every step, `Neumann` adds it to the field
  in every step, `Initial` sets it once at the start.
- **Variable**: the field, e.g. `Displacements`, `Forces`, `Force Densities` or `Temperature`.
- **Coordinate**: the component `x`, `y` or `z` of a vector field; not needed for scalar fields
  such as `Temperature`.
- **Node Set**: a node set of the [Discretization](mesh_input.md); several are joined with `+`.
- **Value**: a number or an expression in `t` (time), `st` (time within the current solver
  step) and the coordinates `x`, `y`, `z`.
- **Step ID**: in a [Multistep Solver](solver/solver.md#Multistep-Solver) the steps the condition is
  active in, one id or a comma-separated list; all steps if not given.

```@eval
using PeriLab, Markdown
Markdown.parse(PeriLab.parameter_table("Boundary Conditions"))
```
