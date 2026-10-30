# Solver

Keys of a `Solver` section (and of each step of a `Multistep Solver`); the solver is chosen by one of the sub-sections below:

```@eval
using PeriLab, Markdown
Markdown.parse(PeriLab.parameter_table("Solver"))
```

## Multistep Solver

Instead of `Solver`, a `Multistep Solver` runs several steps one after another, each with its
own time range, active models and solver. Steps run in ascending `Step ID`; a step starts at
the current time if its `Initial Time` is not given or earlier, and `Additional Time` can be
given instead of `Final Time`. A [boundary condition](../boundary_conditions.md) with a
`Step ID` is active only in these steps.

```yaml
Multistep Solver:
  Step 1:
    Step ID: 1
    Initial Time: 0.0
    Final Time: 2700.0
    Material Models: false
    Thermal Models: true
    Verlet:
      Fixed dt: 1.0
  Step 2:
    Step ID: 2
    Additional Time: 1.0e-4
    Material Models: true
    Thermal Models: false
    Verlet:
      Safety Factor: 0.95
```

## Fallback

A step that fails (e.g. a Static step that does not converge) is discarded. With a `Fallback`
solver, the step is recomputed by it from the last converged state, and the fallback runs to the
final time of the solver step; without one the run is cancelled.

```yaml
Solver:
  Initial Time: 0.0
  Final Time: 1.0e-3
  Static:
    Maximum number of iterations: 50
  Fallback:
    Verlet:
      Safety Factor: 0.9
```

During the fallback each output keeps the time between two outputs, so a fallback with many
short steps does not write every step. `Fallback Output Frequency` and `Fallback Number of
Output Steps` in an [output](../output.md) set the frequency for the fallback steps instead.

## Verlet

```@eval
using PeriLab, Markdown
Markdown.parse(PeriLab.parameter_table("Verlet"))
```

!!! warning "Fixed dt"
    If a fixed time step is defined, the time integration can become unstable.

The Verlet time integration is used as standard solver for dynamic hyperbolic differential equation of motion. It is also used in Peridigm [LittlewoodDJ2023](@cite). The displacements for step $i+1$ are solved as follows

$$\mathbf{u}_{i+1} = \mathbf{u}_{i} + \Delta t\dot{\mathbf{u}}_{i} + \frac12 \Delta t^2\ddot{\mathbf{u}}_{i}$$

with

$$\ddot{\mathbf{u}}_{i} = \frac{\mathbf{F}_i}{\rho}$$

where $\rho$ is the mass density of the point and $\mathbf{F}_i=\mathbf{F}_{external}-\mathbf{F}_{internal}$ for the current time step.

For parabolic time integration as in temperature models the following schema is used

$$\boldsymbol{\tau}_{i+1} =  \boldsymbol{\tau}_i - \Delta t \frac{\mathbf{H}}{\rho C_v}$$

where $\rho$ is the mass density, $C_v$ is the [specific heat capacity](https://en.wikipedia.org/wiki/Specific_heat_capacity) and $\mathbf{H}$ is the heat flux of each point [OterkusS2014b](@cite).

For the time intergration a stable increment has to be determined.

## Static

```@eval
using PeriLab, Markdown
Markdown.parse(PeriLab.parameter_table("Static"))
```

The static solver from [NLsolve.jl](https://github.com/JuliaNLSolvers/NLsolve.jl) has been included. Specifically the [method = :anderson](https://github.com/JuliaNLSolvers/NLsolve.jl#anderson-acceleration) is used.

The solver computes the residual of the internal reaction force densities and the external applied force densities
$$r =  \left[\underline{\mathbf{T}}_{external} + \underline{\mathbf{T}}_{internal}\right] / s_{Residual\,scaling}$$

Right now the default value $m$ of the Anderson acceleration method is chosen.

$s_{Residual\,scaling}$ should be in the range of the Young's modulus and is than scaled by the volume.

$$s_{Residual\,scaling} /= minimum(volume)^2$$

!!! warning "Multiphysics"
    Currently only the mechanical solver is included!

---

**Start Value**

The start value defines the values taken for the first iteration. The default is zero. Two options are possible. The first option defines start values in the mesh file. The name is

    start_values_x
    start_values_y
    start_values_z (optional)

These numbers are stored in the datamanager and are used for the iteration.

The alternative is to define a Linear Start Value as it is given as option in the list above. This list defines values for 3D as

$$start_{val} = [A_{x-min}\,A_{x-max}\,A_{y-min}\,A_{y-max}\,A_{z-min}\,A_{z-max}]$$

and for 2D

$$start_{val} = [A_{x-min}\,A_{x-max}\,A_{y-min}\,A_{y-max}]$$

where $A$ are amplitude values chosen by the user.

With these values the field ''start_values'' is computed as

$$start\_value(x,y,z (optional))= \frac{a(x,y,z(optional))\cdot coordinates+n(x,y,z(optional))}{nsteps}$$

$$a(x,y,z (optional))=\frac{\text{max}(start_{val})-\text{min}(start_{val})}{\text{max}(coordinates)-\text{min}(coordinates)}$$

and

$$n(x,y,z (optional))=\text{max}(start_{val})-a(x,y,z (optional))\text{max}(coordinates)$$

## Linear Static Matrix Based

```@eval
using PeriLab, Markdown
Markdown.parse(PeriLab.parameter_table("Linear Static Matrix Based"))
```

The solver computes the stiffness matrix of the problem [WillbergC2026](@cite). It is solved by

$$\mathbf{u}=\mathbf{K}_{PD}^{-1}\mathbf{F}_{external}$$

The number of steps defines the virtual time step.

### Matrix Update

Uses the previous time step as original configuration. This allows the analysis of geometrically non-linear deformations. It must be updated if damages or additive models are used.

!!! warning "Models"
    Not all models are fully tested yet in this framework.


## Verlet Matrix Based

It is the same solver as the Verlet based solver above. The main difference is, that the matrix style is used [WillbergC2026](@cite).

!!! warning "Models"
    If update is active it is not very efficient, because the creation of new matrix is more costly than the material point approach.

```@eval
using PeriLab, Markdown
Markdown.parse(PeriLab.parameter_table("Verlet Matrix Based"))
```

With `Model Reduction`, see also the [condensation methods](condensation.md):

```@eval
using PeriLab, Markdown
Markdown.parse(PeriLab.parameter_table("Verlet Matrix Based → Model Reduction"))
```


## Computational effort

A single-core performance comparison between the matrix-based solver
and the Verlet integrator for varying discretizations and horizon ratios $m =
\delta/\Delta x$ is shown in the table. All times in seconds. The sparsity pattern is a one-time preprocessing
cost and is excluded from the per-step effort in the last column. For the updated matrix algorithm this leads to a more efficient process. It was computed with v2.01

| Points  |  m | Assembly (s) | Solve (s) | Pattern (s) | Verlet (s) | Effort ratio incl. Pattern (−) | Effort ratio excl. Pattern (−) |
|--------:|---:|-------------:|----------:|------------:|-----------:|-------------------------------:|-------------------------------:|
| 160,801 |  2 |         3.10 |     19.60 |        3.83 |      0.234 |                          113.4 |                           97.0 |
|  40,401 |  2 |         0.80 |      4.23 |        1.00 |      0.040 |                          150.8 |                          125.8 |
|         |  3 |         4.00 |      9.16 |        3.64 |      0.108 |                          155.6 |                          121.9 |
|         |  4 |        12.30 |     19.80 |        8.56 |      0.129 |                          315.2 |                          248.8 |
|  10,201 |  2 |         0.19 |      0.52 |        0.15 |      0.012 |                           71.7 |                           59.2 |
|         |  3 |         0.90 |      1.81 |        1.95 |      0.020 |                          233.0 |                          135.5 |
|         |  4 |         2.98 |      2.31 |        2.50 |      0.071 |                          109.7 |                           74.5 |



## Newmark

```@eval
using PeriLab, Markdown
Markdown.parse(PeriLab.parameter_table("Newmark"))
```

The solver uses the Newmark solver [NewmarkNM1959](@cite).


### Matrix Update

Uses the previous time step as original configuration. This allows the analysis of geometrically non-linear deformations. It must be updated if damages or additive models are used.

!!! warning "Models"
    Not all models are fully tested yet in this framework.
