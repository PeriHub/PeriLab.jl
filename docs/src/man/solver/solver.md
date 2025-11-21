# Solver

| Parameter       | Type  | Optional | Description                                                                    |
| --------------- | ----- | -------- | ------------------------------------------------------------------------------ |
| Material Models | Bool  | Yes      | Activates the time integration for materials and material evaluation           |
| Damage Models   | Bool  | Yes      | Activates the damage evaluation                                                |
| Thermal Models  | Bool  | Yes      | Activates the time integration for thermal models and thermal model evaluation |
| Additive Models | Bool  | Yes      | Activates the additive model evaluation                                        |
| Maximum Damage  | Float | Yes      | Defines the maximum damage in one point                                        |
| Initial Time    | Float | No       | Defines the initial time                                                       |
| Final Time      | Float | No       | Defines the final time                                                         |
| Fixed dt        | Float | Yes      | Defines a fixed time step                                                      |
| Number of Steps | Int   | Yes      | Defines a fixed number of steps                                                |
| Verlet          | Dict  | Yes      | Defines the Verlet solver                                                      |
| Static          | Dict  | Yes      | Defines the Static solver                                                      |

## Verlet

| Parameter         | Type  | Optional | Description                                     |
| ----------------- | ----- | -------- | ----------------------------------------------- |
| Safety Factor     | Float | Yes      | Defines a scaling factor for the time increment |
| Numerical Damping | Float | Yes      | Defines a damping factor                        |

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

| Parameter                    | Type          | Optional | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| ---------------------------- | ------------- | -------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| NLsolve                      | Bool          | Yes      | Place Holder                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| Solution tolerance           | Float         | Yes      | Defines how much change between two iterations of the solution variable is allowed.                                                                                                                                                                                                                                                                                                                                                                                       |
| Residual tolerance           | Float         | Yes      | Defines how much change between two iterations of the maximum residual variable is allowed.                                                                                                                                                                                                                                                                                                                                                                               |
| Maximum number of iterations | Int           | Yes      | Maximum number of iteration of the solver.                                                                                                                                                                                                                                                                                                                                                                                                                                |
| Show solver iteration        | Bool          | Yes      | Shows additional information                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| Residual scaling             | Float         | Yes      | Scales the residual and the variable in same order. Should be in the range of the Young's modulus.                                                                                                                                                                                                                                                                                                                                                                        |
| Solver Type                  | String        | Yes      | not implmented yet                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| m                            | Int           | Yes      | Only for Anderson solver; It does not use Jacobian information or linesearch, but has a history whose size is controlled by the m parameter: $m=0$ corresponds to the simple fixed-point iteration above, and higher values use a larger history size to accelerate the iterations. Higher values of m usually increase the speed of convergence, but increase the storage and computation requirements and might lead to instabilities. $m=15$ is set as standard value. |
| Linear Start Value           | Vector{Float} | Yes      | Defines start and end values of a linear function over the length of the model (detailed explanation in the text)                                                                                                                                                                                                                                                                                                                                                         |

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

| Parameter     | Type | Optional | Description                                  |
| ------------- | ---- | -------- | -------------------------------------------- |
| Update Matrix | Bool | Yes      | Activates the update of the stiffness matrix |

The solver computes the stiffness matrix of the problem. It is solved by

$$ \mathbf{u} = \mathbf{K}_{PD}^{-1}\mathbf{F}_{external}$$

The number of steps defines the virtual time step.

### Update Matrix

Uses the previous time step as original configuration. This allows the analysis of geometrically non-linear deformations. It must be updated if damages or additive models are used.

!!! warning "Models"
Not all models are fully tested yet in this framework.

!!! warning "Models"
Stress computations are not inclueded yet. Please check the issues.

## Verlet Matrix Based

It is the same solver as the Verlet based solver above. The main difference is, that the matrix style is used.

!!! warning "Models"
If update is active it is not very efficient, because the creation of new matrix is more costly than the material point approach.

| Parameter       | Type | Optional | Description                                  |
| --------------- | ---- | -------- | -------------------------------------------- |
| Update Matrix   | Bool | Yes      | Activates the update of the stiffness matrix |
| Model reduction | Dict | Yes      | Defines the model reduction options          |

If model reduction

| Parameter | Type   | Optional | Description                                              |
| --------- | ------ | -------- | -------------------------------------------------------- |
| Type      | String | No       | Defines the type of model redction (Gyan, Craig-Bampton) |
| tbd       |        |          |                                                          |
