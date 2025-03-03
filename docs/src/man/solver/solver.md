# Solver

| Parameter | Type | Optional | Description |
|---|---|---|---|
| Material Models  | Bool | Yes | Activates the time integration for materials and material evaluation |
| Damage Models    | Bool | Yes | Activates the damage evaluation |
| Thermal Models   | Bool | Yes | Activates the time integration for thermal models and thermal model evaluation |
| Additive Models  | Bool | Yes | Activates the additive model evaluation |
| Maximum Damage   | Float | Yes | Defines the maximum damage in one point |
| Initial Time | Float | No | Defines the initial time |
| Final Time | Float | No | Defines the final time |
| Fixed dt    | Float | Yes | Defines a fixed time step |
| Number of Steps   | Int | Yes | Defines a fixed number of steps |
| Verlet | Dict | Yes | Defines the Verlet solver |
| Static | Dict | Yes | Defines the Static solver |

## Verlet

| Parameter | Type | Optional | Description |
|---|---|---|---|
| Safety Factor  | Float | Yes | Defines a scaling factor for the time increment |
| Numerical Damping | Float | Yes | Defines a damping factor |


!!! warning "Fixed dt"
    If a fixed time step is defined, the time integration can become unstable.


The Verlet time integration is used as standard solver for dynamic hyperbolic differential equation of motion. It is also used in Peridigm [LittlewoodDJ2023](@cite). The displacements for step $i+1$ are solved as follows

$$\mathbf{u}_{i+1} = \mathbf{u}_{i} + \Delta t\dot{\mathbf{u}}_{i} + \frac12 \Delta t^2\ddot{\mathbf{u}}_{i}$$

with

$$\ddot{\mathbf{u}}_{i} = \frac{\mathbf{F}_i}{\rho}$$

where $\rho$ is the mass density of the point and $\mathbf{F}_i=\mathbf{F}_{external}-\mathbf{F}_{internal}$ for the current time step.

For parabolic hyperbolic time integration as in temperature models the following schema is used

$$\boldsymbol{\tau}_{i+1} =  \boldsymbol{\tau}_i - \Delta t \frac{\mathbf{H}}{\rho C_v}$$

where $\rho$ is the mass density, $C_v$ is the [specific heat capacity](https://en.wikipedia.org/wiki/Specific_heat_capacity) and $\mathbf{H}$ is the heat flux of each point [OterkusS2014b](@cite).

For the time intergration a stable increment has to be determined.

## Static



| Parameter | Type | Optional | Description |
|---|---|---|---|
| NLsolve  | Bool | Yes | Place Holder |
| Solution tolerance  | Float | Yes | Defines how much change between two iterations of the solution variable is allowed. |
| Residual tolerance  | Float | Yes |  Defines how much change between two iterations of the maximum residual variable is allowed. |
| Maximum number of iterations  | Bool | Yes | Maximum number of iteration of the solver. |
| Show solver iteration  | Bool | Yes | Shows additional information |
| Residual scaling  | Float | Yes | Scales the residual and the variable in same order, e.g. if you have linear elastic problem $u=\frac{Fl}{EA}$ your scaling should be in the order of $EA$ |
| Solver Type  | String | Yes | not implmented yet |
| m | Int | Yes | Only for Anderson solver;  It does not use Jacobian information or linesearch, but has a history whose size is controlled by the m parameter: $m=0$ corresponds to the simple fixed-point iteration above, and higher values use a larger history size to accelerate the iterations. Higher values of m usually increase the speed of convergence, but increase the storage and computation requirements and might lead to instabilities.  |
| Linear Start Value | Vector{Float} | Yes | Defines start and end values of a linear function over the length of the model (detailed explanation in the text) |
The static solver from [NLsolve.jl](https://github.com/JuliaNLSolvers/NLsolve.jl) has been included. Specifically the [method = :anderson](https://github.com/JuliaNLSolvers/NLsolve.jl#anderson-acceleration) is used.

The solver computes the residual of the internal reaction force densities and the external applied force densities
$$r =  \left[\underline{\mathbf{T}}_{external} + \underline{\mathbf{T}}_{internal}\right / s_{Residual\,scaling}$$

Right now the default value $m$ of the Anderson acceleration method is chosen.



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

$$start\_value(x,y,z (optional))= \frac{\text{max}(start_{val})-\text{min}(start_{val})}{\text{max}(coordinates)-\text{min}(coordinates)}coordinates$$
