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

The static solver from [NLsolve.jl](https://github.com/JuliaNLSolvers/NLsolve.jl) has been included. Specifically the [method = :anderson](https://github.com/JuliaNLSolvers/NLsolve.jl#anderson-acceleration) is used.

!!! warning "Multiphysics"
    Currently only the mechanical solver is included!
