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
| Numerical Damping | Float | Yes | Defines a damping factor |
| Verlet | Dict | Yes | Defines the Verlet solver |

## Verlet

| Parameter | Type | Optional | Description |
|---|---|---|---|
| Safety Factor  | Float | Yes | Defines a scaling factor for the time increment |
| Fixed dt    | Float | Yes | Defines a fixed time step |
| Number of Steps   | Int | Yes | Defines a fixed number of steps |

>Note: If a fixed time step is defined, the time integration can become unstable.


The Verlet time integration is used as standard solver for dynamic hyperbolic differential equation of motion. It is also used in Peridigm [LittlewoodDJ2023](@cite). The displacements for step $i+1$ are solved as follows

$$ \mathbf{u}_{i+1} = \mathbf{u}_{i} + \Delta t\dot{\mathbf{u}}_{i} + \frac12 \Delta t^2\ddot{\mathbf{u}}_{i} $$

with 

$$ \ddot{\mathbf{u}}_{i} = \frac{\mathbf{F}_i}{\rho} $$

where $\rho$ is the mass density of the point and $\mathbf{F}_i=\mathbf{F}_{external}-\mathbf{F}_{internal}$ for the current time step.

For parabolic hyperbolic time integration as in temperature models the following schema is used

$$ \boldsymbol{\tau}_{i+1} =  \boldsymbol{\tau}_i - \Delta t \frac{\mathbf{H}}{\rho C_v} $$

where $\rho$ is the mass density, $C_v$ is the [specific heat capacity](https://en.wikipedia.org/wiki/Specific_heat_capacity) and $\mathbf{H}$ is the heat flow of each point [OterkusS2014b](@cite).

For the time intergration a stable increment has to be determined.