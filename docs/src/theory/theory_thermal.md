# Thermal models
The theory is taken from [WillbergC2024b](@cite)
## Thermo-mechanics
To introduce a thermo-mechanical coupling the mechanical strains in \autoref{eq:GreenLagrangeStrains} have to be expanded with the thermal strains:
$$
 \boldsymbol{\varepsilon}=\boldsymbol{\varepsilon}_{mechanical} + \boldsymbol{\varepsilon}_{thermal}.

$$
The thermal strains are defined as
$$
 \boldsymbol{\varepsilon}_{thermal} =- \boldsymbol{\alpha}\tau
$$
with $\tau$ as the temperature increment and $\boldsymbol{\alpha}$ as matrix of the coefficients of thermal expansion. Typically, this matrix is diagonal. The coupled stresses for linear elastic material utilizing [Hook's law](https://en.wikipedia.org/wiki/Hooke%27s_law) is given as
$$
\boldsymbol{\sigma}=\mathbf{C}\cdot\cdot\left(\boldsymbol{\varepsilon}_{mechanical} - \boldsymbol{\alpha}\tau \right)
$$

## Thermal flux
The mechanical response due to temperature changes is included in the Peridynamics model. However, the heat flux must be included as well.  Under the assumption that mechanical deformations do not change the temperature, the thermodynamic equilibrium equation can be studied separately.
$$
\rho C_v\dot{\tau} = \int_{\mathcal{H}}(\underline{h}(\textbf{x},t)\langle\boldsymbol{\xi}\rangle-\underline{h}(\textbf{x}',t)\langle\boldsymbol{\xi}'\rangle)dV_{\textbf{x}}+ S_i
$$
The parameters are $\rho$ the mass density, $C_v$ the specific heat capacity, $\dot{\tau}$ the temperature gradient in time, $dV_{\textbf{x}}$ the volume and $S_i$ the heat sink or heat source.
The heat flux of a bond is defined as
$$
\underline{h}(\textbf{x},t)\langle\boldsymbol{\xi}\rangle = \mathbf{q}^T\mathbf{K}^{-1}(\textbf{x})\boldsymbol{\xi}
$$
with $\mathbf{q}$ as classical heat flux and $\mathbf{K}$ as the shape tensor. It follows
$$
\nabla\cdot\mathbf{q} = \int_{\mathcal{H}}\left[\mathbf{q}(\textbf{x}')^T\mathbf{K}^{-1}(\textbf{x}')+\mathbf{q}(\textbf{x})^T\mathbf{K}^{-1}(\textbf{x})\right]\boldsymbol{\xi}dV_{\textbf{x}}
$$
which can be derived utilizing the spatial gradient of the temperature $\nabla\tau$ as

$$
\mathbf{q} = -\boldsymbol{\lambda}\nabla\tau
$$

$\boldsymbol{\lambda}$ is the $3\times3$ matrix of thermal conductivity. Typically it is a diagonal matrix.
Following \cite{BrighentiR2021} the spatial temperature gradient $\nabla\tau$ can be derived as

$$
\nabla\tau = \mathbf{K}^{-1}\int_{\mathcal{H}}\left[\tau(\mathbf{x}')-\tau(\mathbf{x})\right]\boldsymbol{\xi}\underline{\omega}\langle\boldsymbol{\xi}\rangle dV_{\mathbf{x}}
$$

The numerical solving process is then

$$
\rho C_v \frac{\tau^{t+dt}-\tau^{t}}{dt}=\nabla\mathbf{q} + S_i
$$

$$
\tau^{t+dt} = dt\frac{\nabla\mathbf{q} + S_i}{\rho C_v} + \tau^{t}
$$

## Heat transfer to environment

Following \cite{GuX2019c,OterkusS2014b} the heat volumetric density at the surface for a assigned heat flux normal to the surface $q_{bc}$ is:
$$
S_i = \frac{q_{bc}}{\Delta}
$$
where $\Delta$ can be set to $dx$.
Thereby, $q_{bc}$ is
$$
q_{bc} = \kappa (\tau-\tau_{env})
$$
where $\kappa$ is the heat convection coefficient between solid and environment and $\tau_{env}$ the environmental temperature. For a mesh free model the question arises how the outer surface and the corresponded surface can be identified. For the outer surface the Peridynamics neighborhood $\mathcal{H}$ is utilized. It is assumed that is circle for 2D and a sphere for 3D.
Therefore, the following criteria has to be fulfilled for 2D

$$
V_{2D}=2\pi\delta^2 h \geq \int_{\mathcal{H}}dV
$$
and 3D
$$
V_{3D}=\frac43\pi\delta^3 \geq \int_{\mathcal{H}}dV
$$

Each point which is next to the surface will have less volume represented be the discrete material points. Defining a limit value
$$
f_{limit} \leq V_{specific} =  \frac{\int_{\mathcal{H}}dV}{V_{2D\,or\,3D}}
$$
allows an easy identification of surface nodes $i$ during the printing process. Combining \autoref{eq:newTemp} and \autoref{eq:Heat_transfer_to_environment} allows the calculation of the change in temperature for these nodes $i$ as
$$
\tau_i^{t+dt} = dt\frac{\nabla\mathbf{q}_i + \frac{\kappa (\tau_i^{t}-\tau_{env})}{dx}}{(\rho C_v)_i} + \tau^{t}_i\,.
$$

## Time step
The minimum time step for the explicit time integration of the temperature field to obtain a stable solution is given by
$$
    \Delta t < \text{min}\left(\frac{\left(\rho C_v\right)_i}{\sum_{j=1}^{N}\frac{\text{max}(\text{eig}(\boldsymbol{\lambda}))}{|\mathbf{\xi}_{ij}|}V_j}\right)
$$
with $N$ the number of neighbors of point $i$ [OterkusS2014](@cite).
