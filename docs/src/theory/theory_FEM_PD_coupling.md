# Finite Element - Peridynamics coupling

| Module | Related Model in PeriLab |
|---|---|
| Arlequin_coupling | [Arlequin Method](https://github.com/PeriHub/PeriLab.jl/blob/main/src/FEM/Coupling/Arlequin_coupling.jl) |


To couple PD with FEM the distance between the expected PD deformation and the FEM deformation has to be computed [PernatiiA2023](@cite).

$$\mathbf{z}(\mathbf{u})=\mathbf{z}_0\mathbf{Z}\mathbf{u}$$
Assume the displacement constraints in the overlapping zone between local and non-local domains are described as:
$$\mathbf{z}=d-\sum_{\Omega}\mathbf{N}_i(\boldsymbol{\xi})\mathbf{u}_i$$
with $\mathbf{d}$ the displacement of the PD point and $\sum_{\Omega}\mathbf{N}_i(\boldsymbol{\xi})\mathbf{u}_i$ the displacement within the finite element.

With

$$\mathbf{K}_z=\kappa\begin{bmatrix}
\mathbf{I} & \mathbf{N}_d \\
\mathbf{N}^T_d & \mathbf{N}^T_d\mathbf{N}_d
\end{bmatrix}\begin{bmatrix}
\mathbf{d}  \\
\mathbf{u}
\end{bmatrix}$$

## Arlequin
Following [PernatiiA2023](@cite) for the Arlequng method the equation of motion of the coupled system in discretized form looks as follows:


$$\kappa\begin{bmatrix}
\frac{\alpha}{V_{el}}\mathbf{M}_{FE} &  \\
 & (1-\alpha)\rho_{PD}
\end{bmatrix}\begin{bmatrix}
\ddot{\mathbf{d}}  \\
\ddot{\mathbf{u}}
\end{bmatrix} + \begin{bmatrix}
\frac{\alpha}{V_{el}}\mathbf{K}_{FE} &  \\
 & (1-\alpha)\mathbf{f}_{PD}
\end{bmatrix}\begin{bmatrix}
\mathbf{d}  \\
\mathbf{u}
\end{bmatrix} + \mathbf{K}_z\begin{bmatrix}
\mathbf{d}_0  \\
\mathbf{u}_0
\end{bmatrix}=\begin{bmatrix}
\frac{\alpha}{V_{el}}\mathbf{F}_{FE}  \\
(1-\alpha)\mathbf{b}_{PD}
\end{bmatrix}$$
