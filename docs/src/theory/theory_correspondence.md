# Correspondence Peridynamics

The correspondence formulation is a non-ordinary state-based formulation provided by [SillingSA2007](@cite). It has the goal to apply classical models to Peridynamics.

The non-local deformation gradient is defined as

$$\underline{\mathbf{F}}=\int_{\mathcal{H}}\underline{\omega}\langle \boldsymbol{\xi}\rangle\underline{\mathbf{Y}}\otimes\underline{\mathbf{X}}dV \cdot \underline{\mathbf{K}}^{-1}$$

with the shape tensor as

$$\underline{\mathbf{K}}=\int_{\mathcal{H}}\underline{\omega}\langle \boldsymbol{\xi}\rangle\underline{\mathbf{X}}\otimes\underline{\mathbf{X}}dV$$

Based on this definition strain measures can be created to calculate the Cauchy stresses

$$\boldsymbol{\sigma} = f(\underline{\mathbf{F}}, t, T, ...)$$

To get the force densities the First-Piola Kirchhoff stress tensor has to be calculated by

$$\underline{\mathbf{P}} = \text{det}(\underline{\mathbf{F}})\boldsymbol{\sigma}\underline{\mathbf{F}}$$

and finaly the force density vector can be determined as

$$\underline{\mathbf{T}} = \underline{\omega}\langle \boldsymbol{\xi}\rangle\underline{\mathbf{P}}\underline{\mathbf{K}}^{-1}\mathbf{\xi}$$

The 2D plane strain or plane stress models are represented in the Cauchy stresses by assuming that the strain in the third direction are zero or the stresses.

## Zero-energy modes

For correspondence models, the so called zero-energy modes could occur [TupekMR2014](@cite). These modes are non-physical and lead to unstable or unreasonable solutions. Several stabilization methods were published to overcome this problem [BreitenfeldMS2014](@cite), [ChenH2018](@cite), [LiP2018](@cite), [TupekMR2014b](@cite),[WuCT2014](@cite),[WuCT2015](@cite).

A promising approach implemented as global control in PeriLab was published by Wan et al. in 2019 [WanJ2019](@cite). Instead of a bond-based stabilization method proposed by Silling [SillingSA2017](@cite) Wan et al. developed a state-based stabilization method. As positive side effect this method stabilizes the solution for anisotropic material as well. The corrected force density state $\underline{\mathbf{T}}^C$ with suppression of the zero-energy mode is:

$$\underline{\mathbf{T}}^C=\underline{\mathbf{T}}+\underline{\mathbf{T}}^S.$$

Following Wan et al. \cite{WanJ2019} the suppression force density state $\underline{\mathbf{T}}^S$ is:

$$\underline{\mathbf{T}}^S\langle \boldsymbol{\xi}\rangle = \underline{\omega}\langle\boldsymbol{\xi}\rangle\mathbf{C}_1\underline{\mathbf{z}}.$$

with $\underline{\mathbf{z}}$ as the non-uniform deformation state

$$\underline{\mathbf{z}}\langle \boldsymbol{\xi}\rangle= \underline{\mathbf{Y}}\langle\boldsymbol{\xi} \rangle-\mathbf{F}\boldsymbol{\xi}$$

caused by the zero-energy mode. If the approximated non-local deformation gradient
$\mathbf{F}$ exactly maps each undeformed bond to the deformed configuration  no zero-energy mode occur. In that case the non-uniform deformation state is zero and the corrected force density state $\underline{\mathbf{T}}^C$ is equal to the force density state $\underline{\mathbf{T}}$. The second order tensor $\mathbf{C}_1$ is given as

$$\mathbf{C}_1=\mathbf{C}\cdot\cdot\mathbf{K}^{-1},$$

utilizing the elasticity tensor.
