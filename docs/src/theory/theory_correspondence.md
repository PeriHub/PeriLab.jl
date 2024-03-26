# Correspondence Peridynamics

The correspondence formulation is a non-ordinary state-based formulation provided by [SillingSA2007](@cite). It has the goal to apply classical models to Peridynamics.

The non-local deformation gradient is defined as
$$ \underline{\mathbf{F}}=\int_{\mathcal{H}}\underline{\omega}\langle \boldsymbol{\xi}\rangle\underline{\mathbf{Y}}\otimes\underline{\mathbf{X}}dV \cdot \underline{\mathbf{K}}^{-1}$$

with the shape tensor as

$$ \underline{\mathbf{K}}=\int_{\mathcal{H}}\underline{\omega}\langle \boldsymbol{\xi}\rangle\underline{\mathbf{X}}\otimes\underline{\mathbf{X}}dV $$

Based on this definition strain measures can be created to calculate the Cauchy stresses

$$ \boldsymbol{\sigma} = f(\underline{\mathbf{F}}, t, T, ...) $$

To get the force densities the First-Piola Kirchhoff stress tensor has to be calculated by
$$ \underline{\mathbf{P}} = \text{det}(\underline{\mathbf{F}})\boldsymbol{\sigma}\underline{\mathbf{F}}  $$
and finaly the force density vector can be determined as 
$$ \underline{\mathbf{T}} = \underline{\omega}\langle \boldsymbol{\xi}\rangle\underline{\mathbf{P}}\underline{\mathbf{K}}^{-1}\mathbf{\xi} $$

The 2D plane strain or plane stress models are represented in the Cauchy stresses by assuming that the strain in the third direction are zero or the stresses.