Peridynamics is an integral continuum mechanics formulation. For a pure mechanical description it can be formulated as:

$$\int_{\mathcal{H}}\underline{\mathbf{T}}\langle\mathbf{x},t\rangle-\underline{\mathbf{T}}\langle\mathbf{x}',t\rangle dV- \mathbf{b}=\rho\ddot{\mathbf{u}} $$

The parameters are:

| Parameter|Name |
|---|---|
| $\mathcal{H}$| Neighborhood [-]|
| $V$| Volume [$m^3$]|
| $\mathbf{x}$| Position of point [$m$]|
| $\mathbf{x}'$| Position of neighbor [$m$]|
| $t$| Time [$s$] |
| $\mathbf{b}$| Body force densities [$N/m^3$] |
| $\mathbf{u}$| Displacements [$m$] |
| $\ddot{\mathbf{u}}$| Accelerations [$m/s^2$] |
| $\underline{\mathbf{T}}$| Force density state [$N/m^6$] |
| $\rho$| Mass density [$kg/m^3$]|

To solve this three main types of formulations are usable; bond-based, ordinary state-based and non-ordinary state-based.

| Method | Related Model in PeriLab |
|---|---|
| Bond-based | [Bond-based Elastic](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Physics/Material/BondBased/Bondbased_Elastic.jl) |
| Ordinary state-based | [PD Solid Elastic](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Physics/Material/Material_Models/PD_Solid_Elastic.jl) [PD Solid Plastic](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Physics/Material/Material_Models/PD_Solid_Plastic.jl) |
|Non-ordinary state-based| [Correspondence Elastic](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Physics/Material/Material_Models/Correspondence_Elastic.jl) [Correspondence Plastic](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Physics/Material/Material_Models/Correspondence_Plastic.jl)|


## Bond-based


## Ordinary State-based
> More details can be found here [WillbergC2019](@cite)

For an isotropic Peridynamic solid and small deformations we can define $$\underline{x}=|\underline{\mathbf{X}}|$$ and  $$\underline{y}=|\underline{\mathbf{Y}}|$$
and 
$$\underline{e}=\underline{y}-\underline{x}=|\boldsymbol{\eta}|$$

>$\underline{y}-\underline{x}\neq|\boldsymbol{\eta}|$ for the general case
The force density scalar state can be defined as
$$\underline{t}=|\underline{\mathbf{T}}|$$

The weighted volume is
$$m_V = \int_{\mathcal{H}} \underline{\omega}\langle \boldsymbol{\xi}\rangle \underline{x} \underline{x} dV$$

The dilatation is given as 
$$\theta = \frac{3}{m_V} = \int_{\mathcal{H}}\underline{\omega}\langle \boldsymbol{\xi}\rangle \underline{x} \underline{e}\langle \boldsymbol{\xi}\rangle dV$$

$$ \underline{t} = \frac{\omega\langle \boldsymbol{\xi}\rangle }{m_v}\left[3K \theta \underline{x} + 15G \underline{e}^d  \right] $$

with the decomposition in the devatoring and isotropic part of the strain

$$ \underline{e}^d\langle \boldsymbol{\xi}\rangle = \epsilon_{ij}^d\xi_i\frac{x_j}{|\boldsymbol{\xi}|} $$ 

and 
$$ \underline{e}^i\langle \boldsymbol{\xi}\rangle = \epsilon_{ij}^i\xi_i\frac{x_j}{|\boldsymbol{\xi}|} $$ 

The force density can be determined as

$$ \underline{\mathbf{T}}=\underline{t}\frac{\underline{\mathbf{Y}}}{|\underline{\mathbf{Y}}|} $$

## Correspondence