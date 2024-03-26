# Ordinary state-based Peridynamics
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
For plane stress and plane strain the equations are taken form [BobaruF2016](@cite).

TODO
