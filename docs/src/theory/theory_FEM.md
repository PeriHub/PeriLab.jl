# Finite Element Method
The elements in PeriLab are formulated using the matrix of shape functions $\mathbf{N}$ and the matrix of derivatives $\mathbf{B}$ [Zienkiewicz2013](@cite), [WillbergC2012](@cite).

$$\mathbf{M}=\mathbf{N}^T\mathbf{N}\rho$$

To run a point wise analysis and bring it into the same structure as the PD formulation the lumped mass matrix form is used

$$m_i=\sum\limits_{j=1}^{n} M_{i,j}$$

to get a diagonal mass matrix.
To compute the stiffness matrix the following form is used.

$$\mathbf{K}=\int_V\mathbf{B}^T\boldsymbol{\sigma}(\boldsymbol{\varepsilon})dV$$

The element strain is given as

$$\boldsymbol{\varepsilon}=\mathbf{B}^T\mathbf{u}$$

This formulation allows the flexible integration of material laws. Using the linear elastic material with the elasticity matrix $\mathbf{C}$

$$\mathbf{K}=\int_V \mathbf{B}^T\mathbf{CB}dV$$

| Module | Related Model in PeriLab |
|---|---|
| Lagrange_element | [Lagrange](https://github.com/PeriHub/PeriLab.jl/blob/main/src/FEM/Element_formulation/Lagrange_element.jl) |

## Lagrange functions
Lagrange polynomials can be used to formulated finite elements
[AbramowitzStegun1983](@cite). These polynomials can be defined recursively for a polynomial $p$.

$
L(x) =  \prod\limits_{\begin{smallmatrix}0\le m\le p\\ m\neq j\end{smallmatrix}} \frac{x-x_m}{x_j-x_m}
$

The values $x$ are defined in local coordinated $[-1\,1]$. These shape functions can be defined seperatly for each direction, also with different polynomial orders. In combination these functions are used in the matrix $\mathbf{N}$.

The derivative can be computed recursively as well.

$$L_j'(x)=L_j(x)\sum\limits_{\begin{smallmatrix}i=0\\ i\neq j\end{smallmatrix}}^{p}\frac{1}{x-x_i}$$

The number of nodes per element is depended on the degrees of freedom (dof) $(p+1)^{dof}$.
