# Non-ordinary state based Peridynamics

## Correspondence Peridynamics

The correspondence formulation is a non-ordinary state-based formulation provided by [SillingSA2007](@cite). It has the goal to apply classical models to Peridynamics.

The non-local deformation gradient is defined as

$$\mathbf{F}=\int_{\mathcal{H}}\underline{\omega}\langle \boldsymbol{\xi}\rangle\underline{\mathbf{Y}}\langle \boldsymbol{\xi}\rangle\otimes\underline{\mathbf{X}}\langle \boldsymbol{\xi}\rangle dV \cdot \mathbf{K}^{-1}$$

with the positive definite shape tensor as

$$\mathbf{K}=\int_{\mathcal{H}}\underline{\omega}\langle \boldsymbol{\xi}\rangle\underline{\mathbf{X}}\langle \boldsymbol{\xi}\rangle\otimes\underline{\mathbf{X}}\langle \boldsymbol{\xi}\rangle dV$$

!!! info "Positive definiteness in numerics"
    In numerical applications if bonds break the shape tensor is positive semi definite. $\det\mathbf{K}=0$ can occur and the inversion of the shape tensor won't work.

Based on this definition strain measures can be created to calculate the Cauchy stresses

$$\boldsymbol{\sigma} = f(\mathbf{F}, t, T, ...)$$

To get the force densities the First-Piola Kirchhoff stress tensor has to be calculated by

$$\mathbf{P} = \text{det}(\mathbf{F})\boldsymbol{\sigma}\mathbf{F}^{-T}$$

and finaly the force density vector can be determined as

$$\underline{\mathbf{T}} = \underline{\omega}\langle \boldsymbol{\xi}\rangle\mathbf{P}\mathbf{K}^{-1}\mathbf{\xi}$$

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
$\mathbf{F}$ exactly maps each undeformed bond to the deformed configuration no zero-energy mode occur. In that case the non-uniform deformation state is zero and the corrected force density state $\underline{\mathbf{T}}^C$ is equal to the force density state $\underline{\mathbf{T}}$. The second order tensor $\mathbf{C}_1$ is given as

$$\mathbf{C}_1=\mathbf{C}\cdot\cdot\mathbf{K}^{-1},$$

utilizing the elasticity tensor.

## Bond-associated correspondence Peridynamics

## Matrix based approach

!!! info "Redefinition of shape tensor"
    For the derivation the shape tensor is defined as $\mathbf{D}$.

In the discretized form, the shape tensor for material point $i$ is computed as the weighted sum over all neighbors:
$$\mathbf{D}_i = \sum_{j \in \mathcal{H}_i} \underline{\omega}_{ij} V_j \underline{\mathbf{X}}_{ij} \otimes \underline{\mathbf{X}}_{ij}$$
where $\underline{\omega}_{ij}$ is the influence function, $V_j$ is the volume of material point $j$, and $\underline{\mathbf{X}}_{ij} = \mathbf{x}_j - \mathbf{x}_i$ is the undeformed bond vector state.

The deformation gradient tensor for material point $i$ is calculated using the discretized correspondence relation:

$$\mathbf{F}_i = \mathbf{I} + \left[ \sum_{j \in \mathcal{H}_i} \underline{\omega}_{ij} V_j \underline{\mathbf{U}}_{ij} \otimes \underline{\mathbf{X}}_{ij} \right] \mathbf{D}_i^{-1}$$

where $\mathbf{I}$ is the identity tensor and $\underline{\mathbf{U}}_{ij} = \mathbf{u}_j - \mathbf{u}_i$ is the displacement vector state between material points $i$ and $j$.

Under the assumption of small deformations, the linearized strain tensor is computed from the symmetric part of the displacement gradient:

$$\boldsymbol{\varepsilon}_i = \frac{1}{2}(\mathbf{F}_i + \mathbf{F}_i^T) - \mathbf{I}$$

Substituting leads to:

$$\boldsymbol{\varepsilon}_i = \frac{1}{2} \sum_{k \in \mathcal{H}_i} V_k\underline{\omega}_{ik} \left[\left( \underline{\mathbf{U}}_{ik} \otimes \underline{\mathbf{X}}_{ik}\right) \mathbf{D}_i^{-1} + \mathbf{D}_i^{-T} \left(  \underline{\mathbf{X}}_{ik} \otimes \underline{\mathbf{U}}_{ik}\right)\right]$$

To facilitate the matrix formulation, a strain-displacement operator is introduced. This operator relates the displacement field to the strain tensor through linear relationships.

Two auxiliary vectors are defined for each bond $ik$:

$$\mathbf{B}_{1,ik} = \underline{\mathbf{X}}_{ik}^T \mathbf{D}_i^{-1}$$

and

$$\mathbf{B}_{2,ik} = \mathbf{D}_i^{-T} \underline{\mathbf{X}}_{ik}$$

The strain-displacement operator $\mathbf{B}_{ik}$ is constructed as a third-order tensor with components:

$$B_{ik,mnp} = \frac{1}{2}\left(\delta_{mp} B_{1,ik,n} + \delta_{np} B_{2,ik,m}\right)$$

where $\delta_{ij}$ is the Kronecker delta, and $m,n,p$ are spatial indices. This tensor formulation allows the strain computation to be expressed compactly in index notation:

$$\varepsilon_{i,mn} = \sum_{k \in \mathcal{H}_i} V_k\underline{\omega}_{ik} B_{ik,mno} \underline{U}_{ik,o} = \sum_{k \in \mathcal{H}_i} V_k\underline{\omega}_{ik} B_{ik,mno}(U_{k,o}-U_{i,o})$$

The bond force density vector is computed using the linearized constitutive relation:

$$\begin{aligned}
\underline{\mathbf{T}}_{ij} &= \underline{\omega}_{ij}\boldsymbol{\sigma}_i \mathbf{D}_i^{-1}\underline{\mathbf{X}}_{ij}  \\
&= \underline{\omega}_{ij}\left(\mathbf{C}_i:\boldsymbol{\varepsilon}_i\right) \mathbf{D}_i^{-1}\underline{\mathbf{X}}_{ij}  \\
 &=  \underline{\omega}_{ij}\underline{\mathbf{X}}_{ij}^T \mathbf{D}_i^{-T} (\mathbf{C}_i:\boldsymbol{\varepsilon}_i)^T\\
 &=  \underline{\omega}_{ij}\underline{\mathbf{X}}_{ij}^T \mathbf{D}^{-1}_i (\mathbf{C}_i:\boldsymbol{\varepsilon}_i)
\end{aligned}$$

The stiffness matrix components are derived separating the displacement vector from the strain tensor. To do so, the tensor contraction $\mathbf{C}_i : \mathbf{B}_{ik}$ in index notation becomes:

$$[\mathbf{C} : \mathbf{B}_{ik}]_{mnq} = C_{i,mnop} B_{ik,opq}$$

Combining the strain-displacement relationship with the force calculation from the stiffness matrix components are:

$$\mathbf{K}_{ij} = -\underline{\omega}_{ij} V_j V_k \underline{\omega}_{ik} \mathbf{D}_i^{-1} \underline{\mathbf{X}}_{ij} [\mathbf{C} : \mathbf{B}_{ik}]$$

with the explicit index notation:

$$K_{ij,mo} = -\underline{\omega}_{ij} V_j V_k \underline{\omega}_{ik} \sum_{n,p} [\mathbf{C} : \mathbf{B}_{ik}]_{mno} [D_i^{-1}]_{np} X_{ij,p}$$

### Zero energy mode compensation

The global control [WanJ2019](@cite) is introcuded in matrix form.

The corrected force density state $\underline{\mathbf{T}}^C$ combines the original correspondence force with a stabilization term:

$$\underline{\mathbf{T}}^C=\underline{\mathbf{T}}+\underline{\mathbf{T}}^S$$

where $\underline{\mathbf{T}}$ is the original correspondence force density state and $\underline{\mathbf{T}}^S$ is the suppression force density state. Following Wan et al. [WanJ2019](@cite), the suppression force density state is defined as:

$$\underline{\mathbf{T}}^S\langle \boldsymbol{\xi}\rangle =\underline{\omega}\langle\boldsymbol{\xi}\rangle\mathbf{Z}\underline{\mathbf{z}}$$

where $\mathbf{Z}$ is the zero-energy stiffness tensor and $\underline{\mathbf{z}}$ is the non-uniform deformation state.

The non-uniform deformation state $\underline{\mathbf{z}}$ quantifies the deviation between the actual deformed configuration and the configuration predicted by the correspondence deformation gradient:

$$\underline{\mathbf{z}}\langle \boldsymbol{\xi}\rangle=\underline{\mathbf{Y}}\langle\boldsymbol{\xi} \rangle-\mathbf{F}\boldsymbol{\xi}$$

The second-order zero-energy stiffness tensor $\mathbf{Z}$ is constructed using the elasticity tensor from the constitutive relation:

$$\mathbf{Z}=\mathbf{C}:\mathbf{D}^{-1}$$

where $\mathbf{C}$ is the elasticity tensor and $\mathbf{D}^{-1}$ is the inverse shape tensor.

For the discrete implementation, the non-uniform deformation state for bond $ij$ is expressed in terms of displacement differences:

$$\begin{aligned}\underline{\mathbf{z}}_{ij} &= \underline{\mathbf{Y}}_{ij}-\mathbf{F}\underline{\mathbf{X}}_{ij}\\&= \underline{\mathbf{X}}_{ij}+\underline{\mathbf{U}}_{ij}-\mathbf{F}\underline{\mathbf{X}}_{ij}\\&= \underline{\mathbf{U}}_{ij} - \left( \sum_{k \in \mathcal{H}_i} \underline{\omega}_{ik} V_k \underline{\mathbf{U}}_{ik} \otimes \underline{\mathbf{X}}_{ik}\right) \mathbf{D}_i^{-1}\underline{\mathbf{X}}_{ij}\\&= \underline{\mathbf{U}}_{ij} - \sum_{k \in \mathcal{H}_i} \underline{\omega}_{ik} V_k \underline{\mathbf{U}}_{ik} (\underline{\mathbf{X}}_{ik}^T \mathbf{D}_i^{-1}\underline{\mathbf{X}}_{ij})\\&=\underline{\mathbf{U}}_{ij}\left[\mathbf{I} - \sum_{k \in \mathcal{H}_i} \underline{\omega}_{ik} V_k (\underline{\mathbf{X}}_{ik}^T \mathbf{D}_i^{-1}\underline{\mathbf{X}}_{ij})\right]\end{aligned}$$

$$\begin{aligned}\underline{\mathbf{T}}^S_{ij} &= \mathbf{Z}_i\underline{\mathbf{U}}_{ij}\left[\mathbf{I}- \sum_{k \in \mathcal{H}_i} \underline{\omega}_{ik} V_k (\underline{\mathbf{X}}_{ik} \otimes \mathbf{D}_i^{-1}\underline{\mathbf{X}}_{ij})\right]\\&=(\mathbf{C}_i:\mathbf{D}_i^{-1})\underline{\mathbf{U}}_{ij}\left[\mathbf{I}-\sum_{k \in \mathcal{H}_i} \underline{\omega}_{ik} V_k (\underline{\mathbf{X}}_{ik}^T \mathbf{D}_i^{-1}\underline{\mathbf{X}}_{ij})\right]\end{aligned}$$

The stabilization terms can be directly incorporated into the matrix formulation. The additional stiffness matrix components arising from the zero-energy mode compensation are:

$$\mathbf{K}^S_{ij}  = -\left[\mathbf{I} - \sum_{k \in \mathcal{H}_i} \underline{\omega}_{ik} V_k (\underline{\mathbf{X}}_{ik}^T \mathbf{D}_i^{-1}\underline{\mathbf{X}}_{ij})\right](\mathbf{C}:\mathbf{D}_i^{-1})V_i$$

The corresponding transpose terms are:

$$\mathbf{K}^S_{ji} =  -\left[\mathbf{I} - \sum_{k \in \mathcal{H}_i} \underline{\omega}_{ik} V_k (\underline{\mathbf{X}}_{ik}^T \mathbf{D}_i^{-1}\underline{\mathbf{X}}_{ij})\right](\mathbf{C}:\mathbf{D}_i^{-1}) V_j$$

The diagonal terms ensure equilibrium by summing all off-diagonal contributions:

$$\mathbf{K}^S_{ii} = -\sum_{\substack{\ell \in \mathcal{H}_i \\ \ell \neq i}} \mathbf{K}^S_{i\ell}$$

The advantage of this stabilization approach lies in its seamless integration with the matrix-based formulation. The stabilization stiffness terms $\mathbf{K}^S$ are simply added to the correspondence stiffness matrix $\mathbf{K}$, yielding a stabilized system:

$$(\mathbf{K} + \mathbf{K}^S)\mathbf{u} = \mathbf{F}_{ext}$$
