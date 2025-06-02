## Seminar 4: From bond-based to state-based I (Theory)

**Definition of states**

- A state is not in general a linear function of $\xi$ .
- A state is not in general a continuous function of $\xi$.
- The real Euclidean space $V$ is infinite-dimensional, while the real Euclidean space $\mathcal{L}_2$ (the set of second order tensors) has dimension 9

**Definitions**

*scalar state*

$\underline{a}\langle \boldsymbol{\xi}\rangle$

*vector state*

$\underline{\mathbf{A}}\langle \boldsymbol{\xi}\rangle$

*Shape tensor*

$$\mathbf{K} = \underline{\mathbf{X}}*\underline{\mathbf{X}} = \int_{\mathcal{H}}\underline{\omega}\langle \boldsymbol{\xi}\rangle\underline{\mathbf{X}}\langle \boldsymbol{\xi}\rangle\otimes\underline{\mathbf{X}}\langle \boldsymbol{\xi}\rangle dV$$

- is positive definite

![](../assets/tensors_vs_states.png)

Figure taken from [SillingSA2007](@cite)


## Constitutive Models

![](../assets/pd_states.png)


*Deformation vector state field*

$$\boldsymbol{\xi} = \mathbf{x}'-\mathbf{x}$$

$$\underline{\mathbf{Y}}[\mathbf{x},t]\langle \boldsymbol{\xi}\rangle=\mathbf{y}(\mathbf{x}+\boldsymbol{\xi},t)-\mathbf{y}(\mathbf{x},t)$$

$$\underline{x}=|\underline{\mathbf{X}}|=|\boldsymbol{\xi}|\quad\underline{y}=|\underline{\mathbf{Y}}|$$

$$\underline{e}=\underline{y}-\underline{x}=|\boldsymbol{\eta}|$$

$$\underline{y}-\underline{x}\neq|\boldsymbol{\eta}|$$
$$\underline{t}=|\underline{\mathbf{T}}|$$

**Weighted volume**

$$m_V = \int_{\mathcal{H}} \underline{\omega}\langle \boldsymbol{\xi}\rangle \underline{x} \underline{x} dV$$

**Dilatation**

$$\theta = \frac{3}{m_V} = \int_{\mathcal{H}}\underline{\omega}\langle \boldsymbol{\xi}\rangle \underline{x} \underline{e}\langle \boldsymbol{\xi}\rangle dV$$

$$\underline{t} = \frac{\omega\langle \boldsymbol{\xi}\rangle }{m_v}\left[3K \theta \underline{x} + 15G \underline{e}^d  \right]$$

*Decomposition in the devatoring and isotropic part of the strain*

$$\underline{e}^d\langle \boldsymbol{\xi}\rangle = \epsilon_{ij}^d\xi_i\frac{x_j}{|\boldsymbol{\xi}|}$$


$$\underline{e}^i\langle \boldsymbol{\xi}\rangle = \epsilon_{ij}^i\xi_i\frac{x_j}{|\boldsymbol{\xi}|}$$


The force density can be determined as

$$\underline{\mathbf{T}}=\underline{t}\frac{\underline{\mathbf{Y}}}{|\underline{\mathbf{Y}}|}$$


## Correspondence
$$\underline{\mathbf{F}}=\int_{\mathcal{H}}\underline{\omega}\langle \boldsymbol{\xi}\rangle\underline{\mathbf{Y}}\langle \boldsymbol{\xi}\rangle\otimes\underline{\mathbf{X}}\langle \boldsymbol{\xi}\rangle dV \cdot \underline{\mathbf{K}}^{-1}$$


$$\underline{\mathbf{K}}=\int_{\mathcal{H}}\underline{\omega}\langle \boldsymbol{\xi}\rangle\underline{\mathbf{X}}\langle \boldsymbol{\xi}\rangle\otimes\underline{\mathbf{X}}\langle \boldsymbol{\xi}\rangle dV$$

$$\boldsymbol{\sigma} = f(\underline{\mathbf{F}}, t, T, ...)$$

$$\underline{\mathbf{P}} = \text{det}(\underline{\mathbf{F}})\boldsymbol{\sigma}\underline{\mathbf{F}}$$


$$\underline{\mathbf{T}} = \underline{\omega}\langle \boldsymbol{\xi}\rangle\underline{\mathbf{P}}\underline{\mathbf{K}}^{-1}\mathbf{\xi}$$


## Zero-energy modes

For correspondence models, the so called zero-energy modes could occur [TupekMR2014](@cite). These modes are non-physical and lead to unstable or unreasonable solutions. Several stabilization methods were published to overcome this problem [BreitenfeldMS2014](@cite), [ChenH2018](@cite), [LiP2018](@cite), [TupekMR2014b](@cite),[WuCT2014](@cite),[WuCT2015](@cite).


$$\underline{\mathbf{T}}^C=\underline{\mathbf{T}}+\underline{\mathbf{T}}^S$$


$$\underline{\mathbf{T}}^S\langle \boldsymbol{\xi}\rangle = \underline{\omega}\langle\boldsymbol{\xi}\rangle\mathbf{C}_1\underline{\mathbf{z}}$$

$$\underline{\mathbf{z}}\langle \boldsymbol{\xi}\rangle= \underline{\mathbf{Y}}\langle\boldsymbol{\xi} \rangle-\mathbf{F}\boldsymbol{\xi}$$


$$\mathbf{C}_1=\mathbf{C}\cdot\cdot\mathbf{K}^{-1}$$
