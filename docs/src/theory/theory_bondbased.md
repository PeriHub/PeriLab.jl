# Bond-based Peridynamics

Bond-based Peridynamics is a nonlocal extension of classical continuum mechanics, designed to address discontinuities like cracks within materials. Unlike traditional methods, which use partial differential equations (PDEs) and are based on local interactions, Peridynamics operates on integral equations and accounts for long-range forces.

## Fundamental Concepts

In bond-based Peridynamics, the material is considered as a continuum of particles that interact with each other through bonds. These interactions are not limited to immediate neighbors, allowing the framework to naturally handle the initiation and propagation of cracks.

### Equation of Motion

The equation of motion in Peridynamics is an integral equation, differing from the local PDEs in classical mechanics. For a particle at position $\mathbf{x}$, the equation is:

$$\rho(\mathbf{x}) \ddot{\mathbf{u}}(\mathbf{x}, t) = \int_{\mathcal{H}} \mathbf{f}(\mathbf{x}', \mathbf{x}, t) \, dV + \mathbf{b}(\mathbf{x}, t)$$

where:
$\rho(\mathbf{x})$ is the mass density at $\mathbf{x}$.
$\ddot{\mathbf{u}}(\mathbf{x}, t)$ is the acceleration at point $\mathbf{x}$.
$\mathcal{H}$ represents the horizon around $\mathbf{x}$, within which interactions occur.
$\mathbf{b}(\mathbf{x}, t)$ is the body force term.


### Interaction

The fundamental interaction in bond-based Peridynamics is between pairs of points or particles within a certain horizon distance. The force vector between two points, $x$ and $x'$, is given by:

$$\mathbf{f}(\mathbf{x}', \mathbf{x}) = \underline{\omega}\langle \boldsymbol{\xi} \rangle c \, (\mathbf{u}(\mathbf{x}') - \mathbf{u}(\mathbf{x}))$$

where:
$\mathbf{f}(\mathbf{x}', \mathbf{x})$ is the force vector exerted by the particle at $\mathbf{x}'$ on the particle at $\mathbf{x}$.
$c$ is a bond modulus representing the stiffness of the bond.
$\underline{\omega}\langle \boldsymbol{\xi} \rangle$ is a bond-associated influence function.
$\mathbf{u}(\mathbf{x})$ is the displacement of the particle at $\mathbf{x}$.

### Bond Moduli
From [TrageserJ2020](@cite) we get

| Dimension | Bond Stiffness | Poisson's ratio
|---|---|---|
|plane strain:| $c = \frac{48 E}{\pi 5 \delta^3}$ | fixed $\nu=0.25$|
|plane stress:| $c = \frac{9 E}{\pi \delta^3}$ | fixed $\nu=1/3$|
|3D:| $c = \frac{12 E}{\pi  \delta^4}$ | fixed $\nu=0.25$|

## Unified Bond-based Peridynamics
This model was developed by Guan et al. [GuanJ2024](@cite) and utilizes the strain identification technique (SIT) [GuanJ2023](@cite). It extends the original bond-based theory and extends it. Poisson's ratios which differs from the proposed ones are allowed.
For 3D they defined

$$\mathbf{f}=c\left[\frac{3  (4\nu-1) }{2 (1 + \nu)}\varepsilon_m- \frac{5  (1 - 2\nu) }{2 (1 + \nu)} s  \right]\frac{\boldsymbol{\xi}}{|\boldsymbol{\xi}|}$$
$\varepsilon_m=\frac{\varepsilon}{3}$

where $\varepsilon$ is the bond strain, $s$ is the bond stretch and $c$ is the bond stiffness from the bond based formulation.

For 2D plane stress and plane stress the equations are more complex

$$\mathbf{f}_{2D}=c_2\left[R_as_{ij}+R_b\varepsilon_m\right]\frac{\xi_{ij}}{|\xi_{ij}|}$$

Several parameter have to be computed

$$c_2=\left\{\begin{aligned} &\frac{6}{\pi  \delta  (1 - 2  \nu)(1 + \nu)} & & \text{plane strain}\\
&\frac{6E}{\pi\delta^3(1-\nu)}& & \text{plane stress}
\end{aligned}\right.$$


$$R_a=\left\{\begin{aligned} &2(1 + \nu) \beta I_1 & & \text{plane strain}\\
& \frac{2  (1 - \nu)}{(1 - 2\nu)} \beta  I_1& & \text{plane stress}
\end{aligned}\right.$$



$$R_b=\left\{\begin{aligned} &2 (1 + \nu)  (1 - \beta) I_2 & & \text{plane strain}\\
& \frac{6\nu (1 - \nu)}{(1 - 2 \nu)^2}\beta I_1 +
         2  \frac{1 - \nu}{1 - 2 \nu}\left(1 - \frac{1 + \nu}{1 - 2 \nu} \beta\right) I_2& & \text{plane stress}
\end{aligned}\right.$$

with
$$\beta = \frac{5(1-2\nu)}{2(1+\nu)},$$

$$\begin{aligned}I_1 = n  \sqrt{1 - n^2}&&\text{and}&& I_2 = n  a  \sinh(\sqrt{1 / n^2 - 1})\end{aligned}$$
where $n$, $\xi_{ik}$
and $s_{ij}$ are given as
$$n = \frac{\xi}{\delta},$$
$$\xi_{ik}=(x_{i,k},y_{i,k},z_{i,k}).$$

$$s_{ij}=\varepsilon_x\frac{x_{ij}^2}{|\xi_{ik}|^2} + \varepsilon_z\frac{z_{ij}^2}{|\xi_{ik}|^2},$$


The meaning of $a$ is not specified by Guan et al. [GuanJ2024](@cite) and set to one. As can be seen in the equations $E$ and $\nu$ can be specified independently.
