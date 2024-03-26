# Bond-based Peridynamics

Bond-based Peridynamics is a nonlocal extension of classical continuum mechanics, designed to address discontinuities like cracks within materials. Unlike traditional methods, which use partial differential equations (PDEs) and are based on local interactions, Peridynamics operates on integral equations and accounts for long-range forces.

## Fundamental Concepts

In bond-based Peridynamics, the material is considered as a continuum of particles that interact with each other through bonds. These interactions are not limited to immediate neighbors, allowing the framework to naturally handle the initiation and propagation of cracks.

### Equation of Motion

The equation of motion in Peridynamics is an integral equation, differing from the local PDEs in classical mechanics. For a particle at position $\mathbf{x}$, the equation is:

$$ \rho(\mathbf{x}) \ddot{\mathbf{u}}(\mathbf{x}, t) = \int_{\mathcal{H}} \mathbf{f}(\mathbf{x}', \mathbf{x}, t) \, dV + \mathbf{b}(\mathbf{x}, t) $$

where:
$\rho(\mathbf{x})$ is the mass density at $\mathbf{x}$.
$\ddot{\mathbf{u}}(\mathbf{x}, t)$ is the acceleration at point $\mathbf{x}$.
$\mathcal{H}$ represents the horizon around $\mathbf{x}$, within which interactions occur.
$\mathbf{b}(\mathbf{x}, t)$ is the body force term.


### Interaction

The fundamental interaction in bond-based Peridynamics is between pairs of points or particles within a certain horizon distance. The force vector between two points, $x$ and $x'$, is given by:

$$ \mathbf{f}(\mathbf{x}', \mathbf{x}) = \underline{\omega}\langle \boldsymbol{\xi} \rangle c \, (\mathbf{u}(\mathbf{x}') - \mathbf{u}(\mathbf{x})) $$

where:
$\mathbf{f}(\mathbf{x}', \mathbf{x})$ is the force vector exerted by the particle at $\mathbf{x}'$ on the particle at $\mathbf{x}$.
$c$ is a bond modulus representing the stiffness of the bond.
$\underline{\omega}\langle \boldsymbol{\xi} \rangle$ is a bond-associated influence function.
$\mathbf{u}(\mathbf{x})$ is the displacement of the particle at $\mathbf{x}$.

### Bond Moduli
From [TrageserJ2020](@cite) we get

plane strain: $c = \frac{48 E}{\pi 5 \delta^3}$ with a fixed $\nu=0.25$
plane stress: $c = \frac{9 E}{\pi \delta^3}$ with a fixed $\nu=1/3$
3D: $c = \frac{18 K}{\pi  \delta^4}$ with a fixed $\nu=0.25$

         