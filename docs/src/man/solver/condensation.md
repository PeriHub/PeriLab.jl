## Condensation
Condensation methods can be used to reduce the number of degress of freedom (dof). The condensated dofs should not be updated and therefore no fracture or non-linear behavior should exist in this region. In theory it is possible, but it leads to continuous update of this region and is inefficient.

## Guyan Condensation
The system is partitioned into condensed $l$ and retained $r$ dof
[GuyanRJ1965](@cite):

$$\begin{equation}\begin{bmatrix}
\mathbf{K}_{rr} & \mathbf{K}_{rl} \\
\mathbf{K}_{lr} & \mathbf{K}_{ll}
\end{bmatrix}
\begin{bmatrix} \mathbf{u}_r \\ \mathbf{u}_l \end{bmatrix} =
\begin{bmatrix} \mathbf{F}_r \\ \mathbf{0} \end{bmatrix}
\end{equation}$$

Since $\mathbf{F}_l = \mathbf{0}$, the second block row gives:

$$\begin{equation}
\mathbf{K}_{lr}\mathbf{u}_r + \mathbf{K}_{ll}\mathbf{u}_l = \mathbf{0}
\quad\Rightarrow\quad
\mathbf{u}_l = \mathbf{T}\,\mathbf{u}_r, \qquad
\mathbf{T} = -\mathbf{K}_{ll}^{-1}\mathbf{K}_{lr}
\end{equation}$$

Substituting into the first block row yields:

$$\begin{equation}
\mathbf{K}_{rr}\mathbf{u}_r + \mathbf{K}_{rl}\mathbf{T}\mathbf{u}_r
= \mathbf{F}_r
\end{equation}$$

which defines the condensed system
$\hat{\mathbf{K}}_{rr}\,\mathbf{u}_r = \mathbf{F}_r$ with the
condensed stiffness matrix:

$$\begin{equation}
\hat{\mathbf{K}}_{rr} = \mathbf{K}_{rr} - \mathbf{K}_{rl}\mathbf{K}_{ll}^{-1}\mathbf{K}_{lr}
= \mathbf{K}_{rr} - \mathbf{K}_{rl}\mathbf{T}
\end{equation}$$

Following Guyan [GuyanRJ1965](@cite) for dynamic problems, the mass matrix $\mathbf{M}$, which is diagonal  in the PD discretization, is condensed analogously using the
same transformation $\mathbf{T}$:

$$\begin{equation}
\hat{\mathbf{M}}_{rr} = \mathbf{M}_{rr} + \mathbf{T}^T\mathbf{M}_{ll}\mathbf{T}
\end{equation}$$

where $\mathbf{M}_{rr}$ and $\mathbf{M}_{ll}$ are the diagonal mass
submatrices of the retained and condensed partitions, respectively. The
product $\mathbf{T}^T\mathbf{M}_{ll}\mathbf{T}$ introduces
off-diagonal coupling, so $\hat{\mathbf{M}}_{rr}$ is generally full.
The condensed dynamic system is:

$$\begin{equation}
\hat{\mathbf{M}}_{rr}\ddot{\mathbf{u}}_r +
\hat{\mathbf{K}}_{rr}\,\mathbf{u}_r = \mathbf{F}_r
\end{equation}$$

Several limitations apply. The factorization of $\mathbf{K}_{ll}$
is a one-time preprocessing cost, amortized over all load steps but
potentially significant for very large $\Omega_l$. More importantly,
$\mathbf{T}$ is computed once from the initial stiffness, so
$\Omega_l$ must remain linear elastic and undamaged throughout the
simulation. Finally, inertial effects in $\Omega_l$ are neglected;
high-frequency waves impinging on the $\Omega_l$/$\Omega_r$ interface
are partially reflected rather than correctly transmitted, restricting
the approach to quasi-static or low-frequency dynamic fracture
problems.

### PD--Matrix Coupling within $\Omega_r$
No PD bond may reach into $\Omega_l$, which is enforced by
requiring $\Omega_c$ to be at least one horizon $\delta$ wide around
$\Omega_p$:

$$\mathcal{H}_i \subseteq \Omega_r = \Omega_c \cup \Omega_p
\quad \forall\, i \in \Omega_p$$

If you use model reduction the active part is the point wise method [WillbergC2026](@cite). Fracture can easily be implemented. The equation shows how it works. You have regions $cc$ which includes all matrix parts. In the case of reduced models the active and condensed nodes (.)^c. This region couples in the material point region $pc$ and $cp$. If you want to ''cut'' parts of the matrix you have to delete the $pc$ and $cp$ parts of the material point method $(.)^p$.

```math
\begin{bmatrix} \mathbf{f}_c \\ \mathbf{f}_p \end{bmatrix} =
\underbrace{
\begin{bmatrix}
\hat{\mathbf{K}}_{cc} & \hat{\mathbf{K}}^{c}_{cp}+\hat{\mathbf{K}}^{p}_{cp} \\
\hat{\mathbf{K}}^{c}_{pc}+\hat{\mathbf{K}}^{p}_{pc} & \hat{\mathbf{K}}_{pp}
\end{bmatrix}
\begin{bmatrix} \mathbf{u}_c \\ \mathbf{u}_p \end{bmatrix}
}_{\text{Matrix part}}
+
\underbrace{
\begin{bmatrix}
\mathbf{f}_c^{\mathrm{pd}} \\
\mathbf{f}_p^{\mathrm{pd}}
\end{bmatrix}
}_{\text{Material point part}}
```
where $\hat{\mathbf{K}}^{p}_{cp}=\hat{\mathbf{K}}^{p}_{pc}=\hat{\mathbf{K}}_{pp}=\mathbf{0}$ by construction.
The distance to the reduced nodes is large enough. The figure shows that $2\delta$ should be at least the distance from the fracture.


![](../../assets/coupling_nodes.png)

- **Condensed nodes** (red, $\Omega_l$): Far-field region,
    condensed out. Must remain linear elastic, undamaged, and carry
    no external loads.
- **Matrix nodes** (light blue, $\Omega_c \subset \Omega_r$):
    Solved via the stiffness matrix. Captures load introduction,
    boundary conditions, and correct stiffness distributions. Must
    remain undamaged.
- **PD nodes** (blue, $\Omega_p \subset \Omega_r$):
    Damage region. Bond forces are evaluated via the material point
    formulation at each time step.


![w:600](../../assets/effect_of_coupling_size.png)

## Craig Bampton condensation

Guyan condensation carries the mass of the condensed region through
$\mathbf{T}^T\mathbf{M}_{ll}\mathbf{T}$, but admits no deformation of $\Omega_l$ other
than the static one prescribed by $\mathbf{T}$. The condensed region contributes inertia
without deformation modes of its own, which introduces errors at higher frequencies.
Craig-Bampton condensation [CraigRR1968](@cite) retains the static relation and adds those modes.

The partitioning into condensed $l$ and retained $r$ dof is the one used above. The
displacement field is expressed through two sets of shape functions:

$$\begin{equation}
\begin{bmatrix} \mathbf{u}_r \\ \mathbf{u}_l \end{bmatrix} =
\underbrace{\begin{bmatrix}
\mathbf{I} & \mathbf{0} \\
\boldsymbol{\Phi}_c & \boldsymbol{\Phi}_n
\end{bmatrix}}_{\mathbf{T}}
\begin{bmatrix} \mathbf{u}_r \\ \boldsymbol{\eta} \end{bmatrix}
\end{equation}$$

The constraint modes $\boldsymbol{\Phi}_c = -\mathbf{K}_{ll}^{-1}\mathbf{K}_{lr}$ are
the static response of $\Omega_l$ to a unit displacement of each retained dof and are
identical to the Guyan transformation. The fixed-interface normal modes
$\boldsymbol{\Phi}_n$ follow from the eigenvalue problem of the condensed region with all
retained dof held fixed:

$$\begin{equation}
\mathbf{K}_{ll}\boldsymbol{\Phi}_n = \mathbf{M}_{ll}\boldsymbol{\Phi}_n\boldsymbol{\Lambda},
\qquad \boldsymbol{\Lambda} = \mathrm{diag}(\omega_1^2,\dots,\omega_{n}^2)
\end{equation}$$

Only the $n$ lowest modes are retained; the modal amplitudes $\boldsymbol{\eta}$ are not
displacements. Mass normalisation of the modes,

$$\begin{equation}
\boldsymbol{\Phi}_n^T\mathbf{M}_{ll}\boldsymbol{\Phi}_n = \mathbf{I},
\qquad
\boldsymbol{\Phi}_n^T\mathbf{K}_{ll}\boldsymbol{\Phi}_n = \boldsymbol{\Lambda}
\end{equation}$$

determines the structure of the reduced matrices
$\hat{\mathbf{K}} = \mathbf{T}^T\mathbf{K}\mathbf{T}$ and
$\hat{\mathbf{M}} = \mathbf{T}^T\mathbf{M}\mathbf{T}$:

$$\begin{equation}
\hat{\mathbf{K}} =
\begin{bmatrix}
\hat{\mathbf{K}}_{rr} & \mathbf{0} \\
\mathbf{0} & \boldsymbol{\Lambda}
\end{bmatrix},
\qquad
\hat{\mathbf{M}} =
\begin{bmatrix}
\hat{\mathbf{M}}_{rr} & \boldsymbol{\Phi}_c^T\mathbf{M}_{ll}\boldsymbol{\Phi}_n \\
\boldsymbol{\Phi}_n^T\mathbf{M}_{ll}\boldsymbol{\Phi}_c & \mathbf{I}
\end{bmatrix}
\end{equation}$$

The blocks $\hat{\mathbf{K}}_{rr}$ and $\hat{\mathbf{M}}_{rr}$ are those of the Guyan
reduction; setting $n = 0$ removes the modal block and recovers it exactly. The
off-diagonal blocks of $\hat{\mathbf{K}}$ vanish because
$\mathbf{K}_{ll}\boldsymbol{\Phi}_c = -\mathbf{K}_{lr}$ makes the contributions
$\mathbf{K}_{rl}\boldsymbol{\Phi}_n$ and
$\boldsymbol{\Phi}_c^T\mathbf{K}_{ll}\boldsymbol{\Phi}_n$ cancel: constraint modes and
fixed-interface modes are $\mathbf{K}$-orthogonal. The mass exhibits no such
cancellation, so both parts remain coupled through inertia. The condensed system reads:

$$\begin{equation}
\hat{\mathbf{M}}
\begin{bmatrix} \ddot{\mathbf{u}}_r \\ \ddot{\boldsymbol{\eta}} \end{bmatrix}
+
\hat{\mathbf{K}}
\begin{bmatrix} \mathbf{u}_r \\ \boldsymbol{\eta} \end{bmatrix}
=
\begin{bmatrix} \mathbf{F}_r \\ \mathbf{0} \end{bmatrix}
\end{equation}$$

It has $n_r + n$ dof. Boundary conditions and external loads act on the $n_r$ physical
retained dof; the $n$ modal amplitudes are integrated alongside them. Loads on condensed
dof would require the projection
$\mathbf{F}_\eta = \boldsymbol{\Phi}_n^T\mathbf{F}_l$ and are not supported.

$\Omega_l$ must remain linear elastic and undamaged, since $\boldsymbol{\Phi}_c$ and
$\boldsymbol{\Phi}_n$ are computed once from the initial stiffness. The eigenvalue
problem adds to the preprocessing cost of the factorization of $\mathbf{K}_{ll}$. In
return, waves crossing the $\Omega_l$/$\Omega_r$ interface are transmitted correctly up
to approximately $\omega_n$; above that frequency they are reflected.

## Multi-Level Craig-Bampton (Cascade)

For a condensed region too large to factorize and solve the eigenvalue problem for in one
piece -- both scale with the total number of condensed dof -- the reduction can instead be
done in steps. The condensed material points are ordered as a wavefront, from the
farthest to the closest to the rest of the model (breadth first search over the bonds),
and split into `Number of Subregions` chunks $\Omega_1, \dots, \Omega_K$, every
material point with all its dof. The chunks are layers parallel to the boundary, so each
one is coupled only to its neighbouring layers, and the retained points enter the
reduction only in the last steps. $\Omega_1$ is condensed out of the full system exactly as above, which leaves a
smaller system; $\Omega_2$ is condensed out of that one, and so on. Every step
factorizes only its own subregion together with the fill-in earlier steps left in it.

The static part is exact for any partition: condensing $\Omega_1$ and then $\Omega_2$
gives the same Schur complement as condensing $\Omega_1 \cup \Omega_2$ at once.
`Number of Modes` `= 0` is therefore Guyan condensation of the whole region. How the
region is split does not matter for this, only for the memory: the fill-in of every step
is dense among the points coupled to the condensed part, which the wavefront order keeps
to one layer plus the retained boundary. The distances a reduction relies on -- from the
crack path, from loaded or constrained points -- are those of the reduction blocks to the
rest of the model and are set by the block definition.

The modes are passed on from level to level. Every level condenses the previous level's
modal coordinates together with its subregion, and computes up to `Number of Modes` new
fixed-interface modes of that combined set, with everything outside of it held. Since
the mass is no longer diagonal once a subregion has been condensed, this is a general
eigenvalue problem. Only the modes of the last level remain in the reduced system; with
`Maximum Frequency` set, every level keeps only the modes up to that frequency. The
truncation errors of the levels accumulate, so the cutoff should lie well above the
excited frequency band. For a single
subregion the cascade is identical to the single-level reduction above.
