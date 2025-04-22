## Lecture 1: Bond based truss solution
For a truss everything is a scalar. The norms are given for better illustration.

$i$ - point

$j$ - neighbor

$\xi = x_j-x_i$

$\eta = u_j - u_i$

__Bond force__

$$f_{ij} = c\omega_{ij} s\frac{\xi+\eta}{||\xi + \eta||} = c\omega_{ij} s$$
$$s = \frac{||\xi + \eta|| - ||\xi||}{||\xi||} = \frac{\eta}{\xi}$$


$$\rho u_i = \frac{1}{2}\sum_j (f_{ji}-f_{ij})V_j+b_i$$

point based is solved like this

__Matrix based__

$$\begin{bmatrix}
K_{11} & \cdots & K_{1j} & 0& ...&0\\
\vdots & \ddots & \\
K_{j1} &  & K_{jj} &&&\\
0 & && \\
\vdots & & &&\ddots\\
0& & \cdots&&&K_{nn}
\end{bmatrix}\begin{bmatrix}u_1\\\vdots\\u_j\\u_{j+1}\\\vdots \\u_n\end{bmatrix} = \mathbf{F}_{internal}$$

$K_{ij} = -\frac{c}{|\xi_{ij}|} \omega_{ij}V_j \quad \text{for } i \ne j,\quad K_{ii} = -\sum_{j \ne i} K_{ij}; i=1\,...\,n \quad \text{and} \quad j=i+1\,...\,n_{neighbors,i}+1$

$\mathbf{F}_{internal}$ is the force density in $\left[\frac{N}{m^3}\right]$


1D [MasoumiA2023](@cite) $\rightarrow$ $c = \frac{2E}{A\delta^2}$

For FEM

$$\mathbf{K}=\frac{EA}{L}\begin{bmatrix}1&-1\\-1&1\end{bmatrix}$$


---

**Assumptions**

 $\delta = L,\, n_{neighbors,i}=1,\, \omega_{ij}=1,\, V_i=V=const.$

$V=AL$


!!! info "Volume"
    Must be $AL$, because the sum of the neighborvolume must represent the whole volume of the neighborhood, which is $AL$.



$c=0.5\frac{2E}{AL^2}$

$K_{11}=-K_{12}=-0.5\frac{c}{\xi} V=-0.5\frac{c}{L} V=-\frac{E}{AL^3}AL=-\frac{E}{L^2}$



$$\mathbf{K}=\frac{E}{L^2}\begin{bmatrix}1&-1\\-1&1\end{bmatrix}$$

bring it in local form of the stiffness matrix the forces and not the force densities the stiffness has to be multiplied by $V$.


$^{forces}K_{11}=-^{forces}K_{12}=-\frac{E}{L^2}AL=\frac{EA}{L}$


$$\mathbf{K}=\frac{EA}{L}\begin{bmatrix}1&-1\\-1&1\end{bmatrix}$$




---
## Mesh
```plaintext
header: x y block_id volume
0 0 1 1
1 0 1 1
2 0 1 1
3 0 1 1
4 0 1 1
```

## Yaml
```yaml
PeriLab:
  Discretization:
    Node Sets:
      Node Set 1: 1
      Node Set 2: 5
    Type: "Text File"
    Input Mesh File: "truss.txt"
  Models:
    Material Models:
      Test:
        Material Model: "Bond-based Elastic"
        Symmetry: "isotropic plane stress"
        Young's Modulus: 7000
        Poisson's Ratio: 0.3
  Blocks:
    block_1:
      Block ID: 1
      Material Model: "Test"
      Density: 20
      Horizon: 2
  Boundary Conditions:
    BC_1:
      Variable: "Displacement"
      Node Set: "Node Set 1"
      Coordinate: "x"
      Value: "100*t"
      Type: Dirichlet
    BC_2:
      Variable: "Displacement"
      Coordinate: "x"
      Node Set: "Node Set 2"
      Value: "0.1*t"
      Type: Dirichlet
  Solver:
    Material Models: True
    Initial Time: 0.0
    Final Time: 1.0
    Number of Steps: 20
    Static:
      Show solver iteration: true
      Residual tolerance: 1e-7
      Solution tolerance: 1e-8
      Residual scaling: 7000
      m: 550
      Maximum number of iterations: 100
  Outputs:
    Output1:
      Output Filename: "truss"
      Output File Type: Exodus
      Number of Output Steps: 20
      Output Variables:
        Displacements: True
        Number of Neighbors: True
        Forces: True
```
