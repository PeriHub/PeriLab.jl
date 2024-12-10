# Additive Models

| Additive Model           | Simple |
|------------------------|:----------------:|
| Print Temperature         | ✔️|


## Simple

To realize an additive model all bonds from a point are disconnected from it's neighbors by setting all $\underline{\omega}\langle\boldsymbol{\xi}\rangle = 0$. Within the mesh input an activation time is specified $t_{activate}$. If this time is reached during the simulation process the point is activated.

$
\underline{\omega}\langle\boldsymbol{\xi}\rangle \in \mathcal{H}_{\textbf{x}} =    \left\{\begin{array}{l}
0 \qquad\text{for }t<t_{activate}\\
1  \qquad\text{for }t \geq t_{activate}\\
\end{array}\right.
$

Depended on the process modeled, additional information can be passed to the point or the bonds connected with him.
For this simple printing process a printing temperature is added utilized the heat source $S_i$. In the presented model, the bonding is ideal an no phase or chemical changes occur. However, in principal such models are applicable. Also it must be noted, that due to mechanical or thermo-mechanical loading bonds can be damaged, if a damage model is applied. In that case it won't be activated again.

Within this process the $t_{activate}$ can be user defined. However, to reproduce real processes an interface with the G-code is needed. This interface provides the information when the tool arrives at a specific point and defines the activation time as shown in the figure (taken from [WillbergC2024b](@cite)).



![Virtual printing process](https://onlinelibrary.wiley.com/cms/asset/b5c3f70e-2585-4b01-854b-88febf14741e/adts202400818-fig-0009-m.jpg)
