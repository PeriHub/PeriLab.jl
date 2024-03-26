## Theory manual
The theory manual condenses the implemented features. 

Peridynamics is an integral continuum mechanics formulation. For a pure mechanical description it can be formulated as:

$$\int_{\mathcal{H}}\underline{\mathbf{T}}\langle\mathbf{x},t\rangle-\underline{\mathbf{T}}\langle\mathbf{x}',t\rangle dV- \mathbf{b}=\rho\ddot{\mathbf{u}}$$

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
| [Bond-based](theory/theory_bondbased.md) | [Bond-based Elastic](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Physics/Material/BondBased/Bondbased_Elastic.jl) |
| [Ordinary state-based](theory/theory_ordinary.md) | [PD Solid Elastic](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Physics/Material/Material_Models/PD_Solid_Elastic.jl) [PD Solid Plastic](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Physics/Material/Material_Models/PD_Solid_Plastic.jl) |
|[Non-ordinary state-based](theory/theory_correspondence.md)| [Correspondence Elastic](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Physics/Material/Material_Models/Correspondence_Elastic.jl) [Correspondence Plastic](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Physics/Material/Material_Models/Correspondence_Plastic.jl)|