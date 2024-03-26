# Material Models

## Existing Models
The physics related functions can be found [here](@ref "Physics - Functions").

| Method | Related Model in PeriLab |
|---|---|
| [Bond-based](@ref "Bond-Based Peridynamics") | [Bond-based Elastic](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Physics/Material/BondBased/Bondbased_Elastic.jl) |
| [Ordinary state-based](@ref "Ordinary state-based") | [PD Solid Elastic](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Physics/Material/Material_Models/PD_Solid_Elastic.jl), [PD Solid Plastic](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Physics/Material/Material_Models/PD_Solid_Plastic.jl) |
|[Non-ordinary state-based](@ref "Correspondence")| [Correspondence Elastic](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Physics/Material/Material_Models/Correspondence_Elastic.jl), [Correspondence Plastic](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Physics/Material/Material_Models/Correspondence_Plastic.jl)|


| Material Model             | Bond-based Elastic | PD Solid Elastic | PD Solid Plastic | Correspondence Elastic | Correspondence Plastic |
|----------------------------|-------------------|------------------|------------------|------------------------|------------------------|
| Symmetry                   | ✔️| ✔️ | ✔️ | ✔️     | ✔️     |
| Poisson's/Young's/Bulk/Shear | ✔️| ✔️ | ✔️ | ✔️     | ✔️     |
| Yield Stress               |                   |                  | ✔️|                        | ✔️     |
| Zero Energy Control        |                   |                  |                  | ✔️     | ✔️     |
| C11/C12/.../C66            | (✔️)| (✔️) | (✔️) | (✔️)     | (✔️)     |

## Bond-based Elastic

The theory of the bond-based elastic material can be found [here](@ref "Theory Basics")

Material parameter
tabelle


| Parameter | Unit | Description |
|---|---|---|
|Youngs Modulus | $\left[N/m^2\right]$| [Young's modulus](https://en.wikipedia.org/wiki/Young%27s_modulus) or elasticity modulus
|Shear Modulus | $\left[N/m^2\right]$| [Shear modulus](https://en.wikipedia.org/wiki/Shear_modulus)
|Bulk Modulus | $\left[N/m^2\right]$| [Bulk modulus](https://en.wikipedia.org/wiki/Bulk_modulus) or compression modulus

One of theses parameters have to be defined. 
>Note: In the bond-based formulation the Poisson's ratio is fixed by [0.25 for 2D plane strain and 1/3 for 3D and 2D plane stress](https://link.springer.com/article/10.1007/s42102-019-00021-x), respectively.

## PD Solid Elastic

| Parameter | Unit | Description |
|---|---|---|
|Youngs Modulus | $\left[N/m^2\right]$| [Young's modulus](https://en.wikipedia.org/wiki/Young%27s_modulus) or elasticity modulus
|Shear Modulus | $\left[N/m^2\right]$| [Shear modulus](https://en.wikipedia.org/wiki/Shear_modulus)
|Bulk Modulus | $\left[N/m^2\right]$| [Bulk modulus](https://en.wikipedia.org/wiki/Bulk_modulus) or compression modulus
|Poissons Ratio Modulus | $\left[-\right]$| [Poisson's ratio](https://en.wikipedia.org/wiki/Poisson%27s_ratio)

Two of these parameters have to be defined. The other two are determined automatically and can be used in the model if needed.

## PD Solid Plastic

| Parameter | Unit | Description |
|---|---|---|
|Yield Stress | $\left[N/m^2\right]$| [Yield stress](https://en.wikipedia.org/wiki/Yield_(engineering))


## Correspondence Elastic


| Parameter | Unit | Description |
|---|---|---|
|Youngs Modulus | $\left[N/m^2\right]$| [Young's modulus](https://en.wikipedia.org/wiki/Young%27s_modulus) or elasticity modulus
|Shear Modulus | $\left[N/m^2\right]$| [Shear modulus](https://en.wikipedia.org/wiki/Shear_modulus) or elasticity modulus
|Bulk Modulus | $\left[N/m^2\right]$| [Bulk modulus](https://en.wikipedia.org/wiki/Bulk_modulusseh) or compression modulus
|Poissons Ratio Modulus | $\left[-\right]$| [Poisson's ratio](https://en.wikipedia.org/wiki/Poisson%27s_ratio)
|C11, C12, ..., C66 (optional) | $\left[N/m^2\right]$| [Parameter of the Hook matrix](https://en.wikipedia.org/wiki/Hooke%27s_law#Matrix_representation_(stiffness_tensor))

For Correspondence Elastic you can provide all 27 elastic parameters if you like by adding C11,...,C66. 

>Note: For the time step calculation two of the four isotropic elastic parameter have to be defined.

>Note 2: If you define a field "Angles" for 2D or "Anglesx", "Anglesy" and "Anglesz" for 3D in the mesh file your material will be rotated. This helps to create an arbitrary material orientation.

## Correspondence Plastic

| Parameter | Unit | Description |
|---|---|---|
|Yield Stress | $\left[N/m^2\right]$| [Yield stress](https://en.wikipedia.org/wiki/Yield_(engineering))

## Correspondence UMAT

in developement


## Model merging

In PeriLab you are able to combine models with each other, by simply adding a +. Therefore, modules can be merged and double coding can be avoided.

>Note: If you want to run elastic platic please use Correspondence Elastic + Correspondence Plastic or PD Solid Elastic + PD Solid Plastic

You can find the paramter by opening the routines. By I will add them as soon as possible in the function description and the documentation.

For isotropic material two of these parameters have to be provided Shear Modulus, Poissons Ration, Youngs Modulus or Bulk Modulus.



