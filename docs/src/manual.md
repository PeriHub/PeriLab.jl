# Model defintion
A basic PeriLab model exists two files. A mesh file and a yaml file.

## Mesh file
The mesh file is an ASCII file and has the following structure for 2D.
| header: x | y | block_id | volume |
|:---|:---|:---|:---|
|-1.5| 0 | 1 | 1.e-01 |
|-0.5| 0 | 1 | 1.e-01 |
|0.5 | 0 | 1 | 1.e-01 |
|1.5 | 0 | 1 | 1.e-01 |
|-1.5| 1 | 1 | 1.e-01 |

and 3D 

| header: x | y | z | block_id | volume |
|:---|:---|:---|:---|:---|
|-1.5| 0 | 2 | 1 | 1.e-01 |
|-0.5| 0 | 2 | 1 | 1.e-01 |
|0.5 | 0 | 2 | 1 | 1.e-01 |
|1.5 | 0 | 0 | 1 | 1.e-01 |
|-1.5| 1 | 1 | 1 | 1.e-01 |

The difference is found automatically. If no z occurs PeriLab identifies it as 2D problem and requests a plane stress or plane strain definition.

Additional parameter can be applied in the header. They will be added in the datamanager and can be used in the programm. If you add x,y,z to the parameter a multidimensional field will be created, e.g.
MyVarx, MyVary will be created as field MyVar with 2 degrees of freedom.


**Variable definition**
| Parameter | Datamanager name | Header name 2D | Header name 3D | Type | 
|:---|:---|:---|:---|:---|
|x, y, z (optional) coordinate of the node | Coordinates | x, y | x, y, z | Float64, Int64|
| Definition to which block the node corresponds. Is needed in the Yaml file to define properties | Block_Id | block_id | block_id | Int64|
| Volume the node represents. | Volume | volume | volume | Float64, Int64|


**Optional**
| Parameter | Datamanager name | Header name 2D | Header name 3D | Type | 
|:---|:---|:---|:---|:---|
| Orientation of a Node | Angles | Angles | Anglesx, Anglesy, Anglesz | Float64, Int64|
| Activation time of a node, e.g. used for additive manufacturing to define when the node will be acativated | Activation_Time | Activation_Time | Activation_Time | Float64, Int64|
| Status of the node. If it is false the node is deactivated, but exists. This variable is automatically created if additive models are used and set everywhere to false, if it is not predefined | Active | Active | Active | Bool |

## Yaml file

tbd

## Node set file (optional)


# Existing Models
The physics related functions can be found [here](
https://perihub.github.io/PeriLab.jl/stable/lib/physics_functions/
).
## Material Models

| Method | Related Model in PeriLab |
|---|---|
| [Bond-based](theory/basics.md#Bond-based) | [Bond-based Elastic](https://gitlab.com/dlr-perihub/PeriLab.jl/-/blob/main/src/Physics/Material/BondBased/Bondbased_Elastic.jl) |
| [Ordinary state-based](theory/basics.md#PD_Solind) | [PD Solid Elastic](https://gitlab.com/dlr-perihub/PeriLab.jl/-/blob/main/src/Physics/Material/Material_Models/PD_Solid_Elastic.jl) [PD Solid Plastic](https://gitlab.com/dlr-perihub/PeriLab.jl/-/blob/main/src/Physics/Material/Material_Models/PD_Solid_Plastic.jl) |
|[Non-ordinary state-based](theory/basics.md#Correspondence)| [Correspondence Elastic](https://gitlab.com/dlr-perihub/PeriLab.jl/-/blob/main/src/Physics/Material/Material_Models/Correspondence_Elastic.jl) [Correspondence Plastic](https://gitlab.com/dlr-perihub/PeriLab.jl/-/blob/main/src/Physics/Material/Material_Models/Correspondence_Plastic.jl)|

**Bond-based Elastic**

The theory of the bond-based elastic material can be found [here](theory/basics.md#Bond-based)

Material parameter
tabelle


| Parameter | Unit | Description |
|---|---|---|
|Youngs Modulus | $\left[N/m^2\right]$| [Young's modulus](https://en.wikipedia.org/wiki/Young%27s_modulus) or elasticity modulus
|Shear Modulus | $\left[N/m^2\right]$| [Shear modulus](https://en.wikipedia.org/wiki/Shear_modulus)
|Bulk Modulus | $\left[N/m^2\right]$| [Bulk modulus](https://en.wikipedia.org/wiki/Bulk_modulus) or compression modulus

One of theses parameters have to be defined. 
>Note: In the bond-based formulation the Poisson's ratio is fixed by [0.25 for 2D plane strain and 1/3 for 3D and 2D plane stress](https://link.springer.com/article/10.1007/s42102-019-00021-x), respectively.

**PD Solid Elastic**

| Parameter | Unit | Description |
|---|---|---|
|Youngs Modulus | $\left[N/m^2\right]$| [Young's modulus](https://en.wikipedia.org/wiki/Young%27s_modulus) or elasticity modulus
|Shear Modulus | $\left[N/m^2\right]$| [Shear modulus](https://en.wikipedia.org/wiki/Shear_modulus)
|Bulk Modulus | $\left[N/m^2\right]$| [Bulk modulus](https://en.wikipedia.org/wiki/Bulk_modulus) or compression modulus
|Poissons Ratio Modulus | $\left[-\right]$| [Poisson's ratio](https://en.wikipedia.org/wiki/Poisson%27s_ratio)

Two of these parameters have to be defined. The other two are determined automatically and can be used in the model if needed.

**PD Solid Plastic**

| Parameter | Unit | Description |
|---|---|---|
|Yield Stress | $\left[N/m^2\right]$| [Yield stress](https://en.wikipedia.org/wiki/Yield_(engineering))


**Correspondence Elastic**


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


**Correspondence Plastic**

| Parameter | Unit | Description |
|---|---|---|
|Yield Stress | $\left[N/m^2\right]$| [Yield stress](https://en.wikipedia.org/wiki/Yield_(engineering))

**Correspondence UMAT**

in developement


**Model merging**

In PeriLab you are able to combine models with each other, by simply adding a +. Therefore, modules can be merged and double coding can be avoided.

>Note: If you want to run elastic platic please use Correspondence Elastic + Correspondence Plastic or PD Solid Elastic + PD Solid Plastic

You can find the paramter by opening the routines. By I will add them as soon as possible in the function description and the documentation.

For isotropic material two of these parameters have to be provided Shear Modulus, Poissons Ration, Youngs Modulus or Bulk Modulus.



