# Material Models

## Existing Models
The models related functions can be found [here](@ref "Model Factory - Functions").

| Method | Related Model in PeriLab |
|---|---|
| [Bond-based](@ref "Bond-based Peridynamics") | [Bond-based Elastic](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Models/Material/BondBased/Bondbased_Elastic.jl) |
| [Ordinary state-based](@ref "Ordinary state-based Peridynamics") | [PD Solid Elastic](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Models/Material/Material_Models/Ordinary/PD_Solid_Elastic.jl), [PD Solid Plastic](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Models/Material/Material_Models/Ordinary/PD_Solid_Plastic.jl) |
|[Non-ordinary state-based](@ref "Correspondence Peridynamics")| [Correspondence Elastic](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Models/Material/Material_Models/Correspondence/Correspondence_Elastic.jl), [Correspondence Plastic](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Models/Material/Material_Models/Correspondence/Correspondence_Plastic.jl), [Correspondence UMAT](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Models/Material/Material_Models/Correspondence/Correspondence_UMAT.jl)|


| Material Model             | Bond-based Elastic | PD Solid Elastic | PD Solid Plastic | Correspondence Elastic | Correspondence Plastic |
|----------------------------|-------------------|------------------|------------------|------------------------|------------------------|
| Symmetry                   | ✔️| ✔️ | ✔️ | ✔️     | ✔️     |
| Poisson's/Young's/Bulk/Shear | ✔️| ✔️ | ✔️ | ✔️     | ✔️     |
| Yield Stress               |                   |                  | ✔️|                        | ✔️     |
| Zero Energy Control        |                   |                  |                  | ✔️     | ✔️     |
| C11/C12/.../C66            | (✔️)| (✔️) | (✔️) | (✔️)     | (✔️)     |

| Parameter | Unit | Description |
|---|---|---|
| Density |  $\left[\frac{kg}{m^3}\right]$ | Specific heat capacity of the block
| Horizon |  $[m]$ | Radius of the neighborhood |

## Bond-based Elastic

 The [Bond-based Elastic](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Models/Material/BondBased/Bondbased_Elastic.jl) calculates the linear elastic behavior of a simple bond-based material. The theory of the bond-based elastic material can be found [here](@ref "Bond-based Peridynamics").

| Parameter | Unit | Description |
|---|---|---|
|Youngs Modulus | $\left[N/m^2\right]$| [Young's modulus](https://en.wikipedia.org/wiki/Young%27s_modulus) or elasticity modulus
|Shear Modulus | $\left[N/m^2\right]$| [Shear modulus](https://en.wikipedia.org/wiki/Shear_modulus)
|Bulk Modulus | $\left[N/m^2\right]$| [Bulk modulus](https://en.wikipedia.org/wiki/Bulk_modulus) or compression modulus

One of theses parameters have to be defined.

!!! info "Fixed Poisson's ratio"
    In the bond-based formulation the Poisson's ratio is fixed by [0.25 for 2D plane strain and 1/3 for 3D and 2D plane stress](https://link.springer.com/article/10.1007/s42102-019-00021-x), respectively.

## PD Solid Elastic

The [PD Solid Elastic](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Models/Material/Material_Models/Ordinary/PD_Solid_Elastic.jl) module calculates the isotropic linear elastic material law for a peridynamic solid material. The underlying theory can be found [here](@ref "Ordinary state-based Peridynamics").

| Parameter | Unit | Description |
|---|---|---|
|Youngs Modulus | $\left[N/m^2\right]$| [Young's modulus](https://en.wikipedia.org/wiki/Young%27s_modulus) or elasticity modulus
|Shear Modulus | $\left[N/m^2\right]$| [Shear modulus](https://en.wikipedia.org/wiki/Shear_modulus)
|Bulk Modulus | $\left[N/m^2\right]$| [Bulk modulus](https://en.wikipedia.org/wiki/Bulk_modulus) or compression modulus
|Poissons Ratio Modulus | $\left[-\right]$| [Poisson's ratio](https://en.wikipedia.org/wiki/Poisson%27s_ratio)

Two of these parameters have to be defined. The other two are determined automatically and can be used in the model if needed.

## PD Solid Plastic

The [PD Solid Plastic](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Models/Material/Material_Models/Ordinary/PD_Solid_Plastic.jl) material uses elastic stresses and calculate the plastic part for a peridynamic solid material. Has to be combined with a function, which provides elastic stresses.

| Parameter | Unit | Description |
|---|---|---|
|Yield Stress | $\left[N/m^2\right]$| [Yield stress](https://en.wikipedia.org/wiki/Yield_(engineering))


## Correspondence Elastic
The [Correspondence Elastic](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Models/Material/Material_Models/Correspondence/Correspondence_Elastic.jl) module calculates the fully anisotropic linear elastic material law. The underlying correspondence theory can be found [here](@ref "Correspondence Peridynamics").

| Parameter | Unit | Description |
|---|---|---|
|Youngs Modulus | $\left[N/m^2\right]$| [Young's modulus](https://en.wikipedia.org/wiki/Young%27s_modulus) or elasticity modulus
|Shear Modulus | $\left[N/m^2\right]$| [Shear modulus](https://en.wikipedia.org/wiki/Shear_modulus) or elasticity modulus
|Bulk Modulus | $\left[N/m^2\right]$| [Bulk modulus](https://en.wikipedia.org/wiki/Bulk_modulusseh) or compression modulus
|Poissons Ratio Modulus | $\left[-\right]$| [Poisson's ratio](https://en.wikipedia.org/wiki/Poisson%27s_ratio)
|C11, C12, ..., C66 (optional) | $\left[N/m^2\right]$| [Parameter of the Hook matrix](https://en.wikipedia.org/wiki/Hooke%27s_law#Matrix_representation_(stiffness_tensor))

For Correspondence Elastic you can provide all 27 elastic parameters if you like by adding C11,...,C66.

!!! tip "Isotropic elastic parameter"
    For the time step calculation two of the four isotropic elastic parameter have to be defined.

!!! tip "Material Rotation"
    If you define a field "Angles" for 2D or "Anglesx", "Anglesy" and "Anglesz" for 3D in the mesh file your material will be rotated. This helps to create an arbitrary material orientation.

## Correspondence Plastic
The [Correspondence Plastic](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Models/Material/Material_Models/Correspondence/Correspondence_Plastic.jl) material uses elastic stresses and calculate the plastic part. Has to be combined with a function, which provides elastic stresses.

| Parameter | Unit | Description |
|---|---|---|
|Yield Stress | $\left[N/m^2\right]$| [Yield stress](https://en.wikipedia.org/wiki/Yield_(engineering))

## Correspondence UMAT
The [Correspondence UMAT](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Models/Material/Material_Models/Correspondence_UMAT.jl) can be used to include Abaqus user materials.

!!! warning "Replace ABAQUS Functions"
    No extra Abaqus functions should be called in the user subroutine. For example `INCLUDE 'ABA_PARAM.INC'` needs to be replaced by `implicit real(8) (a-h,o-z)`

You can call the user subroutine by defining path with a compiled Fotran library.
[WillbergC2023](@cite) gives an overview about the interface for Peridigm. In PeriLab all fields in the UMAT interface are supported execpt these integer (NOEL, NPT, LAYER, KSPT, JSTEP, KINC) and float values (PNEWDT, CELENT). In the material module these values defined for the interfase and named as not_supported_int and not_supported_float, respectively.

| Parameter | Type and Range | Description | Optional |
|---|---|---|---|
| File | String | Path and filename of the UMAT, e.g. "./src/Models/Material/UMATs/libusertest.so" | No |
| Number of State Variables | Int $\geq$ 0 | Number of state variables; Defines the size of state variable field Data_Manager.create_constant_node_scalar_field("State Variables", Float64, num_state_vars)| yes |
| Number of Properties | Int $\geq$ 1 | Properties of the material; Needed for the propterty field Data_Manager.create_constant_free_size_field("Properties", Float64, (num_props, 1))
 | No |
| Property_$iID| Float | iID is 1...Number of Properties. It has to be in order and can be utilized in the UMAT. | No|
|UMAT Material Name|String (maximum are 80 characters)| Defines material names defined in the UMAT to destinguash between different areas of the Fortran routine | Yes |
|UMAT name| in development | Should allow the definition of own subroutine name. The standard will be UMAT| in development |
|Predefined Field Names| String separated by spaces $\geq0$ | Define all the fields in the mesh file which should be used as pre-defined values. An increment field is than defined as well. E.g. Predefined Field Names: "Temperature" "Color"; Temperature and Color must exist in the mesh file. **They must be defined as Float or Int in that case**.| Yes|

### Compilation of the UMAT subroutine
In order to compile the UMAT subroutine you need to install gfortran. Have a look at [this](https://fortran-lang.org/en/learn/os_setup/install_gfortran/) page.
After installation you can compile the Fortran subroutine with the following command:

```shell
gfortran -shared -fPIC -o libusermat.so base.f
```

## Model merging

In PeriLab you are able to combine models with each other, by simply adding a +. Therefore, modules can be merged and double coding can be avoided.

!!! tip "Elastic platic"
    If you want to run elastic platic please use Correspondence Elastic + Correspondence Plastic or PD Solid Elastic + PD Solid Plastic

!!! warning "Model order"
    The order is defined by the user. Because the plastic routines need stresses to work, make sure the materials which provide these stresses are before the plastic models.
