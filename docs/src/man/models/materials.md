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

`Density` and `Horizon` are keys of the [block](../../generated/input_sections.md#Blocks). Every material model entry accepts these shared keys:

```@eval
using PeriLab, Markdown
Markdown.parse(PeriLab.parameter_table("Material Models"))
```

## Bond-based Elastic

 The [Bond-based Elastic](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Models/Material/BondBased/Bondbased_Elastic.jl) calculates the linear elastic behavior of a simple bond-based material. The theory of the bond-based elastic material can be found [here](@ref "Bond-based Peridynamics").

It needs one of the shared keys `Young's Modulus`, `Bulk Modulus` and `Shear Modulus`.

!!! info "Fixed Poisson's ratio"
    In the bond-based formulation the Poisson's ratio is fixed by [0.25 for 2D plane strain and 1/3 for 3D and 2D plane stress](https://link.springer.com/article/10.1007/s42102-019-00021-x), respectively.

## Unified Bond-based Elastic

The unified bond-based material after Guan et al. [GuanJ2024](@cite) allows Poisson's ratios
other than the fixed bond-based ones; see the [theory](@ref "Unified Bond-based Peridynamics").
The PeriLab version is simplified: no rigid body rotation and no surface correction. It needs
two of the shared isotropic elastic keys and has no keys of its own.

## PD Solid Elastic

The [PD Solid Elastic](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Models/Material/Material_Models/Ordinary/PD_Solid_Elastic.jl) module calculates the isotropic linear elastic material law for a peridynamic solid material. The underlying theory can be found [here](@ref "Ordinary state-based Peridynamics").

It needs two of the shared keys `Young's Modulus`, `Poisson's Ratio`, `Bulk Modulus` and `Shear Modulus`. The other two are determined automatically and can be used in the model if needed.

## PD Solid Plastic

The [PD Solid Plastic](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Models/Material/Material_Models/Ordinary/PD_Solid_Plastic.jl) material uses elastic stresses and calculate the plastic part for a peridynamic solid material. Has to be combined with a function, which provides elastic stresses.

```@eval
using PeriLab, Markdown
Markdown.parse(PeriLab.parameter_table("PD Solid Plastic"))
```


## Correspondence Elastic
The [Correspondence Elastic](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Models/Material/Material_Models/Correspondence/Correspondence_Elastic.jl) module calculates the fully anisotropic linear elastic material law. The underlying correspondence theory can be found [here](@ref "Correspondence Peridynamics").

It uses the shared elastic keys. With an anisotropic `Symmetry` you can provide the entries C11, ..., C66 of the [Hooke matrix](https://en.wikipedia.org/wiki/Hooke%27s_law#Matrix_representation_(stiffness_tensor)).

!!! tip "Isotropic elastic parameter"
    For the time step calculation two of the four isotropic elastic parameter have to be defined.

!!! tip "Material Rotation"
    If you define a field "Angles" for 2D or "Anglesx", "Anglesy" and "Anglesz" for 3D in the mesh file your material will be rotated. This helps to create an arbitrary material orientation.

## Correspondence Plastic
The [Correspondence Plastic](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Models/Material/Material_Models/Correspondence/Correspondence_Plastic.jl) material uses elastic stresses and calculate the plastic part. Has to be combined with a function, which provides elastic stresses.

```@eval
using PeriLab, Markdown
Markdown.parse(PeriLab.parameter_table("Correspondence Plastic"))
```

## Correspondence UMAT
The [Correspondence UMAT](https://github.com/PeriHub/PeriLab.jl/blob/main/src/Models/Material/Material_Models/Correspondence_UMAT.jl) can be used to include Abaqus user materials.

!!! warning "Replace ABAQUS Functions"
    No extra Abaqus functions should be called in the user subroutine. For example `INCLUDE 'ABA_PARAM.INC'` needs to be replaced by `implicit real(8) (a-h,o-z)`

You can call the user subroutine by defining path with a compiled Fotran library.
[WillbergC2023](@cite) gives an overview about the interface for Peridigm. In PeriLab all fields in the UMAT interface are supported execpt these integer (NOEL, NPT, LAYER, KSPT, JSTEP, KINC) and float values (PNEWDT, CELENT). In the material module these values defined for the interfase and named as not_supported_int and not_supported_float, respectively.

```@eval
using PeriLab, Markdown
Markdown.parse(PeriLab.parameter_table("Correspondence UMAT"))
```

!!! tip "Matrix based"
    If you use the matrix based approach you can use UMATS as well. But you must provide the tangential stiffness matrix. As solver right now you can use Matrix Verlet and Linear Static Matrix Based. The material definition remains equal.

### Compilation of the UMAT subroutine
In order to compile the UMAT subroutine you need to install gfortran. Have a look at [this](https://fortran-lang.org/en/learn/os_setup/install_gfortran/) page.
After installation you can compile the Fortran subroutine with the following command:

```shell
gfortran -shared -fPIC -ffixed-form -ffixed-line-length-132 -fno-range-check -o libusermat.so base.f
```

## Correspondence VUMAT

Like the [Correspondence UMAT](#Correspondence-UMAT), but for an Abaqus VUMAT (explicit)
subroutine compiled into a shared library. Its properties are given as `Property_1` …
`Property_N`.

```@eval
using PeriLab, Markdown
Markdown.parse(PeriLab.parameter_table("Correspondence VUMAT"))
```

## Rigid

A rigid block: no internal forces and no deformation. It has no keys of its own.

## Model merging

In PeriLab you are able to combine models with each other, by simply adding a +. Therefore, modules can be merged and double coding can be avoided.

!!! tip "Elastic platic"
    If you want to run elastic platic please use Correspondence Elastic + Correspondence Plastic or PD Solid Elastic + PD Solid Plastic

!!! warning "Model order"
    The order is defined by the user. Because the plastic routines need stresses to work, make sure the materials which provide these stresses are before the plastic models.
