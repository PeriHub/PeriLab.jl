# Thermal Models
The main models are thermal expansion, thermal flow and heat transfer. `Specific Heat Capacity`, `Density` and `Horizon` are keys of the [block](../../generated/input_sections.md#Blocks). Every thermal model entry accepts these shared keys:

```@eval
using PeriLab, Markdown
Markdown.parse(PeriLab.parameter_table("Thermal Models"))
```
## Thermal Expansion

Calculates the thermal expansion of the material.

```@eval
using PeriLab, Markdown
Markdown.parse(PeriLab.parameter_table("Thermal Expansion"))
```

!!! info "Thermal expansion"
    PeriLab supports currently only isotropic thermal expansion.

## Thermal Flow

```@eval
using PeriLab, Markdown
Markdown.parse(PeriLab.parameter_table("Thermal Flow"))
```
## Heat Transfer

```@eval
using PeriLab, Markdown
Markdown.parse(PeriLab.parameter_table("Heat Transfer"))
```
## HETVAL

User-defined heat generation from an Abaqus HETVAL subroutine, compiled into a shared library
like a [UMAT](@ref "Compilation of the UMAT subroutine"). Its properties are given as
`Property_1` … `Property_N`.

```@eval
using PeriLab, Markdown
Markdown.parse(PeriLab.parameter_table("HETVAL"))
```

## Model merging

In PeriLab you are able to combine models with each other, by simply adding a +. Therefore, modules can be merged and double coding can be avoided. This is necessary if you want to model the heating of a model and its expansion.

!!! tip "Full Model"
    If you want to run a full thermal model Thermal Flow + Heat Transfer + Thermal Expansion.

!!! tip "Model order"
    The order is defined by the user. Therfore, in this example first the flow, than the transfer to the environment and than the expansion will be calculated.
