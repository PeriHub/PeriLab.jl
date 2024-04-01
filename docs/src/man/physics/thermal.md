# Thermal Models
The list shows the three main model, thermal expansion, thermal flow and heat transfer. All can use simular parameters to work. 


| Thermal Model                 | Thermal Expansion | Thermal Flow | Heat Transfer |
|-------------------------------|-------------------|--------------|---------------|
| Type                          | ✔️ | ✔️ | ✔️ |
| Heat Transfer Coefficient     | ✔️ | ✔️ | ✔️ |
| Environmental Temperature     | ✔️ | ✔️ | ✔️ |
| Thermal Conductivity          |                   |              | ✔️ |
| Thermal Expansion Coefficient | ✔️ |              |               |
| Thermal Conductivity Print Bed|                   |              | ✔️ |
| Print Bed Temperature         |                   |              | ✔️ |

There are block defined properties for needed for the thermal models.

| Parameter | Unit | Description |
|---|---|---|
| Specific Heat Capacity |  $\left[\frac{J}{kgK}\right]$ | Specific heat capacity of the block  |
| Density |  $\left[\frac{kg}{m^3}\right]$ | Specific heat capacity of the block  
| Horizon |  $[m]$ | Radius of the neighborhood |
## Thermal Expansion

Calculates the thermal expansion of the material. 

| Parameter | Unit | Description |
|---|---|---|
|Thermal Expansion Coefficient | $\left[1/K\right]$| [Thermal expansion](https://en.wikipedia.org/wiki/Thermal_expansion) can be a $3\times3$ matrix. |

>Note: PeriLab supports currently only isotropic thermal expansion. 

## Thermal Flow

| Parameter | Unit | Description |
|---|---|---|
| Thermal Conductivity |  $\left[\frac{W}{mK}\right]$  | |
## Heat Transfer

| Parameter | Unit | Description |
|---|---|---|
| Heat Transfer Coefficient |  $\left[\frac{W}{m^2K}\right]$ | Coefficient describing the heat transfer between a solid and a gas or liquid |
## Model merging

In PeriLab you are able to combine models with each other, by simply adding a +. Therefore, modules can be merged and double coding can be avoided. This is necessary if you want to model the heating of a model and its expansion.

>Note: If you want to run a full thermal model Thermal Flow + Heat Transfer + Thermal Expansion.

>Note 2: The order is defined by the user. Therfore, in this example first the flow, than the transfer to the environment and than the expansion will be calculated. 