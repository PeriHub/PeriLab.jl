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


## Thermal Expansion

Calculates the thermal expansion of the material. 

| Parameter | Unit | Description |
|---|---|---|
|Thermal Expansion Coefficient | $\left[1/K\right]$| [Thermal expansion](https://en.wikipedia.org/wiki/Thermal_expansion) can be a $3\times3$ matrix. |

>Note: PeriLab supports currently only isotropic thermal expansion. 

## Thermal Flow
| Parameter | Unit | Description |
|---|---|---|
|  |   |
## Heat Transfer
| Parameter | Unit | Description |
|---|---|---|
|  |   |
## Model merging

In PeriLab you are able to combine models with each other, by simply adding a +. Therefore, modules can be merged and double coding can be avoided. This is necessary if you want to model the heating of a model and its expansion.

>Note: If you want to run a full thermal model Thermal Flow + Heat Transfer + Thermal Expansion.

>Note 2: The order is defined by the user. Therfore, in this example first the flow, than the transfer to the environment and than the expansion will be calculated. 