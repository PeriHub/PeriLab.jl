# Functions

## Multithreading support



## Index
```@index
Pages = ["functions.md"]
```

```@meta
CurrentModule = PeriLab.Data_manager
```
## Data_manager
```@docs
set_glob_to_loc
get_property
get_rank
set_nset
create_constant_bond_field
get_local_nodes
create_node_field
get_properties
get_max_rank
set_num_controller
set_num_responder
create_constant_node_field
```

```@meta
CurrentModule = PeriLab.IO
```
## IO
```@docs
Read_Mesh.Geometry.bond_geometry
```

```@meta
CurrentModule = PeriLab.Solver.Physics
```
## Physics
```@docs
Material.Set_modules.Material_template.material_name
Material.Set_modules.Correspondence.material_name
Material.Set_modules.Material_template.compute_force
Material.Set_modules.Correspondence.Correspondence_Elastic.compute_stresses
Material.Set_modules.Correspondence.Correspondence_Elastic.correspondence_name
Additive.Set_modules.Additive_template.additive_name
Additive.Set_modules.simple_additive.additive_name
Additive.Set_modules.simple_additive.compute_additive
Additive.Set_modules.Additive_template.compute_additive
Thermal.Set_modules.Heat_transfer.compute_thermal_model
Thermal.Set_modules.Thermal_Flow.compute_thermal_model
Thermal.Set_modules.Thermal_expansion.compute_thermal_model
Thermal.Set_modules.Thermal_template.compute_thermal_model
Thermal.Set_modules.Thermal_expansion.thermal_model_name
Thermal.Set_modules.Heat_transfer.thermal_model_name
Thermal.Set_modules.Thermal_Flow.thermal_model_name
Thermal.Set_modules.Thermal_template.thermal_model_name
```

```@meta
CurrentModule = PeriLab.Solver.Damage
```
## Damage
```@docs
Set_modules.Damage_template.compute_damage_pre_calculation
Set_modules.Critical_stretch.compute_damage
Set_modules.Damage_template.compute_damage
Set_modules.Critical_Energy_Model.compute_damage
Set_modules.Critical_Energy_Model.damage_name
Set_modules.Critical_stretch.damage_name
Set_modules.Damage_template.damage_name
```
