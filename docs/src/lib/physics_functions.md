# Physics - Functions

## Index
```@index
Pages = ["physics_functions.md"]
```

```@meta
CurrentModule = PeriLab.Solver.Physics
```
## Physics
```@docs
compute_models
compute_damage_pre_calculation
get_block_model_definition
init_material_model_fields
init_damage_model_fields
init_models
init_thermal_model_fields
init_additive_model_fields
init_pre_calculation
read_properties
set_heat_capacity
```

```@meta
CurrentModule = PeriLab.Solver.Physics.Additive
```
## Additive
```@docs
Set_modules.Additive_template.additive_name
Set_modules.Additive_template.compute_additive_model
```

```@meta
CurrentModule = PeriLab.Solver.Damage
```
## Damage
```@docs
Set_modules.Damage_template.damage_name
Set_modules.Damage_template.compute_damage
Set_modules.Damage_template.compute_damage_pre_calculation
damage_index
set_bond_damage
init_interface_crit_values
init_aniso_crit_values
Set_modules.Critical_Energy_Model.get_quad_horizon
```

```@meta
CurrentModule = PeriLab.Solver.Physics.Material
```
## Material
```@docs
Set_modules.Material_template.init_material_model
Set_modules.Material_template.material_name
Set_modules.Material_template.compute_forces
determine_isotropic_parameter
check_material_symmetry
distribute_force_densities
get_all_elastic_moduli
get_Hooke_matrix
distribute_forces
matrix_to_voigt
voigt_to_matrix
check_symmetry
get_symmmetry
Set_modules.PD_Solid_Elastic.elastic
Set_modules.PD_Solid_Elastic.calculate_symmetry_params
Set_modules.PD_Solid_Elastic.Ordinary.compute_weighted_volume
Set_modules.PD_Solid_Elastic.Ordinary.compute_dilatation
Set_modules.Correspondence.zero_energy_mode_compensation
Set_modules.Correspondence.calculate_bond_force
Set_modules.Correspondence.rotate
Set_modules.Correspondence.rotate_second_order_tensor
Set_modules.Correspondence.Correspondence_Elastic.compute_stresses
Set_modules.Correspondence.Correspondence_Elastic.correspondence_name
Set_modules.Correspondence.global_zero_energy_control.control_name
Set_modules.Correspondence.global_zero_energy_control.compute_control
Set_modules.Correspondence.global_zero_energy_control.get_zero_energy_mode_force
Set_modules.Correspondence.global_zero_energy_control.create_zero_energy_mode_stiffness
Set_modules.Correspondence.global_zero_energy_control.rotate_fourth_order_tensor
```

```@meta
CurrentModule = PeriLab.Solver.Physics.Thermal
```
## Thermal
```@docs
Set_modules.Thermal_template.thermal_model_name
Set_modules.Thermal_template.compute_thermal_model
distribute_heat_flows
Thermal.Set_modules.Heat_transfer.calculate_specific_volume
Thermal.Set_modules.Thermal_expansion.thermal_deformation
Thermal.Set_modules.Thermal_Flow.compute_heat_flow_state_correspondence
Thermal.Set_modules.Thermal_Flow.compute_heat_flow_state_bond_based
```