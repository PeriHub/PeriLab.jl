# Models - Functions

## Index
```@index
Pages = ["models_functions.md"]
```

```@meta
CurrentModule = PeriLab.Solver.Model_Factory
```
## Models
```@docs
compute_models
get_block_model_definition
init_models
read_properties
set_heat_capacity
```

```@meta
CurrentModule = PeriLab.Solver.Model_Factory.Additive
```
## Additive
```@docs
Set_modules.Additive_template.additive_name
Set_modules.Additive_template.compute_model
```

```@meta
CurrentModule = PeriLab.Solver.Model_Factory.Damage
```
## Damage
```@docs
Set_modules.Damage_template.damage_name
Set_modules.Damage_template.compute_model
damage_index
set_bond_damage
init_interface_crit_values
init_aniso_crit_values
Set_modules.Critical_Energy_Model.get_quad_horizon
```

```@meta
CurrentModule = PeriLab.Solver.Model_Factory.Material
```
## Material
```@docs
Set_modules.Material_template.init_model
Set_modules.Material_template.material_name
Set_modules.Material_template.compute_model
determine_isotropic_parameter
check_material_symmetry
distribute_force_densities
get_all_elastic_moduli
get_Hooke_matrix
distribute_forces
matrix_to_voigt
voigt_to_matrix
check_symmetry
get_symmetry
get_von_mises_stress
Set_modules.PD_Solid_Elastic.elastic
Set_modules.PD_Solid_Elastic.Ordinary.calculate_symmetry_params
Set_modules.PD_Solid_Elastic.Ordinary.compute_weighted_volume
Set_modules.PD_Solid_Elastic.Ordinary.compute_dilatation
Set_modules.Correspondence.zero_energy_mode_compensation
Set_modules.Correspondence.calculate_bond_force
Set_modules.Correspondence.Set_modules.Correspondence_Elastic.compute_stresses
Set_modules.Correspondence.Set_modules.Correspondence_Elastic.correspondence_name
Set_modules.Correspondence.Global_zero_energy_control.control_name
Set_modules.Correspondence.Global_zero_energy_control.compute_control
Set_modules.Correspondence.Global_zero_energy_control.get_zero_energy_mode_force
Set_modules.Correspondence.Global_zero_energy_control.create_zero_energy_mode_stiffness
Set_modules.Correspondence.Global_zero_energy_control.rotate_fourth_order_tensor
```

```@meta
CurrentModule = PeriLab.Solver.Model_Factory.Thermal
```
## Thermal
```@docs
Set_modules.Thermal_template.thermal_model_name
Set_modules.Thermal_template.compute_model
distribute_heat_flows
Thermal.Set_modules.Heat_transfer.calculate_specific_volume
Thermal.Set_modules.Thermal_expansion.thermal_deformation
Thermal.Set_modules.Thermal_Flow.compute_heat_flow_state_correspondence
Thermal.Set_modules.Thermal_Flow.compute_heat_flow_state_bond_based
```
