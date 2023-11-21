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
CurrentModule = PeriLab.Solver
```
## Solver
```@docs
init
get_blockNodes
set_density
set_horizon
solver
synchronise_field
write_results
```

<!-- ```@meta
CurrentModule = PeriLab.Solver.Verlet
```
## Verlet
```@docs
compute_thermodynamic_critical_time_step
compute_mechanical_critical_time_step
test_timestep
compute_crititical_time_step
init_solver
get_integration_steps
run_solver
``` -->

```@meta
CurrentModule = PeriLab.Solver.Boundary_conditions
```
## Boundary_conditions
```@docs
check_valid_bcs
init_BCs
boundary_condition
apply_bc
clean_up
eval_bc
```

```@meta
CurrentModule = PeriLab.IO
```
## IO
```@docs
merge_exodus_files
open_result_file
close_result_file
close_result_files
delete_files
get_file_size
clearNP1
get_results_mapping
initialize_data
init_write_results
read_input_file
write_results
get_global_values
find_global_core_value!
show_block_summary
Read_Mesh.Geometry.bond_geometry
calculate_nodelist
calculate_block
global_value_sum
global_value_max
global_value_min
global_value_avg
```

```@meta
CurrentModule = PeriLab.IO.Write_Exodus_Results
```
## Write_Exodus_Results
```@docs
create_result_file
paraview_specifics
get_paraview_coordinates
get_block_nodes
init_results_in_exodus
write_step_and_time
write_nodal_results_in_exodus
write_global_results_in_exodus
merge_exodus_file
```

```@meta
CurrentModule = PeriLab.IO.Write_CSV_Results
```
## Write_Exodus_Results
```@docs
create_result_file
write_global_results_in_csv
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
