# IO - Functions

## Index
```@index
Pages = ["io_functions.md"]
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
init_orientations
init_write_results
write_results
get_global_values
find_global_core_value!
show_block_summary
calculate_nodelist
calculate_block
global_value_sum
global_value_max
global_value_min
global_value_avg
```

```@meta
CurrentModule = PeriLab.IO
```
## Read_Mesh
```@docs
init_data
create_and_distribute_bond_norm
get_local_element_topology
get_local_overlap_map
local_nodes_from_dict
distribute_neighborhoodlist_to_cores
get_local_neighbors
get_bond_geometry
define_nsets
distribution_to_cores
check_mesh_elements
read_mesh
area_of_polygon
set_dof
load_and_evaluate_mesh
create_neighborhoodlist
get_number_of_neighbornodes
node_distribution
_init_overlap_map_
create_overlap_map
create_distribution
create_distribution_node_based
create_distribution_neighbor_based
neighbors
bond_intersects_disc
bond_intersect_infinite_plane
bond_intersect_rectangle_plane
apply_bond_filters
disk_filter
rectangular_plane_filter
glob_to_loc
check_for_duplicate_in_dataframe
check_types_in_dataframe
```

```@meta
CurrentModule = PeriLab.IO.Geometry
```
## Geometry
```@docs
bond_geometry
shape_tensor
compute_deformation_gradient
compute_strain
rotation_tensor
```

```@meta
CurrentModule = PeriLab.IO
```
## Read_Input_Deck
```@docs
read_input
read_input_file
```

```@meta
CurrentModule = PeriLab.IO
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
CurrentModule = PeriLab.IO
```
## Write_CSV_Results
```@docs
write_global_results_in_csv
```

```@meta
CurrentModule = PeriLab.Solver.Helpers
```
## Helpers
```@docs
find_indices
find_active
get_active_update_nodes
find_files_with_ending
check_inf_or_nan
matrix_style
progress_bar
get_fourth_order
find_inverse_bond_id
rotate
rotate_second_order_tensor
```

```@meta
CurrentModule = PeriLab.Solver.Model_Factory.Pre_calculation
```
## Pre_calculation
```@docs
compute
init_pre_calculation
Bond_Deformation.compute
Deformation_Gradient.compute
Shape_Tensor.compute
```

```@meta
CurrentModule = PeriLab.IO.Parameter_Handling
```
## parameter_handling
```@docs
validate_yaml
validate_structure_recursive
get_all_keys
get_density
get_heat_capacity
get_horizon
get_values
get_number_of_blocks
get_block_models
get_computes_names
get_output_variables
get_computes
get_mesh_name
get_bond_filters
get_node_sets
check_for_duplicates
get_output_filenames
get_output_type
get_flush_file
get_write_after_damage
get_output_fieldnames
get_outputs
get_output_frequency
get_model_parameter
get_models_option
get_solver_name
get_initial_time
get_final_time
get_safety_factor
get_fixed_dt
get_numerical_damping
get_max_damage
get_solver_options
get_header
```

```@meta
CurrentModule = PeriLab.Solver.Model_Factory.Material.Set_modules
```
## Set_modules
```@docs
find_jl_files
find_module_files
include_files
create_module_specifics
```

```@meta
CurrentModule = PeriLab.Logging_module
```
## Logging_module
```@docs
set_result_files
print_table
progress_filter
init_logging
```
