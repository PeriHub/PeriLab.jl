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
read_input_file
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
CurrentModule = PeriLab.IO.Read_Mesh
```
## Read_Mesh
```@docs
init_data
get_local_overlap_map
local_nodes_from_dict
distribute_neighborhoodlist_to_cores
get_local_neighbors
get_bond_geometry
define_nsets
distribution_to_cores
check_mesh_elements
read_mesh_from_txt
read_mesh_from_exodus
set_dof
load_and_evaluate_mesh
create_neighborhoodlist
get_number_of_neighbornodes
node_distribution
_init_overlap_map_
create_overlap_map
create_base_chunk
neighbors
bondIntersectsDisk
bondIntersectInfinitePlane
bondIntersectRectanglePlane
apply_bond_filters
disk_filter
rectangular_plane_filter
glob_to_loc
```

```@meta
CurrentModule = PeriLab.IO.Read_Mesh.Geometry
```
## Geometry
```@docs
bond_geometry
shape_tensor
deformation_gradient
strain
rotation_tensor
angle_between_vectors
```

```@meta
CurrentModule = PeriLab.IO.Read_Input_Deck
```
## Read_Input_Deck
```@docs
read_input
read_input_file
```

```@meta
CurrentModule = PeriLab.IO
```
## MPI
```@docs
send_single_value_from_vector
synch_responder_to_controller
synch_controller_to_responder
synch_controller_bonds_to_responder
split_vector
synch_controller_bonds_to_responder_flattened
send_vectors_to_cores
send_vector_from_root_to_core_i
send_value
get_vector
recv_vector_from_root
find_and_set_core_value_min
find_and_set_core_value_max
find_and_set_core_value_sum
find_and_set_core_value_avg
gather_values
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
## Write_CSV_Results
```@docs
create_result_file
write_global_results_in_csv
```

```@meta
CurrentModule = PeriLab.Solver
```
## Helpers
```@docs
find_indices
find_active
get_header
find_files_with_ending
check_inf_or_nan
matrix_style
progress_bar
get_fourth_order
find_inverse_bond_id
```

```@meta
CurrentModule = PeriLab.Solver.Physics.Pre_calculation
```
## Pre_calculation
```@docs
compute
init_pre_calculation
Bond_Deformation.compute
Bond_Deformation_Gradient.compute
Bond_Shape_Tensor.compute
Deformation_Gradient.compute
Shape_Tensor.compute
```

```@meta
CurrentModule = PeriLab.IO
```
## parameter_handling
```@docs
validate_yaml
validate_structure_recursive
get_all_keys
get_density
get_heatcapacity
get_horizon
get_values
get_number_of_blocks
get_block_models
get_computes_names
get_output_variables
get_computes
get_node_set
get_mesh_name
get_topology_name
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
get_physics_option
get_solver_name
get_initial_time
get_final_time
get_safety_factor
get_fixed_dt
get_numerical_damping
get_max_damage
get_solver_options
```

```@meta
CurrentModule = PeriLab.Solver.Physics.Material.Set_modules
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
progress_filter
init_logging
```