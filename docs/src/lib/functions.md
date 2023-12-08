# Functions

## Index
```@index
Pages = ["functions.md"]
```

```@meta
CurrentModule = PeriLab.Data_manager
```
## Data_manager
```@docs
get_comm
set_comm
check_property
create_bond_field
create_constant_bond_field
create_constant_node_field
create_field
create_node_field
get_all_field_keys
get_block_list
get_crit_values_matrix
get_dof
get_field
get_field_type
get_inverse_nlist
get_local_nodes
get_nlist
get_nnodes
get_NP1_to_N_Dict
get_nnsets
get_nsets
get_num_responder
get_overlap_map
get_synch_fields
get_physics_options
get_properties
get_property
get_rank
get_max_rank
loc_to_glob
init_property
rotation_data
set_block_list
set_crit_values_matrix
set_distribution
set_dof
set_glob_to_loc
set_inverse_nlist
set_nnodes
set_num_controller
set_nnsets
set_nset
set_num_responder
set_overlap_map
set_physics_options
set_property
set_properties
set_rank
set_max_rank
set_synch
switch_NP1_to_N
synch_manager
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

```@meta
CurrentModule = PeriLab.Solver.Verlet
```
## Verlet
```@docs
compute_thermodynamic_critical_time_step
get_cs_denominator
compute_mechanical_critical_time_step
test_timestep
compute_crititical_time_step
init_solver
get_integration_steps
run_solver
```

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