# Functions

## Multithreading support

By default, selected operations in DataFrames.jl automatically use multiple threads
when available. Multi-threading is task-based and implemented using the `@spawn`
macro from Julia Base. Tasks are therefore scheduled on the `:default` threadpool.
Functions that take user-defined functions and may run it in parallel
accept a `threads` keyword argument which allows disabling multithreading
when the provided function requires serial execution or is not thread-safe.

This is a list of operations that currently make use of multi-threading:
- `DataFrame` constructor with `copycols=true`; also recursively all functions
  that call this constructor, e.g. `copy`.
- `getindex` when multiple columns are selected.
- `groupby` (both when hashing is required and when fast path using `DataAPI.refpool`
  is used).
- `*join` functions for composing output data frame (but currently not for finding
  matching rows in joined data frames).
- `combine`, `select[!]`, and `transform[!]` on `GroupedDataFrame` when either of the conditions below is met:
  * multiple transformations are performed (each transformation is spawned in a separate task)
  * a transformation produces one row per group and the passed transformation
    is a custom function (i.e. not for standard reductions, which use
    optimized single-threaded methods).
- `dropmissing` when the provided data frame has more than 1 column and `view=false` 
  (subsetting of individual columns is spawned in separate tasks).

In general at least Julia 1.4 is required to ensure that multi-threading is used
and the Julia process must be started with more than one thread. Some operations
turn on multi-threading only if enough rows in the processed data frame are present
(the exact threshold when multi-threading is enabled is considered to be undefined
and might change in the future).

Except for the list above, where multi-threading is used automatically,
all functions provided by DataFrames.jl that update a data frame are not thread safe.
This means that while they can be called from any thread, the caller is responsible
for ensuring that a given `DataFrame` object is never modified by one thread while
others are using it (either for reading or writing). Using the same `DataFrame`
at the same time from different threads is safe as long as it is not modified.

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
set_nmasters
set_nslaves
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
