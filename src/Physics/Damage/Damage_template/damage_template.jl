# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Damage_template
export compute_damage
export compute_damage_pre_calculation
export damage_name
"""
    damage_name()

Gives the damage name. It is needed for comparison with the yaml input deck.

# Arguments

# Returns
- `name::String`: The name of the damage.

Example:
```julia
println(damage_name())
"Damage Template"
```
"""
function damage_name()
  return "Damage Template"
end

"""
c ompute_damage(datamanager, nodes, damage_parameter, block, time, dt)

Calculates the damage criterion of each bond. This template has to be copied, the file renamed and edited by the user to create a new material. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `damage_parameter::Dict(String, Any)`: Dictionary with material parameter.
- `block::Int64`: Block number
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
# Returns
- `datamanager::Data_manager`: Datamanager.
Example:
```julia
```
"""
function compute_damage(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, damage_parameter::Dict, block::Int64, time::Float64, dt::Float64)
  @info "Please write a damage model name in damage_name()."
  @info "You can call your routine within the yaml file."
  @info "Fill the compute_damage(datamanager, nodes, damage_parameter, block, time, dt) function."
  @info "The datamanager and damage_parameter holds all you need to solve your problem on material level."
  @info "Add own files and refer to them. If a module does not exist. Add it to the project or contact the developer."
  return datamanager
end

"""
    compute_damage_pre_calculation(datamanager, nodes, synchronise_field, time, dt)

Pre computes values needed for the damage calculation.Calculates the damage criterion of each bond. This template has to be copied, the file renamed and edited by the user to create a new damages. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `block::Int64`: Block number
- `synchronise_field`: Synchronise function to distribute parameter through cores.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
# Returns
- `datamanager::Data_manager`: Datamanager.
Example:
```julia
```
"""
function compute_damage_pre_calculation(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, block::Int64, synchronise_field, time::Float64, dt::Float64)
  @info "Please write a damage model name in damage_name()."
  @info "You can call your routine within the yaml file."
  @info "Fill the compute_damage_pre_calculation(datamanager, nodes, damage_parameter, time, dt) function if needed."
  @info "The datamanager and damage_parameter holds all you need to solve your problem on material level."
  @info "Add own files and refer to them. If a module does not exist. Add it to the project or contact the developer."
  return datamanager
end
end