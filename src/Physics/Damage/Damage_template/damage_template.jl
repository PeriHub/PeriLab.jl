# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Damage_template
export compute_damage
export damage_name
"""
   damage_name()

   Gives the damage name. It is needed for comparison with the yaml input deck.

   Parameters:

   Returns:
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
   compute_damage(datamanager, nodes, damage_parameter, time, dt)

   Calculates the damage criterion of each bond. This template has to be copied, the file renamed and edited by the user to create a new material. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

   Parameters:
        - `datamanager::Data_manager`: Datamanager.
        - `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
        - `damage_parameter::Dict(String, Any)`: Dictionary with material parameter.
        - `time::Float64`: The current time.
        - `dt::Float64`: The current time step.
   Returns:
        - - `datamanager::Data_manager`: Datamanager.
   Example:
   ```julia
     ```
   """
function compute_damage(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, damage_parameter::Dict, time::Float64, dt::Float64)
  @info "Please write a damage model name in damage_name()."
  @info "You can call your routine within the yaml file."
  @info "Fill the compute_force(datamanager, nodes, damage_parameter, time, dt) function."
  @info "The datamanager and damage_parameter holds all you need to solve your problem on material level."
  @info "add own files and refer to them. If a module does not exist. Add it to the project or contact the developer."
  return datamanager
end

end