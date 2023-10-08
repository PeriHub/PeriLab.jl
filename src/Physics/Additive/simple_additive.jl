# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module simple_additive
export compute_additive
export additive_name
"""
   additive_name()

   Gives the additive name. It is needed for comparison with the yaml input deck.

   Parameters:

   Returns:
   - `name::String`: The name of the additive model.

   Example:
   ```julia
   println(additive_name())
   "additive Template"
   ```
   """
function additive_name()
  return "Simple"
end
"""
   compute_additive(datamanager, nodes, additive_parameter, time, dt)

   Calculates the force densities of the additive. This template has to be copied, the file renamed and edited by the user to create a new additive. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

   Parameters:
        - `datamanager::Data_manager`: Datamanager.
        - `nodes::Vector{Int64}`: List of block nodes.
        - `additive parameter::Dict(String, Any)`: Dictionary with additive parameter.
        - `time::Float32`: The current time.
        - `dt::Float32`: The current time step.
   Returns:
        - - `datamanager::Data_manager`: Datamanager.
   Example:
   ```julia
     ```
   """
function compute_additive(datamanager, nodes, additive_parameter, time, dt)
  flux = datamanager.get_field("Temperature Flux", "NP1")
  activation_time = datamanager_get_field("Activation Time")
  bond_damage = datamanager_get_field("Bond Damage", "NP1")
  active = datamanager.get_field("Active")
  heatCapacity = additive_parameter["Heat Capacity"]
  density = additive_parameter["Density"]
  printTemperature = additive_parameter["Print Temperature"]

  for iID in nodes
    if time - dt <= activation_time[iID] < time
      active[iID] = true
      bond_damage[iID][:, :] = 1.0
      flux[iID] = -printTemperature * heatCapacity * density / dt
    end
  end

  return datamanager
end

end