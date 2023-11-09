# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Thermal_template
export compute_thermal_model
export thermal_model_name
"""
   thermal_model_name()

   Gives the thermal model name. It is needed for comparison with the yaml input deck.

   Parameters:

   Returns:
   - `name::String`: The name of the thermal flow model.

   Example:
   ```julia
   println(flow_name())
   "Thermal Template"
   ```
   """
function thermal_model_name()
  return "Thermal Template"
end
"""
   compute_thermal_model(datamanager, nodes, thermal_parameter, time, dt)

   Calculates the thermal behavior of the material. This template has to be copied, the file renamed and edited by the user to create a new flow. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

   Parameters:
        - `datamanager::Data_manager`: Datamanager.
        - `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
        - `flow parameter::Dict(String, Any)`: Dictionary with flow parameter.
        - `time::Float64`: The current time.
        - `dt::Float64`: The current time step.
   Returns:
        - - `datamanager::Data_manager`: Datamanager.
   Example:
   ```julia
     ```
   """
function compute_thermal_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, thermal_parameter::Dict, time::Float64, dt::Float64)
  @info "Please write a thermal model name in thermal_name()."
  @info "You can call your routine within the yaml file."
  @info "Fill the compute_thermal_model(datamanager, nodes, thermal_parameter, time, dt) function."
  @info "The datamanager and thermal_parameter holds all you need to solve your problem on thermal flow level."
  @info "Add own files and refer to them. If a module does not exist. Add it to the project or contact the developer."
  return datamanager
end

end