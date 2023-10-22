# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Material_template
export compute_force
export material_name
"""
   material_name()

   Gives the material name. It is needed for comparison with the yaml input deck.

   Parameters:

   Returns:
   - `name::String`: The name of the material.

   Example:
   ```julia
   println(material_name())
   "Material Template"
   ```
   """
function material_name()
  return "Material Template"
end
"""
   compute_force(datamanager, nodes, material_parameter, time, dt)

   Calculates the force densities of the material. This template has to be copied, the file renamed and edited by the user to create a new material. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

   Parameters:
        - `datamanager::Data_manager`: Datamanager.
        - `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
        - `material_parameter::Dict(String, Any)`: Dictionary with material parameter.
        - `time::Float64`: The current time.
        - `dt::Float64`: The current time step.
   Returns:
        - - `datamanager::Data_manager`: Datamanager.
   Example:
   ```julia
     ```
   """
function compute_force(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, material_parameter::Dict, time::Float64, dt::Float64)
  @info "Please write a material name in material_name()."
  @info "You can call your routine within the yaml file."
  @info "Fill the compute_force(datamanager, nodes, material_parameter, time, dt) function."
  @info "The datamanager and material_parameter holds all you need to solve your problem on material level."
  @info "add own files and refer to them. If a module does not exist. Add it to the project or contact the developer."
  return datamanager
end

end