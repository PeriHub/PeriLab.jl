# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Pre_calculation_template
export compute_Pre_calculation
export Pre_calculation_name
"""
   pre_calculation_name()

   Gives the pre_calculation name. It is needed for comparison with the yaml input deck.

   Parameters:

   Returns:
   - `name::String`: The name of the Pre_calculation.

   Example:
   ```julia
   println(pre_calculation_name())
   "Pre_calculation Template"
   ```
   """
function pre_calculation_name()
  return "pre_calculation Template"
end
"""
   pre_calculation(datamanager, nodes, Pre_calculation_parameter, time, dt)

   Calculates the force densities of the material. This template has to be copied, the file renamed and edited by the user to create a new material. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

   Parameters:
        - `datamanager::Data_manager`: Datamanager.
        - `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
        - `Pre_calculation_parameter::Dict(String, Any)`: Dictionary with material parameter.
        - `time::Float32`: The current time.
        - `dt::Float32`: The current time step.
   Returns:
        - - `datamanager::Data_manager`: Datamanager.
   Example:
   ```julia
     ```
   """

function pre_calculation(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, Pre_calculation_parameter::Dict, time::Float32, dt::Float32)
  @info "Please write a possible precalculation routines in pre_calculation_name()."
  @info "You can call your routine within the yaml file."
  @info "Fill the compute_force(datamanager, nodes, Pre_calculation_parameter, time, dt) function."
  @info "The datamanger and Pre_calculation_parameter holds all you need to solve your problem on material level."
  @info "add own files and refer to them. If a module does not exist. Add it to the project or contact the developer."
  return datamanager
end

end