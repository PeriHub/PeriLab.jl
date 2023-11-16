# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Correspondence_Elastic
include("../material_basis.jl")
export compute_stresses
export correspondence_name
"""
   correspondence_name()

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
function correspondence_name()
  return "Correspondence Elastic"
end
"""
   compute_stresses(datamanager, nodes, material_parameter, time, dt)

   Calculates the force densities of the material. This template has to be copied, the file renamed and edited by the user to create a new material. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

   Parameters:
        - `datamanager::Data_manager`: Datamanager.
        - `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
        - `dof::Int64`: Degrees of freedom
        - `material_parameter::Dict(String, Any)`: Dictionary with material parameter.
        - `time::Float64`: The current time.
        - `dt::Float64`: The current time step.
        - `strainInc::SubArray`: Strain increment.
        - `stressN::SubArray`: Stress of step N.
        - `stressNP1::SubArray`: Stress of step N+1.
   Returns:
        - `datamanager::Data_manager`: Datamanager.
        - `stressNP1::SubArray`: updated stresses
   Example:
   ```julia
     ```
   """
function compute_stresses(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, dof::Int64, material_parameter::Dict, time::Float64, dt::Float64, strainInc::SubArray, stressN::SubArray, stressNP1::SubArray)

  hookeMatrix = get_Hooke_matrix(material_parameter, material_parameter["Symmetry"], dof)

  for iID in nodes
    stressNP1[iID, :, :] = voigt_to_matrix(hookeMatrix * matrix_to_voigt(strainInc[iID, :, :])) + stressN[iID, :, :]
  end

  return stressNP1, datamanager
end




end