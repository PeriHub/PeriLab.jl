# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Correspondence_Plastic
using LinearAlgebra
export fe_support
export init_material_model
export material_name
export compute_forces

"""
  fe_support()

Gives the information if the material supports the FEM part of PeriLab

# Arguments

# Returns
- bool: true - for FEM support; false - for no FEM support

Example:
```julia
println(fe_support())
false
```
"""
function fe_support()
  return false
end

"""
  init_material_model(datamanager::Module)

Initializes the material model.

# Arguments
  - `datamanager::Data_manager`: Datamanager.

# Returns
  - `datamanager::Data_manager`: Datamanager.
"""
function init_material_model(datamanager::Module)

  datamanager.create_node_field("von Mises Stresses", Float64, "Matrix", 1)
  return datamanager
end
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
  return "Plastic"
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

  """
  dof::Int64 = datamanager.get_dof()
  stress = datamanager.get_field("Cauchy Stresses", "NP1")
  von_Mises_stress = datamanager.get_field("von Mises Stresses", "NP1")
  yield_stress::Float64 = material_parameter["Yield Stress"]
  sphericalStressNP1::Float64 = 0
  deviatoricStressNP1::Matrix{Float64} = zeros(dof, dof)
  for iID in nodes
    sphericalStressNP1 = sum(stress[iID][i][i] for i in 1:dof) / 3
    deviatoricStressNP1 = stress[iID] - sphericalStressNP1 .* I(dof)
  end
  stressPlasticNP1 = 0 #->field

  von_Mises_stress[iID] = sqrt(3.0 / 2.0 * sum(deviatoricStressNP1[:, :]))

  #reducedYieldStress = CORRESPONDENCE::FLAWFUNCTION(isFlaw,yieldStress, flawMagnitude, flawSize, modelCoord, flawLocationX, flawLocationY, flawLocationZ, 1);
  reduced_yield_stress = yield_stress
  if von_Mises_stress[iID] > reduced_yield_stress
    println()
  end

 """
  return datamanager
end



end