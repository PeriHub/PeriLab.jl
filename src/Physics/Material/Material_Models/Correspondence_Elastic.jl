# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Correspondence_Elastic
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
  return "Correspondence_Elastic"
end
"""
   compute_force(datamanager, nodes, material_parameter, time, dt)

   Calculates the force densities of the material. This template has to be copied, the file renamed and edited by the user to create a new material. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

   Parameters:
        - `datamanager::Data_manager`: Datamanager.
        - `nodes::Vector{Int64}`: List of block nodes.
        - `material_parameter::Dict(String, Any)`: Dictionary with material parameter.
        - `time::Float32`: The current time.
        - `dt::Float32`: The current time step.
   Returns:
        - - `datamanager::Data_manager`: Datamanager.
   Example:
   ```julia
     ```
   """
function compute_force(datamanager, nodes, material_parameter, time, dt)
  rotation::Bool = false
  if "Angles" in datamanager.get_all_field_keys()
    rotation = true
    angles = datamanager.get_field("Angles")
  end
  dof = datamanager.get_dof()
  strainN = datamanager.get_field("Strain", "N")
  strainNP1 = datamanager.get_field("Strain", "NP1")
  stressN, stressNP1 = datamanager.create_field("Cauchy Stress")
  bond_forceN, bond_forceNP1 = datamanager.create_bond_field("Bond Force")
  force_densities = datamanager.get_field("Force Densities", "NP1")

  hookMatrix = get_Hook_matrix(material_parameter, dof)
  bond_force = elastic(nodes, hookMatrix, angles, rotation, bond_force)

  force_densities[:] = distribute_forces(nodes, nlist, bond_force, volume, force_densities)

  return datamanager
end


function elastic(nodes, hookMatrix, angles, rotation, bond_force)
  for node in nodes
    if rotation

    end
    MATRICES::tensorRotation(angles, sigmaNP1, false, sigmaNP1)

    MATRICES::tensorRotation(angles, sigmaNP1, false, sigmaNP1)
    if rotation

    end

  end
  return bond_force
end

end