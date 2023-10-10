# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Critical_Energy_Correpondence_Model
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
   "Critical Energy Correspondence"
   ```
   """
function damage_name()
  return "Critical Energy Correspondence"
end
"""
   compute_damage(datamanager, nodes, damage_parameter, time, dt)

   Calculates the force densities of the material. This template has to be copied, the file renamed and edited by the user to create a new material. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

   Parameters:
        - `datamanager::Data_manager`: Datamanager.
        - `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
        - `damage_parameter::Dict(String, Any)`: Dictionary with material parameter.
        - `time::Float32`: The current time.
        - `dt::Float32`: The current time step.
   Returns:
        - - `datamanager::Data_manager`: Datamanager.
   Example:
   ```julia
     ```
   """
function compute_damage(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, damage_parameter::Dict, time::Float32, dt::Float32)
  # synch of bond force missing
  dof::Int64 = datamanager.get_dof()
  update_list = datamanager.get_field("Update List")
  horizon = datamanager.get_field("Horizon")
  bond_damage = datamanager.get_field("Bond Damage", "NP1")
  deformed_bondN = datamanager.get_field("Deformed Bond", "N")
  deformed_bondNP1 = datamanager.get_field("Deformed Bond", "NP1")
  bond_force = datamanager.get_field("Bond Forces")


  # parameter
  critical_Energy = damage_parameter["Critical Energy"]
  tension::Bool = false
  if haskey(damage_parameter, "Only Tension")
    tension = damage_parameter["Only Tension"]
  end
  nneighbors = datamanager.get_field("Number of Neighbors")
  bond_energy::Float32 = 0.0

  for iID in nodes
    update_list[iID] = false
    for jID in 1:nneighbors[iID]
      if tension
        dist = deformed_bondNP1[iID][jID, end] - deformed_bondN[iID][jID, end]
        if dist < 0
          continue
        end
      end

      bond_energy = 0.5 * sum(deformed_bondNP1[iID][jID, 1:end-1] * bond_force[iID][jID, 1:end])
      if critical_Energy < bond_energy / get_quad_horizon(horizon[iID], dof)
        bond_damage[iID][jID] = 0.0
        update_list[iID] = true
      end
    end
    # pre calculation of shape Tensor and defGrad
    # possibility that only an update of the existing is needed which saves time -> change list
    return datamanager
  end














  return datamanager
end

function get_quad_horizon(horizon::Float32, dof::Int64)
  if dof == 2
    thickness::Float32 = 1 # is a placeholder
    return Float32(3 / (pi * horizon^3 * thickness))
  end
  return Float32(4 / (pi * horizon^4))
end
end