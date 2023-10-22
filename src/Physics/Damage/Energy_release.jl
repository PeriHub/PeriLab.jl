# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Critical_Energy_Model
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
   "Critical Energy"
   ```
   """
function damage_name()
  return "Critical Energy"
end
"""
   compute_damage(datamanager, nodes, damage_parameter, time, dt)

   Calculates the elastic energy of each bond and compares it to a critical one. If it is exceeded, the bond damage value is set to zero.
   [WillbergC2019](@cite), [FosterJT2011](@cite)

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
  dof::Int64 = datamanager.get_dof()
  nlist = datamanager.get_nlist()
  datamanager.synch_field("Force Densities")
  update_list = datamanager.get_field("Update List")
  horizon = datamanager.get_field("Horizon")
  bond_damage = datamanager.get_field("Bond Damage", "NP1")
  bondGeometry = datamanager.get_field("Bond Geometry")
  forceDensities = datamanager.get_field("Force Densities", "NP1")
  deformed_bond = datamanager.get_field("Deformed Bond Geometry", "NP1")
  critical_Energy = damage_parameter["Critical Energy"]
  tension::Bool = false
  bond_energy::Float64 = 0.0
  dist::Float64 = 0.0
  if haskey(damage_parameter, "Only Tension")
    tension = damage_parameter["Only Tension"]
  end
  nneighbors = datamanager.get_field("Number of Neighbors")

  for iID in nodes
    for jID in 1:nneighbors[iID]
      if tension
        dist = deformed_bond[iID][jID, end] - bondGeometry[iID][jID, end]
        if dist < 0
          continue
        end
      end
      # 0.25 ?!

      # kraft muss projizizert werden und synchronisiert
      bond_energy = 0.5 * sum((forceDensities[iID] - forceDensities[nlist[iID][jID]]) .* deformed_bond[iID][jID, 1:end] .* deformed_bond[iID][jID, 1:end] ./ deformed_bond[iID][jID, end])
      if critical_Energy < bond_energy / get_quad_horizon(horizon[iID], dof)
        bond_damage[iID][jID] = 0.0
        update_list[iID] = true
      end
    end
  end
  return datamanager
end

function get_quad_horizon(horizon::Float64, dof::Int64)
  if dof == 2
    thickness::Float64 = 1 # is a placeholder
    return Float64(3 / (pi * horizon^3 * thickness))
  end
  return Float64(4 / (pi * horizon^4))
end
end