# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Critical_Energy_Model
include("../Material/Material_Factory.jl")
using .Material
using LinearAlgebra
export compute_damage
export compute_damage_pre_calculation
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
   compute_damage(datamanager, nodes, damage_parameter, block, time, dt)

   Calculates the elastic energy of each bond and compares it to a critical one. If it is exceeded, the bond damage value is set to zero.
   [WillbergC2019](@cite), [FosterJT2011](@cite)

   Parameters:
        - `datamanager::Data_manager`: Datamanager.
        - `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
        - `damage_parameter::Dict(String, Any)`: Dictionary with material parameter.
        - `block::Int64`: Block number.
        - `time::Float64`: The current time.
        - `dt::Float64`: The current time step.
   Returns:
        - - `datamanager::Data_manager`: Datamanager.
   Example:
   ```julia
     ```
   """
function compute_damage(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, damage_parameter::Dict, block::Int64, time::Float64, dt::Float64)
  dof::Int64 = datamanager.get_dof()
  nlist = datamanager.get_nlist()
  blockList = datamanager.get_block_list()
  blockIds = datamanager.get_field("Block_Id")
  update_list = datamanager.get_field("Update List")
  horizon = datamanager.get_field("Horizon")
  bond_damage = datamanager.get_field("Bond Damage", "NP1")
  bondGeometry = datamanager.get_field("Bond Geometry")
  forceDensities = datamanager.get_field("Force Densities", "NP1")
  deformed_bond = datamanager.get_field("Deformed Bond Geometry", "NP1")
  critical_Energy = damage_parameter["Critical Value"]
  tension::Bool = false
  interBlockDamage::Bool = false
  bond_energy::Float64 = 0.0
  dist::Float64 = 0.0
  projected_force = zeros(Float64, dof)
  if haskey(damage_parameter, "Only Tension")
    tension = damage_parameter["Only Tension"]
  end
  if haskey(damage_parameter, "Interblock Damage")
    interBlockDamage = damage_parameter["Interblock Damage"]
  end
  if interBlockDamage
    inter_critical_Energy = datamanager.get_crit_values_matrix()
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

      projected_force = dot((forceDensities[iID, :] - forceDensities[nlist[iID][jID], :]), deformed_bond[iID][jID, 1:dof]) / (deformed_bond[iID][jID, end] * deformed_bond[iID][jID, end]) .* deformed_bond[iID][jID, 1:dof]

      bond_energy = 0.5 * sum(projected_force[1:dof] .* deformed_bond[iID][jID, 1:dof])
      crit_energy = critical_Energy
      if interBlockDamage
        crit_energy = inter_critical_Energy[blockIds[iID], blockIds[nlist[iID][jID]], block]
      end
      if crit_energy < bond_energy / get_quad_horizon(horizon[iID], dof)
        bond_damage[iID][jID] = 0.0
        update_list[iID] = true
      end
    end
  end
  return datamanager
end

function compute_damage_pre_calculation(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, block::Int64, synchronise_field, time::Float64, dt::Float64)

  datamanager = Material.compute_forces(datamanager, nodes, datamanager.get_properties(block, "Material Model"), time, dt)
  datamanager = Material.distribute_force_densities(datamanager, nodes)

  synchronise_field(datamanager.get_comm(), datamanager.get_synch_fields(), datamanager.get_overlap_map(), datamanager.get_field, "Force DensitiesNP1", "download_from_cores")
  synchronise_field(datamanager.get_comm(), datamanager.get_synch_fields(), datamanager.get_overlap_map(), datamanager.get_field, "Force DensitiesNP1", "upload_to_cores")
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