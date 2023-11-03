# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Critical_stretch
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
   "Critical Stretch"
   ```
   """
function damage_name()
    return "Critical Stretch"
end
"""
   compute_damage(datamanager, nodes, damage_parameter, block, time, dt)

   Calculates the stretch of each bond and compares it to a critical one. If it is exceeded, the bond damage value is set to zero.

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

    nlist = datamanager.get_nlist()
    bondDamageNP1 = datamanager.get_field("Bond Damage", "NP1")
    bondGeom = datamanager.get_field("Bond Geometry")
    deformed_bond = datamanager.get_field("Deformed Bond Geometry", "NP1")
    nneighbors = datamanager.get_field("Number of Neighbors")
    blockIds = datamanager.get_field("Block_Id")
    cricital_stretch = damage_parameter["Critical Value"]
    interBlockDamage::Bool = false
    if haskey(damage_parameter, "Interblock Damage")
        interBlockDamage = damage_parameter["Interblock Damage"]
    end
    if interBlockDamage
        inter_critical_stretch = datamanager.get_crit_values_matrix()
    end
    for iID in nodes
        for jID in nneighbors[iID]
            stretch = (deformed_bond[iID][jID, end] - bondGeom[iID][jID, end]) / bondGeom[iID][jID, end]
            crit_stretch = cricital_stretch
            if interBlockDamage
                crit_stretch = inter_critical_stretch[blockIds[iID], blockIds[nlist[iID][jID]], block]
            end
            if stretch > crit_stretch
                bondDamageNP1[iID][jID] = 0
            end
        end
    end
    return datamanager
end
function compute_damage_pre_calculation(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, block::Int64, synchronise_field, time::Float64, dt::Float64)
    return datamanager
end
end