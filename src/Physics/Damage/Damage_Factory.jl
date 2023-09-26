# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Damage
include("./Critical_Stretch.jl")
using .Critical_stretch
export compute_damage

function compute_damage(datamanager, nodes, damage_parameter, time, dt)
    bondDamageN = datamanager.get_field("Bond Damage", "N")
    bondDamageNP1 = datamanager.get_field("Bond Damage", "NP1")
    bondDamageNP1 = copy(bondDamageN)
    if damage_parameter["Damage Model"] == "Test"
        return testing_damage(datamanager, time)
    end
    if damage_parameter["Damage Model"] == Critical_stretch.damage_name()
        datamanager = Critical_stretch.compute_damage(datamanager, nodes, damage_parameter, time, dt)
    end
    datamanager = damage_index(datamanager, nodes)
    return datamanager
end
"""
  damage_index(nodes,datamananager)

    Function calculates the damage index related to the neighborhood volume for a set of corresponding nodes. 
    The damage index is defined as damaged volume in relation the neighborhood volume.
    damageIndex = sum_i (brokenBonds_i * volume_i) / volumeNeighborhood

    Parameters:
    - `nodes::Vector{Int64}`: corresponding nodes to this model
    - `datamanager::Data_manager`: all model data

"""
function damage_index(datamananager, nodes)
    nlist = datamananager.get_nlist()
    volume = datamananager.get_field("Volume")
    bondDamageNP1 = datamananager.get_field("Bond Damage", "NP1")
    damage = datamananager.get_field("Damage", "NP1")
    for iID in nodes
        undamaged_volume = sum(volume[nlist[iID][:]])
        totalDamage = sum((1 .- bondDamageNP1[iID][:]) .* volume[nlist[iID][:]])
        damage[iID] = totalDamage / undamaged_volume
    end
    return datamananager

end
end