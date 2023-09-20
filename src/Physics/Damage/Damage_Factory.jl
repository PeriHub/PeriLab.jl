module Damage
export get_damage

function get_damage(nodes, params, datamanager)
    if damage["Damage Model"] == "Test"
        return testing_damage(datamanager, time)
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
function damage_index(nodes, datamananager)
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