# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Damage
include("../../Core/Module_inclusion/set_Modules.jl")
using .Set_modules
global module_list = Set_modules.find_module_files(@__DIR__, "damage_name")
Set_modules.include_files(module_list)

export compute_damage

function compute_damage(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, model_param::Dict, time::Float64, dt::Float64)
    bondDamageN = datamanager.get_field("Bond Damage", "N")
    bondDamageNP1 = datamanager.get_field("Bond Damage", "NP1")
    bondDamageNP1 = copy(bondDamageN)
    update_list = datamanager.get_field("Update List")
    update_list .= true

    specifics = Dict{String,String}("Call Function" => "compute_damage", "Name" => "damage_name")
    datamanager = Set_modules.create_module_specifics(model_param["Damage Model"], module_list, specifics, (datamanager, nodes, model_param, time, dt))
    if isnothing(datamanager)
        @error "No damage model of name " * model_param["Damage Model"] * " exists."
    end

    datamanager = damage_index(datamanager, nodes)
    return datamanager
end
"""
  damage_index(datamananager,::Union{SubArray, Vector{Int64})

    Function calculates the damage index related to the neighborhood volume for a set of corresponding nodes. 
    The damage index is defined as damaged volume in relation the neighborhood volume.
    damageIndex = sum_i (brokenBonds_i * volume_i) / volumeNeighborhood

    Parameters:
    - `datamanager::Data_manager`: all model data
    - `nodes::Union{SubArray, Vector{Int64}}`: corresponding nodes to this model

"""
function damage_index(datamananager::Module, nodes::Union{SubArray,Vector{Int64}})
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