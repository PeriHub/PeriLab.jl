# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Damage
include("../../Core/Module_inclusion/set_Modules.jl")
using .Set_modules
using TimerOutputs
global module_list = Set_modules.find_module_files(@__DIR__, "damage_name")
Set_modules.include_files(module_list)

export compute_damage
export compute_damage_pre_calculation
export init_interface_crit_values

function compute_damage(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, model_param::Dict, block::Int64, time::Float64, dt::Float64)

    specifics = Dict{String,String}("Call Function" => "compute_damage", "Name" => "damage_name")
    datamanager = Set_modules.create_module_specifics(model_param["Damage Model"], module_list, specifics, (datamanager, nodes, model_param, block, time, dt))
    if isnothing(datamanager)
        @error "No damage model of name " * model_param["Damage Model"] * " exists."
    end
    datamanager = damage_index(datamanager, nodes)
    return datamanager
end

function compute_damage_pre_calculation(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, block::Int64, model_param::Dict, synchronise_field, time::Float64, dt::Float64)
    specifics = Dict{String,String}("Call Function" => "compute_damage_pre_calculation", "Name" => "damage_name")
    datamanager = Set_modules.create_module_specifics(model_param["Damage Model"], module_list, specifics, (datamanager, nodes, block, synchronise_field, time, dt))
    if isnothing(datamanager)
        @error "No damage model of name " * model_param["Damage Model"] * " exists."
    end
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
    bond_damageNP1 = datamananager.get_field("Bond Damage", "NP1")
    damage = datamananager.get_field("Damage", "NP1")
    for iID in nodes
        undamaged_volume = sum(volume[nlist[iID][:]])
        totalDamage = sum((1 .- bond_damageNP1[iID][:]) .* volume[nlist[iID][:]])
        if damage[iID] < totalDamage / undamaged_volume
            damage[iID] = totalDamage / undamaged_volume
        end
    end
    return datamananager

end

function set_bond_damage(datamanager::Module, nodes::Union{SubArray,Vector{Int64}})
    bond_damageN = datamanager.get_field("Bond Damage", "N")
    bond_damageNP1 = datamanager.get_field("Bond Damage", "NP1")
    for iID in nodes
        bond_damageNP1[iID][:] = bond_damageN[iID][:]
    end
    return datamanager
end

function init_interface_crit_values(datamanager::Module, params::Dict)
    blockList = datamanager.get_block_list()

    inter_critical_Value = zeros(Float64, (length(blockList), length(blockList), length(blockList)))
    for block_id in blockList
        if !haskey(params["Blocks"]["block_$block_id"], "Damage Model")
            continue
        end
        damageName = params["Blocks"]["block_$block_id"]["Damage Model"]
        damage_parameter = params["Physics"]["Damage Models"][damageName]
        if !haskey(damage_parameter, "Interblock Damage")
            continue
        end
        if !damage_parameter["Interblock Damage"]
            continue
        end
        critical_value = damage_parameter["Critical Value"]
        for block_iId in blockList
            for block_jId in blockList
                critValueName = "Interblock Critical Value $(block_iId)_$block_jId"
                if haskey(damage_parameter, critValueName)
                    inter_critical_Value[block_iId, block_jId, block_id] = damage_parameter[critValueName]
                else
                    inter_critical_Value[block_iId, block_jId, block_id] = critical_value
                end
            end
        end
    end
    datamanager.set_crit_values_matrix(inter_critical_Value)

    return datamanager
end
end