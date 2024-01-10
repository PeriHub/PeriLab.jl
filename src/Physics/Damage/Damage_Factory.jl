# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Damage
include("../../Core/Module_inclusion/set_Modules.jl")
using .Set_modules
global module_list = Set_modules.find_module_files(@__DIR__, "damage_name")
Set_modules.include_files(module_list)

export compute_damage
export compute_damage_pre_calculation
export init_interface_crit_values
export init_damage
"""
    compute_damage(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, model_param::Dict, block::Int64, time::Float64, dt::Float64)

Computes the damage model

# Arguments
- `datamanager::Module`: The datamanager
- `nodes::Union{SubArray,Vector{Int64}}`: The nodes
- `model_param::Dict`: The model parameters
- `block::Int64`: The block
- `time::Float64`: The current time
- `dt::Float64`: The time step
# Returns
- `datamanager::Module`: The datamanager
"""
function compute_damage(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, model_param::Dict, block::Int64, time::Float64, dt::Float64)

    specifics = Dict{String,String}("Call Function" => "compute_damage", "Name" => "damage_name")
    datamanager = Set_modules.create_module_specifics(model_param["Damage Model"], module_list, specifics, (datamanager, nodes, model_param, block, time, dt))
    if isnothing(datamanager)
        @error "No damage model of name " * model_param["Damage Model"] * " exists."
        return nothing
    end
    datamanager = damage_index(datamanager, nodes)
    return datamanager
end

"""
    compute_damage_pre_calculation(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, block::Int64, synchronise_field, time::Float64, dt::Float64)

Compute the pre calculation for the damage.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `block::Int64`: Block number
- `synchronise_field`: Synchronise function to distribute parameter through cores.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function compute_damage_pre_calculation(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, block::Int64, model_param::Dict, synchronise_field, time::Float64, dt::Float64)
    specifics = Dict{String,String}("Call Function" => "compute_damage_pre_calculation", "Name" => "damage_name")
    datamanager = Set_modules.create_module_specifics(model_param["Damage Model"], module_list, specifics, (datamanager, nodes, block, synchronise_field, time, dt))
    if isnothing(datamanager)
        @error "No damage model of name " * model_param["Damage Model"] * " exists."
        return nothing
    end
    return datamanager
end

"""
    damage_index(datamananager,::Union{SubArray, Vector{Int64})

Function calculates the damage index related to the neighborhood volume for a set of corresponding nodes. 
The damage index is defined as damaged volume in relation the neighborhood volume.
damageIndex = sum_i (brokenBonds_i * volume_i) / volumeNeighborhood

# Arguments
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

"""
    set_bond_damage(datamanager::Module, nodes::Union{SubArray,Vector{Int64}})

Set the bond damage field to the bond damage field

# Arguments
- `datamanager::Module`: The datamanager
- `nodes::Union{SubArray,Vector{Int64}}`: The nodes
# Returns
- `datamanager::Module`: The datamanager
"""
function set_bond_damage(datamanager::Module, nodes::Union{SubArray,Vector{Int64}})
    bond_damageN = datamanager.get_field("Bond Damage", "N")
    bond_damageNP1 = datamanager.get_field("Bond Damage", "NP1")
    for iID in nodes
        bond_damageNP1[iID][:] = bond_damageN[iID][:]
    end
    return datamanager
end

"""
    init_interface_crit_values(datamanager::Module, params::Dict, block_id::Int64)

Initialize the critical values

# Arguments
- `datamanager::Module`: The datamanager
- `params::Dict`: The parameters
- `block_id::Int64`: current block
# Returns
- `datamanager::Module`: The datamanager
"""
function init_interface_crit_values(datamanager::Module, damage_parameter::Dict, block_id::Int64)
    max_block_id = maximum(datamanager.get_block_list())
    inter_critical_value = zeros(Float64, (max_block_id, max_block_id, max_block_id))
    if !haskey(damage_parameter, "Interblock Damage")
        return datamanager
    end
    critical_value = damage_parameter["Critical Value"]
    for block_iId in 1:datamanager.get_block_list()
        for block_jId in 1:datamanager.get_block_list()
            critical_value_name = "Interblock Critical Value $(block_iId)_$block_jId"
            if haskey(damage_parameter["Interblock Damage"], critical_value_name)
                inter_critical_value[block_iId, block_jId, block_id] = damage_parameter["Interblock Damage"][critical_value_name]
            else
                inter_critical_value[block_iId, block_jId, block_id] = critical_value
            end
        end
    end
    datamanager.set_crit_values_matrix(inter_critical_value)
    return datamanager
end

"""
    init_aniso_crit_values(datamanager::Module, params::Dict, block_id::Int64)

Initialize the anisotropic critical values

# Arguments
- `datamanager::Module`: The datamanager
- `params::Dict`: The parameters
- `block_id::Int64`: current block
# Returns
- `datamanager::Module`: The datamanager
"""
function init_aniso_crit_values(datamanager::Module, damage_parameter::Dict, block_id::Int64)
    aniso_crit::Dict{Int64,Any} = Dict()

    crit_0 = damage_parameter["Critical Value"]
    crit_90 = damage_parameter["Critical Value"]
    if !haskey(damage_parameter, "Anisotropic Damage")
        return datamanager
    end
    crit_0 = damage_parameter["Anisotropic Damage"]["Critical Value X"]
    crit_90 = damage_parameter["Anisotropic Damage"]["Critical Value Y"]
    aniso_crit[block_id] = [crit_0, crit_90]

    datamanager.set_aniso_crit_values(aniso_crit)
    return datamanager
end

function init_damage_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, block::Int64)

    model_param = datamanager.get_properties(block, "Damage Model")
    specifics = Dict{String,String}("Call Function" => "init_damage_model", "Name" => "damage_name")
    if !haskey(model_param, "Damage Model")
        return datamanager
    end

    datamanager = Set_modules.create_module_specifics(model_param["Damage Model"], module_list, specifics, (datamanager, nodes, model_param, block))
    if isnothing(datamanager)
        @error "No damage model of name " * model_param["Damage Model"] * " exists."
        return nothing
    end

    datamanager = Damage.init_interface_crit_values(datamanager, model_param, block)
    datamanager = Damage.init_aniso_crit_values(datamanager, model_param, block)
    return datamanager
end
end