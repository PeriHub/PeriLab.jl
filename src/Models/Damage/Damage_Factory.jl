# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Damage
include("../../Core/Module_inclusion/set_Modules.jl")
global module_list = find_module_files(@__DIR__, "damage_name")
include_files(module_list)
using TimerOutputs
using LoopVectorization
using .....Helpers: find_inverse_bond_id
export fields_for_local_synchronization
export compute_model
export init_interface_crit_values
export init_model
export init_fields

"""
    init_fields(datamanager::Module)

Initialize damage model fields

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `params::Dict`: Parameters.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function init_fields(datamanager::Module)
    dof = datamanager.get_dof()
    datamanager.create_node_field("Damage", Float64, 1)

    anisotropic_damage = false

    nlist = datamanager.get_nlist()
    inverse_nlist = datamanager.set_inverse_nlist(find_inverse_bond_id(nlist))
    return datamanager
end

"""
    compute_model(datamanager::Module, nodes::AbstractVector{Int64}, model_param::Dict, block::Int64, time::Float64, dt::Float64,to::TimerOutput,)

Computes the damage model

# Arguments
- `datamanager::Module`: The datamanager
- `nodes::AbstractVector{Int64}`: The nodes
- `model_param::Dict`: The model parameters
- `block::Int64`: The block
- `time::Float64`: The current time
- `dt::Float64`: The time step
# Returns
- `datamanager::Module`: The datamanager
"""
function compute_model(datamanager::Module,
                       nodes::AbstractVector{Int64},
                       model_param::Dict,
                       block::Int64,
                       time::Float64,
                       dt::Float64,
                       to::TimerOutput)
    mod = datamanager.get_model_module(model_param["Damage Model"])

    @timeit to "run model" datamanager=mod.compute_model(datamanager, nodes, model_param,
                                                         block, time, dt)

    if isnothing(datamanager.get_filtered_nlist())
        @timeit to "compute index" return damage_index(datamanager, nodes)
    end

    @timeit to "compute index" return damage_index(datamanager, nodes,
                                                   datamanager.get_filtered_nlist())
end

"""
    fields_for_local_synchronization(datamanager, model, block)

Defines all synchronization fields for local synchronization

# Arguments
- `datamanager::Module`: datamanager.
- `model::String`: Model class.
- `block::Int64`: block ID
# Returns
- `datamanager::Module`: Datamanager.
"""
function fields_for_local_synchronization(datamanager, model, block)
    model_param = datamanager.get_properties(block, "Damage Model")
    mod = datamanager.get_model_module(model_param["Damage Model"])

    return mod.fields_for_local_synchronization(datamanager, model)
end

"""
    damage_index(datamanager,::Union{SubArray, Vector{Int64})

Function calculates the damage index related to the neighborhood volume for a set of corresponding nodes.
The damage index is defined as damaged volume in relation the neighborhood volume.
damageIndex = sum_i (brokenBonds_i * volume_i) / volumeNeighborhood

# Arguments
- `datamanager::Data_manager`: all model data
- `nodes::AbstractVector{Int64}`: corresponding nodes to this model
"""
function damage_index(datamanager::Module,
                      nodes::AbstractVector{Int64},
                      nlist_filtered_ids::Vector{Vector{Int64}})
    bond_damageNP1 = datamanager.get_bond_damage("NP1")
    for iID in nodes
        bond_damageNP1[iID][nlist_filtered_ids[iID]] .= 1
    end
    return damage_index(datamanager, nodes)
end

function damage_index(datamanager::Module, nodes::AbstractVector{Int64})
    nlist = datamanager.get_nlist()::Vector{Vector{Int64}}
    volume = datamanager.get_field("Volume")::AbstractVector{Float64}
    bond_damageNP1 = datamanager.get_bond_damage("NP1")::Vector{Vector{Float64}}
    damage = datamanager.get_damage("NP1")::AbstractVector{Float64}
    compute_index(damage, nodes, volume, nlist, bond_damageNP1)
    return datamanager
end

function compute_index(damage::AbstractVector{Float64},
                       nodes::AbstractVector{Int64},
                       volume::AbstractVector{Float64},
                       nlist::Vector{Vector{Int64}},
                       bond_damage::Vector{Vector{Float64}})::Nothing
    @inbounds @fastmath for iID in nodes
        undamaged_volume = 0.0  # More explicit than zero(Float64)
        totalDamage = 0.0

        # Cache the vectors to help type inference
        neighbors = nlist[iID]
        bonds = bond_damage[iID]

        @inbounds @fastmath for j in eachindex(neighbors)
            jID = neighbors[j]
            vol_j = volume[jID]
            undamaged_volume += vol_j
            totalDamage += (1.0 - bonds[j]) * vol_j
        end

        # More explicit type handling
        current_damage = damage[iID]
        threshold = current_damage * undamaged_volume
        if threshold < totalDamage
            damage[iID] = totalDamage / undamaged_volume
        end
    end
    return nothing
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
function init_interface_crit_values(datamanager::Module,
                                    damage_parameter::Dict,
                                    block_id::Int64)
    if !haskey(damage_parameter, "Interblock Damage")
        return datamanager
    end
    max_block_id = maximum(datamanager.get_block_id_list())
    inter_critical_value = datamanager.get_crit_values_matrix()
    if inter_critical_value == fill(-1, (1, 1, 1))
        inter_critical_value = fill(Float64(damage_parameter["Critical Value"]),
                                    (max_block_id, max_block_id, max_block_id))
    end
    for block_iId in 1:max_block_id
        for block_jId in 1:max_block_id
            critical_value_name = "Interblock Critical Value $(block_iId)_$block_jId"
            if haskey(damage_parameter["Interblock Damage"], critical_value_name)
                if damage_parameter["Interblock Damage"][critical_value_name] isa Number
                    inter_critical_value[block_iId, block_jId,
                                         block_id] = damage_parameter["Interblock Damage"][critical_value_name]
                end
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
function init_aniso_crit_values(datamanager::Module,
                                damage_parameter::Dict,
                                block_id::Int64)
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

    # if !haskey(damage_parameter, "Anisotropic Damage")
    #     return datamanager
    # end
    # aniso_crit[block_id] =  damage_parameter["Anisotropic Damage"]
    # datamanager.set_aniso_crit_values(aniso_crit)

    return datamanager
end

"""
    init_model(datamanager::Module, nodes::AbstractVector{Int64}, block::Int64)

Initialize the damage models.

# Arguments
- `datamanager::Module`: The data manager module where the degradation model will be initialized.
- `nodes::AbstractVector{Int64}`: Nodes for the degradation model.
- `block::Int64`: Block identifier for the degradation model.

# Returns
- `datamanager`: The modified data manager module with the initialized degradation model.

# Example
```julia
datamanager = init_model(my_data_manager, [1, 2, 3], 1)

"""
function init_model(datamanager::Module, nodes::AbstractVector{Int64},
                    block::Int64)
    model_param = datamanager.get_properties(block, "Damage Model")
    # if haskey(model_param, "Anisotropic Damage")
    #     datamanager.create_bond_field("Bond Damage Anisotropic", Float64, datamanager.get_dof(), 1)
    # end
    mod = create_module_specifics(model_param["Damage Model"],
                                  module_list,
                                  "damage_name")

    if isnothing(mod)
        @error "No damage model of name " * model_param["Damage Model"] * " exists."
        return nothing
    end
    datamanager.set_model_module(model_param["Damage Model"], mod)
    datamanager = mod.init_model(datamanager, nodes, model_param, block)
    datamanager = mod.fields_for_local_synchronization(datamanager, "Damage Model")
    datamanager = Damage.init_interface_crit_values(datamanager, model_param, block)
    datamanager = Damage.init_aniso_crit_values(datamanager, model_param, block)
    return datamanager
end
end
