# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Damage
include("../../Core/Module_inclusion/set_Modules.jl")
using .Set_modules
global module_list = Set_modules.find_module_files(@__DIR__, "damage_name")
Set_modules.include_files(module_list)
include("../../Support/helpers.jl")
using TimerOutputs
using .Helpers: find_inverse_bond_id

export compute_model
export compute_damage_pre_calculation
export init_interface_crit_values
export init_model
export init_fields
export synch_field

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

    nlist = datamanager.get_field("Neighborhoodlist")
    inverse_nlist = datamanager.set_inverse_nlist(find_inverse_bond_id(nlist))
    return datamanager
end

"""
    compute_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, model_param::Dict, block::Int64, time::Float64, dt::Float64,to::TimerOutput,)

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
function compute_model(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    model_param::Dict,
    block::Int64,
    time::Float64,
    dt::Float64,
    to::TimerOutput,
)

    mod = datamanager.get_model_module(model_param["Damage Model"])
    datamanager = mod.compute_model(datamanager, nodes, model_param, block, time, dt)

    if isnothing(datamanager.get_filtered_nlist())
        return damage_index(datamanager, nodes)
    end

    return damage_index(datamanager, nodes, datamanager.get_filtered_nlist())
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
function compute_damage_pre_calculation(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    block::Int64,
    model_param::Dict,
    time::Float64,
    dt::Float64,
)
    mod = datamanager.get_model_module(model_param["Damage Model"])
    return mod.compute_damage_pre_calculation(datamanager, nodes, block, time, dt)
end

"""
    synch_field(datamanager::Module, damage_model::String, synchronise_field)

Field for synchronisation.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `damage_model::String`: The damage model
- `synchronise_field`: Synchronise function to distribute parameter through cores.
"""
function synch_field(datamanager::Module, damage_model::String, synchronise_field)
    mod = datamanager.get_model_module(damage_model)
    return mod.synch_field(datamanager, synchronise_field)
end

"""
    damage_index(datamanager,::Union{SubArray, Vector{Int64})

Function calculates the damage index related to the neighborhood volume for a set of corresponding nodes.
The damage index is defined as damaged volume in relation the neighborhood volume.
damageIndex = sum_i (brokenBonds_i * volume_i) / volumeNeighborhood

# Arguments
- `datamanager::Data_manager`: all model data
- `nodes::Union{SubArray, Vector{Int64}}`: corresponding nodes to this model
"""
function damage_index(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    nlist_filtered_ids::SubArray,
)
    nlist = datamanager.get_nlist()
    volume = datamanager.get_field("Volume")
    bond_damageNP1 = datamanager.get_bond_damage("NP1")
    damage = datamanager.get_damage("NP1")

    for iID in nodes
        undamaged_volume = sum(volume[nlist[iID]])
        bond_damageNP1[iID][nlist_filtered_ids[iID]] .= 1
        totalDamage = sum((1 .- bond_damageNP1[iID]) .* volume[nlist[iID]])
        if damage[iID] < totalDamage / undamaged_volume
            damage[iID] = totalDamage / undamaged_volume
        end
    end

    return datamanager

end

function damage_index(datamanager::Module, nodes::Union{SubArray,Vector{Int64}})
    nlist = datamanager.get_nlist()
    volume = datamanager.get_field("Volume")
    bond_damageNP1 = datamanager.get_bond_damage("NP1")
    damage = datamanager.get_damage("NP1")
    for iID in nodes
        undamaged_volume = sum(volume[nlist[iID]])
        totalDamage = sum((1 .- bond_damageNP1[iID]) .* volume[nlist[iID]])
        if damage[iID] < totalDamage / undamaged_volume
            damage[iID] = totalDamage / undamaged_volume
        end
    end
    return datamanager

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
    bond_damageN = datamanager.get_bond_damage("N")
    bond_damageNP1 = datamanager.get_bond_damage("NP1")
    aniso_bond_damageN = datamanager.get_field("Bond Damage Anisotropic", "N", false)
    aniso_bond_damageNP1 = datamanager.get_field("Bond Damage Anisotropic", "NP1", false)
    for iID in nodes
        bond_damageNP1[iID] .= bond_damageN[iID]
        if !isnothing(aniso_bond_damageN)
            aniso_bond_damageNP1[iID] .= aniso_bond_damageN[iID]
        end
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
function init_interface_crit_values(
    datamanager::Module,
    damage_parameter::Dict,
    block_id::Int64,
)
    if !haskey(damage_parameter, "Interblock Damage")
        return datamanager
    end
    max_block_id = maximum(datamanager.get_block_list())
    inter_critical_value = datamanager.get_crit_values_matrix()
    if inter_critical_value == fill(-1, (1, 1, 1))
        inter_critical_value = fill(
            Float64(damage_parameter["Critical Value"]),
            (max_block_id, max_block_id, max_block_id),
        )
    end
    for block_iId = 1:max_block_id
        for block_jId = 1:max_block_id
            critical_value_name = "Interblock Critical Value $(block_iId)_$block_jId"
            if haskey(damage_parameter["Interblock Damage"], critical_value_name)
                inter_critical_value[block_iId, block_jId, block_id] =
                    damage_parameter["Interblock Damage"][critical_value_name]
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
function init_aniso_crit_values(
    datamanager::Module,
    damage_parameter::Dict,
    block_id::Int64,
)
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
    init_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, block::Int64)

Initialize the damage models.

# Arguments
- `datamanager::Module`: The data manager module where the corrosion model will be initialized.
- `nodes::Union{SubArray,Vector{Int64}}`: Nodes for the corrosion model.
- `block::Int64`: Block identifier for the corrosion model.

# Returns
- `datamanager`: The modified data manager module with the initialized corrosion model.

# Example
```julia
datamanager = init_model(my_data_manager, [1, 2, 3], 1)

"""
function init_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, block::Int64)
    model_param = datamanager.get_properties(block, "Damage Model")
    if haskey(model_param, "Anisotropic Damage")
        datamanager.create_bond_field("Bond Damage Anisotropic", Float64, dof, 1)
    end
    mod = Set_modules.create_module_specifics(
        model_param["Damage Model"],
        module_list,
        "damage_name",
    )

    if isnothing(mod)
        @error "No damage model of name " * model_param["Damage Model"] * " exists."
        return nothing
    end
    datamanager.set_model_module(model_param["Damage Model"], mod)
    datamanager = mod.init_model(datamanager, nodes, model_param, block)
    datamanager.set_damage_models(model_param["Damage Model"])
    datamanager = Damage.init_interface_crit_values(datamanager, model_param, block)
    datamanager = Damage.init_aniso_crit_values(datamanager, model_param, block)
    return datamanager
end
end
