# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Damage

using TimerOutputs: @timeit
using ....Data_Manager
using ...Solver_Manager: find_module_files, create_module_specifics
global module_list = find_module_files(@__DIR__, "damage_name")
for mod in module_list
    include(mod["File"])
end

using LoopVectorization
using .....Helpers: find_inverse_bond_id
export fields_for_local_synchronization
export compute_model
export init_interface_crit_values
export init_model
export init_fields

"""
    init_fields()

Initialize damage model fields

# Arguments
- `params::Dict`: Parameters.
"""
function init_fields()
    dof = Data_Manager.get_dof()
    Data_Manager.create_node_scalar_field("Damage", Float64)

    anisotropic_damage = false

    nlist = Data_Manager.get_nlist()
    inverse_nlist = Data_Manager.set_inverse_nlist(find_inverse_bond_id(nlist))
end

"""
    compute_model(nodes::AbstractVector{Int64}, model_param::Dict, block::Int64, time::Float64, dt::Float64)

Computes the damage model

# Arguments
- `nodes::AbstractVector{Int64}`: The nodes
- `model_param::Dict`: The model parameters
- `block::Int64`: The block
- `time::Float64`: The current time
- `dt::Float64`: The time step
"""
function compute_model(nodes::AbstractVector{Int64},
                       model_param::Dict,
                       block::Int64,
                       time::Float64,
                       dt::Float64)
    mod = Data_Manager.get_model_module(model_param["Damage Model"])

    mod.compute_model(nodes, model_param,
                      block, time, dt)

    if isnothing(Data_Manager.get_filtered_nlist())
        @timeit "compute index" return damage_index(nodes)
    end

    @timeit "compute index" return damage_index(nodes,
                                                Data_Manager.get_filtered_nlist())
end

"""
    fields_for_local_synchronization(model, block)

Defines all synchronization fields for local synchronization

# Arguments
- `model::String`: Model class.
- `block::Int64`: block ID
"""
function fields_for_local_synchronization(model, block)
    model_param = Data_Manager.get_properties(block, "Damage Model")
    mod = Data_Manager.get_model_module(model_param["Damage Model"])

    mod.fields_for_local_synchronization(model)
end

"""
    damage_index(::Union{SubArray, Vector{Int64})

Function calculates the damage index related to the neighborhood volume for a set of corresponding nodes.
The damage index is defined as damaged volume in relation the neighborhood volume.
damageIndex = sum_i (brokenBonds_i * volume_i) / volumeNeighborhood

# Arguments
- `nodes::AbstractVector{Int64}`: corresponding nodes to this model
"""
function damage_index(nodes::AbstractVector{Int64},
                      nlist_filtered_ids::BondScalarState{Int64})
    bond_damageNP1 = Data_Manager.get_bond_damage("NP1")
    for iID in nodes
        bond_damageNP1[iID][nlist_filtered_ids[iID]] .= 1
    end
    return damage_index(nodes)
end

function damage_index(nodes::AbstractVector{Int64})
    nlist = Data_Manager.get_nlist()::BondScalarState{Int64}
    volume = Data_Manager.get_field("Volume")::NodeScalarField{Float64}
    bond_damageNP1 = Data_Manager.get_bond_damage("NP1")::BondScalarState{Float64}
    damage = Data_Manager.get_damage("NP1")::NodeScalarField{Float64}
    compute_index(damage, nodes, volume, nlist, bond_damageNP1)
end

function compute_index(damage::NodeScalarField{Float64},
                       nodes::AbstractVector{Int64},
                       volume::NodeScalarField{Float64},
                       nlist::BondScalarState{Int64},
                       bond_damage::BondScalarState{Float64})::Nothing
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
    init_interface_crit_values(params::Dict, block_id::Int64)

Initialize the critical values

# Arguments
- `params::Dict`: The parameters
- `block_id::Int64`: current block
"""
function init_interface_crit_values(damage_parameter::Dict,
                                    block_id::Int64)
    if !haskey(damage_parameter, "Interblock Damage")
        return
    end
    max_block_id = maximum(Data_Manager.get_block_id_list())
    inter_critical_value = Data_Manager.get_crit_values_matrix()
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
    Data_Manager.set_crit_values_matrix(inter_critical_value)
end

"""
    init_aniso_crit_values(params::Dict, block_id::Int64)

Initialize the anisotropic critical values

# Arguments
- `params::Dict`: The parameters
- `block_id::Int64`: current block
"""
function init_aniso_crit_values(damage_parameter::Dict,
                                block_id::Int64, dof::Int64)
    aniso_crit::Dict{Int64,Any} = Dict()

    if !haskey(damage_parameter, "Anisotropic Damage")
        @warn "has no anisotropic damage"
        return
    end

    crit_x = damage_parameter["Anisotropic Damage"]["Critical Value X"]
    crit_y = damage_parameter["Anisotropic Damage"]["Critical Value Y"]
    if dof == 2
        aniso_crit[block_id] = [crit_x, crit_y]
    else
        crit_z = get(damage_parameter["Anisotropic Damage"], "Critical Value Z", crit_y)
        aniso_crit[block_id] = [crit_x, crit_y, crit_z]
    end
    Data_Manager.set_aniso_crit_values(aniso_crit)

    # if !haskey(damage_parameter, "Anisotropic Damage")
    #     return
    # end
    # aniso_crit[block_id] =  damage_parameter["Anisotropic Damage"]
    # Data_Manager.set_aniso_crit_values(aniso_crit)

    return
end

"""
    init_model(nodes::AbstractVector{Int64}, block::Int64)

Initialize the damage models.

# Arguments
- `nodes::AbstractVector{Int64}`: Nodes for the degradation model.
- `block::Int64`: Block identifier for the degradation model.

# Example
```julia
init_model(my_data_manager, [1, 2, 3], 1)
```
"""
function init_model(nodes::AbstractVector{Int64},
                    block::Int64)
    model_param = Data_Manager.get_properties(block, "Damage Model")
    # if haskey(model_param, "Anisotropic Damage")
    #     Data_Manager.create_bond_vector_state("Bond Damage Anisotropic", Float64, Data_Manager.get_dof(), 1)
    # end
    mod = create_module_specifics(model_param["Damage Model"],
                                  module_list,
                                  @__MODULE__,
                                  "damage_name")

    if isnothing(mod)
        @error "No damage model of name " * model_param["Damage Model"] * " exists."
        return nothing
    end
    Data_Manager.set_model_module(model_param["Damage Model"], mod)
    mod.init_model(nodes, model_param, block)
    mod.fields_for_local_synchronization("Damage Model")
    Damage.init_interface_crit_values(model_param, block)
    Damage.init_aniso_crit_values(model_param, block, Data_Manager.get_dof())
end
end
