# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Critical_Stretch
using .....Data_Manager
using .....Helpers: sub_in_place!, div_in_place!

export compute_model
export damage_name
export init_model
export fields_for_local_synchronization

"""
    damage_name()

Gives the damage name. It is needed for comparison with the yaml input deck.

# Arguments

# Returns
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
    compute_model(nodes, damage_parameter, block, time, dt)

Calculates the stretch of each bond and compares it to a critical one. If it is exceeded, the bond damage value is set to zero.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `damage_parameter::Dict(String, Any)`: Dictionary with material parameter.
- `block::Int64`: Block number.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
Example:
```julia
```
"""
function compute_model(nodes::AbstractVector{Int64},
                       damage_parameter::Dict,
                       block::Int64,
                       time::Float64,
                       dt::Float64)
    nlist::BondScalarState{Int64} = Data_Manager.get_nlist()
    bond_damageNP1::BondScalarState{Float64} = Data_Manager.get_bond_damage("NP1")
    update_list::NodeScalarField{Bool} = Data_Manager.get_field("Update")
    undeformed_bond_length::BondScalarState{Float64} = Data_Manager.get_field("Bond Length")
    deformed_bond_length::BondScalarState{Float64} = Data_Manager.get_field("Deformed Bond Length",
                                                                            "NP1")
    block_ids::NodeScalarField{Int64} = Data_Manager.get_field("Block_Id")
    temp::BondScalarState{Float64} = Data_Manager.get_field("Temporary Bond Field")
    critical_field = Data_Manager.has_key("Critical_Value")
    if critical_field
        critical_stretch = Data_Manager.get_field("Critical_Value")
    else
        critical_stretch = damage_parameter["Critical Value"]
    end
    tension::Bool = get(damage_parameter, "Only Tension", true)
    inter_block_damage::Bool = haskey(damage_parameter, "Interblock Damage")
    if inter_block_damage
        inter_critical_stretch::Array{Float64,3} = Data_Manager.get_crit_values_matrix()
    end

    sub_in_place!(temp, deformed_bond_length, undeformed_bond_length)
    div_in_place!(temp, temp, undeformed_bond_length)

    stretch::Float64 = 0.0
    crit_stretch::Float64 = 0.0

    for iID in nodes
        @fastmath @inbounds @simd for jID in eachindex(nlist[iID])
            # stretch = (deformed_bond_length[iID][jID] - undeformed_bond_length[iID][jID]) /
            #           undeformed_bond_length[iID][jID]

            if critical_field
                crit_stretch = critical_stretch[iID]
            else
                crit_stretch = inter_block_damage ?
                               inter_critical_stretch[block_ids[iID],
                block_ids[nlist[iID][jID]],
                block] : critical_stretch
            end

            stretch = tension ? temp[iID][jID] : abs(temp[iID][jID])
            if stretch > crit_stretch
                bond_damageNP1[iID][jID] = 0.0
                update_list[iID] = true
            end
        end
    end
end

"""
    fields_for_local_synchronization(model::String)

Returns a user developer defined local synchronization. This happens before each model.



# Arguments

"""
function fields_for_local_synchronization(model::String)
end

function init_model(odes::AbstractVector{Int64},
                    damage_parameter::Dict,
                    block::Int64)
end
end
