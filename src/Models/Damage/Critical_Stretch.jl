# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Critical_stretch
include("../../Support/Helpers.jl")
using .Helpers: sub_in_place!, div_in_place!
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
    compute_model(datamanager, nodes, damage_parameter, block, time, dt)

Calculates the stretch of each bond and compares it to a critical one. If it is exceeded, the bond damage value is set to zero.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `damage_parameter::Dict(String, Any)`: Dictionary with material parameter.
- `block::Int64`: Block number.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
# Returns
- `datamanager::Data_manager`: Datamanager.
Example:
```julia
```
"""
function compute_model(datamanager::Module,
                       nodes::AbstractVector{Int64},
                       damage_parameter::Dict,
                       block::Int64,
                       time::Float64,
                       dt::Float64)
    nlist = datamanager.get_nlist()
    bond_damageNP1 = datamanager.get_bond_damage("NP1")
    update_list = datamanager.get_field("Update")
    undeformed_bond_length = datamanager.get_field("Bond Length")
    deformed_bond_length = datamanager.get_field("Deformed Bond Length", "NP1")
    block_ids = datamanager.get_field("Block_Id")
    temp = datamanager.get_field("Temporary Bond Field")
    critical_field = datamanager.has_key("Critical_Value")
    if critical_field
        critical_stretch = datamanager.get_field("Critical_Value")
    else
        critical_stretch = damage_parameter["Critical Value"]
    end
    tension::Bool = get(damage_parameter, "Only Tension", true)
    inter_block_damage::Bool = haskey(damage_parameter, "Interblock Damage")
    if inter_block_damage
        inter_critical_stretch::Array{Float64,3} = datamanager.get_crit_values_matrix()
    end

    sub_in_place!(temp, deformed_bond_length, undeformed_bond_length)
    div_in_place!(temp, temp, undeformed_bond_length)

    if !critical_field
        if !any(any(x -> x > critical_stretch, tension ? vec : abs.(vec)) for vec in temp)
            # return if no stretch value is larger than critical_stretch
            return datamanager
        end
    end
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
    return datamanager
end

"""
    fields_for_local_synchronization(datamanager::Module, model::String)

Returns a user developer defined local synchronization. This happens before each model.



# Arguments

"""
function fields_for_local_synchronization(datamanager::Module, model::String)
    return datamanager
end

function init_model(datamanager::Module,
                    nodes::AbstractVector{Int64},
                    damage_parameter::Dict,
                    block::Int64)
    return datamanager
end
end
