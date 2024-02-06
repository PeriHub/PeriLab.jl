# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Critical_stretch
include("../Pre_calculation/Pre_Calculation_Factory.jl")
using .Pre_calculation
export compute_damage
export compute_damage_pre_calculation
export damage_name
export init_damage_model
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
    compute_damage(datamanager, nodes, damage_parameter, block, time, dt)

Calculates the stretch of each bond and compares it to a critical one. If it is exceeded, the bond damage value is set to zero.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
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
function compute_damage(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, damage_parameter::Dict, block::Int64, time::Float64, dt::Float64)

    nlist = datamanager.get_nlist()
    bond_damageNP1 = datamanager.get_bond_damage("NP1")
    update_list = datamanager.get_field("Update List")
    undeformed_bond = datamanager.get_field("Bond Geometry")
    deformed_bond = datamanager.get_field("Deformed Bond Geometry", "NP1")
    nneighbors = datamanager.get_field("Number of Neighbors")
    block_ids = datamanager.get_field("Block_Id")
    cricital_stretch = damage_parameter["Critical Value"]
    tension::Bool = get(damage_parameter, "Only Tension", true)
    inter_block_damage::Bool = haskey(damage_parameter, "Interblock Damage")
    if inter_block_damage
        inter_critical_stretch::Array{Float64,3} = datamanager.get_crit_values_matrix()
    end
    for iID in nodes
        for jID in nneighbors[iID]
            stretch = (deformed_bond[iID][jID, end] - undeformed_bond[iID][jID, end]) / undeformed_bond[iID][jID, end]

            crit_stretch = inter_block_damage ? inter_critical_stretch[block_ids[iID], block_ids[nlist[iID][jID]], block] : cricital_stretch

            if !tension
                stretch = abs(stretch)
            end
            if stretch > crit_stretch
                bond_damageNP1[iID][jID, :] .= 0.0
                update_list[iID] = true
            end
        end
    end
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
function compute_damage_pre_calculation(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, block::Int64, synchronise_field, time::Float64, dt::Float64)
    return datamanager
end

function init_damage_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, damage_parameter::Dict, block::Int64)
    return datamanager
end
end