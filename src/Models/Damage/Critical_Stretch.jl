# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Critical_Stretch
using .....Data_Manager
using .....Geometry: compute_stretch!
export compute_model
export init_model
export fields_for_local_synchronization

using ......ParameterSpec: @params, register_damage, value
@params struct CriticalStretchParams
    only_tension::Bool = opt("Only Tension"; default = true)
end
__init__() = register_damage("Critical Stretch", CriticalStretchParams)

"""
    compute_model(nodes, p, damage, block, time, dt)

Calculates the stretch of each bond and compares it to a critical one. If it is exceeded, the bond damage value is set to zero.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `p::CriticalStretchParams`: The model parameters.
- `damage::BlockDamage`: The typed damage model of the block (`damage.base`: Critical Value, Interblock Damage).
- `block::Int64`: Block number.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
"""
function compute_model(nodes::AbstractVector{Int64}, p::CriticalStretchParams, damage,
                       block::Int64, time::Float64, dt::Float64)
    return _critical_stretch!(nodes, damage.base.critical_value, p.only_tension,
                              damage.base.interblock_damage !== nothing, block)
end

# function barrier: `critical_value` has a concrete type (Constant or Table1D) here
function _critical_stretch!(nodes::AbstractVector{Int64}, critical_value, tension::Bool,
                            inter_block_damage::Bool, block::Int64)
    nlist::BondScalarState{Int64} = Data_Manager.get_nlist()
    bond_damageNP1::BondScalarState{Float64} = Data_Manager.get_bond_damage("NP1")
    update_list::NodeScalarField{Bool} = Data_Manager.get_field("Update")
    undeformed_bond_length::BondScalarState{Float64} = Data_Manager.get_field("Bond Length")
    deformed_bond_length::BondScalarState{Float64} = Data_Manager.get_field("Deformed Bond Length",
                                                                            "NP1")
    stretch::BondScalarState{Float64} = Data_Manager.get_field("Bond Stretch")
    block_ids::NodeScalarField{Int64} = Data_Manager.get_field("Block_Id")

    critical_field = Data_Manager.has_key("Critical_Value")
    if critical_field
        critical_stretch = Data_Manager.get_field("Critical_Value")
    end
    if inter_block_damage
        inter_critical_stretch::Array{Float64,3} = Data_Manager.get_crit_values_matrix()
    end
    compute_stretch!(stretch, deformed_bond_length, undeformed_bond_length)
    # TBD all points; for MPI it must be only the nodes

    stretch_check::Float64 = 0.0
    crit_stretch::Float64 = 0.0

    for iID in nodes
        @fastmath @inbounds @simd for jID in eachindex(nlist[iID])
            if critical_field
                crit_stretch = critical_stretch[iID]
            else
                crit_stretch = inter_block_damage ?
                               inter_critical_stretch[block_ids[iID],
                                                      block_ids[nlist[iID][jID]],
                                                      block] : value(critical_value, iID)
            end

            stretch_check = tension ? stretch[iID][jID] : abs(stretch[iID][jID])
            if stretch_check > crit_stretch
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

"""
    init_model(nodes, p, damage, block)

Creates the bond stretch field.
"""
function init_model(nodes::AbstractVector{Int64}, p::CriticalStretchParams, damage,
                    block::Int64)
    Data_Manager.create_constant_bond_scalar_state("Bond Stretch", Float64)
end

end
