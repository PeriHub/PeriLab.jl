# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Critical_Stretch
using LinearAlgebra: norm
using .....Data_Manager
using .....Geometry: compute_stretch!
using .....Helpers: with_dof
using ..Damage: critical_value_lookup, COMPRESSION_MODES, transverse
export compute_model
export init_model
export fields_for_local_synchronization

using ......ParameterSpec: @params, register_damage
@params "Breaks a bond once its stretch exceeds the critical value" struct CriticalStretchParams
    "How a shortened bond counts: full (its absolute stretch), none (skipped) or shear (its transverse displacement over its length)"
    compression::String = opt("Compression"; default = "none", allowed = COMPRESSION_MODES)
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
    crit_of = critical_value_lookup(damage.base.critical_value,
                                    damage.base.interblock_damage !== nothing, block;
                                    bonding_exponent = damage.base.interlayer_bonding ? 1 :
                                                       0)
    return _critical_stretch!(nodes, p.compression, crit_of)
end

# function barrier: `crit_of` has a concrete closure type here, the dof becomes `D`
function _critical_stretch!(nodes::AbstractVector{Int64}, compression::String,
                            crit_of::F) where {F}
    with_dof(D -> _critical_stretch!(nodes, compression, crit_of, D))
end

function _critical_stretch!(nodes::AbstractVector{Int64}, compression::String,
                            crit_of::F, ::Val{D}) where {F,D}
    nlist::BondScalarState{Int64} = Data_Manager.get_nlist()
    bond_damageNP1::BondScalarState{Float64} = Data_Manager.get_bond_damage("NP1")
    update_list::NodeScalarField{Bool} = Data_Manager.get_field("Update")
    undeformed_bond_length::BondScalarState{Float64} = Data_Manager.get_field("Bond Length")
    deformed_bond_length::BondScalarState{Float64} = Data_Manager.get_field("Deformed Bond Length",
                                                                            "NP1")
    stretch::BondScalarState{Float64} = Data_Manager.get_field("Bond Stretch")
    undeformed_bond = Data_Manager.get_field("Bond Geometry")::BondVectorState{Float64,D}
    deformed_bond = Data_Manager.get_field("Deformed Bond Geometry",
                                           "NP1")::BondVectorState{Float64,D}
    compute_stretch!(stretch, deformed_bond_length, undeformed_bond_length)
    # TBD all points; for MPI it must be only the nodes
    full = compression == "full"
    shear = compression == "shear"

    for iID in nodes
        @inbounds for (jID, neighborID) in enumerate(nlist[iID])
            s = stretch[iID][jID]
            stretch_check = if s >= 0 || full
                abs(s)
            elseif shear
                deformed = deformed_bond[iID][jID]
                n = deformed / deformed_bond_length[iID][jID]
                norm(transverse(deformed - undeformed_bond[iID][jID], n)) /
                undeformed_bond_length[iID][jID]
            else
                s
            end
            if stretch_check > crit_of(iID, jID, neighborID)
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
