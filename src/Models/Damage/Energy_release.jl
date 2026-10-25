# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Critical_Energy

using StaticArrays
using LinearAlgebra: mul!, dot
using Polyester: @batch

using ......Data_Manager
using ......Helpers: MINBATCH,
                     rotate,
                     fastdot,
                     div_in_place!,
                     mul_in_place!,
                     with_dof

export compute_model
export init_model
export fields_for_local_synchronization

using ......ParameterSpec: @params, register_damage
using ..Damage: critical_value_lookup, init_inverse_bond_id
@params "Breaks a bond once its energy exceeds the critical energy release rate" struct CriticalEnergyParams
    "Breaks bonds only under tension, not under compression"
    only_tension::Bool = opt("Only Tension"; default = true)
    "Out-of-plane thickness of a 2D model, used to convert the energy release rate to a bond energy"
    thickness::Float64 = opt("Thickness"; default = 1.0, min = 0, quantity = :length)
end
__init__() = register_damage("Critical Energy", CriticalEnergyParams)

"""
    compute_model(nodes, p, damage, block, time, dt)

Calculates the elastic energy of each bond and compares it to a critical one. If it is exceeded, the bond damage value is set to zero.
[WillbergC2019](@cite), [FosterJT2011](@cite)

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `p::CriticalEnergyParams`: The model parameters.
- `damage::BlockDamage`: The typed damage model of the block (`damage.base`: Critical Value, Interblock Damage).
- `block::Int64`: Block number.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
"""
function compute_model(nodes::AbstractVector{Int64}, p::CriticalEnergyParams, damage,
                       block::Int64, time::Float64, dt::Float64)
    crit_of = critical_value_lookup(damage.base.critical_value,
                                    damage.base.interblock_damage !== nothing, block)
    return _critical_energy!(nodes, p.only_tension, crit_of)
end

# function barrier: `crit_of` has a concrete closure type here, the dof becomes `D`
function _critical_energy!(nodes::AbstractVector{Int64}, tension::Bool,
                           crit_of::F) where {F}
    with_dof(D -> _critical_energy!(nodes, tension, crit_of, D))
end

function _critical_energy!(nodes::AbstractVector{Int64}, tension::Bool, crit_of::F,
                           ::Val{D}) where {F,D}
    nlist = Data_Manager.get_nlist()::BondScalarState{Int64}
    update_list = Data_Manager.get_field("Update")::NodeScalarField{Bool}
    bond_damage = Data_Manager.get_bond_damage("NP1")::BondScalarState{Float64}
    bond_energy = Data_Manager.get_field("Bond Energy")::BondScalarState{Float64}
    undeformed_bond = Data_Manager.get_field("Bond Geometry")::BondVectorState{Float64,D}
    undeformed_bond_length = Data_Manager.get_field("Bond Length")::BondScalarState{Float64}
    bond_forces = Data_Manager.get_field("Bond Forces")::BondVectorState{Float64,D}
    deformed_bond = Data_Manager.get_field("Deformed Bond Geometry",
                                           "NP1")::BondVectorState{Float64,D}
    deformed_bond_length = Data_Manager.get_field("Deformed Bond Length",
                                                  "NP1")::BondScalarState{Float64}
    quad_horizons = Data_Manager.get_field("Quad Horizon")::NodeScalarField{Float64}
    inverse_bond_id = Data_Manager.get_field("Inverse Bond ID")::BondScalarState{Int64}

    @batch minbatch=MINBATCH for iID in nodes
        # the bond views of the node, built once instead of per bond
        neighbors = nlist[iID]
        deformed = deformed_bond[iID]
        undeformed = undeformed_bond[iID]
        deformed_length = deformed_bond_length[iID]
        undeformed_length = undeformed_bond_length[iID]
        inverse_ids = inverse_bond_id[iID]
        forces = bond_forces[iID]
        energy = bond_energy[iID]
        damage = bond_damage[iID]
        quad_horizon = quad_horizons[iID]
        @inbounds @fastmath for jID in eachindex(neighbors)
            relative_displacement = deformed[jID] - undeformed[jID]
            norm_displacement = dot(relative_displacement, relative_displacement)
            if norm_displacement == 0 ||
               (tension && deformed_length[jID] - undeformed_length[jID] < 0)
                continue
            end
            # the bond may not exist at the neighbor, due to different horizons
            inverse_neighborID = inverse_ids[jID]
            inverse_neighborID == 0 && continue
            neighborID = neighbors[jID]

            force_difference = forces[jID] - bond_forces[neighborID][inverse_neighborID]
            product = abs(dot(force_difference, relative_displacement))
            relative_displacement = abs.(relative_displacement)
            projected = (product / norm_displacement) * relative_displacement
            bond_energy_j = 0.25 * dot(projected, relative_displacement)
            energy[jID] = bond_energy_j
            # TODO: broken bonds are flagged again every step (their zero energy forces keep
            # the energy above the limit), so their nodes are recomputed in every update pass.
            # A `damage[jID] > 0 &&` guard saves ~10% for correspondence, but changes bond
            # associated results: there a node also depends on the weighted volume of its
            # neighbors, whose changes are not flagged.
            if bond_energy_j > crit_of(iID, neighborID) * quad_horizon
                damage[jID] = 0.0
                update_list[iID] = true
            end
        end
    end
end

"""
    fields_for_local_synchronization(model::String)

Returns a user developer defined local synchronization. This happens before each model.

# Arguments
- `model::String`: The model to be synchronized.
"""
function fields_for_local_synchronization(model::String)
    download_from_cores = false
    upload_to_cores = true
    Data_Manager.set_local_synch(model, "Bond Forces", download_from_cores, upload_to_cores)
end

"""
    get_quad_horizon(horizon::Float64, dof::Int64)

Get the quadric of the horizon.

# Arguments
- `horizon::Float64`: The horizon of the block.
- `dof::Int64`: The degree of freedom.
- `thickness::Float64`: The thickness of the block.
# Returns
- `quad_horizon::Float64`: The quadric of the horizon.
"""
function get_quad_horizon(horizon::Float64, dof::Int64, thickness::Float64,
                          mesh_scaling_sum::Float64)
    avg_horizon = horizon / 3 * mesh_scaling_sum
    if dof == 2
        return Float64(3 / (pi * avg_horizon^3 * thickness))
    end
    #TODO: Use average horizon
    return Float64(4 / (pi * avg_horizon^4))
end

"""
    init_model(nodes, p, damage, block)

Creates the fields of the model and the quadric horizons of the block nodes.
"""
function init_model(nodes::AbstractVector{Int64}, p::CriticalEnergyParams, damage,
                    block::Int64)
    dof = Data_Manager.get_dof()
    quad_horizons = Data_Manager.create_constant_node_scalar_field("Quad Horizon", Float64)
    horizon = Data_Manager.get_field("Horizon")
    thickness::Float64 = p.thickness
    mesh_scaling = Data_Manager.get_horizon_mesh_scaling()
    bond_energy = Data_Manager.create_constant_bond_scalar_state("Bond Energy",
                                                                 Float64)
    init_inverse_bond_id()
    for iID in nodes
        quad_horizons[iID] = get_quad_horizon(horizon[iID], dof, thickness,
                                              sum(inv, mesh_scaling))
    end
end
end
