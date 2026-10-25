# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Critical_Energy_Aniso
using LinearAlgebra: dot
using StaticArrays

using ......Data_Manager
using ......PeriLabExceptions: @abort
using ......Helpers: node_tensor, with_dof
using ..Damage: init_inverse_bond_id

export compute_model
export init_model
export fields_for_local_synchronization

using ......ParameterSpec: @params, register_damage
@params "Critical energy release rate with different critical values along the material axes (Anisotropic Damage)" struct CriticalEnergyAnisotropicParams
    "Breaks bonds only under tension, not under compression"
    only_tension::Bool = opt("Only Tension"; default = false)
    "Out-of-plane thickness of a 2D model, used to convert the energy release rate to a bond energy"
    thickness::Float64 = opt("Thickness"; default = 1.0, min = 0, quantity = :length)
end
__init__() = register_damage("Critical Energy Anisotropic", CriticalEnergyAnisotropicParams)
"""
    compute_model(nodes, p, damage, block, time, dt)

Calculates the elastic energy of each bond and compares it to a critical one. If it is exceeded, the bond damage value is set to zero.
[WillbergC2019](@cite), [FosterJT2011](@cite)

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `p::CriticalEnergyAnisotropicParams`: The model parameters.
- `damage::BlockDamage`: The typed damage model of the block (`damage.base`: Critical Value, Interblock Damage, Anisotropic Damage).
- `block::Int64`: Block number.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
"""
function compute_model(nodes::AbstractVector{Int64}, p::CriticalEnergyAnisotropicParams,
                       damage, block::Int64, time::Float64, dt::Float64)
    return _critical_energy_aniso!(nodes, p.only_tension)
end

function _critical_energy_aniso!(nodes::AbstractVector{Int64}, tension::Bool)
    # function barrier on the dof: the fields below get concrete types
    with_dof(D -> _critical_energy_aniso!(nodes, tension, D))
end

function _critical_energy_aniso!(nodes::AbstractVector{Int64}, tension::Bool,
                                 ::Val{D}) where {D}
    nlist = Data_Manager.get_nlist()::BondScalarState{Int64}
    block_ids = Data_Manager.get_field("Block_Id")::NodeScalarField{Int64}
    update_list = Data_Manager.get_field("Update")::NodeScalarField{Bool}
    bond_damage = Data_Manager.get_bond_damage("NP1")::BondScalarState{Float64}
    undeformed_bond = Data_Manager.get_field("Bond Geometry")::BondVectorState{Float64,D}
    undeformed_bond_length = Data_Manager.get_field("Bond Length")::BondScalarState{Float64}
    bond_forces = Data_Manager.get_field("Bond Forces")::BondVectorState{Float64,D}
    deformed_bond = Data_Manager.get_field("Deformed Bond Geometry",
                                           "NP1")::BondVectorState{Float64,D}
    deformed_bond_length = Data_Manager.get_field("Deformed Bond Length",
                                                  "NP1")::BondScalarState{Float64}
    quad_horizons = Data_Manager.get_field("Quad Horizon")::NodeScalarField{Float64}
    inverse_bond_id = Data_Manager.get_field("Inverse Bond ID")::BondScalarState{Int64}
    rotation_tensor = Data_Manager.get_field("Rotation Tensor")::NodeTensorField{Float64}
    aniso_crit_values = Data_Manager.get_aniso_crit_values()::Dict{Int64,Any}
    rotation = Data_Manager.get_rotation()

    @inbounds for iID in nodes
        # the bond views of the node, built once instead of per bond
        neighbors = nlist[iID]
        deformed = deformed_bond[iID]
        undeformed = undeformed_bond[iID]
        deformed_length = deformed_bond_length[iID]
        undeformed_length = undeformed_bond_length[iID]
        inverse_ids = inverse_bond_id[iID]
        forces = bond_forces[iID]
        damage = bond_damage[iID]
        R = node_tensor(rotation_tensor, iID, Val(D))
        critical = SVector{D,Float64}(aniso_crit_values[block_ids[iID]]::Vector{Float64}) *
                   quad_horizons[iID]
        @fastmath for jID in eachindex(neighbors)
            relative_displacement = deformed[jID] - undeformed[jID]
            norm_displacement = dot(relative_displacement, relative_displacement)
            if norm_displacement == 0 ||
               (tension && deformed_length[jID] - undeformed_length[jID] < 0)
                continue
            end
            # the bond may not exist at the neighbor, due to different horizons
            inverse_neighborID = inverse_ids[jID]
            inverse_neighborID == 0 && continue

            force_difference = forces[jID] - bond_forces[neighbors[jID]][inverse_neighborID]
            product = dot(force_difference, relative_displacement)
            projected_force = (product / norm_displacement) * relative_displacement
            bond_energy = 0.25 * sum(abs.(projected_force) .* abs.(relative_displacement))

            rotated_bond = rotation ? R * deformed[jID] : deformed[jID]
            bond_norm = abs.(rotated_bond) / deformed_length[jID]
            exceeded = bond_energy * bond_norm .> critical
            damage[jID] = max(damage[jID] - sum(bond_norm .* exceeded), 0.0)
            any(exceeded) && (update_list[iID] = true)
        end
    end
end

"""
    fields_for_local_synchronization(model::String)

Returns a user developer defined local synchronization. This happens before each model.



# Arguments

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
function get_quad_horizon(horizon::Float64, dof::Int64, thickness::Float64)
    if dof == 2
        return Float64(3 / (pi * horizon^3 * thickness))
    end
    return Float64(4 / (pi * horizon^4))
end

"""
    init_model(nodes, p, damage, block)

Creates the fields of the model and the quadric horizons of the block nodes; needs
`Anisotropic Damage` and the Angles field.
"""
function init_model(nodes::AbstractVector{Int64}, p::CriticalEnergyAnisotropicParams,
                    damage, block::Int64)
    if damage.base.anisotropic_damage === nothing
        @abort "Critical Energy Anisotropic requires Anisotropic Damage."
    end
    dof = Data_Manager.get_dof()
    quad_horizon = Data_Manager.create_constant_node_scalar_field("Quad Horizon", Float64)
    init_inverse_bond_id()
    horizon = Data_Manager.get_field("Horizon")
    thickness::Float64 = p.thickness
    for iID in nodes
        quad_horizon[iID] = get_quad_horizon(horizon[iID], dof, thickness)
    end

    if !Data_Manager.get_rotation()
        @abort "Anisotropic damage requires Angles field"
    end
end
end
