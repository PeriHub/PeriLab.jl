# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Critical_Energy

using StaticArrays
using LinearAlgebra: mul!, dot

using ......Data_Manager
using ......Helpers:
                     abs!,
                     rotate,
                     fastdot,
                     sub_in_place!,
                     div_in_place!,
                     mul_in_place!

export compute_model
export damage_name
export init_model
export fields_for_local_synchronization

using ......ParameterSpec: @params, register_damage, value
@params struct CriticalEnergyParams
    only_tension::Bool = opt("Only Tension"; default = true)
    thickness::Float64 = opt("Thickness"; default = 1.0, min = 0, quantity = :length)
end
__init__() = register_damage("Critical Energy", CriticalEnergyParams)

"""
    damage_name()

Gives the damage name. It is needed for comparison with the yaml input deck.

# Returns
- `name::String`: The name of the damage.

Example:
```julia
println(damage_name())
"Critical Energy"
```
"""
function damage_name()
    return "Critical Energy"
end

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
    return _critical_energy!(nodes, damage.base.critical_value, p.only_tension,
                             damage.base.interblock_damage !== nothing, block)
end

# function barrier: `critical_value` has a concrete type (Constant or Table1D) here
function _critical_energy!(nodes::AbstractVector{Int64}, critical_value, tension::Bool,
                           inter_block_damage::Bool, block::Int64)
    dof = Data_Manager.get_dof()
    nlist::BondScalarState{Int64} = Data_Manager.get_nlist()
    block_ids::NodeScalarField{Int64} = Data_Manager.get_field("Block_Id")
    update_list::NodeScalarField{Bool} = Data_Manager.get_field("Update")
    bond_damage::BondScalarState{Float64} = Data_Manager.get_bond_damage("NP1")

    bond_energy::BondScalarState{Float64} = Data_Manager.get_field("Bond Energy")
    undeformed_bond::BondVectorState{Float64} = Data_Manager.get_field("Bond Geometry")
    undeformed_bond_length::BondScalarState{Float64} = Data_Manager.get_field("Bond Length")
    bond_forces::BondVectorState{Float64} = Data_Manager.get_field("Bond Forces")
    deformed_bond::BondVectorState{Float64} = Data_Manager.get_field("Deformed Bond Geometry",
                                                                     "NP1")
    deformed_bond_length::BondScalarState{Float64} = Data_Manager.get_field("Deformed Bond Length",
                                                                            "NP1")
    bond_displacements::BondVectorState{Float64} = Data_Manager.get_field("Bond Displacements")
    critical_field = has_key("Critical_Value")
    if critical_field
        critical_energy = Data_Manager.get_field("Critical_Value")
    end
    critical_energy_value::Float64 = 0.0
    quad_horizons::NodeScalarField{Float64} = Data_Manager.get_field("Quad Horizon")
    inverse_nlist::Vector{Dict{Int64,Int64}} = Data_Manager.get_inverse_nlist()

    if inter_block_damage
        inter_critical_energy::Array{Float64,3} = Data_Manager.get_crit_values_matrix()
    end

    norm_displacement::Float64 = 0.0
    product::Float64 = 0.0

    temp_vector::Vector{Float64} = zeros(Float64, dof)

    sub_in_place!(bond_displacements, deformed_bond, undeformed_bond)

    for iID in nodes
        @fastmath @inbounds for jID in eachindex(nlist[iID])
            relative_displacement::Vector{Float64} = bond_displacements[iID][jID]
            norm_displacement = dot(relative_displacement, relative_displacement)
            if norm_displacement == 0 || (tension &&
                deformed_bond_length[iID][jID] - undeformed_bond_length[iID][jID] < 0)
                continue
            end

            neighborID::Int64 = nlist[iID][jID]
            if !haskey(inverse_nlist[neighborID], iID)
                continue
            end
            inverse_neighborID::Int64 = inverse_nlist[neighborID][iID]

            neighbor_block_id::Int64 = block_ids[neighborID]

            # check if the bond also exist at other node, due to different horizons
            # try
            #     neighbor_bond_force .= bond_forces[neighborID][inverse_nlist[neighborID][iID]]
            # catch e
            #     # Handle the case when the key doesn't exist
            # end

            bond_force::Vector{Float64} = bond_forces[iID][jID]
            neighbor_bond_force::Vector{Float64} = bond_forces[neighborID][inverse_neighborID]
            temp_vector .= bond_force .- neighbor_bond_force

            product = abs(dot(temp_vector, relative_displacement))
            abs!(relative_displacement)
            mul!(temp_vector, product / norm_displacement, relative_displacement)
            product = dot(temp_vector, relative_displacement)
            bond_energy[iID][jID] = 0.25 * product
            if critical_field
                critical_energy_value = critical_energy[iID]
            elseif inter_block_damage
                critical_energy_value = inter_critical_energy[block_ids[iID],
                neighbor_block_id, block]
            else
                critical_energy_value = value(critical_value, iID)
            end

            product = critical_energy_value * quad_horizons[iID]
            if bond_energy[iID][jID] > product
                bond_damage[iID][jID] = 0.0
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
    Data_Manager.create_constant_bond_vector_state("Bond Displacements", Float64, dof)
    horizon = Data_Manager.get_field("Horizon")
    thickness::Float64 = p.thickness
    mesh_scaling = Data_Manager.get_horizon_mesh_scaling()
    bond_energy = Data_Manager.create_constant_bond_scalar_state("Bond Energy",
                                                                 Float64)
    for iID in nodes
        quad_horizons[iID] = get_quad_horizon(horizon[iID], dof, thickness,
                                              sum(inv, mesh_scaling))
    end
end
end
