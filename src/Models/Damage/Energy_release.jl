# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Critical_Energy
using ......Helpers:
                     rotate,
                     fastdot,
                     sub_in_place!,
                     div_in_place!,
                     mul_in_place!,
                     interpol_data,
                     is_dependent
using LinearAlgebra: mul!, dot
using StaticArrays
export compute_model
export damage_name
export init_model
export fields_for_local_synchronization
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
    compute_model(datamanager, nodes, damage_parameter, block, time, dt)

Calculates the elastic energy of each bond and compares it to a critical one. If it is exceeded, the bond damage value is set to zero.
[WillbergC2019](@cite), [FosterJT2011](@cite)

# Arguments
- `datamanager::Data_Manager`: Datamanager.
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `damage_parameter::Dict(String, Any)`: Dictionary with material parameter.
- `block::Int64`: Block number.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
# Returns
- `datamanager::Data_Manager`: Datamanager.
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
    dof = datamanager.get_dof()
    nlist = datamanager.get_nlist()
    block_ids = datamanager.get_field("Block_Id")
    update_list = datamanager.get_field("Update")
    bond_damage = datamanager.get_bond_damage("NP1")

    undeformed_bond = datamanager.get_field("Bond Geometry")
    undeformed_bond_length = datamanager.get_field("Bond Length")
    bond_forces = datamanager.get_field("Bond Forces")
    deformed_bond = datamanager.get_field("Deformed Bond Geometry", "NP1")
    deformed_bond_length = datamanager.get_field("Deformed Bond Length", "NP1")
    bond_displacements = datamanager.get_field("Bond Displacements")
    critical_field = datamanager.has_key("Critical_Value")
    critical_energy = critical_field ? datamanager.get_field("Critical_Value") :
                      damage_parameter["Critical Value"]
    critical_energy_value = 0.0
    quad_horizons = datamanager.get_field("Quad Horizon")
    inverse_nlist = datamanager.get_inverse_nlist()

    dependend_value,
    dependent_field = is_dependent("Critical Value", damage_parameter,
                                   datamanager)

    tension::Bool = get(damage_parameter, "Only Tension", false)
    inter_block_damage::Bool = haskey(damage_parameter, "Interblock Damage")
    if inter_block_damage
        inter_critical_energy::Array{Float64,3} = datamanager.get_crit_values_matrix()
    end

    bond_energy::Float64 = 0.0
    norm_displacement::Float64 = 0.0
    product::Float64 = 0.0

    temp_vector::Vector{Float64} = zeros(Float64, dof)

    sub_in_place!(bond_displacements, deformed_bond, undeformed_bond)
    warning_flag = true

    for iID in nodes
        block_id::Int64 = block_ids[iID]
        quad_horizon::Float64 = quad_horizons[iID]
        neighbors::Vector{Int64} = nlist[iID]
        @fastmath @inbounds for jID in eachindex(neighbors)
            relative_displacement::Vector{Float64} = bond_displacements[iID][jID]
            norm_displacement = dot(relative_displacement, relative_displacement)
            if norm_displacement == 0 || (tension &&
                deformed_bond_length[iID][jID] - undeformed_bond_length[iID][jID] < 0)
                continue
            end

            neighborID::Int64 = neighbors[jID]
            neighbor_block_id::Int64 = block_ids[neighborID]

            # check if the bond also exist at other node, due to different horizons
            # try
            #     neighbor_bond_force .= bond_forces[neighborID][inverse_nlist[neighborID][iID]]
            # catch e
            #     # Handle the case when the key doesn't exist
            # end

            bond_force::Vector{Float64} = bond_forces[iID][jID]
            neighbor_bond_force::Vector{Float64} = bond_forces[neighborID][inverse_nlist[neighborID][iID]]
            temp_vector .= bond_force .- neighbor_bond_force

            product = dot(temp_vector, relative_displacement)
            mul!(temp_vector, product / norm_displacement, relative_displacement)
            product = fastdot(temp_vector, relative_displacement, true)
            bond_energy = 0.25 * product

            if critical_field
                critical_energy_value = critical_energy[iID]
            elseif inter_block_damage
                critical_energy_value = inter_critical_energy[block_id, neighbor_block_id, block]

                # param_name = "Interblock Critical Value " * string(block_ids[iID]) * "_" *
                #              string(block_ids[neighborID])

                # dependend_value,
                # dependent_field = is_dependent(param_name, damage_parameter, datamanager)
                # if dependend_value
                #     critical_energy_value = interpol_data(dependent_field[iID],
                #                                           damage_parameter[param_name]["Data"],
                #                                           warning_flag)
                # end
            elseif dependend_value
                critical_energy_value = interpol_data(dependent_field[iID],
                                                      damage_parameter["Critical Value"]["Data"],
                                                      warning_flag)
            else
                critical_energy_value = critical_energy
            end

            product = critical_energy_value * quad_horizon
            if bond_energy > product
                bond_damage[iID][jID] = 0.0
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
    download_from_cores = false
    upload_to_cores = true
    datamanager.set_local_synch(model, "Bond Forces", download_from_cores, upload_to_cores)
    return datamanager
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
    #TODO: Use average horizon
    if dof == 2
        return Float64(3 / (pi * horizon^3 * thickness))
    end
    return Float64(4 / (pi * horizon^4))
end

function init_model(datamanager::Module,
                    nodes::AbstractVector{Int64},
                    damage_parameter::Dict,
                    block::Int64)
    dof = datamanager.get_dof()
    quad_horizon = datamanager.create_constant_node_field("Quad Horizon", Float64, 1)
    datamanager.create_constant_bond_field("Bond Displacements", Float64, dof)
    horizon = datamanager.get_field("Horizon")
    thickness::Float64 = get(damage_parameter, "Thickness", 1)
    for iID in nodes
        quad_horizon[iID] = get_quad_horizon(horizon[iID], dof, thickness)
    end

    return datamanager
end
end
