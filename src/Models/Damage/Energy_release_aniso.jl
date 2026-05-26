# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Critical_Energy_Aniso
using TimerOutputs: @timeit
using LinearAlgebra: mul!, dot
using StaticArrays

using ......Data_Manager
using ......Helpers:
                     rotate,
                     fastdot,
                     sub_in_place!,
                     div_in_place!,
                     mul_in_place!,
                     interpol_data,
                     is_dependent

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
    return "Critical Energy Anisotropic"
end

"""
    compute_model(nodes, damage_parameter, block, time, dt)

Calculates the elastic energy of each bond and compares it to a critical one. If it is exceeded, the bond damage value is set to zero.
[WillbergC2019](@cite), [FosterJT2011](@cite)

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
    @timeit "init fields" begin
        dof::Int64 = Data_Manager.get_dof()
        nlist::BondScalarState{Int64} = Data_Manager.get_nlist()
        block_ids::NodeScalarField{Int64} = Data_Manager.get_field("Block_Id")
        update_list::NodeScalarField{Bool} = Data_Manager.get_field("Update")
        bond_damage::BondScalarState{Float64} = Data_Manager.get_bond_damage("NP1")
        undeformed_bond::BondVectorState{Float64} = Data_Manager.get_field("Bond Geometry")
        undeformed_bond_length::BondScalarState{Float64} = Data_Manager.get_field("Bond Length")
        bond_forces::BondVectorState{Float64} = Data_Manager.get_field("Bond Forces")
        deformed_bond::BondVectorState{Float64} = Data_Manager.get_field("Deformed Bond Geometry",
                                                                         "NP1")
        deformed_bond_length::BondScalarState{Float64} = Data_Manager.get_field("Deformed Bond Length",
                                                                                "NP1")
        bond_displacements::BondVectorState{Float64} = Data_Manager.get_field("Bond Displacements")
        critical_field::Bool = Data_Manager.has_key("Critical_Value")
        critical_energy = critical_field ? Data_Manager.get_field("Critical_Value") :
                          damage_parameter["Critical Value"]
        quad_horizons::NodeScalarField{Float64} = Data_Manager.get_field("Quad Horizon")
        inverse_nlist::Vector{Dict{Int64,Int64}} = Data_Manager.get_inverse_nlist()
        rotation_tensor::NodeTensorField{Float64} = Data_Manager.get_field("Rotation Tensor")
        aniso_crit_values::Dict{Int64,Any} = Data_Manager.get_aniso_crit_values()
    end

    @timeit "init params" begin
        dependend_value, dependent_field = is_dependent("Critical Value", damage_parameter)
        rotation = Data_Manager.get_rotation()
        tension = get(damage_parameter, "Only Tension", false)
        inter_block_damage = haskey(damage_parameter, "Interblock Damage")
        if inter_block_damage
            inter_critical_energy = Data_Manager.get_crit_values_matrix()
        end
        warning_flag = true
        critical_energy_value = 0.0
    end

    @timeit "preallocate" begin
        rotated_bond = zeros(Float64, dof)
        bond_norm_all = zeros(Float64, dof)
        condition = zeros(Bool, dof)
        temp = zeros(Float64, dof)
        temp2 = zeros(Float64, dof)
        deformed_buf = zeros(Float64, dof)
        relative_disp_buf = zeros(Float64, dof)
        bond_force_buf = zeros(Float64, dof)
        neighbor_force_buf = zeros(Float64, dof)
        bond_force_delta = zeros(Float64, dof)
        projected_force = zeros(Float64, dof)
        bond_energy = 0.0
        norm_displacement = 0.0
        product = 0.0
    end

    @timeit "sub_in_place" sub_in_place!(bond_displacements, deformed_bond, undeformed_bond)

    for iID in nodes
        bond_displacement = bond_displacements[iID]   # Vector — kein view
        bond_force_vec = bond_forces[iID]           # Vector — kein view
        quad_horizon = quad_horizons[iID]         # Vector — kein view
        rotation_temp = view(rotation_tensor,iID,:,:)  # 2D-Slice — view korrekt
        crit_vals = aniso_crit_values[block_ids[iID]]  # Vector — kein view

        @timeit "jID loop" begin
            @fastmath @inbounds for (jID, neighborID) in enumerate(nlist[iID])
                @timeit "relative disp" relative_disp_buf.=bond_displacement[jID]

                @timeit "norm disp" norm_displacement=dot(relative_disp_buf,
                                                          relative_disp_buf)

                @timeit "tension check" begin
                    if norm_displacement == 0 || (tension &&
                        deformed_bond_length[iID][jID] -
                        undeformed_bond_length[iID][jID] < 0)
                        continue
                    end
                end

                @timeit "check_keys" begin
                    if !haskey(inverse_nlist[neighborID], iID)
                        continue
                    end
                end

                @timeit "neighbor bond force" begin
                    inID = inverse_nlist[neighborID][iID]
                    neighbor_force_buf .= bond_forces[neighborID][inID]
                    bond_force_buf .= bond_force_vec[jID]
                    bond_force_delta .= bond_force_buf .- neighbor_force_buf
                    product = dot(bond_force_delta, relative_disp_buf)
                    mul!(projected_force, product / norm_displacement, relative_disp_buf)
                    product = fastdot(projected_force, relative_disp_buf, true)
                    bond_energy = 0.25 * product
                end

                @timeit "critical energy" begin
                    if critical_field
                        critical_energy_value = critical_energy[iID]
                    elseif inter_block_damage
                        critical_energy_value = inter_critical_energy[block_ids[iID],
                        block_ids[neighborID],
                        block]
                        param_name = "Interblock Critical Value " *
                                     string(block_ids[iID]) * "_" *
                                     string(block_ids[neighborID])
                        dependend_value,
                        dependent_field = is_dependent(param_name,
                                                       damage_parameter)
                        if dependend_value
                            critical_energy_value = interpol_data(dependent_field[iID],
                                                                  damage_parameter[param_name]["Data"],
                                                                  warning_flag)
                        end
                    elseif dependend_value
                        critical_energy_value = interpol_data(dependent_field[iID],
                                                              damage_parameter["Critical Value"]["Data"],
                                                              warning_flag)
                    else
                        critical_energy_value = critical_energy
                    end
                end

                @timeit "rotated bond" begin
                    deformed_buf .= deformed_bond[iID][jID]
                    if !rotation
                        rotated_bond .= deformed_buf
                    else
                        mul!(rotated_bond, rotation_temp, deformed_buf)
                    end
                end

                len = deformed_bond_length[iID][jID]::Float64
                div_in_place!(bond_norm_all, rotated_bond, len, true)

                mul!(temp, bond_energy, bond_norm_all)

                mul!(temp2, crit_vals, quad_horizon)

                condition .= temp .> temp2

                mul_in_place!(temp, bond_norm_all, condition)

                update_bond_damage!(bond_damage, iID, jID, temp)

                @inbounds for c in condition
                    if c
                        update_list[iID] = true
                        break
                    end
                end
            end  # jID loop
        end  # jID loop timeit
    end  # iID loop
end

@inline function update_bond_damage!(bond_damage::BondScalarState{Float64}, iID::Int64,
                                     jID::Int64, temp::Vector{Float64})
    bond_damage[iID][jID] -= sum(temp)
    if bond_damage[iID][jID] < 0
        bond_damage[iID][jID] = 0.0
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

function init_model(nodes::AbstractVector{Int64},
                    damage_parameter::Dict,
                    block::Int64)
    dof = Data_Manager.get_dof()
    quad_horizon = Data_Manager.create_constant_node_scalar_field("Quad Horizon", Float64)
    Data_Manager.create_constant_bond_vector_state("Bond Displacements", Float64, dof)
    horizon = Data_Manager.get_field("Horizon")
    thickness::Float64 = get(damage_parameter, "Thickness", 1)
    for iID in nodes
        quad_horizon[iID] = get_quad_horizon(horizon[iID], dof, thickness)
    end

    if haskey(damage_parameter, "Anisotropic Damage")
        rotation::Bool = Data_Manager.get_rotation()
        if !rotation
            @error "Anisotropic damage requires Angles field"
            return nothing
        end
    end
end
end
