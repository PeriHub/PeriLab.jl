# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Critical_Energy_Model
include("../Material/Material_Factory.jl")
include("../../Support/geometry.jl")
include("../Pre_calculation/Pre_Calculation_Factory.jl")
using .Material
using .Geometry
using .Pre_calculation
using LinearAlgebra
export compute_damage
export compute_damage_pre_calculation
export damage_name

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
    calculate_energy_components(total_energy::Float64, angle_degrees::Float64)

Calculates the horizontal and vertical components of an elastic energy release.

# Arguments
- `total_energy::Float64`: Total elastic energy.
- `angle_degrees::Float64`: Angle in degrees.

# Returns
- `energy_horizontal::Float64`: Horizontal elastic energy.
- `energy_vertical::Float64`: Vertical elastic energy.
"""
function calculate_energy_components(total_energy::Float64, angle_degrees::Float64)
    # Convert angle from degrees to radians
    angle_radians = deg2rad(angle_degrees)

    # Calculate horizontal and vertical components
    energy_horizontal = total_energy * cos(angle_radians)
    energy_vertical = total_energy * sin(angle_radians)

    return abs(energy_horizontal), abs(energy_vertical)
end

"""
    compute_damage(datamanager, nodes, damage_parameter, block, time, dt)

Calculates the elastic energy of each bond and compares it to a critical one. If it is exceeded, the bond damage value is set to zero.
[WillbergC2019](@cite), [FosterJT2011](@cite)

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
    dof::Int64 = datamanager.get_dof()
    nlist = datamanager.get_nlist()
    block_list = datamanager.get_block_list()
    block_ids = datamanager.get_field("Block_Id")
    update_list = datamanager.get_field("Update List")
    horizon = datamanager.get_field("Horizon")
    bond_damage = datamanager.get_field("Bond Damage", "NP1")
    bond_damage_aniso = datamanager.get_field("Bond Damage Anisotropic")
    undeformed_bond = datamanager.get_field("Bond Geometry")
    bond_forces = datamanager.get_field("Bond Forces")
    deformed_bond = datamanager.get_field("Deformed Bond Geometry", "NP1")
    critical_energy = damage_parameter["Critical Value"]
    quad_horizon = datamanager.get_field("Quad Horizon")
    inverse_nlist = datamanager.get_inverse_nlist()
    # for anisotropic damage models
    rotation::Bool, angles = datamanager.rotation_data()

    tension::Bool = get(damage_parameter, "Only Tension", true)
    inter_block_damage::Bool = haskey(damage_parameter, "Interblock Damage")
    if inter_block_damage
        inter_critical_energy::Array{Float64,3} = datamanager.get_crit_values_matrix()
    end
    aniso_damage::Bool = haskey(damage_parameter, "Anisotropic Damage")
    if aniso_damage
        aniso_crit_values = datamanager.get_aniso_crit_values()
        bond_norm::Float64 = 0.0
    end


    bond_energy::Float64 = 0.0
    norm_displacement::Float64 = 0.0
    x_vector::Vector{Float64} = zeros(Float64, dof)
    x_vector[1] = 1.0

    neighbor_bond_force::Vector{Float64} = zeros(Float64, dof)
    projected_force::Vector{Float64} = zeros(Float64, dof)

    for iID in nodes
        @views relative_displacement_vector = deformed_bond[iID][:, 1:dof] .- undeformed_bond[iID][:, 1:dof]
        if aniso_damage
            rotation_tensor = Geometry.rotation_tensor(angles[iID, :])
        end

        for (jID, neighborID) in enumerate(nlist[iID])
            norm_displacement = norm(relative_displacement_vector[jID, :])

            if norm_displacement == 0
                continue
            end
            if tension && deformed_bond[iID][jID, end] - undeformed_bond[iID][jID, end] < 0
                continue
            end

            # check if the bond also exist at other node, due to different horizons
            if haskey(inverse_nlist[neighborID], iID)
                @views neighbor_bond_force .= bond_forces[neighborID][inverse_nlist[neighborID][iID], :]
            else
                fill!(neighbor_bond_force, 0.0)
            end

            @views projected_force .= dot(bond_forces[iID][jID, :] - neighbor_bond_force, relative_displacement_vector[jID, :]) / (norm_displacement * norm_displacement) .* relative_displacement_vector[jID, :]

            bond_energy = 0.25 * dot(abs.(projected_force), abs.(relative_displacement_vector[jID, :]))
            if bond_energy < 0
                @error "Bond energy smaller zero"
                return nothing
            end

            crit_energy = inter_block_damage ? inter_critical_energy[block_ids[iID], block_ids[neighborID], block] : critical_energy


            if aniso_damage
                rotated_bond = rotation_tensor[1:dof, 1:dof]' * deformed_bond[iID][jID, 1:dof]
                for i in 1:dof
                    if bond_damage_aniso[iID][jID, i] == 0 || rotated_bond[i] == 0
                        continue
                    end
                    @views bond_norm = abs(rotated_bond[i]) / deformed_bond[iID][jID, end]
                    if bond_energy / quad_horizon[iID] * bond_norm > aniso_crit_values[block_ids[iID]][i]
                        @inbounds bond_damage[iID][jID] -= bond_norm # TODO: check if this is correct
                        if bond_damage[iID][jID] < 0
                            bond_damage[iID][jID] = 0
                        end
                        bond_damage_aniso[iID][jID, i] = 0
                        @inbounds update_list[iID] = true
                    end
                end
            else
                if (bond_energy / quad_horizon[iID]) > crit_energy
                    @inbounds bond_damage[iID][jID] = 0.0
                    @inbounds update_list[iID] = true
                end
            end
            if 1 < bond_damage[iID][jID] || bond_damage[iID][jID] < 0
                @error "Bond damage out of bounds"
                return nothing
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

    synchronise_field(datamanager.get_comm(), datamanager.get_synch_fields(), datamanager.get_overlap_map(), datamanager.get_field, "Bond Forces", "upload_to_cores")

    return datamanager
end

"""
    get_quad_horizon(horizon::Float64, dof::Int64)

Get the quadric of the horizon.

# Arguments
- `horizon::Float64`: The horizon of the block.
- `dof::Int64`: The degree of freedom.
# Returns
- `quad_horizon::Float64`: The quadric of the horizon.
"""
function get_quad_horizon(horizon::Float64, dof::Int64)
    if dof == 2
        thickness::Float64 = 1 # is a placeholder
        return Float64(3 / (pi * horizon^3 * thickness))
    end
    return Float64(4 / (pi * horizon^4))
end
function init_damage_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, damage_parameter::Dict, block::Int64)
    quad_horizon = datamanager.create_constant_node_field("Quad Horizon", Float64, 1)
    horizon = datamanager.get_field("Horizon")
    dof = datamanager.get_dof()
    for iID in nodes
        quad_horizon[iID] = get_quad_horizon(horizon[iID], dof)
    end
    return datamanager
end
end