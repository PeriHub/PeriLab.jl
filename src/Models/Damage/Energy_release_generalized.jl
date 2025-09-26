# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Energy_release_generalized
include("../Material/Material_Factory.jl")
include("../../Support/Geometry.jl")
include("../../Support/Helpers.jl")
using .Material
using .Geometry
using .Helpers:
                rotate,
                fastdot,
                sub_in_place!,
                div_in_place!,
                mul_in_place!,
                interpol_data,
                is_dependent
using LinearAlgebra
using StaticArrays
export compute_model
export damage_name
export init_model
export fields_for_local_synchronization
function fields_for_local_synchronization(datamanager::Module, model::String)
    download_from_cores = false
    upload_to_cores = true
    datamanager.set_local_synch(model, "Bond Forces", download_from_cores, upload_to_cores)
    return datamanager
end
"""
	damage_name()

Gives the damage name. It is needed for comparison with the yaml input deck.

# Returns
- `name::String`: The name of the damage.

Example:
```julia
println(damage_name())
"Generalized Critical Energy"
```
"""
function damage_name()
    return "Generalized Critical Energy"
end

"""
	compute_model(datamanager, nodes, damage_parameter, block, time, dt)

Calculates the elastic energy of each bond and compares it to a critical one. If it is exceeded, the bond damage value is set to zero.
Enhanced with mixed-mode capability using directional energy optimization.
[WillbergC2019](@cite), [FosterJT2011](@cite)

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `damage_parameter::Dict(String, Any)`: Dictionary with material parameter.
- `block::Int64`: Block number.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
# Returns
- `datamanager::Data_manager`: Datamanager.
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

    GIC = damage_parameter["Critical Value"]
    GIIC = damage_parameter["GIIC"]
    quad_horizons = datamanager.get_field("Quad Horizon")
    inverse_nlist = datamanager.get_inverse_nlist()

    # for anisotropic damage models
    rotation::Bool = datamanager.get_rotation()

    tension::Bool = get(damage_parameter, "Only Tension", false)

    norm_displacement::Float64 = 0.0
    product::Float64 = 0.0
    x_vector::Vector{Float64} = @SVector zeros(Float64, dof)
    x_vector[1] = 1.0

    bond_force::Vector{Float64} = zeros(Float64, dof)
    neighbor_bond_force::Vector{Float64} = zeros(Float64, dof)
    bond_force_delta::Vector{Float64} = zeros(Float64, dof)
    projected_force::Vector{Float64} = zeros(Float64, dof)
    relative_displacement::Vector{Float64} = zeros(Float64, dof)

    sub_in_place!(bond_displacements, deformed_bond, undeformed_bond)
    warning_flag = true

    for iID in nodes
        bond_displacement = bond_displacements[iID]
        bond_force_vec = bond_forces[iID]
        quad_horizon = quad_horizons[iID]

        @fastmath @inbounds for jID in eachindex(nlist[iID])
            relative_displacement = bond_displacement[jID]
            norm_displacement = fastdot(relative_displacement, relative_displacement)
            if norm_displacement == 0 || (tension &&
                deformed_bond_length[iID][jID] - undeformed_bond_length[iID][jID] < 0)
                continue
            end

            bond_force .= bond_force_vec[jID]
            bond_force_delta .= bond_force .- neighbor_bond_force

            product = fastdot(bond_force_delta, relative_displacement)
            mul!(projected_force, product / norm_displacement, relative_displacement)
            product = fastdot(projected_force, relative_displacement, true)
            bond_energy_GIC = 0.25 * product

            product = GIC * quad_horizon
            if bond_energy_GIC > product
                bond_damage[iID][jID] = 0.0
                update_list[iID] = true
                # continue
            end
            if dof == 2
                orth_relative_displacement = [
                    relative_displacement[1],
                    -relative_displacement[2]
                ]
                product = fastdot(bond_force_delta, orth_relative_displacement)
                mul!(projected_force, product / norm_displacement,
                     orth_relative_displacement)
                product = fastdot(projected_force, orth_relative_displacement, true)
                bond_energy_GIIC = 0.25 * product

                product = GIIC * quad_horizon

                if bond_energy_GIIC > abs(product)
                    bond_damage[iID][jID] = 0.0
                    update_list[iID] = true
                end
            end
        end
    end
    return datamanager
end

#function eval_damage!(bond_damage::Vector{Vector{Float64}}, update_list::Vector{Bool} critical_energy_value::Float64)
#
#    product = fastdot(bond_force_delta, relative_displacement)
#    mul!(projected_force, product / norm_displacement, relative_displacement)
#    product = fastdot(projected_force, relative_displacement, true)
#    bond_energy_GIC = 0.25 * product
#
#    if bond_energy > critical_energy_value
#        bond_damage[iID][jID] = 0.0
#        update_list[iID] = true
#    end
#end

"""
	init_model(datamanager, nodes, damage_parameter, block)

Initialize the damage model with enhanced mixed-mode support.
"""
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

    # Enhanced parameter validation for mixed-mode
    mixed_mode_analysis = get(damage_parameter, "Mixed Mode Analysis", false)

    if mixed_mode_analysis
        # Check for required parameters
        has_GIIC = haskey(damage_parameter, "GIIC")
        has_mixed_total = haskey(damage_parameter, "G_total_mixed")

        if !has_GIIC && !has_mixed_total
            @error "Mixed Mode Analysis enabled but neither GIIC nor G_total_mixed provided for block $block"
        elseif !has_GIIC && has_mixed_total
            @info "GIIC will be extracted from G_total_mixed for block $block"
        end
    else
        # Standard validation
        if !haskey(damage_parameter, "GIIC")
            @warn "GIIC parameter not found in damage model $(damage_name()) for block $block. Using GIC as default."
            damage_parameter["GIIC"] = damage_parameter["Critical Value"]
        end
    end

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
end
