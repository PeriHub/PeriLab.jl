# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using ...Data_Manager
using ...Helpers: find_active_nodes, get_active_update_nodes, add_in_place!, invert
using StaticArrays: MMatrix, SMatrix
using ..Material_Basis:
                        get_strain, get_Hooke_matrix,
                        compute_deviatoric_and_spherical_stresses

"""
    get_forces_from_force_density()

Computes the forces from the force densities.
"""
function get_forces_from_force_density()
    force_density = Data_Manager.get_field("Force Densities", "NP1")
    forces = Data_Manager.get_field("Forces", "NP1")
    volume = Data_Manager.get_field("Volume")
    forces = force_density .* volume
end

"""
    get_partial_stresses(nodes::Vector{Int64})

Computes the partial stresses.

# Arguments
- `nodes::Vector{Int64}`: List of block nodes.
"""
function get_partial_stresses(nodes::AbstractVector{Int64})
    bond_forces = Data_Manager.get_field("Bond Forces")
    undeformed_bond = Data_Manager.get_field("Bond Geometry")
    volume = Data_Manager.get_field("Volume")
    stress = Data_Manager.get_field("Cauchy Stress", "NP1")

    for iID in nodes
        str = @view stress[iID, :, :]
        add_in_place!(str, bond_forces[iID], undeformed_bond[iID], volume[iID])
    end
end

"""
    calculate_von_mises_stress(nodes::AbstractVector{Int64})

Calculate the von Mises stress.

# Arguments
- `nodes::AbstractVector{Int64}`: The nodes.
"""
function calculate_von_mises_stress(nodes::AbstractVector{Int64})
    dof = Data_Manager.get_dof()
    stress_NP1 = Data_Manager.get_field("Cauchy Stress", "NP1")
    von_Mises_stress = Data_Manager.get_field("von Mises Stress", "NP1")

    if dof == 2
        for iID in nodes
            von_Mises_2d!(von_Mises_stress, stress_NP1, iID)
        end
    else  # dof == 3
        for iID in nodes
            von_Mises_3d!(von_Mises_stress, stress_NP1, iID)
        end
    end
end

@inline function von_Mises_2d!(von_Mises_stress::Vector{Float64},
                               stress_NP1::Array{Float64,3}, iID::Int64)
    @inbounds begin
        sigma11 = stress_NP1[iID, 1, 1]
        sigma22 = stress_NP1[iID, 2, 2]
        sigma12 = stress_NP1[iID, 1, 2]

        vm = sigma11 * sigma11 + sigma22 * sigma22 - sigma11 * sigma22 +
             3 * sigma12 * sigma12
        von_Mises_stress[iID] = sqrt(vm)
    end
end

@inline function von_Mises_3d!(von_Mises_stress::Vector{Float64},
                               stress_NP1::Array{Float64,3}, iID::Int64)
    @inbounds begin
        sigma11 = stress_NP1[iID, 1, 1]
        sigma22 = stress_NP1[iID, 2, 2]
        sigma33 = stress_NP1[iID, 3, 3]
        sigma12 = stress_NP1[iID, 1, 2]
        sigma13 = stress_NP1[iID, 1, 3]
        sigma23 = stress_NP1[iID, 2, 3]

        vm = sigma11 * sigma11 + sigma22 * sigma22 + sigma33 * sigma33 -
             sigma11 * sigma22 - sigma11 * sigma33 - sigma22 * sigma33 +
             3 * (sigma12 * sigma12 + sigma13 * sigma13 + sigma23 * sigma23)

        von_Mises_stress[iID] = sqrt(vm)
    end
end

"""
    calculate_strain(nodes::AbstractVector{Int64}, hooke_matrix::Matrix{Float64})

Calculate the von Mises stress.

# Arguments
- `nodes::AbstractVector{Int64}`: The nodes.
- `hooke_matrix::Matrix{Float64}`: The hooke matrix.
"""
function calculate_strain(nodes::AbstractVector{Int64},
                          hooke_matrix::Union{Matrix{Float64},MMatrix,SMatrix})
    stress_NP1 = Data_Manager.get_field("Cauchy Stress", "NP1")
    strain = Data_Manager.get_field("Strain", "NP1")
    for iID in nodes
        strain[iID, :, :] = get_strain(stress_NP1[iID, :, :], hooke_matrix)
    end
end

"""
    calculate_stresses(block_nodes::Dict{Int64,Vector{Int64}}, options::Dict{String, Any})

Computes the stresses.

# Arguments
- `block_nodes::Dict{Int64,Vector{Int64}}`: List of block nodes.
- `options::Dict{String, Any}`: List of options.
"""
function calculate_stresses(block_nodes::Dict{Int64,Vector{Int64}},
                            options::Dict{String,Any})
    active_list = Data_Manager.get_field("Active")
    for block in eachindex(block_nodes)
        correspondence = occursin("Correspondence",
                                  Data_Manager.get_properties(block, "Material Model")["Material Model"])
        if options["Calculate Cauchy"] |
           options["Calculate von Mises stress"] |
           options["Calculate Strain"] | correspondence
            active_nodes = Data_Manager.get_field("Active Nodes")
            active_nodes = find_active_nodes(active_list, active_nodes, block_nodes[block])
            if !correspondence
                get_partial_stresses(active_nodes)
            end
        end
        if options["Calculate von Mises stress"] | correspondence
            calculate_von_mises_stress(active_nodes)
        end
        if options["Calculate Strain"] && !correspondence
            material_parameter = Data_Manager.get_properties(block, "Material Model")
            hookeMatrix = get_Hooke_matrix(material_parameter,
                                           material_parameter["Symmetry"],
                                           Data_Manager.get_dof())
            calculate_strain(active_nodes,
                             invert(hookeMatrix))
        end
    end
end
