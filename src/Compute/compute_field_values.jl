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

    for iID in nodes
        @views von_Mises_stress[iID] = von_Mises(von_Mises_stress[iID],
                                                 stress_NP1[iID, :, :], dof)
    end
end
function von_Mises(von_Mises_stress, stress, dof)
    # faster than putting the first loop in the second one
    von_Mises_stress = 0
    for i in 1:dof
        von_Mises_stress += stress[i, i] * stress[i, i]
    end
    for i in 1:dof
        for j in (i + 1):dof
            von_Mises_stress += -stress[i, i] * stress[j, j] +
                                3 * stress[i, j] * stress[i, j]
        end
    end
    von_Mises_stress = sqrt(von_Mises_stress)
    return von_Mises_stress
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
                             inv(hookeMatrix))
        end
    end
end
