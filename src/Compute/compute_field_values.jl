# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using ...Helpers: find_active_nodes, get_active_update_nodes, add_in_place!, invert
using StaticArrays: MMatrix, SMatrix
using ..Material_Basis:
                        get_strain, get_Hooke_matrix,
                        compute_deviatoric_and_spherical_stresses
"""
    get_forces_from_force_density(datamanager::Module)

Computes the forces from the force densities.

# Arguments
- `datamanager::Data_Manager`: Datamanager.
# Returns
- `datamanager::Data_Manager`: Datamanager.
"""
function get_forces_from_force_density(datamanager::Module)
    force_density = datamanager.get_field("Force Densities", "NP1")
    forces = datamanager.get_field("Forces", "NP1")
    volume = datamanager.get_field("Volume")
    forces = force_density .* volume
    return datamanager
end

"""
    get_partial_stresses(datamanager::Module, nodes::Vector{Int64})

Computes the partial stresses.

# Arguments
- `datamanager::Data_Manager`: Datamanager.
- `nodes::Vector{Int64}`: List of block nodes.
# Returns
- `datamanager::Data_Manager`: Datamanager.
"""
function get_partial_stresses(datamanager::Module, nodes::AbstractVector{Int64})
    bond_forces = datamanager.get_field("Bond Forces")
    undeformed_bond = datamanager.get_field("Bond Geometry")
    volume = datamanager.get_field("Volume")
    stress = datamanager.get_field("Cauchy Stress", "NP1")

    for iID in nodes
        str = @view stress[iID, :, :]
        add_in_place!(str, bond_forces[iID], undeformed_bond[iID], volume[iID])
    end
    return datamanager
end

"""
    calculate_von_mises_stress(datamanager::Module, nodes::AbstractVector{Int64})

Calculate the von Mises stress.

# Arguments
- `datamanager::Data_Manager`: Datamanager.
- `nodes::AbstractVector{Int64}`: The nodes.
# Returns
- `datamanager::Data_Manager`: Datamanager.
"""
function calculate_von_mises_stress(datamanager::Module,
                                    nodes::AbstractVector{Int64})
    dof = datamanager.get_dof()
    stress_NP1 = datamanager.get_field("Cauchy Stress", "NP1")
    von_Mises_stress = datamanager.get_field("von Mises Stress", "NP1")

    for iID in nodes
        @views von_Mises_stress[iID] = von_Mises(von_Mises_stress[iID],
                                                 stress_NP1[iID, :, :], dof)
    end

    return datamanager
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
    calculate_strain(datamanager::Module, nodes::AbstractVector{Int64}, hooke_matrix::Matrix{Float64})

Calculate the von Mises stress.

# Arguments
- `datamanager::Data_Manager`: Datamanager.
- `nodes::AbstractVector{Int64}`: The nodes.
- `hooke_matrix::Matrix{Float64}`: The hooke matrix.
# Returns
- `datamanager::Data_Manager`: Datamanager.
"""
function calculate_strain(datamanager::Module,
                          nodes::AbstractVector{Int64},
                          hooke_matrix::Union{Matrix{Float64},MMatrix,SMatrix})
    stress_NP1 = datamanager.get_field("Cauchy Stress", "NP1")
    strain = datamanager.get_field("Strain", "NP1")
    for iID in nodes
        strain[iID, :, :] = get_strain(stress_NP1[iID, :, :], hooke_matrix)
    end
    return datamanager
end

"""
    calculate_stresses(datamanager::Module, block_nodes::Dict{Int64,Vector{Int64}}, options::Dict{String, Any})

Computes the stresses.

# Arguments
- `datamanager::Data_Manager`: Datamanager.
- `block_nodes::Dict{Int64,Vector{Int64}}`: List of block nodes.
- `options::Dict{String, Any}`: List of options.
# Returns
- `datamanager::Data_Manager`: Datamanager.
"""
function calculate_stresses(datamanager::Module,
                            block_nodes::Dict{Int64,Vector{Int64}},
                            options::Dict{String,Any})
    active_list = datamanager.get_field("Active")
    for block in eachindex(block_nodes)
        correspondence = occursin("Correspondence",
                                  datamanager.get_properties(block, "Material Model")["Material Model"])
        if options["Calculate Cauchy"] |
           options["Calculate von Mises stress"] |
           options["Calculate Strain"] | correspondence
            active_nodes = datamanager.get_field("Active Nodes")
            active_nodes = find_active_nodes(active_list, active_nodes, block_nodes[block])
            if !correspondence
                datamanager = get_partial_stresses(datamanager, active_nodes)
            end
        end
        if options["Calculate von Mises stress"] | correspondence
            datamanager = calculate_von_mises_stress(datamanager, active_nodes)
        end
        if options["Calculate Strain"] && !correspondence
            material_parameter = datamanager.get_properties(block, "Material Model")
            hookeMatrix = get_Hooke_matrix(datamanager,
                                           material_parameter,
                                           material_parameter["Symmetry"],
                                           datamanager.get_dof())
            datamanager = calculate_strain(datamanager,
                                           active_nodes,
                                           invert(hookeMatrix,
                                                  "Hooke matrix not invertable"))
        end
    end

    return datamanager
end
