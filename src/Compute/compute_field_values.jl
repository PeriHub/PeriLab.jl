# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using .Helpers: get_active_update_nodes
using StaticArrays: MMatrix, SMatrix
using .Helpers: invert
include("../Models/Material/material_basis.jl")

"""
    get_forces_from_force_density(datamanager::Module)

Computes the forces from the force densities.

# Arguments
- `datamanager::Data_manager`: Datamanager.
# Returns
- `datamanager::Data_manager`: Datamanager.
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
- `datamanager::Data_manager`: Datamanager.
- `nodes::Vector{Int64}`: List of block nodes.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function get_partial_stresses(datamanager::Module, nodes::Union{SubArray,Vector{Int64}})
    bond_forces = datamanager.get_field("Bond Forces")
    undeformed_bond = datamanager.get_field("Bond Geometry")
    volume = datamanager.get_field("Volume")
    stress = datamanager.get_field("Cauchy Stress", "NP1")
    dof = datamanager.get_dof()
    for iID in nodes
        stress[iID, :, :] .+= bond_forces[iID]' * undeformed_bond[iID] .* volume[iID]
    end
    return datamanager
end

"""
    calculate_von_mises_stress(datamanager::Module, nodes::Union{SubArray,Vector{Int64}})

Calculate the von Mises stress.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: The nodes.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function calculate_von_mises_stress(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
)
    dof = datamanager.get_dof()
    stress_NP1 = datamanager.get_field("Cauchy Stress", "NP1")
    von_Mises_stress = datamanager.get_field("von Mises Stress", "NP1")
    for iID in nodes
        von_Mises_stress[iID], spherical_stress_NP1, deviatoric_stress_NP1 =
            get_von_mises_stress(von_Mises_stress[iID], dof, stress_NP1[iID, :, :])
    end
    return datamanager
end

"""
    calculate_strain(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, hooke_matrix::Matrix{Float64})

Calculate the von Mises stress.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: The nodes.
- `hooke_matrix::Matrix{Float64}`: The hooke matrix.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function calculate_strain(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    hooke_matrix::Union{Matrix{Float64},MMatrix,SMatrix},
)
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
- `datamanager::Data_manager`: Datamanager.
- `block_nodes::Dict{Int64,Vector{Int64}}`: List of block nodes.
- `options::Dict{String, Any}`: List of options.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function calculate_stresses(
    datamanager::Module,
    block_nodes::Dict{Int64,Vector{Int64}},
    options::Dict{String,Any},
)

    for block in eachindex(block_nodes)
        if !occursin(
            "Correspondence",
            datamanager.get_properties(block, "Material Model")["Material Model"],
        )
            if options["Calculate Cauchy"] |
               options["Calculate von Mises"] |
               options["Calculate Strain"]
                active = datamanager.get_field("Active")
                update_list = datamanager.get_field("Update List")
                active_nodes, update_nodes =
                    get_active_update_nodes(active, update_list, block_nodes, block)
                datamanager = get_partial_stresses(datamanager, active_nodes)
            end
            if options["Calculate von Mises"]
                datamanager = calculate_von_mises_stress(datamanager, active_nodes)
            end
            if options["Calculate Strain"]
                material_parameter = datamanager.get_properties(block, "Material Model")
                hookeMatrix = get_Hooke_matrix(
                    material_parameter,
                    material_parameter["Symmetry"],
                    datamanager.get_dof(),
                )
                datamanager = calculate_strain(
                    datamanager,
                    active_nodes,
                    invert(hookeMatrix, "Hook matrix not invertable"),
                )
            end
        end
    end

    return datamanager
end
