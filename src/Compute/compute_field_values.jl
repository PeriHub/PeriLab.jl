# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

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
    forces[:] .= force_density .* volume
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
function get_partial_stresses(datamanager::Module)
    nnodes = datamanager.get_nnodes()
    bond_forces = datamanager.get_field("Bond Forces", "NP1")
    undeformed_bond = datamanager.get_field("Bond Geometry")
    volume = datamanager.get_field("Volume")
    stress = datamanager.get_field("Cauchy Stress", "NP1")
    dof = datamanager.get_dof()
    for iID in 1:nnodes
        for jID in eachindex(bond_forces[iID][:, 1])
            for i in 1:dof
                for j in 1:dof
                    stress[iID, i, j] += bond_forces[iID][jID, i] * undeformed_bond[iID][jID, j] * volume[iID]
                end
            end
        end
    end
    return datamanager
end

