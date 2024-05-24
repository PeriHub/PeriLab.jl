# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Bond_Deformation_Gradient
include("../Material/Material_Models/Bond_Associated_Correspondence.jl")
using .Bond_Associated_Correspondence: find_local_neighbors
include("../../Support/geometry.jl")
using .Geometry: calculate_bond_length, bond_associated_deformation_gradient
export compute

"""
    compute(datamanager, nodes)

Compute the bond deformation gradient.

# Arguments
- `datamanager`: Datamanager.
- `nodes`: List of nodes.
# Returns
- `datamanager`: Datamanager.
"""


function compute(datamanager::Module, nodes::Union{SubArray,Vector{Int64}})

    dof = datamanager.get_dof()
    nlist = datamanager.get_nlist()
    volume = datamanager.get_field("Volume")
    omega = datamanager.get_field("Influence Function")
    bond_damage = datamanager.get_bond_damage("NP1")
    coordinates = datamanager.get_field("Coordinates")
    deformed_coordinates = datamanager.get_field("Deformed Coordinates", "NP1")

    inverse_bond_shape_tensor = datamanager.get_field("Inverse Bond Associated Shape Tensor")
    bond_deformation_gradient = datamanager.get_field("Bond Associated Deformation Gradient")


    blocks = datamanager.get_field("Block_Id")
    # TODO optimize out. should not be done in every step. Init is enough
    horizon = datamanager.get_field("Horizon")

    for iID in nodes
        bond_horizon = datamanager.get_property(blocks[iID], "Material Model", "Bond Horizon")

        # TODO optimize out. should not be done in every step. Init is enough
        if isnothing(bond_horizon)
            bond_horizon = horizon[iID]
        end
        for (jID, nID) in enumerate(nlist[iID])
            neighbor_nlist = find_local_neighbors(nID, coordinates, nlist[iID], bond_horizon)

            undeformed_bond, distances = calculate_bond_length(nID, coordinates, neighbor_nlist)
            deformed_bond, distances = calculate_bond_length(nID, deformed_coordinates, neighbor_nlist)
            # TODO Bond damage is not correct and must be adapted
            # indices are not needed for that
            indices = vcat(1:length(neighbor_nlist))

            deformation_gradient[iID][jID, :, :] = bond_associated_deformation_gradient(dof, volume[neighbor_nlist], omega[neighbor_nlist], bond_damage[iID][indices], undeformed_bond, deformed_bond, bond_deformation_gradient[iID][jID, :, :])

            deformation_gradient[iID][jID, :, :] *= inverse_bond_shape_tensor[iID][jID, :, :]
        end

    end
    return datamanager
end


end