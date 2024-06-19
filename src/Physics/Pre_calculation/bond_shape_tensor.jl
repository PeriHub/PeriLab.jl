# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Bond_Shape_Tensor
include("../Material/Material_Models/Bond_Associated_Correspondence.jl")

using .Bond_Associated_Correspondence: find_local_neighbors
include("../../Support/geometry.jl")
using .Geometry: calculate_bond_length, bond_associated_shape_tensor

export compute

"""
    compute(datamanager, nodes)

Compute the bond shape tensor.

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
    bond_shape_tensor = datamanager.get_field("Bond Associated Shape Tensor")
    inverse_bond_shape_tensor = datamanager.get_field("Inverse Bond Associated Shape Tensor")
    blocks = datamanager.get_field("Block_Id")


    for iID in nodes
        bond_horizon = datamanager.get_property(blocks[iID], "Material Model", "Bond Horizon")

        # TODO bond damage
        for (jID, nID) in enumerate(nlist[iID])

            neighbor_nlist = find_local_neighbors(nID, coordinates, nlist[iID], bond_horizon)
            undeformed_bond, distances = calculate_bond_length(nID, coordinates, neighbor_nlist)

            # TODO Bond damage and Omega are not correct and must be adapted
            # indices are not needed for that
            indices = vcat(1:length(neighbor_nlist))


            bond_shape_tensor[iID][jID, :, :], inverse_bond_shape_tensor[iID][jID, :, :] = bond_associated_shape_tensor(dof, volume[neighbor_nlist], omega[iID][indices], bond_damage[iID][indices], undeformed_bond, bond_shape_tensor[iID][jID, :, :], inverse_bond_shape_tensor[iID][jID, :, :])
        end

    end
    return datamanager
end
end