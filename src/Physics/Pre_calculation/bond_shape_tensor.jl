# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Bond_Shape_Tensor
include("../Material/Material_Models/Bond_Associated_Correspondence_Material.jl")
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
undeformed_bond = datamanager.get_field("Bond Geometry")
bond_shape_tensor = datamanager.get_field("Bond Associated Shape Tensor")
inverse_bond_shape_tensor = datamanager.get_field("Inverse Bond Associated Shape Tensor")


blocks = datamanager.get_field("Block_Id")
# TODO optimize out. should not be done in every step. Init is enough
horizon = datamanager.get_field("Horizon")

for iID in nodes
    bond_horizon = datamanager.get_properties(blocks[iID], "Material Model", "Bond Horizon")
    
    # TODO optimize out. should not be done in every step. Init is enough
    if isnothing(bond_horizon)
        bond_horizon = horizon[iID]
    end
    for (jID, nID) in enumerate(nlist[iID])
        neighbor_nlist = find_local_neighbors(neighbor_coordinate[nID], coordinates[nlist[nlist.!=nID],:], nlist[iID], bond_horizon)
        indices = vcat(1:jID-1, jID+1:length(nlist[iID]))
        shape_tensor[iID][jID,:,:], inverse_shape_tensor[iID][jID,:,:] = Geometry.bond_associated_shape_tensor(dof, volume[neighbor_nlist], omega[neighbor_nlist], bond_damage[iID][indices], undeformed_bond[iID][indices], bond_shape_tensor[iID][jID, :, :], inverse_bond_shape_tensor[iID][jID, :, :])      
      end
    
end
return datamanager
end
end

function shape_tensor(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, nlist, volume, omega, bond_damage, undeformed_bond, shape_tensor, inverse_shape_tensor)
    shape_tensor .= 0
    for iID in nodes
        shape_tensor[iID, :, :] = calculate_shape_tensor(shape_tensor[iID, :, :], dof, volume[nlist[iID]], omega[iID], bond_damage[iID], undeformed_bond[iID])
        try
            inverse_shape_tensor[iID, :, :] = inv(shape_tensor[iID, :, :])
        catch ex
            @error "Shape Tensor is singular and cannot be inverted $(ex).\n - Check if your mesh is 3D, but has only one layer of nodes\n - Check number of damaged bonds."
            return nothing, nothing
        end
    end

    return shape_tensor, inverse_shape_tensor
end

function calculate_shape_tensor(shape_tensor::Matrix{Float64}, dof::Int64, volume, omega, bond_damage, undeformed_bond)

for i in 1:dof
    for j in 1:dof
        shape_tensor[i, j] = sum(bond_damage .* undeformed_bond[:, i] .* undeformed_bond[:, j] .* volume .* omega)
    end
end

return shape_tensor
end