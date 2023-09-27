# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Geometry
using LinearAlgebra
export bond_geometry
export shape_tensor
"""
    bond_geometry(nodes::Vector{Int64}, dof::Int64, nlist, coor, bondgeom)

Calculate bond geometries between nodes based on their coordinates.

# Arguments
- `nodes::Vector{Int64}`: A vector of integers representing node IDs.
- `dof::Int64`: An integer representing the degrees of freedom.
- `nlist`: A data structure (e.g., a list or array) representing neighboring node IDs for each node.
- `coor`: A matrix representing the coordinates of each node.
- `bondgeom`: A preallocated array or data structure to store bond geometries.

# Output
- `bondgeom`: An updated `bondgeom` array with calculated bond geometries.

# Description
This function calculates bond geometries between nodes. For each node in `nodes`, it computes the bond vector between the node and its neighboring nodes based on their coordinates. It also calculates the distance (magnitude) of each bond vector.

If the distance of any bond vector is found to be zero, indicating identical point coordinates, an error is raised.

# Example
```julia
nodes = [1, 2, 3]
dof = 2
nlist = [[2, 3], [1, 3], [1, 2]]
coor = [0.0 0.0; 1.0 0.0; 0.0 1.0]
bondgeom = zeros(Float64, length(nodes), length(nlist[1]), dof + 1)

bond_geometry(nodes, dof, nlist, coor, bondgeom)
"""
function bond_geometry(nodes::Vector{Int64}, dof::Int64, nlist, coor, bondgeom)
    for iID in nodes
        for jID in eachindex(nlist[iID])
            bondgeom[iID][jID, 1:dof] = coor[nlist[iID][jID], :] - coor[iID, :]
            bondgeom[iID][jID, dof+1] = norm(bondgeom[iID][jID, 1:dof])
            if bondgeom[iID][jID, dof+1] == 0
                @error "Identical point coordinates with no distance $iID, $jID"
            end
        end
    end
    return bondgeom
end
"""
    shape_tensor(nodes::Vector{Int64}, dof::Int64, nlist, volume, omega, bondDamage, bondGeometry, shapeTensor, invShapeTensor)

Calculate the shape tensor and its inverse for a set of nodes in a computational mechanics context.

# Arguments
- `nodes::Vector{Int64}`: A vector of integers representing node IDs.
- `dof::Int64`: An integer representing the degrees of freedom.
- `nlist`: A data structure (e.g., a list or array) representing neighboring node IDs for each node.
- `volume`: A vector or array containing volume information for each node.
- `omega`: A vector or array containing omega information for each node.
- `bondDamage`: A data structure representing bond damage for each node.
- `bondGeometry`: A data structure representing bond geometries for each node.
- `shapeTensor`: A preallocated 3D array to store the shape tensors for each node.
- `invShapeTensor`: A preallocated 3D array to store the inverse shape tensors for each node.

# Output
- `shapeTensor`: An updated `shapeTensor` array with calculated shape tensors.
- `invShapeTensor`: An updated `invShapeTensor` array with calculated inverse shape tensors.

# Description
This function calculates the shape tensor and its inverse for a set of nodes in a computational mechanics context. The shape tensor is a key quantity used in continuum mechanics to describe material deformation. It is calculated based on bond damage, bond geometries, volume, and omega information for each node.

For each node in `nodes`, the function iterates through degrees of freedom (`dof`) and computes elements of the shape tensor. The inverse of the shape tensor is also calculated and stored in `invShapeTensor`.

# Example
```julia
nodes = [1, 2, 3]
dof = 3
nlist = [[2, 3], [1, 3], [1, 2]]
volume = [0.1, 0.2, 0.3]
omega = [0.5, 0.4, 0.6]
bondDamage = zeros(Float32, length(nodes), length(nlist[1]))
bondGeometry = rand(Float32, length(nodes), length(nlist[1]), dof)
shapeTensor = zeros(Float32, length(nodes), dof, dof)
invShapeTensor = zeros(Float32, length(nodes), dof, dof)

shape_tensor(nodes, dof, nlist, volume, omega, bondDamage, bondGeometry, shapeTensor, invShapeTensor)
"""

function shape_tensor(nodes::Vector{Int64}, dof::Int64, nlist, volume, omega, bondDamage, bondGeometry, shapeTensor, invShapeTensor)

    for iID in nodes
        shapeTensor[iID, :, :] = zeros(Float32, dof, dof)
        for i in 1:dof
            for j in 1:dof
                shapeTensor[iID, i, j] = sum(bondDamage[iID][:] .* bondGeometry[iID][:, i] .* bondGeometry[iID][:, j] .* volume[nlist[iID][:]] .* omega[iID][:])
            end
        end
        invShapeTensor[iID, :, :] = inv(shapeTensor[iID, :, :])
    end

    return shapeTensor, invShapeTensor
end

end