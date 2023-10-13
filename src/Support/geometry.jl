# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Geometry
using LinearAlgebra
using Rotations
export bond_geometry
export shape_tensor
"""
     bond_geometry(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, nlist, coor, bondgeom)

 Calculate bond geometries between nodes based on their coordinates.

# # Arguments
 - `nodes::Union{SubArray,Vector{Int64}}`: A vector of integers representing node IDs.
 - `dof::Int64`: An integer representing the degrees of freedom.
 - `nlist`: A data structure (e.g., a list or array) representing neighboring node IDs for each node.
 - `coor`: A matrix representing the coordinates of each node.
 - `bondgeom`: A preallocated array or data structure to store bond geometries.

# # Output
 - `bondgeom`: An updated `bondgeom` array with calculated bond geometries.

# # Description
 This function calculates bond geometries between nodes. For each node in `nodes`, it computes the bond vector between the node and its neighboring nodes based on their coordinates. It also calculates the distance (magnitude) of each bond vector.

 If the distance of any bond vector is found to be zero, indicating identical point coordinates, an error is raised.

# # Example
 ```julia
 nodes = [1, 2, 3]
 dof = 2
 nlist = [[2, 3], [1, 3], [1, 2]]
 coor = [0.0 0.0; 1.0 0.0; 0.0 1.0]
 bondgeom = zeros(Float64, length(nodes), length(nlist[1]), dof + 1)

 bond_geometry(nodes, dof, nlist, coor, bondgeom)
"""
function bond_geometry(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, nlist, coor, bondgeom)
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
    shape_tensor(nodes::Union{SubArray, Vector{Int64}}, dof::Int64, nlist, volume, omega, bondDamage, bondGeometry, shapeTensor, invShapeTensor)

Calculate the shape tensor and its inverse for a set of nodes in a computational mechanics context.

# Arguments
- `nodes::Union{SubArray, Vector{Int64}}`: A vector of integers representing node IDs.
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
bondDamage = zeros(Float64, length(nodes), length(nlist[1]))
bondGeometry = rand(Float64, length(nodes), length(nlist[1]), dof)
shapeTensor = zeros(Float64, length(nodes), dof, dof)
invShapeTensor = zeros(Float64, length(nodes), dof, dof)

shape_tensor(nodes, dof, nlist, volume, omega, bondDamage, bondGeometry, shapeTensor, invShapeTensor)
"""

function shape_tensor(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, nlist, volume, omega, bondDamage, bondGeometry, shapeTensor, invShapeTensor)

    for iID in nodes
        shapeTensor[iID, :, :] = zeros(Float64, dof, dof)
        for i in 1:dof
            for j in 1:dof
                shapeTensor[iID, i, j] = sum(bondDamage[iID][:] .* bondGeometry[iID][:, i] .* bondGeometry[iID][:, j] .* volume[nlist[iID][:]] .* omega[iID][:])
            end

        end
        try
            invShapeTensor[iID, :, :] = inv(shapeTensor[iID, :, :])
        catch ex
            @error "Shape Tensor is singular and cannot be inverted $(ex).\n - Check if your mesh is 3D, but has only one layer of nodes\n - Check number of damaged bonds."
        end
    end

    return shapeTensor, invShapeTensor
end

"""
    deformation_gradient(nodes::Union{SubArray, Vector{Int64}}, dof::Int64, nlist, volume, omega, bondDamage, bondGeometry, deformed_bond, invShapeTensor, defGrad)

Calculate the deformation gradient tensor for a set of nodes in a computational mechanics context.

# Arguments
- `nodes::Union{SubArray, Vector{Int64}}`: A vector of integers representing node IDs.
- `dof::Int64`: An integer representing the degrees of freedom.
- `nlist`: A data structure (e.g., a list or array) representing neighboring node IDs for each node.
- `volume`: A vector or array containing volume information for each node.
- `omega`: A vector or array containing omega information for each node.
- `bondDamage`: A data structure representing bond damage for each node.
- `bondGeometry`: A data structure representing bond geometries for each node.
- `deformed_bond`: A data structure representing deformed bond properties for each node.
- `invShapeTensor`: A data structure representing the inverse shape tensors for each node.
- `defGrad`: A preallocated 3D array to store the deformation gradient tensors for each node.

# Output
- `defGrad`: An updated `defGrad` array with calculated deformation gradient tensors.

# Description
This function calculates the deformation gradient tensor for a set of nodes in a computational mechanics context. The deformation gradient tensor characterizes the deformation of a material.

For each node in `nodes`, the function iterates through degrees of freedom (`dof`) and computes elements of the deformation gradient tensor based on bond damage, deformed bond properties, bond geometries, volume, and omega information. The calculated deformation gradient tensor is stored in `defGrad`.

# Example
```julia
nodes = [1, 2, 3]
dof = 3
nlist = [[2, 3], [1, 3], [1, 2]]
volume = [0.1, 0.2, 0.3]
omega = [0.5, 0.4, 0.6]
bondDamage = zeros(Float64, length(nodes), length(nlist[1]))
bondGeometry = rand(Float64, length(nodes), length(nlist[1]), dof)
deformed_bond = rand(Float64, length(nodes), length(nlist[1]), dof)
invShapeTensor = rand(Float64, length(nodes), dof, dof)
defGrad = zeros(Float64, length(nodes), dof, dof)

deformation_gradient(nodes, dof, nlist, volume, omega, bondDamage, bondGeometry, deformed_bond, invShapeTensor, defGrad)
"""

function deformation_gradient(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, nlist, volume, omega, bondDamage, deformed_bond::Union{SubArray,Vector{Matrix{Float64}}}, bondGeometry, invShapeTensor, defGrad)
    for iID in nodes
        defGrad[iID, :, :] = zeros(Float64, dof, dof)
        for i in 1:dof
            for j in 1:dof
                defGrad[iID, i, j] = sum(bondDamage[iID][:] .* deformed_bond[iID][:, i] .* bondGeometry[iID][:, j] .* volume[nlist[iID][:]] .* omega[iID][:])
            end
        end
        defGrad[iID, :, :] *= invShapeTensor[iID, :, :]
    end

    return defGrad

end
"""
    function strain(nodes::Union{SubArray, Vector{Int64}}, defGrad, strainInc)

Calculate strain increments for specified nodes based on deformation gradients.

## Arguments
-  `nodes::Union{SubArray, Vector{Int64}}`: List of nodes
- `defGrad`: Deformation gradient at current time step (2D or 3D array).

## Returns

- Updated `strainInc` array containing strain increments.

This function iterates over the specified nodes and computes strain increments using the given deformation gradients.

"""

function strain(nodes::Union{SubArray,Vector{Int64}}, defGrad, strain)
    # https://en.wikipedia.org/wiki/Strain_(mechanics)
    # First equation gives Strain increment as shown
    for iID in nodes
        strain[iID, :, :] = 0.5 * (transpose(defGrad[iID, :, :]) * defGrad[iID, :, :] - I)
    end
    return strain
end

"""
    function rotation_tensor(angles::Vector{Float64})

Creates the rotation tensor for 2D or 3D applications. Uses Rotations.jl package.

## Arguments
-  `angles::Vector{Float64}`: Vector of angles definede in degrees of length one or three


## Returns

- Rotation tensor

"""

function rotation_tensor(angles::Vector{Float64})
    if length(angles) == 3
        return RotXYZ(angles[1] / 180 * pi, angles[2] / 180 * pi, angles[3] / 180 * pi)
    end
    return RotXYZ(0, 0, angles[1] / 180 * pi)
end

end