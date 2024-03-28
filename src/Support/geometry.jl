# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Geometry
using LinearAlgebra
using Rotations
export bond_geometry
export shape_tensor

"""
     bond_geometry(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, nlist, coor, undeformed_bond, undeformed_bond_length)

Calculate bond geometries between nodes based on their coordinates.

# Arguments
 - `nodes::Union{SubArray,Vector{Int64}}`: A vector of integers representing node IDs.
 - `dof::Int64`: An integer representing the degrees of freedom.
 - `nlist`: A data structure (e.g., a list or array) representing neighboring node IDs for each node.
 - `coor`: A matrix representing the coordinates of each node.
 - `undeformed_bond`: A preallocated array or data structure to store bond geometries.
 - `undeformed_bond_length`: A preallocated array or data structure to store bond distances.

# Output
 - `undeformed_bond`: An updated `undeformed_bond` array with calculated bond geometries.

# Description
 This function calculates bond geometries between nodes. For each node in `nodes`, it computes the bond vector between the node and its neighboring nodes based on their coordinates. It also calculates the distance (magnitude) of each bond vector.

 If the distance of any bond vector is found to be zero, indicating identical point coordinates, an error is raised.

# Example
 ```julia
 nodes = [1, 2, 3]
 dof = 2
 nlist = [[2, 3], [1, 3], [1, 2]]
 coor = [0.0 0.0; 1.0 0.0; 0.0 1.0]
 undeformed_bond = zeros(Float64, length(nodes), length(nlist[1]), dof + 1)

 undeformed_bond(nodes, dof, nlist, coor, undeformed_bond)
"""
function bond_geometry(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, nlist, coor, undeformed_bond, undeformed_bond_length)

    for iID in nodes

        # Calculate bond vector and distance
        bond_vectors = coor[nlist[iID], :] .- coor[iID, :]'
        distances = sqrt.(sum(bond_vectors .^ 2, dims=2))

        # Check for identical point coordinates
        if any(distances .== 0)
            @error "Identical point coordinates with no distance $iID"
            return nothing
        end

        undeformed_bond[iID] .= bond_vectors
        undeformed_bond_length[iID] .= distances

    end
    return undeformed_bond, undeformed_bond_length
end

"""
    shape_tensor(nodes::Union{SubArray, Vector{Int64}}, dof::Int64, nlist, volume, omega, bond_damage, undeformed_bond, shapeTensor, inverse_shape_tensor)

Calculate the shape tensor and its inverse for a set of nodes in a computational mechanics context.

# Arguments
- `nodes::Union{SubArray, Vector{Int64}}`: A vector of integers representing node IDs.
- `dof::Int64`: An integer representing the degrees of freedom.
- `nlist`: A data structure (e.g., a list or array) representing neighboring node IDs for each node.
- `volume`: A vector or array containing volume information for each node.
- `omega`: A vector or array containing omega information for each node.
- `bond_damage`: A data structure representing bond damage for each node.
- `undeformed_bond`: A data structure representing bond geometries for each node.
- `shapeTensor`: A preallocated 3D array to store the shape tensors for each node.
- `inverse_shape_tensor`: A preallocated 3D array to store the inverse shape tensors for each node.

# Output
- `shapeTensor`: An updated `shapeTensor` array with calculated shape tensors.
- `inverse_shape_tensor`: An updated `inverse_shape_tensor` array with calculated inverse shape tensors.

# Description
This function calculates the shape tensor and its inverse for a set of nodes in a computational mechanics context. The shape tensor is a key quantity used in continuum mechanics to describe material deformation. It is calculated based on bond damage, bond geometries, volume, and omega information for each node.

For each node in `nodes`, the function iterates through degrees of freedom (`dof`) and computes elements of the shape tensor. The inverse of the shape tensor is also calculated and stored in `inverse_shape_tensor`.

# Example
```julia
nodes = [1, 2, 3]
dof = 3
nlist = [[2, 3], [1, 3], [1, 2]]
volume = [0.1, 0.2, 0.3]
omega = [0.5, 0.4, 0.6]
bond_damage = zeros(Float64, length(nodes), length(nlist[1]))
undeformed_bond = rand(Float64, length(nodes), length(nlist[1]), dof)
shapeTensor = zeros(Float64, length(nodes), dof, dof)
inverse_shape_tensor = zeros(Float64, length(nodes), dof, dof)

shape_tensor(nodes, dof, nlist, volume, omega, bond_damage, undeformed_bond, shapeTensor, inverse_shape_tensor)
"""
function shape_tensor(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, nlist, volume, omega, bond_damage, undeformed_bond, shapeTensor, inverse_shape_tensor)
    shapeTensor .= 0
    for iID in nodes
        for i in 1:dof
            for j in 1:dof
                shapeTensor[iID, i, j] = sum(bond_damage[iID] .* undeformed_bond[iID][:, i] .* undeformed_bond[iID][:, j] .* volume[nlist[iID]] .* omega[iID])
            end

        end
        try
            inverse_shape_tensor[iID, :, :] = inv(shapeTensor[iID, :, :])
        catch ex
            @error "Shape Tensor is singular and cannot be inverted $(ex).\n - Check if your mesh is 3D, but has only one layer of nodes\n - Check number of damaged bonds."
            return nothing, nothing
        end
    end

    return shapeTensor, inverse_shape_tensor
end

"""
    deformation_gradient(nodes::Union{SubArray, Vector{Int64}}, dof::Int64, nlist, volume, omega, bond_damage, undeformed_bond, deformed_bond, inverse_shape_tensor, deformation_gradient)

Calculate the deformation gradient tensor for a set of nodes in a computational mechanics context.

# Arguments
- `nodes::Union{SubArray, Vector{Int64}}`: A vector of integers representing node IDs.
- `dof::Int64`: An integer representing the degrees of freedom.
- `nlist`: A data structure (e.g., a list or array) representing neighboring node IDs for each node.
- `volume`: A vector or array containing volume information for each node.
- `omega`: A vector or array containing omega information for each node.
- `bond_damage`: A data structure representing bond damage for each node.
- `undeformed_bond`: A data structure representing bond geometries for each node.
- `deformed_bond`: A data structure representing deformed bond properties for each node.
- `inverse_shape_tensor`: A data structure representing the inverse shape tensors for each node.
- `deformation_gradient`: A preallocated 3D array to store the deformation gradient tensors for each node.

# Output
- `deformation_gradient`: An updated `deformation_gradient` array with calculated deformation gradient tensors.

# Description
This function calculates the deformation gradient tensor for a set of nodes in a computational mechanics context. The deformation gradient tensor characterizes the deformation of a material.

For each node in `nodes`, the function iterates through degrees of freedom (`dof`) and computes elements of the deformation gradient tensor based on bond damage, deformed bond properties, bond geometries, volume, and omega information. The calculated deformation gradient tensor is stored in `deformation_gradient`.

# Example
```julia
nodes = [1, 2, 3]
dof = 3
nlist = [[2, 3], [1, 3], [1, 2]]
volume = [0.1, 0.2, 0.3]
omega = [0.5, 0.4, 0.6]
bond_damage = zeros(Float64, length(nodes), length(nlist[1]))
undeformed_bond = rand(Float64, length(nodes), length(nlist[1]), dof)
deformed_bond = rand(Float64, length(nodes), length(nlist[1]), dof)
inverse_shape_tensor = rand(Float64, length(nodes), dof, dof)
deformation_gradient = zeros(Float64, length(nodes), dof, dof)

deformation_gradient(nodes, dof, nlist, volume, omega, bond_damage, undeformed_bond, deformed_bond, inverse_shape_tensor, deformation_gradient)
"""
function deformation_gradient(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, nlist::SubArray, volume::SubArray, omega::SubArray, bond_damage::SubArray, deformed_bond::Union{SubArray,Vector{Matrix{Float64}}}, undeformed_bond::SubArray, inverse_shape_tensor::SubArray, deformation_gradient::SubArray)
    deformation_gradient .= 0
    for iID in nodes
        for i in 1:dof
            for j in 1:dof
                deformation_gradient[iID, i, j] = sum(bond_damage[iID] .* deformed_bond[iID][:, i] .* undeformed_bond[iID][:, j] .* volume[nlist[iID]] .* omega[iID])
            end
        end
        deformation_gradient[iID, :, :] *= inverse_shape_tensor[iID, :, :]
    end

    return deformation_gradient
end

"""
    function strain(nodes::Union{SubArray, Vector{Int64}}, deformation_gradient, strain)

Calculate strains for specified nodes based on deformation gradients.

## Arguments
-  `nodes::Union{SubArray, Vector{Int64}}`: List of nodes
- `deformation_gradient`: Deformation gradient at current time step (2D or 3D array).

## Returns

- Updated `strain` array containing strains.

This function iterates over the specified nodes and computes strain at each node using the given deformation gradients.

"""
function strain(nodes::Union{SubArray,Vector{Int64}}, deformation_gradient::SubArray, strain::SubArray)
    # https://en.wikipedia.org/wiki/Strain_(mechanics)
    for iID in nodes
        strain[iID, :, :] = 0.5 * (transpose(deformation_gradient[iID, :, :]) * deformation_gradient[iID, :, :] - I)
    end
    return strain
end

"""
    function rotation_tensor(angles::Vector{Float64})

Creates the rotation tensor for 2D or 3D applications. Uses Rotations.jl package.

# Arguments
-  `angles::Vector{Float64}`: Vector of angles definede in degrees of length one or three

# Returns
- Rotation tensor

"""
function rotation_tensor(angles::Union{Vector{Float64},Vector{Int64}})
    if length(angles) == 3
        return RotXYZ(angles[1] / 180 * pi, angles[2] / 180 * pi, angles[3] / 180 * pi)
    end
    # return RotXYZ(0, 0, angles[1] / 180 * pi)
    return RotXY(0, angles[1] / 180 * pi)
end

end