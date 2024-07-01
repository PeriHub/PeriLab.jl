# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Geometry
using LinearAlgebra
using TensorOperations
using Combinatorics: levicivita
using Rotations
export bond_geometry
export shape_tensor

"""
     bond_geometry(nodes::Union{SubArray,Vector{Int64}}, nlist, coor, undeformed_bond, undeformed_bond_length)

Calculate bond geometries between nodes based on their coordinates.

# Arguments
 - `nodes::Union{SubArray,Vector{Int64}}`: A vector of integers representing node IDs.
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
function bond_geometry(nodes::Union{SubArray,Vector{Int64}}, nlist::Union{SubArray,Vector{Int64}}, coor::Union{SubArray,Matrix{Float64},Matrix{Int64}}, undeformed_bond, undeformed_bond_length)

    for iID in nodes
        undeformed_bond[iID], undeformed_bond_length[iID] = calculate_bond_length(iID, coor, nlist[iID])
        if any(undeformed_bond_length[iID] .== 0)
            @error "Identical point coordinates with no distance $iID"
            return nothing
        end
    end
    return undeformed_bond, undeformed_bond_length
end

function calculate_bond_length(iID::Int64, coor::Union{SubArray,Matrix{Float64},Matrix{Int64}}, nlist::Vector{Int64})

    bond_vectors = coor[nlist, :] .- coor[iID, :]'
    # distances = sqrt.(sum(bond_vectors .^ 2, dims=2))[:]

    # Check for identical point coordinates
    return bond_vectors, norm.(eachrow(bond_vectors))
end

"""
    shape_tensor(nodes::Union{SubArray, Vector{Int64}}, dof::Int64, nlist, volume, omega, bond_damage, undeformed_bond, shape_tensor, inverse_shape_tensor)

Calculate the shape tensor and its inverse for a set of nodes in a computational mechanics context.

# Arguments
- `nodes::Union{SubArray, Vector{Int64}}`: A vector of integers representing node IDs.
- `dof::Int64`: An integer representing the degrees of freedom.
- `nlist`: A data structure (e.g., a list or array) representing neighboring node IDs for each node.
- `volume`: A vector or array containing volume information for each node.
- `omega`: A vector or array containing omega information for each node.
- `bond_damage`: A data structure representing bond damage for each node.
- `undeformed_bond`: A data structure representing bond geometries for each node.
- `shape_tensor`: A preallocated 3D array to store the shape tensors for each node.
- `inverse_shape_tensor`: A preallocated 3D array to store the inverse shape tensors for each node.

# Output
- `shape_tensor`: An updated `shape_tensor` array with calculated shape tensors.
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
shape_tensor = zeros(Float64, length(nodes), dof, dof)
inverse_shape_tensor = zeros(Float64, length(nodes), dof, dof)

shape_tensor(nodes, dof, nlist, volume, omega, bond_damage, undeformed_bond, shape_tensor, inverse_shape_tensor)
"""
function shape_tensor(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, nlist, volume, omega, bond_damage, undeformed_bond, shape_tensor, inverse_shape_tensor)

    for iID in nodes
        shape_tensor[iID, :, :] = calculate_shape_tensor(dof, volume[nlist[iID]], omega[iID], bond_damage[iID], undeformed_bond[iID])
        try
            inverse_shape_tensor[iID, :, :] = inv(shape_tensor[iID, :, :])
        catch ex
            @error "Shape Tensor is singular and cannot be inverted $(ex).\n - Check if your mesh is 3D, but has only one layer of nodes\n - Check number of damaged bonds."
            return nothing, nothing
        end
    end

    return shape_tensor, inverse_shape_tensor
end


function calculate_shape_tensor(dof::Int64, volume, omega, bond_damage, undeformed_bond)
    shape_tensor = zeros(Float64, dof, dof)
    # Compute the element-wise product once
    weighted_bond_damage = bond_damage .* volume .* omega

    @views Threads.@threads for i in 1:dof
        for j in 1:dof
            shape_tensor[i, j] = sum(weighted_bond_damage .* undeformed_bond[:, i] .* undeformed_bond[:, j])
        end
    end

    return shape_tensor
end


function bond_associated_shape_tensor(dof::Int64, volume, omega, bond_damage, undeformed_bond, shape_tensor, inverse_shape_tensor)
    # bond geometries -> zwischen nachbar und seinen nachbarn
    shape_tensor[:, :] = calculate_shape_tensor(dof, volume, omega, bond_damage, undeformed_bond)
    try
        inverse_shape_tensor[:, :] = inv(shape_tensor[:, :])
    catch ex
        @error "Shape Tensor is singular and cannot be inverted $(ex).\n - Check if your mesh is 3D, but has only one layer of nodes\n - Check number of damaged bonds."
        return nothing, nothing
    end

    return shape_tensor, inverse_shape_tensor

end



function bond_associated_deformation_gradient(dof::Int64, volume, omega, bond_damage, undeformed_bond, deformed_bond, deformation_gradient)
    # bond deformation -> zwischen nachbar und seinen nachbarn
    return calculate_deformation_gradient(deformation_gradient, dof, bond_damage, deformed_bond, undeformed_bond, volume, omega)

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
        deformation_gradient[iID, :, :] = calculate_deformation_gradient(deformation_gradient[iID, :, :], dof, bond_damage[iID], deformed_bond[iID], undeformed_bond[iID], volume[nlist[iID]], omega[iID])
        deformation_gradient[iID, :, :] *= inverse_shape_tensor[iID, :, :]
    end

    return deformation_gradient
end


"""
compute_left_stretch_tensor(nodes::Union{SubArray,Vector{Int64}}, deformation_gradient::Union{SubArray,Array{Float64}}, left_stretch_tensor::Union{SubArray,Array{Float64}})

Calculate bond geometries between nodes based on their coordinates.

# Arguments
 - `nodes::Union{SubArray,Vector{Int64}}`: A vector of integers representing node IDs.
 - `deformation_gradient::Union{SubArray,Array{Float64}}`: Deformation gradient.
 - `left_stretch_tensor::Union{SubArray,Array{Float64}}`: Left stretch tensor.

# Output
 - `left_stretch_tensor`: Left stretch tensor.

# Description
 This function computes the left stretch tensor ``\\mathbf{U]``. The deformation gradient ``\\mathbf{F]`` can be splited in a rotation part and left or right stretch tensor.
``\\mathbf{F}=\\mathbf{VR}=\\mathbf{RU}``
To get the left stretch tensor, you can make use of the orthogonal character of your rotation tensor.
``\\mathbf{RR}^T=\\mathbf{I}``
Therefore,
``\\mathbf{U}=\\sqrt{\\mathbf{F}\\mathbf{F}^T}=\\sqrt{\\mathbf{UR}\\mathbf{R}^T\\mathbf{U}^T}=\\sqrt{\\mathbf{UIU}^T}=\\sqrt{\\mathbf{U}^2}``
 
"""

function compute_left_stretch_tensor(nodes::Union{SubArray,Vector{Int64}}, deformation_gradient::Union{SubArray,Array{Float64}}, left_stretch_tensor::Union{SubArray,Array{Float64}})
    for iID in nodes
        left_stretch_tensor[iID, :, :] = deformation_gradient[iID, :, :] * transpose(deformation_gradient[iID, :, :])
    end
    return sqrt.(left_stretch_tensor)
end

function calculate_deformation_gradient(deformation_gradient, dof::Int64, bond_damage, deformed_bond, undeformed_bond, volume::Union{Vector{Int64},Vector{Float64}}, omega)
    for i in 1:dof
        for j in 1:dof
            deformation_gradient[i, j] = sum(bond_damage .* deformed_bond[:, i] .* undeformed_bond[:, j] .* volume .* omega)
        end
    end
    return deformation_gradient
end

function compute_weighted_deformation_gradient(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, nlist, volume, gradient_weight, displacement, velocity, deformation_gradient, deformation_gradient_dot)

    for iID in nodes
        deformation_gradient[iID, :, :] = Matrix{Float64}(I(dof))
        deformation_gradient_dot[iID, :, :] .= 0

        disp_state = displacement[nlist[iID], :] .- displacement[iID, :]'
        velocity_state = velocity[nlist[iID], :] .- velocity[iID, :]'
        for (jID, nID) in enumerate(nlist[iID])
            deformation_gradient[iID, :, :] += disp_state[jID, :] * transpose(gradient_weight[jID, :]) .* volume[nID]
            deformation_gradient_dot[iID, :, :] += velocity_state[jID, :] * transpose(gradient_weight[jID, :]) .* volume[nID]
        end

    end
    return deformation_gradient, deformation_gradient_dot
end


function no_name_yet(nodes, datamanager)
    dof = datamanager.get_dof()
    nlist = datamanager.get_nlist()
    volume = datamanager.get_field("Volume")
    gradient_weight = datamanager.get_field("Gradient Weight")
    displacement = datamanager.get_field("Displacements", "NP1")
    velocity = datamanager.get_field("Velocity", "NP1")
    deformation_gradient = datamanager.get_field("Deformation Gradient", "NP1")
    deformation_gradient_dot = datamanager.get_field("Deformation Gradient Dot", "NP1")
    left_stretch_tensorN = datamanager.get_field("Left Stretch Tensor", "N") # 
    left_stretch_tensorNP1 = datamanager.get_field("Left Stretch Tensor", "NP1")
    unrotated_rate_of_deformation = datamanager.get_field("Unrotated Rate Of Deformation")
    rot_tensorN = datamanager.get_field("Rotation Tensor", "N")
    rot_tensorNP1 = datamanager.get_field("Rotation Tensor", "NP1")
    deformation_gradient, deformation_gradient_dot = compute_weighted_deformation_gradient(nodes, dof, nlist, volume, gradient_weight, displacement, velocity, deformation_gradient, deformation_gradient_dot)

    left_stretch_tensorNP1 = compute_left_stretch_tensor(nodes, deformation_gradient, left_stretch_tensorNP1)

    #unrotated_rate_of_deformation, rotTensorNP1 = deformation_gradient_decomposition(nodes, dof, deformation_gradient, deformation_gradient_dot, left_stretch_tensorN, left_stretch_tensorNP1, rot_tensorN, rot_tensorNP1, unrotated_rate_of_deformation)
end
function no_name_2(nodes, datamanager)
    dof = datamanager.get_dof()
    nlist = datamanager.get_nlist()
    volume = datamanager.get_field("Volume")
    gradient_weight = datamanager.get_field("Gradient Weight")
    displacement = datamanager.get_field("Displacements", "NP1")
    velocity = datamanager.get_field("Velocity", "NP1")
    deformation_gradient = datamanager.get_field("Deformation Gradient", "NP1")
    deformation_gradient_dot = datamanager.get_field("Deformation Gradient Dot")
    left_stretch_tensor = datamanager.get_field("Left Stretch Tensor") # 
    unrotated_rate_of_deformation = datamanager.get_field("Unrotated Rate Of Deformation")
    rot_tensor = datamanager.get_field("Rotation Tensor", "NP1")
    deformation_gradient, deformation_gradient_dot = compute_weighted_deformation_gradient(nodes, dof, nlist, volume, gradient_weight, displacement, velocity, deformation_gradient, deformation_gradient_dot)
    left_stretch_tensor = compute_left_stretch_tensor(nodes, deformation_gradient, left_stretch_tensor)

end

function deformation_gradient_decomposition(nodes, deformation_gradient, deformation_gradient_dot, left_stretch_tensor, rot_tensor, unrotated_rate_of_deformation)


    for iID in nodes
        # Follow the algorithm at page 315
        # Compute rate-of-deformation tensor, D = 1/2 * (L + Lt) -> from Eq (4)
        velocity_gradient = deformation_gradient_dot[iID, :, :] * inv(deformation_gradient[iID, :, :])
        rate_of_deformation = 0.5 .* (transpose(velocity_gradient) + velocity_gradient)
        rot_tensor[iID, :, :] = inv(left_stretch_tensor[iID, :, :]) * deformation_gradient[iID, :, :]
        unrotated_rate_of_deformation[iID, :, :] = transpose(rot_tensor[iID, :, :]) * rate_of_deformation * rot_tensor[iID, :, :]
    end
    return unrotated_rate_of_deformation, rot_tensor
end
"""
function deformation_gradient_decomposition(nodes, dof, deformation_gradient, deformation_gradient_dot, left_stretch_tensorN, left_stretch_tensorNP1, rot_tensorN, rot_tensorNP1, unrotated_rate_of_deformation)
    # D. P. Flanagan, L. M. Taylor, Stress integration with finite rotations; https://doi.org/10.1016/0045-7825(87)90065-X
    z = zeros(Float64, dof)
    w = zeros(Float64, dof)
    omega_vector = zeros(Float64, dof)
    omega_matrix = zeros(Float64, dof, dof)
    #scalarTemp = (1.0 - bondDamage) * omega * neighborVolume;
    #
    #*(FdotFirstTerm)   += scalarTemp * velStateX * undeformedBondX;
    #*(FdotFirstTerm+1) += scalarTemp * velStateX * undeformedBondY;
    #*(FdotFirstTerm+2) += scalarTemp * velStateX * undeformedBondZ;
    #*(FdotFirstTerm+3) += scalarTemp * velStateY * undeformedBondX;
    #*(FdotFirstTerm+4) += scalarTemp * velStateY * undeformedBondY;
    #*(FdotFirstTerm+5) += scalarTemp * velStateY * undeformedBondZ;
    #*(FdotFirstTerm+6) += scalarTemp * velStateZ * undeformedBondX;
    #*(FdotFirstTerm+7) += scalarTemp * velStateZ * undeformedBondY;
    #*(FdotFirstTerm+8) += scalarTemp * velStateZ * undeformedBondZ;
    for iID in nodes
        # Follow the algorithm at page 315
        # Compute rate-of-deformation tensor, D = 1/2 * (L + Lt) -> from Eq (4)
        velocity_gradient = deformation_gradient_dot[iID, :, :] * inv(deformation_gradient[iID, :, :])
        rate_of_deformation = 0.5 .* (transpose(velocity_gradient) + velocity_gradient)
        # Compute spin tensor, W = 1/2 * (L - Lt) -> from Eq (4)
        spin = 0.5 .* (transpose(velocity_gradient) - velocity_gradient)
        @tensor begin # Find the vector w_i = -1/2 * \epsilon_{ijk} * W_{jk} (T&F Eq. 11)
            z[i] = levicivita([i, j, k]) * rate_of_deformation[iID, j, m] * left_stretch_tensorN[iID, m, k]
            w[i] = -0.5 * levicivita([i, j, k]) * spin[iID, j, k]
        end
        #Find omega vector, i.e. omega = w +  (trace(V) I - V)^(-1) * z (T&F Eq. 12)
        omega_vector = w - (sum(diag(left_stretch_tensorN[iID, :, :])) * Matrix{Float64}(I, dof, dof) - left_stretch_tensorN[iID, :, :]) * z
        @tensor begin # Find the vector w_i = -1/2 * epsilon_{ijk} * W_{jk} (T&F Eq. 11)
            omega_matrix[i, j] = levicivita([i, k, j]) * omega_vector[k]
        end
        rot_tensorNP1[iID, :, :] = inv((Matrix{Float64}(I, dof, dof) - 0.5 * dt .* omega_matrix)) * (Matrix{Float64}(I, dof, dof) + 0.5 * dt .* omega_matrix) * rot_tensorN[iID, :, :]

        left_stretch_tensorNP1[iID, :, :] = left_stretch_tensorN[iID, :, :] + dt .* ((rate_of_deformation + spin) * left_stretch_tensorN[iID, :, :] - left_stretch_tensorN[iID, :, :] * omega_matrix)
        #L*V-V*Omegamatrix -> rate of stretch
        #vi+1 = vi+dt rateofstretch
        #unroted rate of def = RNP1^T*rateofdef*RNP1
        unrotated_rate_of_deformation[iID, :, :] = rot_tensorNP1[iID, :, :] * rate_of_deformation * transpose(rot_tensorNP1[iID, :, :])
    end

end
"""

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
    return Angle2d(angles[1] / 180 * pi)
end

end