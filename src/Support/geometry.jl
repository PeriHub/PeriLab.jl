# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Geometry
using LinearAlgebra
using TensorOperations
using Combinatorics: levicivita
using Rotations
include("helpers.jl")
using .Helpers: invert
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
function bond_geometry(
    nodes::Union{SubArray,Vector{Int64}},
    nlist::Union{SubArray,Vector{Int64}},
    coor::Union{SubArray,Matrix{Float64},Matrix{Int64}},
    undeformed_bond,
    undeformed_bond_length,
)

    for iID in nodes
        undeformed_bond[iID], undeformed_bond_length[iID] =
            calculate_bond_length(iID, coor, nlist[iID])
    end
    if any(any.(x -> x == 0, [lengths for lengths in undeformed_bond_length[nodes]]))
        @error "Identical point coordinates with no distance"
        return nothing
    end
    return undeformed_bond, undeformed_bond_length
end

function calculate_bond_length(
    iID::Int64,
    coor::Union{SubArray,Matrix{Float64},Matrix{Int64}},
    nlist::Vector{Int64},
)

    @views bond_vectors = coor[nlist, :] .- coor[iID, :]'
    # distances = sqrt.(sum(bond_vectors .^ 2, dims=2))[:]

    # Check for identical point coordinates
    return bond_vectors, norm.(eachrow(bond_vectors))
end

"""
    shape_tensor(nodes::Union{SubArray, Vector{Int64}}, nlist, volume, omega, bond_damage, undeformed_bond, shape_tensor, inverse_shape_tensor)

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

shape_tensor(nodes, nlist, volume, omega, bond_damage, undeformed_bond, shape_tensor, inverse_shape_tensor)
"""
function shape_tensor(
    nodes::Union{SubArray,Vector{Int64}},
    nlist,
    volume,
    omega,
    bond_damage,
    undeformed_bond,
    shape_tensor,
    inverse_shape_tensor,
)

    for iID in nodes

        shape_tensor[iID, :, :] = calculate_shape_tensor(
            volume[nlist[iID]],
            omega[iID],
            bond_damage[iID],
            undeformed_bond[iID],
        )
        inverse_shape_tensor[iID, :, :] .= invert(
            shape_tensor[iID, :, :],
            "Shape Tensor is singular and cannot be inverted).\n - Check if your mesh is 3D, but has only one layer of nodes\n - Check number of damaged bonds.",
        )
    end

    return shape_tensor, inverse_shape_tensor
end


function calculate_shape_tensor(volume, omega, bond_damage, undeformed_bond)

    return (bond_damage .* volume .* omega .* undeformed_bond)' * undeformed_bond

end


function bond_associated_deformation_gradient(
    dof::Int64,
    volume,
    omega,
    bond_damage,
    undeformed_bond,
    deformed_bond,
)

    return calculate_deformation_gradient(
        bond_damage,
        deformed_bond,
        undeformed_bond,
        volume,
        omega,
    )

end

"""
    compute_deformation_gradient(nodes::Union{SubArray, Vector{Int64}}, nlist, volume, omega, bond_damage, undeformed_bond, deformed_bond, inverse_shape_tensor, deformation_gradient)

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

compute_deformation_gradient(nodes, nlist, volume, omega, bond_damage, undeformed_bond, deformed_bond, inverse_shape_tensor, deformation_gradient)
"""
function compute_deformation_gradient(
    nodes::Union{SubArray,Vector{Int64}},
    nlist::SubArray,
    volume::SubArray,
    omega::SubArray,
    bond_damage::SubArray,
    deformed_bond::Union{SubArray,Vector{Matrix{Float64}}},
    undeformed_bond::SubArray,
    inverse_shape_tensor::SubArray,
    deformation_gradient::SubArray,
)
    for iID in nodes
        deformation_gradient[iID, :, :] = calculate_deformation_gradient(
            bond_damage[iID],
            deformed_bond[iID],
            undeformed_bond[iID],
            volume[nlist[iID]],
            omega[iID],
        )

        @views deformation_gradient[iID, :, :] *= inverse_shape_tensor[iID, :, :]
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

function compute_left_stretch_tensor(deformation_gradient::Matrix{Float64})
    return sqrt(deformation_gradient * deformation_gradient')
end

function calculate_deformation_gradient(
    bond_damage,
    deformed_bond,
    undeformed_bond,
    volume::Union{Vector{Int64},Vector{Float64}},
    omega,
)
    return (bond_damage .* volume .* omega .* deformed_bond)' * undeformed_bond

end


function compute_weighted_deformation_gradient(
    nodes::Union{SubArray,Vector{Int64}},
    dof::Int64,
    nlist,
    volume,
    gradient_weight,
    displacement,
    deformation_gradient,
)

    for iID in nodes
        deformation_gradient[iID, :, :] = @views I(dof) +
               (displacement[nlist[iID], :]' .- displacement[iID, :]) *
               (gradient_weight[iID][:, :] .* volume[nlist[iID]])
    end
    return deformation_gradient

end


function deformation_gradient_decomposition(
    nodes::Union{Base.OneTo{Int64},Vector{Int64}},
    deformation_gradient,
    rot_tensor,
)

    for iID in nodes
        # EQ (15) in https://doi.org/10.1016/j.ijsolstr.2008.10.029
        rot_tensor[iID, :, :] =
            invert(
                compute_left_stretch_tensor(deformation_gradient[iID, :, :]),
                "Inversion of left stretch tensor fails in function deformation_gradient_decomposition.",
            ) * deformation_gradient[iID, :, :]
    end
    return rot_tensor
end



"""
    function compute_strain(nodes::Union{Base.OneTo{Int64},Vector{Int64}, SubArray}, deformation_gradient, strain)

Calculate strains for specified nodes based on deformation gradients.

## Arguments
-  `nodes::Union{SubArray, Vector{Int64}}`: List of nodes
- `deformation_gradient`: Deformation gradient at current time step (2D or 3D array).

## Returns

- Updated `strain` array containing strains.

This function iterates over the specified nodes and computes strain at each node using the given deformation gradients.

"""
function compute_strain(
    nodes::Union{Base.OneTo{Int64},Vector{Int64},SubArray},
    deformation_gradient::Union{Array{Float64,3},SubArray},
    strain::Union{Array{Float64,3},SubArray},
)

    # https://en.wikipedia.org/wiki/Strain_(mechanics)
    for iID in nodes
        strain[iID, :, :] = @views 0.5 .* (
            transpose(deformation_gradient[iID, :, :]) * deformation_gradient[iID, :, :] - I
        )
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


function compute_bond_level_rotation_tensor(
    nodes,
    nlist,
    ba_deformation_gradient,
    ba_rotation_tensor,
)
    # all deformation gradients, etc. are in NP1. The increment is calculated outside this routine.
    for iID in nodes
        ba_rotation_tensor[iID][:, :, :] = deformation_gradient_decomposition(
            eachindex(nlist[iID]),
            ba_deformation_gradient[iID][:, :, :],
            ba_rotation_tensor[iID][:, :, :],
        )
    end
    return ba_rotation_tensor
end

function compute_bond_level_deformation_gradient(
    nodes,
    nlist,
    dof,
    bond_geometry,
    bond_length,
    bond_deformation,
    deformation_gradient,
    ba_deformation_gradient,
)
    for iID in nodes
        for (jID, nID) in enumerate(nlist[iID])
            mean_deformation_gradient =
                0.5 .* (deformation_gradient[iID, :, :] + deformation_gradient[nID, :, :])
            for i = 1:dof
                scalar_temp::Float64 =
                    mean_deformation_gradient[i, :]' * bond_geometry[iID][jID, :]
                ba_deformation_gradient[iID][jID, i, :] =
                    mean_deformation_gradient[i, :] +
                    (bond_deformation[iID][jID, i] - scalar_temp) .*
                    bond_geometry[iID][jID, :] /
                    (bond_length[iID][jID] * bond_length[iID][jID])

                # scalarTemp = *(meanDefGrad+0) * undeformedBondX + *(meanDefGrad+1) * undeformedBondY + *(meanDefGrad+2) * undeformedBondZ;

                # *(defGrad+0) = *(meanDefGrad+0) + (defStateX - scalarTemp) * undeformedBondX/undeformedBondLengthSq;
                # *(defGrad+1) = *(meanDefGrad+1) + (defStateX - scalarTemp) * undeformedBondY/undeformedBondLengthSq;
                # *(defGrad+2) = *(meanDefGrad+2) + (defStateX - scalarTemp) * undeformedBondZ/undeformedBondLengthSq;
            end
        end
    end
    return ba_deformation_gradient
end

function compute_bond_level_unrotated_rate_of_deformation_and_rotation_tensor_from_Peridigm()
    mean_deformation_gradient = zeros(Float64, dof, dof)
    mean_deformation_gradient_rate = zeros(Float64, dof, dof)
    velocity_state = zeros(Float64, dof)
    omega = zeros(Float64, dof)
    for iID in nodes
        for (jID, nID) in enumerate(nlist[iID])
            velocity_state = velocity[nID, :] - velocity[iID, :]
            mean_deformation_gradient =
                0.5 .* (deformation_gradient[iID, :, :] + deformation_gradient[nID, :, :])
            mean_deformation_gradient_rate =
                0.5 .*
                (deformation_gradient_dot[iID, :, :] + deformation_gradient_dot[nID, :, :])



            for i = 1:dof
                scalar_temp =
                    sum(mean_deformation_gradient[i, :] .* bond_geometry[iID][jID, i])
                scalar_temp_dot =
                    sum(mean_deformation_gradient_rate[i, :] .* bond_geometry[iID][jID, i])
                deformation_gradient[iID, i, :] =
                    mean_deformation_gradient[i, :] .+
                    (bond_deformation[iID][jID, i] - scalar_temp) .*
                    bond_geometry[iID][jID, i] / bond_length[iID][jID]
                deformation_gradient_dot[iID, i, :] =
                    mean_deformation_gradient_rate[i, :] .+
                    (velocity_state[i] - scalar_temp) .* bond_geometry[iID][jID, i] /
                    bond_length[iID][jID]
            end
            jacobian[iID] = determinant(deformation_gradient[iID, :, :])

            # Compute the Eulerian velocity gradient L = Fdot * Finv
            eulerian_vel_grad =
                deformation_gradient_dot[iID, :, :] * invert(
                    deformation_gradient[iID, :, :],
                    "Deformation gradient is singular in compute_bond_level_unrotated_rate_of_deformation_and_rotation_tensor",
                )

            # Compute rate-of-deformation tensor, D = 1/2 * (L + Lt)
            rate_of_deformation = 0.5 .* (eulerian_vel_grad + eulerian_vel_grad')
            # Compute spin tensor, W = 1/2 * (L - Lt)
            spin = 0.5 .* (eulerian_vel_grad - eulerian_vel_grad')
            ###########################################################################
            # Following Flanagan & Taylor (T&F)
            #
            # Find the vector z_i = epsilon_{ikj} * D_{jm} * V_{mk} (T&F Eq. 13)
            # and find w_i = -1/2 * epsilon_{ijk} * W_{jk} (T&F Eq. 11)
            # where epsilon_{ikj} is the alternator tensor.
            #
            # Components below copied from computer algebra solution to the expansion
            # above
            ###########################################################################
            z .= 0
            w .= 0
            for i = 1:dof
                for k = 1:dof
                    for k = 1:dof
                        z[i] +=
                            levicivita([i, k, j]) .* rate_of_deformation[iID, j, m] *
                            left_stretch_tensorN[iID, m, k]
                        w[i] += -0.5 * levicivita([i, j, k]) .* spin[iID, j, k]
                    end
                end
            end

            # Find trace(V)

            traceV = (sum(diag(left_stretch_tensorN[iID, :, :])))
            # Compute (trace(V) * I - V) store in temp
            temp = trace .* I(dof) - left_stretch_tensorN[iID, :, :]
            # Find omega vector, i.e. \omega = w +  (trace(V) I - V)^(-1) * z (T&F Eq. 12)
            omega = w + invert(temp) * z
            omega_tensor .= 0
            for i = 1:dof
                for k = 1:dof
                    for k = 1:dof
                        omega_tensor[i, j] += levicivita([i, k, j]) * omega[k]
                    end
                end
            end
            ## Increment R with (T&F Eq. 36 and 44) as opposed to solving (T&F 39) this
            ## is desirable for accuracy in implicit solves and has no effect on
            ## explicit solves (other than a slight decrease in speed).
            ##
            ##  Compute Q with (T&F Eq. 44)
            ##
            ##  Omega^2 = w_i * w_i (T&F Eq. 42)
            OmegaSq = omegaX * omegaX + omegaY * omegaY + omegaZ * omegaZ
            # Omega = sqrt{OmegaSq}
            Omega = sqrt(OmegaSq)

            # Avoid a potential divide-by-zero
            """
            if(OmegaSq > 1.e-30){

              // Compute Q = I + sin( dt * Omega ) * OmegaTensor / Omega - (1. - cos(dt * Omega)) * omegaTensor^2 / OmegaSq
              //           = I + scaleFactor1 * OmegaTensor + scaleFactor2 * OmegaTensorSq
              scaleFactor1 = sin(dt*Omega) / Omega;
              scaleFactor2 = -(1.0 - cos(dt*Omega)) / OmegaSq;
              MatrixMultiply(false, false, 1.0, OmegaTensor, OmegaTensor, OmegaTensorSq);
              *(QMatrix)   = 1.0 + scaleFactor1 * *(OmegaTensor)   + scaleFactor2 * *(OmegaTensorSq)   ;
              *(QMatrix+1) =       scaleFactor1 * *(OmegaTensor+1) + scaleFactor2 * *(OmegaTensorSq+1) ;
              *(QMatrix+2) =       scaleFactor1 * *(OmegaTensor+2) + scaleFactor2 * *(OmegaTensorSq+2) ;
              *(QMatrix+3) =       scaleFactor1 * *(OmegaTensor+3) + scaleFactor2 * *(OmegaTensorSq+3) ;
              *(QMatrix+4) = 1.0 + scaleFactor1 * *(OmegaTensor+4) + scaleFactor2 * *(OmegaTensorSq+4) ;
              *(QMatrix+5) =       scaleFactor1 * *(OmegaTensor+5) + scaleFactor2 * *(OmegaTensorSq+5) ;
              *(QMatrix+6) =       scaleFactor1 * *(OmegaTensor+6) + scaleFactor2 * *(OmegaTensorSq+6) ;
              *(QMatrix+7) =       scaleFactor1 * *(OmegaTensor+7) + scaleFactor2 * *(OmegaTensorSq+7) ;
              *(QMatrix+8) = 1.0 + scaleFactor1 * *(OmegaTensor+8) + scaleFactor2 * *(OmegaTensorSq+8) ;

            } else {
              *(QMatrix)   = 1.0 ; *(QMatrix+1) = 0.0 ; *(QMatrix+2) = 0.0 ;
              *(QMatrix+3) = 0.0 ; *(QMatrix+4) = 1.0 ; *(QMatrix+5) = 0.0 ;
              *(QMatrix+6) = 0.0 ; *(QMatrix+7) = 0.0 ; *(QMatrix+8) = 1.0 ;
            };

            // Compute R_STEP_NP1 = QMatrix * R_STEP_N (T&F Eq. 36)
            MatrixMultiply(false, false, 1.0, QMatrix, rotTensorN, rotTensorNP1);

            // Compute rate of stretch, Vdot = L*V - V*Omega
            // First tempA = L*V,
            MatrixMultiply(false, false, 1.0, eulerianVelGrad, leftStretchN, tempA);

            // tempB = V*Omega
            MatrixMultiply(false, false, 1.0, leftStretchN, OmegaTensor, tempB);

            //Vdot = tempA - tempB
            for(int i=0 ; i<9 ; ++i)
              *(rateOfStretch+i) = *(tempA+i) - *(tempB+i);

            //V_STEP_NP1 = V_STEP_N + dt*Vdot
            for(int i=0 ; i<9 ; ++i)
              *(leftStretchNP1+i) = *(leftStretchN+i) + dt * *(rateOfStretch+i);

            // Compute the unrotated rate-of-deformation, d, i.e., temp = D * R
            MatrixMultiply(false, false, 1.0, rateOfDef, rotTensorNP1, temp);

            // d = Rt * temp
            MatrixMultiply(true, false, 1.0, rotTensorNP1, temp, unrotRateOfDef);

            // Store back in element-wise format
            *leftStretchXXNP1 = *(leftStretchNP1+0); *leftStretchXYNP1 = *(leftStretchNP1+1); *leftStretchXZNP1 = *(leftStretchNP1+2);
            *leftStretchYXNP1 = *(leftStretchNP1+3); *leftStretchYYNP1 = *(leftStretchNP1+4); *leftStretchYZNP1 = *(leftStretchNP1+5);
            *leftStretchZXNP1 = *(leftStretchNP1+6); *leftStretchZYNP1 = *(leftStretchNP1+7); *leftStretchZZNP1 = *(leftStretchNP1+8);
            *rotTensorXXNP1 = *(rotTensorNP1+0); *rotTensorXYNP1 = *(rotTensorNP1+1); *rotTensorXZNP1 = *(rotTensorNP1+2);
            *rotTensorYXNP1 = *(rotTensorNP1+3); *rotTensorYYNP1 = *(rotTensorNP1+4); *rotTensorYZNP1 = *(rotTensorNP1+5);
            *rotTensorZXNP1 = *(rotTensorNP1+6); *rotTensorZYNP1 = *(rotTensorNP1+7); *rotTensorZZNP1 = *(rotTensorNP1+8);
            *unrotRateOfDefXX = *(unrotRateOfDef+0); *unrotRateOfDefXY = *(unrotRateOfDef+1); *unrotRateOfDefXZ = *(unrotRateOfDef+2);
            *unrotRateOfDefYX = *(unrotRateOfDef+3); *unrotRateOfDefYY = *(unrotRateOfDef+4); *unrotRateOfDefYZ = *(unrotRateOfDef+5);
            *unrotRateOfDefZX = *(unrotRateOfDef+6); *unrotRateOfDefZY = *(unrotRateOfDef+7); *unrotRateOfDefZZ = *(unrotRateOfDef+8);
          }
          else{
            *leftStretchXXNP1 = *leftStretchXXN; *leftStretchXYNP1 = *leftStretchXYN; *leftStretchXZNP1 = *leftStretchXZN;
            *leftStretchYXNP1 = *leftStretchYXN; *leftStretchYYNP1 = *leftStretchYYN; *leftStretchYZNP1 = *leftStretchYZN;
            *leftStretchZXNP1 = *leftStretchZXN; *leftStretchZYNP1 = *leftStretchZYN; *leftStretchZZNP1 = *leftStretchZZN;
            *rotTensorXXNP1 = *rotTensorXXN; *rotTensorXYNP1 = *rotTensorXYN; *rotTensorXZNP1 = *rotTensorXZN;
            *rotTensorYXNP1 = *rotTensorYXN; *rotTensorYYNP1 = *rotTensorYYN; *rotTensorYZNP1 = *rotTensorYZN;
            *rotTensorZXNP1 = *rotTensorZXN; *rotTensorZYNP1 = *rotTensorZYN; *rotTensorZZNP1 = *rotTensorZZN;
            *unrotRateOfDefXX = 0.0; *unrotRateOfDefXY = 0.0; *unrotRateOfDefXZ = 0.0;
            *unrotRateOfDefYX = 0.0; *unrotRateOfDefYY = 0.0; *unrotRateOfDefYZ = 0.0;
            *unrotRateOfDefZX = 0.0; *unrotRateOfDefZY = 0.0; *unrotRateOfDefZZ = 0.0;
          }
        }
      }
      return returnCode;
    }


    TODO calculate def_Grad and synch def_Grad
    """
        end
    end
end
end
