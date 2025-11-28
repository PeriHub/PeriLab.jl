# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

using LinearAlgebra

"""
Structure to store material properties for elastic materials
"""
struct MaterialProperties
    E::Float64      # Young's modulus [Pa]
    ν::Float64      # Poisson's ratio [-]
    G::Float64      # Shear modulus [Pa]
    K_bulk::Float64 # Bulk modulus [Pa]
end

"""
Compute the shape tensor for a given node in peridynamics
This function calculates the shape tensor used in bond-based peridynamics
for relating deformation to bond forces.

Arguments:
- shape_tensor: 3D array to store computed shape tensors [node, i, j]
- volume: Vector of bond volumes
- bond_damage: Vector of bond damage parameters (0=broken, 1=intact)
- undeformed_bond: Vector of undeformed bond vectors
- iID: Current node ID
"""
function compute_shape_tensor!(shape_tensor::Array{Float64,3},
                               volume::AbstractVector{Float64},
                               bond_damage::Vector{Float64},
                               undeformed_bond::Vector{Vector{Float64}},
                               iID::Int64)
    @inbounds @fastmath for m in axes(shape_tensor, 2), n in axes(shape_tensor, 3)
        Cmn = zero(eltype(shape_tensor))
        @inbounds @fastmath for k in axes(undeformed_bond, 1)
            # Sum contributions from all bonds: C_mn = Σ ξ_m * ξ_n * V * μ
            Cmn += undeformed_bond[k][m] *
                   undeformed_bond[k][n] *
                   volume[k] *
                   bond_damage[k]
        end
        shape_tensor[iID, m, n] = Cmn
    end
end

"""
Compute the inverse of the shape tensor for a given node
Uses LinearAlgebra.inv() to compute matrix inverse

Arguments:
- shape_tensor_inv: 3D array to store inverse shape tensors
- shape_tensor: 3D array containing computed shape tensors
- iID: Current node ID
"""
function compute_shape_tensor_inverse!(shape_tensor_inv::Array{Float64,3},
                                       shape_tensor::Array{Float64,3},
                                       iID::Int64)
    shape_tensor_inv[iID, :, :] = inv(shape_tensor[iID, :, :])
end

"""
Constructor for MaterialProperties from dictionary with shear and bulk moduli
Converts from G and K to full set of elastic constants

Arguments:
- material: Dictionary containing "Shear Modulus" and "Bulk Modulus"

Returns:
- MaterialProperties struct with E, ν, G, K_bulk
"""
function MaterialProperties(material::Dict{String,Float64})
    G = material["Shear Modulus"]
    K = material["Bulk Modulus"]

    # Convert to Young's modulus: E = 9KG/(3K + G)
    E = 9 * K * G / (3 * K + G)

    # Convert to Poisson's ratio: ν = (3K - 2G)/(2(3K + G))
    ν = (3 * K - 2 * G) / (2 * (3 * K + G))

    return MaterialProperties(E, ν, G, K)
end

"""
Compute 2D plane stress elasticity matrix in Voigt notation
Returns the 3x3 elasticity matrix for plane stress conditions

Arguments:
- material: Dictionary containing "Shear Modulus" and "Bulk Modulus"

Returns:
- C: 3x3 elasticity matrix in Voigt notation [C11 C12 0; C12 C22 0; 0 0 C66]
"""
function elasticity_matrix_2d_plane_stress(material)
    G = material["Shear Modulus"]
    K = material["Bulk Modulus"]

    # Convert to engineering constants
    E = 9 * K * G / (3 * K + G)
    ν = (3 * K - 2 * G) / (2 * (3 * K + G))

    println("Young's Modulus E = $E Pa")
    println("Poisson's Ratio ν = $ν")

    # Plane stress elasticity matrix
    factor = E / (1 - ν^2)
    C = [factor factor*ν 0.0
         factor*ν factor 0.0
         0.0 0.0 factor*(1-ν)/2]

    return C
end

"""
Convert 2D Voigt elasticity matrix to 4th order tensor representation
Transforms from Voigt notation (3x3) to full tensor notation (2x2x2x2)

In 2D Voigt notation:
- Index 1: xx component
- Index 2: yy component
- Index 3: xy component (engineering shear strain γ_xy = 2ε_xy)

Arguments:
- C_voigt: 3x3 elasticity matrix in Voigt notation

Returns:
- T4: 4th order elasticity tensor (2x2x2x2)
"""
function voigt_to_tensor4_2d(C_voigt::AbstractMatrix{Float64})
    T4 = zeros(2, 2, 2, 2)

    # Standard 2D Voigt notation mapping:
    # Voigt[1,1] = C_1111  (σ_11 ← ε_11)
    # Voigt[2,2] = C_2222  (σ_22 ← ε_22)
    # Voigt[3,3] = C_1212  (σ_12 ← γ_12 = 2ε_12)
    # Voigt[1,2] = C_1122  (σ_11 ← ε_22)
    # Voigt[2,1] = C_2211  (σ_22 ← ε_11)

    # Main diagonal terms (normal-normal couplings):
    T4[1, 1, 1, 1] = C_voigt[1, 1]  # C_1111
    T4[2, 2, 2, 2] = C_voigt[2, 2]  # C_2222

    # Normal-normal cross terms:
    T4[1, 1, 2, 2] = C_voigt[1, 2]  # C_1122
    T4[2, 2, 1, 1] = C_voigt[2, 1]  # C_2211

    # Pure shear terms:
    T4[1, 2, 1, 2] = C_voigt[3, 3]  # C_1212
    T4[1, 2, 2, 1] = C_voigt[3, 3]  # C_1221 (symmetry)
    T4[2, 1, 1, 2] = C_voigt[3, 3]  # C_2112 (symmetry)
    T4[2, 1, 2, 1] = C_voigt[3, 3]  # C_2121 (symmetry)

    # Normal-shear coupling terms:
    T4[1, 1, 1, 2] = C_voigt[1, 3]  # C_1112 (σ_11 ← γ_12)
    T4[1, 1, 2, 1] = C_voigt[1, 3]  # C_1121 (symmetry)
    T4[2, 2, 1, 2] = C_voigt[2, 3]  # C_2212 (σ_22 ← γ_12)
    T4[2, 2, 2, 1] = C_voigt[2, 3]  # C_2221 (symmetry)

    # Shear-normal coupling terms:
    T4[1, 2, 1, 1] = C_voigt[3, 1]  # C_1211 (σ_12 ← ε_11)
    T4[2, 1, 1, 1] = C_voigt[3, 1]  # C_2111 (symmetry)
    T4[1, 2, 2, 2] = C_voigt[3, 2]  # C_1222 (σ_12 ← ε_22)
    T4[2, 1, 2, 2] = C_voigt[3, 2]  # C_2122 (symmetry)

    return T4
end
