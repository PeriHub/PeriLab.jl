# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
using LinearAlgebra
using FastGaussQuadrature
include("Element_formulation/lagrange_element.jl")
using .Lagrange_element
function compute_strain(nodes::Union{SubArray,Vector{Int64}}, topology, B, u, strain)

    for iID in nodes
        strain[iID, :] = dot(B[iID, :, :], u[topology[iID, :]])

    end
    return strain
end
"""
    get_weights_and_integration_points(dof::Int64, p::Vector{Int64})

Compute integration points and weights using Gauss-Legendre quadrature.

# Arguments
- `dof::Int64`: The number of degrees of freedom.
- `p::Vector{Int64}`: A vector containing the polynomial degrees for each degree of freedom.

# Returns
- `x::Matrix{Float64}`: Matrix of integration points. Each row corresponds to a degree of freedom, and columns contain the integration points.
- `w::Matrix{Float64}`: Matrix of integration weights. Each row corresponds to a degree of freedom, and columns contain the integration weights.

# Example
```julia
dof = 3
p = [2, 3, 1]
x, w = get_weights_and_integration_points(dof, p)
"""
function get_weights_and_integration_points(dof::Int64, p::Vector{Int64})

    x::Matrix{Float64} = zeros(Float64, maximum(p) + 1, dof)
    w::Matrix{Float64} = zeros(Float64, maximum(p) + 1, dof)
    for idof in 1:dof
        x[1:p[idof]+1, idof], w[1:p[idof]+1, idof] = gausslegendre(p[idof] + 1)
    end
    return w, x
end
"""
function set_weights(w::Matrix{Float64})
    Determines the weights at the integration points
"""
function set_weights(w::Matrix{Float64})

    return weigths
end


function get_Jacobian(B::SubArray, coordinates::SubArray)
    return det(dot(B, coordinates))
end

"""
    create_element_matrices()

Compute the matrix of shape functions (N) and its derivative (B) for 2D and 3D.
N^TN*\rho give than the mass matrix and B^TCB the stiffness matrix [WillbergC2013](@cite)

# Arguments
- `dof::Int64`: The number of degrees of freedom.
- `p::Vector{Int64}`: A vector containing the polynomial degrees for each degree of freedom.

# Returns
- N = [N_x 0   0   ...
       0   N_y 0   ...
       0   0   N_z ...]

- B = [B_x 0   0   ...
       0   B_y 0   ...
       0   0   B_z ...]
       B_y B_x 0   ...]
       B_z 0 B_x   ...]
       0   B_z B_y ...]
# Example

"""


function create_element_matrices(dof::Int64, p::Vector{Int64})
    weights, integration_points = get_weights_and_integration_points(dof, p)
    return Lagrange_element.create_element_matrices(dof, p, weights, integration_points)
end