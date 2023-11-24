# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
using LinearAlgebra
using FastGaussQuadrature
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

    x::Matrix{Float64} = zeros(Float64, dof, maximum(p) + 1)
    w::Matrix{Float64} = zeros(Float64, dof, maximum(p) + 1)
    for idof in 1:dof
        x[idof, 1:p[idof]+1], w[idof, 1:p[idof]+1] = gausslegendre(p[idof] + 1)
    end
    return x, w
end

function get_Jacobian(B::SubArray, coordinates::SubArray)
    return det(dot(B, coordinates))
end

"""
    define_lagrangian_grid_space(dof::Int64, p::Vector{Int64})

Define a Lagrangian grid space for shape functions.

# Arguments
- `dof::Int64`: The number of degrees of freedom.
- `p::Vector{Int64}`: A vector containing the polynomial degrees for each degree of freedom.

# Returns
- `xi::Matrix{Float64}`: Matrix representing the Lagrangian grid space. Each row corresponds to a degree of freedom, and columns contain the grid points.

# Example
```julia
dof = 3
p = [2, 3, 1]
xi = define_lagrangian_grid_space(dof, p)
"""
function define_lagarangian_grid_space(dof::Int64, p::Vector{Int64})
    xi::Matrix{Float64} = zeros(Float64, dof, maximum(p) + 1)
    len::Float64 = 0.0
    for idof in 1:dof
        len = 2.0 / p[idof]
        for ip in 1:p[idof]+1
            xi[idof, ip] = -1 + (ip - 1) * len
        end
    end
    return xi
end
