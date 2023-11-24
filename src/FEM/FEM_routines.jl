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
    return w, x
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
    if minimum(p) == 0
        @error "p order for lagarangian grid space must be at least p = 1 and not zero"
        return nothing
    end
    for idof in 1:dof
        len = 2.0 / p[idof]
        for ip in 1:p[idof]+1
            xi[idof, ip] = -1 + (ip - 1) * len
        end
    end
    return xi
end


function get_recursive_lagrange_shape_functions(xi::Vector{Float64}, element_coordinate::Union{Float64,Int64}, p::Int64)
    N::Vector{Float64} = zeros(Float64, p + 1)
    for ip::Int64 in 1:p+1
        N[ip] = 1
        for jp::Int64 in 1:p+1
            if ip == jp
                continue
            end
            N[ip] *= (element_coordinate - xi[jp]) / (xi[ip] - xi[jp])
        end
    end
    return N
end
