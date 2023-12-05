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
- `w::Matrix{Float64}`: Matrix of integration weights. Each row corresponds to a degree of freedom, and columns contain the integration weights.
- `x::Matrix{Float64}`: Matrix of integration points. Each row corresponds to a degree of freedom, and columns contain the integration points.

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
    get_multi_dimensional_integration_point_data(dof::Int64, p::Vector{Int64}, value::Matrix{Float64})

Restructure integration point information for multi-dimensional problems.

# Arguments
- `dof::Int64`: Degree of freedom, only 2 and 3 are supported.
- `p::Vector{Int64}`: Vector containing the number of integration points in each dimension.
- `value::Matrix{Float64}`: Matrix containing the coordinates of the reference element.

# Returns
- `integration_point::Matrix{Float64}`: Matrix of integration point data.

# Example
```julia
dof = 2
p = [2, 3]
value = rand(3, dof)
result = get_multi_dimensional_integration_points(dof, p, value)
"""
function get_multi_dimensional_integration_point_data(dof::Int64, p::Vector{Int64}, value::Matrix{Float64})
    count::Int64 = 0
    integration_point::Matrix{Float64} = zeros(prod(p .+ 1), dof)
    if dof == 2
        for pid in 1:p[1]+1
            for pjd in 1:p[2]+1
                count += 1
                integration_point[count, 1] = value[pid, 1]
                integration_point[count, 2] = value[pjd, 2]
            end
        end
    elseif dof == 3
        for pid in 1:p[1]+1
            for pjd in 1:p[2]+1
                for pkd in 1:p[3]+1
                    count += 1
                    integration_point[count, 1] = value[pid, 1]
                    integration_point[count, 2] = value[pjd, 2]
                    integration_point[count, 3] = value[pkd, 3]
                end
            end
        end
    else
        @error "degree of freedom = $dof is not supported, only 2 and 3."
        return nothing
    end
    return integration_point
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

function create_element_matrices(dof::Int64, p::Vector{Int64}, create_matrices)
    weights, integration_points = get_weights_and_integration_points(dof, p)
    ip_coordinates = get_multi_dimensional_integration_point_data(dof, p, integration_points)
    if isnothing(ip_coordinates)
        return nothing, nothing
    end
    ip_weights = get_multi_dimensional_integration_point_data(dof, p, weights)
    return create_matrices(dof, p, ip_weights, ip_coordinates)
end