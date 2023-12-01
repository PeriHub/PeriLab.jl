# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Lagrange_element



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

"""
    get_recursive_lagrange_shape_functions(xi::Vector{Float64}, element_coordinate::Union{Float64,Int64}, p::Int64)

Generate the values of a recursive Lagrange shape function at a value.

# Arguments
- `xi`: Vector of coordinates for the interpolation points.
- `value`: Value for which shape functions are computed. Between xi[1] and xi[end].
- `p`: Degree of the Lagrange polynomial.

# Returns
An array of shape functions.

# Example
```julia
xi = [0.0, 1.0, 2.0]
element_coordinate = 1.5
p = 2
N = get_recursive_lagrange_shape_functions(xi, value, p)
"""

function get_recursive_lagrange_shape_functions(xi::Vector{Float64}, value::Union{Float64,Int64}, p::Int64)
    N::Vector{Float64} = zeros(Float64, p + 1)
    for ip::Int64 in 1:p+1
        N[ip] = 1
        for jp::Int64 in 1:p+1
            if ip == jp
                continue
            end
            N[ip] *= (value - xi[jp]) / (xi[ip] - xi[jp])
        end
    end
    return N
end

function get_recursive_lagrange_shape_functions_derivative(xi::Vector{Float64}, value::Union{Float64,Int64}, p::Int64)
    # https://en.wikipedia.org/wiki/Lagrange_polynomial#Derivation[6]
    # sympy calculated
    B::Vector{Float64} = zeros(Float64, p + 1)
    for ip in 1:p+1
        for jp in 1:p+1
            if ip == jp
                continue
            end
            temp = 1 / (xi[ip] - xi[jp])
            for mp in 1:p+1
                if mp != ip && mp != jp
                    temp *= (value - xi[mp]) / (xi[ip] - xi[mp])
                end
            end
            B[ip] += temp
        end
    end
    return B
end
function get_2D_matrices(p::Vector{Int64}, ip_weights::Matrix{Float64}, ip_coordinates::Matrix{Float64}, xi::Matrix{Float64}, N::Matrix{Float64}, B::Matrix{Float64})
    Nxi::Vector{Float64} = zeros(p[1] + 1)
    Neta::Vector{Float64} = zeros(p[2] + 1)

    Bxi::Vector{Float64} = zeros(p[1] + 1)
    Beta::Vector{Float64} = zeros(p[2] + 1)

    for (point_id, (ip_coordinate, ip_weight)) in enumerate(zip(eachrow(ip_coordinates), eachrow(ip_weights)))

        for jID in 1:p[2]+1
            Neta = get_recursive_lagrange_shape_functions(xi[2, :], ip_coordinate[2], p[2])[jID] * ip_weight[2]
            Beta = get_recursive_lagrange_shape_functions_derivative(xi[2, :], ip_coordinate[2], p[2])[jID] * ip_weights[2]
            for iID in 1:p[1]+1
                offset = (jID - 1) * (p[1] + 1)
                offset_B = ((jID - 1) * (p[1] + 1)) * dof
                Nxi = get_recursive_lagrange_shape_functions(xi[1, :], ip_coordinate[1], p[1])[iID] * ip_weight[1]
                Bxi = get_recursive_lagrange_shape_functions_derivative(xi[1, :], ip_coordinate[1], p[1])[iID] * ip_weight[1]
                N[point_id, iID+offset] = Nxi[iID] * Neta[jID]
                B[point_id, iID+offset_B] = Bxi[iID] * Neta[jID]

            end
        end
    end
    return N, B
end

function get_3D_matrices(p::Vector{Int64}, ip_weights::Matrix{Float64}, ip_coordinates::Matrix{Float64}, xi::Matrix{Float64}, N::Matrix{Float64}, B::Matrix{Float64})
    Nxi::Vector{Float64} = zeros(p[1] + 1)
    Neta::Vector{Float64} = zeros(p[2] + 1)
    Npsi::Vector{Float64} = zeros(p[3] + 1)
    Bxi::Vector{Float64} = zeros(p[1] + 1)
    Beta::Vector{Float64} = zeros(p[2] + 1)
    Bpsi::Vector{Float64} = zeros(p[3] + 1)
    for (point_id, (ip_coordinate, ip_weights)) in enumerate(zip(eachrow(ip_coordinates), eachrow(ip_weights)))
        for kID in 1:p[3]+1
            Npsi = get_recursive_lagrange_shape_functions(xi[3, :], ip_coordinate[3, 3], p[3])[kID] * ip_weight[3]
            Bpsi = get_recursive_lagrange_shape_functions_derivative(xi[3, :], ip_coordinate[3, 3], p[3])[kID] * ip_weights[3]
            for jID in 1:p[2]+1
                Neta = get_recursive_lagrange_shape_functions(xi[2, :], ip_coordinate[2, 2], p[2])[jID] * ip_weight[2]
                Beta = get_recursive_lagrange_shape_functions_derivative(xi[2, :], ip_coordinate[2, 2], p[2])[jID] * ip_weights[2]
                for iID in 1:p[1]+1
                    offset = (jID - 1) * (p[1] + 1) + (kID - 1) * (p[2] + 1) * (p[1] + 1)
                    offset_B = ((jID - 1) * (p[1] + 1) + (kID - 1) * (p[2] + 1) * (p[1] + 1)) * dof
                    Nxi = get_recursive_lagrange_shape_functions(xi[1, :], ip_coordinate[1, 1], p[1])[iID] * ip_weights[1]
                    Bxi = get_recursive_lagrange_shape_functions_derivative(xi[1, :], ip_coordinate[1, 1], p[1])[iID] * ip_weights[1]
                    N[point_id, iID+offset] = Nxi[iID] * Neta[jID] * Npsi[kID]
                    B[point_id, iID+offset_B] = Bxi[iID] * Neta[jID] * Npsi[kID]
                    B[point_id, iID+offset_B+1] = Nxi[iID] * Beta[jID] * Npsi[kID]
                    B[point_id, iID+offset_B+2] = Nxi[iID] * Neta[jID] * Bpsi[kID]
                end
            end
        end
    end
    return N, B
end

"""
    create_element_matrices(dof::Int64, p::Vector{Int64}, ip_weights::Matrix{Float64}, ip_coordinates::Matrix{Float64})
    
    creates the element matrices for each integration point. This is needed, because if 
    the stresses are not constant at all integration points or the coordinates are rotated these 
    fields are needed. However, for each element type this has to be done ones.       

    # Arguments

    Return
    N, B
    Size of N and B is for each integration point number of nodes times dof for N and for Beta
    number of nodes times dof + shear part (1 for 2D and 3 for 3D)
    For other elements the number of nodes are not equal to the number of integration points in
    a fully integrated element

"""
function create_element_matrices(dof::Int64, p::Vector{Int64}, ip_weights::Matrix{Float64}, ip_coordinates::Matrix{Float64})
    if dof > 3 || dof < 2
        @error "Not support degree of freedom for the finite element matrix creation"
        return nothing
    end

    N::Matrix{Float64} = zeros(Float64, prod(p .+ 1), prod(p .+ 1))
    B::Matrix{Float64} = zeros(Float64, prod(p .+ 1), prod(p .+ 1) * dof)
    xi = define_lagarangian_grid_space(dof, p)

    if dof == 2
        return get_2D_matrices(p, ip_weights, ip_coordinates, xi, N, B)
    end
    return get_3D_matrices(p, ip_weights, ip_coordinates, xi, N, B)
end
end