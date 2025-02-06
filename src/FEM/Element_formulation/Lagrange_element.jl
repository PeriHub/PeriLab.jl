# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Lagrange_element

export element_name
export init_element
#export compute_element

function element_name()
    return "Lagrange"
end

"""
    init_element(datamanager::Module, elements::Union{SubArray,Vector{Int64}}, element_params::Dict, p::Vector{Int64})

Init the Lagrange element of a given polynomial degree. This degree can be different for each direction

# Arguments
- `datamanager::Module`: Datamanager.
- `elements::Union{SubArray,Vector{Int64}}`: listed element numbers
- `element_params::Dict`: Element specific data.
- `p::Vector{Int64}`: A vector containing the polynomial degrees for each degree of freedom.

# Returns
- `Datamanager::Module`: Datamanager containing additional fields needed for the Lagrange element.
"""

function init_element(
    datamanager::Module,
    elements::Union{SubArray,Vector{Int64}},
    element_params::Dict,
    p::Vector{Int64},
)

    return datamanager
end
"""
    define_lagrangian_grid_space(dof::Int64, p::Vector{Int64})

Define a Lagrangian grid space for shape functions. The grid size is defined by the polynomial degree and not the number of integration points.

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
function define_lagrangian_grid_space(dof::Int64, p::Vector{Int64})
    xi::Matrix{Float64} = zeros(Float64, dof, maximum(p) + 1)
    len::Float64 = 0.0
    if minimum(p) == 0
        @error "p order for lagarangian grid space must be at least p = 1 and not zero"
        return nothing
    end
    for idof = 1:dof
        len = 2.0 / p[idof]
        for ip = 1:p[idof]+1
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

function get_recursive_lagrange_shape_functions(
    xi::Vector{Float64},
    value::Union{Float64,Int64},
    p::Int64,
)
    N::Vector{Float64} = zeros(Float64, p + 1)
    for ip::Int64 = 1:p+1
        N[ip] = 1
        for jp::Int64 = 1:p+1
            if ip == jp
                continue
            end
            N[ip] *= (value - xi[jp]) / (xi[ip] - xi[jp])
        end
    end
    ## TODO check if this is valid, if higher order elements are in place
    return N
end

function get_recursive_lagrange_shape_functions_derivative(
    xi::Vector{Float64},
    value::Union{Float64,Int64},
    p::Int64,
)
    # https://en.wikipedia.org/wiki/Lagrange_polynomial#Derivation[6]
    # sympy calculated
    B::Vector{Float64} = zeros(Float64, p + 1)
    for ip = 1:p+1
        for jp = 1:p+1
            if ip == jp
                continue
            end
            temp = 1 / (xi[ip] - xi[jp])
            for mp = 1:p+1
                if mp != ip && mp != jp
                    temp *= (value - xi[mp]) / (xi[ip] - xi[mp])
                end
            end
            B[ip] += temp
        end
    end
    ## TODO check if this is valid, if higher order elements are in place
    return B
end
function get_2D_matrices(
    p::Vector{Int64},
    ip_weights::Matrix{Float64},
    ip_coordinates::Matrix{Float64},
    xi::Matrix{Float64},
    N::Array{Float64},
    B::Array{Float64},
)
    Nxi::Vector{Float64} = zeros(p[1] + 1)
    Neta::Vector{Float64} = zeros(p[2] + 1)

    Bxi::Vector{Float64} = zeros(p[1] + 1)
    Beta::Vector{Float64} = zeros(p[2] + 1)

    for (point_id, (ip_coordinate, ip_weight)) in
        enumerate(zip(eachrow(ip_coordinates), eachrow(ip_weights)))

        for jID = 1:p[2]+1
            Neta[jID] =
                get_recursive_lagrange_shape_functions(xi[2, :], ip_coordinate[2], p[2])[jID] *
                ip_weight[2]
            Beta[jID] =
                get_recursive_lagrange_shape_functions_derivative(
                    xi[2, :],
                    ip_coordinate[2],
                    p[2],
                )[jID] * ip_weights[2]
            for iID = 1:p[1]+1
                pos = ((jID - 1) * (p[1] + 1) + iID - 1) * 2

                Nxi[iID] =
                    get_recursive_lagrange_shape_functions(
                        xi[1, :],
                        ip_coordinate[1],
                        p[1],
                    )[iID] * ip_weight[1]
                Bxi[iID] =
                    get_recursive_lagrange_shape_functions_derivative(
                        xi[1, :],
                        ip_coordinate[1],
                        p[1],
                    )[iID] * ip_weight[1]
                N[point_id, pos+1, 1] = Nxi[iID] * Neta[jID]
                N[point_id, pos+2, 2] = Nxi[iID] * Neta[jID]

                B[point_id, pos+1, 1] = Bxi[iID] * Neta[jID]
                B[point_id, pos+2, 2] = Nxi[iID] * Beta[jID]
                B[point_id, pos+1, 3] = Nxi[iID] * Beta[jID]
                B[point_id, pos+2, 3] = Bxi[iID] * Neta[jID]
            end
        end
    end
    return N, B
end

function get_3D_matrices(
    p::Vector{Int64},
    ip_weights::Matrix{Float64},
    ip_coordinates::Matrix{Float64},
    xi::Matrix{Float64},
    N::Array{Float64},
    B::Array{Float64},
)
    Nxi::Vector{Float64} = zeros(p[1] + 1)
    Neta::Vector{Float64} = zeros(p[2] + 1)
    Npsi::Vector{Float64} = zeros(p[3] + 1)
    Bxi::Vector{Float64} = zeros(p[1] + 1)
    Beta::Vector{Float64} = zeros(p[2] + 1)
    Bpsi::Vector{Float64} = zeros(p[3] + 1)
    for (point_id, (ip_coordinate, ip_weight)) in
        enumerate(zip(eachrow(ip_coordinates), eachrow(ip_weights)))
        for kID = 1:p[3]+1
            Npsi[kID] =
                get_recursive_lagrange_shape_functions(xi[3, :], ip_coordinate[3], p[3])[kID] *
                ip_weight[3]
            Bpsi[kID] =
                get_recursive_lagrange_shape_functions_derivative(
                    xi[3, :],
                    ip_coordinate[3],
                    p[3],
                )[kID] * ip_weights[3]
            for jID = 1:p[2]+1
                Neta[jID] =
                    get_recursive_lagrange_shape_functions(
                        xi[2, :],
                        ip_coordinate[2],
                        p[2],
                    )[jID] * ip_weight[2]
                Beta[jID] =
                    get_recursive_lagrange_shape_functions_derivative(
                        xi[2, :],
                        ip_coordinate[2],
                        p[2],
                    )[jID] * ip_weights[2]
                for iID = 1:p[1]+1

                    pos =
                        (
                            (jID - 1) * (p[1] + 1) +
                            (kID - 1) * (p[2] + 1) * (p[1] + 1) +
                            iID - 1
                        ) * 3
                    Nxi[iID] =
                        get_recursive_lagrange_shape_functions(
                            xi[1, :],
                            ip_coordinate[1],
                            p[1],
                        )[iID] * ip_weights[1]
                    Bxi[iID] =
                        get_recursive_lagrange_shape_functions_derivative(
                            xi[1, :],
                            ip_coordinate[1],
                            p[1],
                        )[iID] * ip_weights[1]
                    N[point_id, pos+2, 2] = Nxi[iID] * Neta[jID] * Npsi[kID]
                    N[point_id, pos+1, 1] = Nxi[iID] * Neta[jID] * Npsi[kID]
                    N[point_id, pos+3, 3] = Nxi[iID] * Neta[jID] * Npsi[kID]

                    B[point_id, pos+1, 1] = Bxi[iID] * Neta[jID] * Npsi[kID]
                    B[point_id, pos+2, 2] = Nxi[iID] * Beta[jID] * Npsi[kID]
                    B[point_id, pos+3, 3] = Nxi[iID] * Neta[jID] * Bpsi[kID]
                    B[point_id, pos+2, 4] = Nxi[iID] * Neta[jID] * Bpsi[kID]
                    B[point_id, pos+3, 4] = Nxi[iID] * Beta[jID] * Npsi[kID]
                    B[point_id, pos+1, 5] = Nxi[iID] * Neta[jID] * Bpsi[kID]
                    B[point_id, pos+3, 5] = Bxi[iID] * Neta[jID] * Npsi[kID]
                    B[point_id, pos+1, 6] = Nxi[iID] * Beta[jID] * Npsi[kID]
                    B[point_id, pos+2, 6] = Bxi[iID] * Neta[jID] * Npsi[kID]


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
function create_element_matrices(
    dof::Int64,
    num_int::Vector{Int64},
    p::Vector{Int64},
    ip_weights::Matrix{Float64},
    ip_coordinates::Matrix{Float64},
)
    if dof > 3 || dof < 2
        @error "Not support degree of freedom for the finite element matrix creation"
        return nothing
    end

    N::Array{Float64} = zeros(Float64, prod(num_int), prod(p .+ 1) * dof, dof)
    B::Array{Float64} = zeros(Float64, prod(num_int), prod(p .+ 1) * dof, 3 * dof - 3)
    xi = define_lagrangian_grid_space(dof, p)

    if dof == 2
        return get_2D_matrices(p, ip_weights, ip_coordinates, xi, N, B)
    end
    return get_3D_matrices(p, ip_weights, ip_coordinates, xi, N, B)
end
end
