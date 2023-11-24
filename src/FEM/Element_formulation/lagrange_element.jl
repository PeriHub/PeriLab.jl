# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

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


function get_recursive_lagrange_shape_functions_derivative()

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
