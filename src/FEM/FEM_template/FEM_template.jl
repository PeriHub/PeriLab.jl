# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module FEM_template

using ..Data_Manager
using .......PeriLabExceptions: @abort

"""
  element_name()

Gives the element name. It is compared with `Element Type` of the `FEM` section of
the input deck.

# Arguments

# Returns
- `name::String`: The name of the element model.

Example:
```julia
println(element_name())
"element Template"
```
"""
function element_name()
    return "element Template"
end

"""
    init_element(elements, fem, p)

Initializes the element formulation. This template has to be copied, the file renamed and edited by the user to create a new element. Additional files can be called from here using include and `import .any_module` or `using .any_module`.

# Arguments
- `elements::AbstractVector{Int64}`: List of elements.
- `fem::FEMParams`: The `FEM` section of the input deck.
- `p::Vector{Int64}`: The polynomial degree in each direction.
"""
function init_element(elements::AbstractVector{Int64},
                      fem,
                      p::Vector{Int64})
end

"""
    create_element_matrices(dof, num_int, p, ip_weights, ip_coordinates)

Creates the shape function matrix N and the strain-displacement matrix B at the
integration points of the element.

# Arguments
- `dof::Int64`: The degrees of freedom (2 or 3).
- `num_int::Vector{Int64}`: The number of integration points in each direction.
- `p::Vector{Int64}`: The polynomial degree in each direction.
- `ip_weights::Matrix{Float64}`: The weights of the integration points.
- `ip_coordinates::Matrix{Float64}`: The coordinates of the integration points.
# Returns
- `N`, `B`: see `Lagrange_element.create_element_matrices`.
"""
function create_element_matrices(dof::Int64,
                                 num_int::Vector{Int64},
                                 p::Vector{Int64},
                                 ip_weights::Matrix{Float64},
                                 ip_coordinates::Matrix{Float64})
    @abort "element Template: create_element_matrices is not implemented."
end

end
