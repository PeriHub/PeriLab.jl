# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

export create_bond_field
export create_constant_free_size_field
export create_constant_bond_field
export create_constant_node_field
export create_constant_element_field
export create_node_field

"""
    create_bond_field(name::String, vartype::Type, dof::Int64)

Creates a bond field with the given name, data type, and degree of freedom.

# Arguments
- `name::String`: The name of the bond field.
- `vartype::Type`: The data type of the bond field.
- `dof::Int64`: The degrees of freedom per bond.
- `VectorOrMatrix::String` (optional) - Vector or Materix; Default is vector

# Returns
- `bond_field::Field`: The created bond field for the current time step.
- `bond_field_np1::Field`: The created bond field for the next time step.

Example:
```julia
create_bond_field("stress", Float64, 6)  # creates a stress bond field with 6 degrees of freedom
```
"""

function create_bond_field(name::String,
                           vartype::Type{T},
                           dof::Int64,
                           default_value::Number = 0;
                           VectorOrMatrix::String = "Vector") where {T<:Union{Int64,Float64,
                                                                              Bool}}
    set_NP1_to_N(name, vartype)
    return create_field(name * "N", vartype, "Bond_Field", dof, T(default_value),
                        VectorOrMatrix),
           create_field(name * "NP1", vartype, "Bond_Field", dof, T(default_value),
                        VectorOrMatrix)
end

"""
    create_constant_bond_field(name::String, vartype::Type, dof::Int64, default_value::Union{Int64,Float64,Bool}=0))

Creates a constant bond field with the given name, data type, and degree of freedom.

# Arguments
- `name::String`: The name of the bond field.
- `vartype::Type`: The data type of the bond field.
- `dof::Int64`: The degrees of freedom per bond.
-  default_value::Union{Int64,Float64,Bool}=0) (optional) - filled with zero or false

# Returns
- `constant_bond_field::Field`: The created constant bond field.

Example:
```julia
create_constant_bond_field("density", Float64, 1)  # creates a density constant bond field
```
"""

function create_constant_bond_field(name::String,
                                    vartype::Type{T},
                                    dof::Int64,
                                    default_value::Number = 0;
                                    VectorOrMatrix::String = "Vector") where {T<:Union{Int64,
                                                                                       Float64,
                                                                                       Bool}}
    return create_field(name, vartype, "Bond_Field", dof, T(default_value), VectorOrMatrix)
end

function create_constant_free_size_field(name::String,
                                         vartype::Type{T},
                                         dof::NTuple{N,Int64},
                                         default_value::Number = 0) where {T<:Union{Int64,
                                                                                    Float64,
                                                                                    Bool},
                                                                           N}
    return create_field(name, vartype, "Free_Size_Field", dof, T(default_value))
end

function create_free_size_field(name::String,
                                vartype::Type{T},
                                dof::Tuple,
                                default_value::Number = 0) where {T<:Union{Int64,
                                                                           Float64,
                                                                           Bool}}
    set_NP1_to_N(name, vartype)
    return create_field(name * "N", vartype, "Free_Size_Field", dof, T(default_value)),
           create_field(name * "NP1", vartype, "Free_Size_Field", dof, T(default_value))
end

"""
    create_constant_node_field(name::String, vartype::Type, dof::Int64)

Creates a constant node field with the given name, data type, and degree of freedom.

# Arguments
- `name::String`: The name of the node field.
- `vartype::Type`: The data type of the node field.
- `dof::Int64`: The degrees of freedom per node.
- `VectorOrMatrix::String` (optional) - Vector or Materix; Default is vector

# Returns
- `constant_node_field::Field`: The created constant node field.

Example:
```julia
create_constant_node_field("temperature", Float64, 1)  # creates a temperature constant node field
```
"""

function create_constant_node_field(name::String,
                                    vartype::Type{T},
                                    dof::Int64,
                                    default_value::Number = 0;
                                    VectorOrMatrix::String = "Vector") where {T<:Union{Int64,
                                                                                       Float64,
                                                                                       Bool}}
    return create_field(name, vartype, "Node_Field", dof, T(default_value), VectorOrMatrix)
end

"""
    create_constant_element_field(name::String, vartype::Type, dof::Int64)

Creates a constant element field with the given name, data type, and degree of freedom.

# Arguments
- `name::String`: The name of the element field.
- `vartype::Type`: The data type of the element field.
- `dof::Int64`: The degrees of freedom per element.
- `VectorOrMatrix::String` (optional) - Vector or Materix; Default is vector

# Returns
- `constant_element_field::Field`: The created constant element field.

Example:
```julia
create_constant_element_field("temperature", Float64, 1)  # creates a temperature constant element field
```
"""
function create_constant_element_field(name::String,
                                       vartype::Type{T},
                                       dof::Int64,
                                       default_value::Number = 0) where {T<:Union{Int64,
                                                                                  Float64,
                                                                                  Bool}}
    return create_field(name, vartype, "Element_Field", dof, T(default_value))
end

"""
    create_node_field(name::String, vartype::Type, dof::Int64)

Creates a node field with the given name, data type, and degree of freedom.

# Arguments
- `name::String`: The name of the node field.
- `type::Type`: The data type of the node field.
- `dof::Int64`: The degree of freedom of each node.
- `VectorOrMatrix::String` (optional) - Vector or Materix; Default is vector
# Returns
- `node_field::Field`: The created node field for the current time step.
- `node_field_np1::Field`: The created node field for the next time step.

Example:
```julia
create_node_field("displacement", Float64, 3)  # creates a displacement node field with 3 degrees of freedom
```
"""

function create_node_field(name::String,
                           vartype::Type{T},
                           dof::Int64,
                           default_value::Number = 0;
                           VectorOrMatrix::String = "Vector") where {T<:Union{Int64,Float64,
                                                                              Bool}}
    set_NP1_to_N(name, vartype)
    return create_field(name * "N", vartype, "Node_Field", dof, T(default_value),
                        VectorOrMatrix),
           create_field(name * "NP1", vartype, "Node_Field", dof, T(default_value),
                        VectorOrMatrix)
end

"""
    create_field(name::String, vartype::Type, bondNode::String, dof::Int64, default_value::Any=0)

Create a field with the given `name` for the specified `vartype`. If the field already exists, return the existing field. If the field does not exist, create a new field with the specified characteristics.

# Arguments
- `name::String`: The name of the field.
- `vartype::Type`: The data type of the field.
- `dof::Int64`: The degrees of freedom per node.
- `default_value::Any`: The default value of the field.

# Returns
The field with the given `name` and specified characteristics.
"""
function create_field(name::String,
                      vartype::Type{T},
                      bond_or_node::String,
                      dof::Q,
                      value::T,
                      VectorOrMatrix::String = "Vector") where {T<:Union{Int64,Float64,
                                                                         Bool},
                                                                Q<:Union{Int64,
                                                                         Tuple{Vararg{Int64}}}}
    if has_key(name)
        if size(_get_field(name), 1) != data["nnodes"]
            @warn "Field $name exists already with different size. Predefined field is returned"
        end
        return _get_field(name)
    end

    if bond_or_node == "Node_Field"
        if dof == 1
            fields[vartype][name] = fill(value, data["nnodes"])
        else
            if VectorOrMatrix == "Matrix"
                fields[vartype][name] = fill(value, (data["nnodes"], dof, dof))
            else
                fields[vartype][name] = fill(value, data["nnodes"], dof)
                # fields[vartype][name] = [fill(value,dof) for j=1:data["nnodes"]]
            end
        end
    elseif bond_or_node == "Bond_Field"
        nBonds = _get_field("Number of Neighbors")
        if dof == 1
            fields[vartype][name] = [fill(value, n) for n in nBonds]
        else
            if VectorOrMatrix == "Matrix"
                fields[vartype][name] = [fill(value, (n, dof, dof)) for n in nBonds]
            else
                # fields[vartype][name] = [fill(value, (n, dof)) for n in nBonds]
                fields[vartype][name] = [[fill(value, dof) for j in 1:n] for n in nBonds]
            end
        end
    elseif bond_or_node == "Element_Field"
        nElements = _get_field("Number of Element Neighbors")
        if dof == 1
            fields[vartype][name] = [fill(value, n) for n in nElements]
        else
            fields[vartype][name] = [fill(value, (n, dof)) for n in nElements]
        end
    elseif bond_or_node == "Free_Size_Field"
        fields[vartype][name] = Array{vartype}(zeros(dof))
    end
    get_function = () -> fields[vartype][name]
    data["field_types"][name] = vartype
    data["field_array_type"][name] = Dict("Type" => VectorOrMatrix, "Dof" => dof,
                                          "get_function" => get_function)
    data["field_names"] = Vector{String}(collect(keys(data["field_types"])))
    return get_function()
end
