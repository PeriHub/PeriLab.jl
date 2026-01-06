# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
export get_all_field_keys

export create_constant_node_scalar_field
export create_node_scalar_field
export create_constant_node_vector_field
export create_node_vector_field
export create_constant_node_tensor_field
export create_node_tensor_field

export create_constant_bond_scalar_state
export create_bond_scalar_state
export create_constant_bond_vector_state
export create_bond_vector_state
export create_constant_bond_tensor_state
export create_bond_tensor_state

export create_constant_element_vector_field

export create_constant_free_size_field
export create_free_size_field

export
       NodeScalarField,
       NodeVectorField,
       NodeTensorField,
       BondScalarState,
       BondVectorState,
       BondTensorState

"""
	get_field(name::String, time::String)

Returns the field with the given name and time.

# Arguments
- `name::String`: The name of the field.
- `time::String`: The time of the field.
# Returns
- `field::Field`: The field with the given name and time.
"""
function get_field(name::String, time::String = "Constant")
    @assert !occursin("NP1", name)
    if time == "Constant"
        return _get_field(name)
    elseif time == "N"
        try
            return _get_field(data["NP1_to_N"][name].N)
        catch
            @error "Field ''" *
                   name *
                   "'' does not exist. Check if it is initialized as constant."
            return nothing
        end
    elseif time == "NP1"
        try
            return _get_field(data["NP1_to_N"][name].NP1)
        catch
            @error "Field ''" *
                   name *
                   "'' does not exist. Check if it is initialized as constant."
            return nothing
        end
    else
        @error "Time $time is not supported. Use 'constant', 'N', or 'NP1'"
        return nothing
    end
end

"""
	get_field_if_exists(name::String, time::String)

Returns the field with the given name if it exists.

# Arguments
- `name::String`: The name of the field.
- `time::String`: The time of the field.
# Returns
- `field::Field`: The field with the given name and time.
"""
function get_field_if_exists(name::String, time::String = "Constant")
    return has_key(name) ? get_field(name, time) : nothing
end

"""
	_get_field(name::String)

Returns the field with the given name.

# Arguments
- `name::String`: The name of the field.
# Returns
- `field::Field`: The field with the given name.
"""
function _get_field(name::String)::Union{Array,Nothing}
    try
        return fieldmanager.fields[name].data
    catch
        @error "Field ''" *
               name *
               "'' does not exist. \n - Check if it is initialized as non-constant. \n - Check if the model is not activated in the solver options, e.g. Pre Calculation Models: False"
        return nothing
    end
end

"""
	get_all_field_keys()

Returns a list of all field keys.
"""
function get_all_field_keys()::Vector{String}
    return data["field_names"]
end

const NodeScalarField = Vector{T} where {T<:Union{Int64,Float64,Bool}}
const NodeVectorField = Matrix{T} where {T<:Union{Int64,Float64,Bool}}
const NodeTensorField = Array{T,N} where {T<:Union{Int64,Float64,Bool},N}

const BondScalarState = Vector{Vector{T}} where {T<:Union{Int64,Float64,Bool}}
const BondVectorState = Vector{Vector{Vector{T}}} where {T<:Union{Int64,Float64,Bool}}
const BondTensorState = Vector{Array{T,N}} where {T<:Union{Int64,Float64,Bool},N}

function create_constant_node_scalar_field(name::String, vartype::Type{T};
                                           default_value::Number = 0) where {T<:Union{Int64,
                                                                                      Float64,
                                                                                      Bool}}
    return _create_node_scalar_field(name, T(default_value))
end

function create_node_scalar_field(name::String, vartype::Type{T};
                                  default_value::Number = 0) where {T<:Union{Int64,Float64,
                                                                             Bool}}
    set_NP1_to_N(name, vartype)

    return _create_node_scalar_field(name * "N", T(default_value)),
           _create_node_scalar_field(name * "NP1", T(default_value))
end

function create_constant_node_vector_field(name::String, vartype::Type{T}, dof::Int64;
                                           default_value::Number = 0) where {T<:Union{Int64,
                                                                                      Float64,
                                                                                      Bool}}
    return _create_node_vector_field(name, T(default_value), dof)
end

function create_node_vector_field(name::String, vartype::Type{T}, dof::Int64;
                                  default_value::Number = 0) where {T<:Union{Int64,Float64,
                                                                             Bool}}
    set_NP1_to_N(name, vartype)

    return _create_node_vector_field(name * "N", T(default_value), dof),
           _create_node_vector_field(name * "NP1", T(default_value), dof)
end

function create_constant_node_tensor_field(name::String, vartype::Type{T}, dof::Int64;
                                           default_value::Number = 0) where {T<:Union{Int64,
                                                                                      Float64,
                                                                                      Bool}}
    return _create_node_tensor_field(name, T(default_value), dof)
end

function create_node_tensor_field(name::String, vartype::Type{T}, dof::Int64;
                                  default_value::Number = 0) where {T<:Union{Int64,Float64,
                                                                             Bool}}
    set_NP1_to_N(name, vartype)

    return _create_node_tensor_field(name * "N", T(default_value), dof),
           _create_node_tensor_field(name * "NP1", T(default_value), dof)
end

function create_constant_bond_scalar_state(name::String, vartype::Type{T};
                                           default_value::Number = 0) where {T<:Union{Int64,
                                                                                      Float64,
                                                                                      Bool}}
    return _create_bond_scalar_state(name, T(default_value))
end

function create_bond_scalar_state(name::String, vartype::Type{T};
                                  default_value::Number = 0) where {T<:Union{Int64,Float64,
                                                                             Bool}}
    set_NP1_to_N(name, vartype)

    return _create_bond_scalar_state(name * "N", T(default_value)),
           _create_bond_scalar_state(name * "NP1", T(default_value))
end

function create_constant_bond_vector_state(name::String, vartype::Type{T}, dof::Int64;
                                           default_value::Number = 0) where {T<:Union{Int64,
                                                                                      Float64,
                                                                                      Bool}}
    return _create_bond_vector_state(name, T(default_value), dof)
end

function create_bond_vector_state(name::String, vartype::Type{T}, dof::Int64;
                                  default_value::Number = 0) where {T<:Union{Int64,Float64,
                                                                             Bool}}
    set_NP1_to_N(name, vartype)

    return _create_bond_vector_state(name * "N", T(default_value), dof),
           _create_bond_vector_state(name * "NP1", T(default_value), dof)
end

function create_constant_bond_tensor_state(name::String, vartype::Type{T}, dof::Int64;
                                           default_value::Number = 0) where {T<:Union{Int64,
                                                                                      Float64,
                                                                                      Bool}}
    return _create_bond_tensor_state(name, T(default_value), dof)
end

function create_bond_tensor_state(name::String, vartype::Type{T}, dof::Int64;
                                  default_value::Number = 0) where {T<:Union{Int64,Float64,
                                                                             Bool}}
    set_NP1_to_N(name, vartype)

    return _create_bond_tensor_state(name * "N", T(default_value), dof),
           _create_bond_tensor_state(name * "NP1", T(default_value), dof)
end

function create_constant_element_vector_field(name::String, vartype::Type{T}, dof::Int64;
                                              default_value::Number = 0) where {T<:Union{Int64,
                                                                                         Float64,
                                                                                         Bool}}
    return _create_element_vector_field(name, T(default_value), dof)
end

function create_constant_free_size_field(name::String, vartype::Type{T}, dof::Tuple;
                                         default_value::Number = 0) where {T<:Union{Int64,
                                                                                    Float64,
                                                                                    Bool}}
    return _create_free_size_field(name, T(default_value), dof)
end
function create_free_size_field(name::String, vartype::Type{T}, dof::Tuple;
                                default_value::Number = 0) where {T<:Union{Int64,Float64,
                                                                           Bool}}
    set_NP1_to_N(name, vartype)

    return _create_free_size_field(name * "N", T(default_value), dof),
           _create_free_size_field(name * "NP1", T(default_value), dof)
end

function create_field!(name::String, _data, field_type::String)
    if has_key(name)
        if size(_get_field(name), 1) != data["nnodes"]
            @warn "Field $name exists already with different size. Predefined field is returned"
        end
        return _get_field(name)
    end

    T = eltype(_data)
    D = typeof(_data)
    fieldmanager.fields[name] = DataField{T,D}(name, _data)

    data["field_types"][name] = Dict("type" => field_type, "vartype" => T)

    data["field_names"] = Vector{String}(collect(keys(data["field_types"])))

    return _data
end

function _create_node_scalar_field(name::String,
                                   value::T) where {T<:Union{Int64,Float64,Bool}}
    return create_field!(name, fill(value, data["nnodes"]), "NodeScalarField")
end
function _create_node_vector_field(name::String, value::T,
                                   dof::Int64) where {T<:Union{Int64,Float64,Bool}}
    @assert dof > 1
    return create_field!(name,
                         fill(value, data["nnodes"], dof), "NodeVectorField")
end
function _create_node_tensor_field(name::String, value::T,
                                   dof::Int64) where {T<:Union{Int64,Float64,Bool}}
    return create_field!(name,
                         fill(value, (data["nnodes"], dof, dof)), "NodeTensorField")
end

function _create_bond_scalar_state(name::String,
                                   value::T) where {T<:Union{Int64,Float64,Bool}}
    nBonds = _get_field("Number of Neighbors")
    return create_field!(name,
                         [fill(value, n) for n in nBonds], "BondScalarState")
end
function _create_bond_vector_state(name::String, value::T,
                                   dof::Int64) where {T<:Union{Int64,Float64,Bool}}
    @assert dof > 1
    nBonds = _get_field("Number of Neighbors")
    return create_field!(name,
                         [[fill(value, dof) for _ in 1:n] for n in nBonds],
                         "BondVectorState")
end
function _create_bond_tensor_state(name::String, value::T,
                                   dof::Int64) where {T<:Union{Int64,Float64,Bool}}
    nBonds = _get_field("Number of Neighbors")
    return create_field!(name,
                         [fill(value, (n, dof, dof)) for n in nBonds], "BondTensorState")
end

function _create_element_scalar_field(name::String,
                                      value::T) where {T<:Union{Int64,Float64,Bool}}
    nElements = _get_field("Number of Element Neighbors")
    return create_field!(name, [fill(value, n) for n in nElements], "ElementScalarField")
end
function _create_element_vector_field(name::String, value::T,
                                      dof::Int64) where {T<:Union{Int64,Float64,Bool}}
    @assert dof > 1
    nElements = _get_field("Number of Element Neighbors")
    return create_field!(name, [fill(value, (n, dof)) for n in nElements],
                         "ElementVectorField")
end

function _create_free_size_field(name::String, value::T,
                                 dof::Tuple) where {T<:Union{Int64,Float64,Bool}}
    return create_field!(name, Array{T}(zeros(dof)), "FreeSizeField")
end
