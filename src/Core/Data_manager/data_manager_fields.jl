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
export get_bond_ptr

export create_constant_free_size_field
export create_free_size_field

export
       NodeScalarField,
       NodeVectorField,
       NodeTensorField,
       NodeField,
       BondScalarState,
       BondVectorState,
       BondVectorField,
       BondTensorState,
       BondState,
       BondField,
       NestedNeighborhoodList,
       FreeSizeField

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
    if time == "Constant"
        return _get_field(name)
    elseif time == "N"
        try
            return _get_field(data["NP1_to_N"][name].N)
        catch
            @abort "Field ''" *
                   name *
                   "'' does not exist. Check if it is initialized as constant."
        end
    elseif time == "NP1"
        try
            return _get_field(data["NP1_to_N"][name].NP1)
        catch
            @abort "Field ''" *
                   name *
                   "'' does not exist. Check if it is initialized as constant."
        end
    else
        @abort "Time $time is not supported. Use 'constant', 'N', or 'NP1'"
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
function _get_field(name::String)::Union{AbstractArray,Nothing}
    try
        return fieldmanager.fields[name].data
    catch
        @abort "Field ''" *
               name *
               "'' does not exist. \n - Check if it is initialized as non-constant. \n - Check if the model is not activated in the solver options, e.g. Pre Calculation Models: False"
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
const NodeTensorField = Array{T,3} where {T<:Union{Int64,Float64,Bool}}
const NodeField{T} = Union{NodeScalarField{T},NodeVectorField{T},NodeTensorField{T}}

# One node of a CSR bond field: a contiguous view into the flat data.
const BondView{E} = SubArray{E,1,Vector{E},Tuple{UnitRange{Int64}},true}
# CSR bond field, flat data of element type `E` sharing the offsets `get_bond_ptr()`.
const BondField{E} = VectorOfArrays{E,1,0,Vector{E},Vector{Int64},Vector{Tuple{}},
                                    BondView{E}}

const BondScalarState = BondField{T} where {T<:Union{Int64,Float64,Bool}}
# Nested neighborhood lists that cannot share the CSR offsets of the bond fields: the global
# list before it is distributed to the ranks (still edited in place by bond filters and FE
# topology) and the filtered bond index lists of the contact bond filter. The local
# neighborhood list on each rank is a `BondScalarState{Int64}`.
const NestedNeighborhoodList = Vector{Vector{Int64}}
# `D` is the dof; code that only knows `dof` at run time reaches the kernels through a
# function barrier (an argument annotated `BondVectorState{T}` specializes on `D`).
const BondVectorState{T,D} = BondField{SVector{D,T}} where {T<:Union{Int64,Float64,Bool},D}
# Vector state of a 2D or 3D model. Use it for locals that hold a vector state fetched with
# `get_field`: the union lets the compiler split kernel calls into static ones, while the
# open `BondVectorState{T}` makes every call dynamic (and allocating).
const BondVectorField{T} = Union{BondVectorState{T,2},BondVectorState{T,3}}
# `D` is the dof, `L = D * D`
const BondTensorState{T,D,L} = BondField{SMatrix{D,D,T,L}} where {T<:Union{Int64,Float64,
                                                                           Bool},
                                                                  D,L}
const BondState{T} = Union{BondScalarState{T},BondVectorState{T},BondTensorState{T}}

const FreeSizeField = Array{T,N} where {T<:Union{Int64,Float64,Bool},N}

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

    return _create_node_scalar_field(data["NP1_to_N"][name].N, T(default_value)),
           _create_node_scalar_field(data["NP1_to_N"][name].NP1, T(default_value))
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

    return _create_node_vector_field(data["NP1_to_N"][name].N, T(default_value), dof),
           _create_node_vector_field(data["NP1_to_N"][name].NP1, T(default_value), dof)
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

    return _create_node_tensor_field(data["NP1_to_N"][name].N, T(default_value), dof),
           _create_node_tensor_field(data["NP1_to_N"][name].NP1, T(default_value), dof)
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

    return _create_bond_scalar_state(data["NP1_to_N"][name].N, T(default_value)),
           _create_bond_scalar_state(data["NP1_to_N"][name].NP1, T(default_value))
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

    return _create_bond_vector_state(data["NP1_to_N"][name].N, T(default_value), dof),
           _create_bond_vector_state(data["NP1_to_N"][name].NP1, T(default_value), dof)
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

    return _create_bond_tensor_state(data["NP1_to_N"][name].N, T(default_value), dof),
           _create_bond_tensor_state(data["NP1_to_N"][name].NP1, T(default_value), dof)
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

    return _create_free_size_field(data["NP1_to_N"][name].N, T(default_value), dof),
           _create_free_size_field(data["NP1_to_N"][name].NP1, T(default_value), dof)
end

function create_field!(name::String, _data, field_type::String)
    if has_key(name)
        if size(_get_field(name), 1) != data["nnodes"]
            @warn "Field $name exists already with different size. Predefined field is returned"
        end
        return _get_field(name)
    end

    T = _vartype(_data)
    D = typeof(_data)
    fieldmanager.fields[name] = DataField{T,D}(name, _data)

    data["field_types"][name] = Dict("type" => field_type, "vartype" => T)

    data["field_names"] = Vector{String}(collect(keys(data["field_types"])))

    return _data
end

_vartype(field) = eltype(field)
# scalar element type, also for SVector / SMatrix bond values
_vartype(field::VectorOfArrays) = eltype(eltype(flatview(field)))

"""
	get_bond_ptr()

Returns the CSR offsets shared by all bond fields: the bonds of node `i` are the flat
indices `ptr[i]:(ptr[i + 1] - 1)`. Built from "Number of Neighbors" on first use. The bond
structure is fixed afterwards, a changed "Number of Neighbors" aborts.
"""
function get_bond_ptr()::Vector{Int64}
    ptr::Vector{Int64} = data["bond_ptr"]
    nbonds = _get_field("Number of Neighbors")
    if isempty(ptr)
        resize!(ptr, length(nbonds) + 1)
        ptr[1] = 1
        @inbounds for i in eachindex(nbonds)
            ptr[i + 1] = ptr[i] + nbonds[i]
        end
        return ptr
    end
    if length(ptr) != length(nbonds) + 1 ||
       any(i -> ptr[i + 1] - ptr[i] != nbonds[i], eachindex(nbonds))
        @abort "\"Number of Neighbors\" changed after the first bond field was created. The bond structure must not change."
    end
    return ptr
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
    ptr = get_bond_ptr()
    return create_field!(name, VectorOfVectors(fill(value, ptr[end] - 1), ptr),
                         "BondScalarState")
end
function _create_bond_vector_state(name::String, value::T,
                                   dof::Int64) where {T<:Union{Int64,Float64,Bool}}
    @assert dof > 1
    return _create_bond_vector_state(name, value, Val(dof))
end
function _create_bond_vector_state(name::String, value::T,
                                   ::Val{D}) where {T<:Union{Int64,Float64,Bool},D}
    ptr = get_bond_ptr()
    data = fill(SVector{D,T}(ntuple(_ -> value, Val(D))), ptr[end] - 1)
    return create_field!(name, VectorOfVectors(data, ptr), "BondVectorState")
end
function _create_bond_tensor_state(name::String, value::T,
                                   dof::Int64) where {T<:Union{Int64,Float64,Bool}}
    return _create_bond_tensor_state(name, value, Val(dof))
end
function _create_bond_tensor_state(name::String, value::T,
                                   ::Val{D}) where {T<:Union{Int64,Float64,Bool},D}
    ptr = get_bond_ptr()
    data = fill(SMatrix{D,D,T,D * D}(ntuple(_ -> value, Val(D * D))), ptr[end] - 1)
    return create_field!(name, VectorOfVectors(data, ptr), "BondTensorState")
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
