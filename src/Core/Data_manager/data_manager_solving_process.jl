# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

export switch_NP1_to_N
export set_NP1_to_N
export switch_bonds!
function switch_bonds!(field_N::Vector{Vector{T}},
                       field_NP1::Vector{Vector{T}}) where {T<:Union{Int64,Float64}}
    for fieldID in eachindex(field_NP1)
        copyto!(field_N[fieldID], field_NP1[fieldID])
    end
end

"""
    set_NP1_to_N(name::String, type::Type)

Sets the NP1_to_N dataset

# Arguments
- `name`::String: The name of the field.
- `type`::Type The field type
"""
function set_NP1_to_N(name::String, ::Type{T}) where {T}
    println(name)
    if !has_key(name)
        data["NP1_to_N"][name] = NP1_to_N(name * "N", name * "NP1", zero(T))
    end
end

function swap_n_np1!(entry::NP1_to_N{T}) where {T}
    tmp = entry.N
    entry.N = entry.NP1
    entry.NP1 = tmp
    return nothing
end

"""
    switch_NP1_to_N()

Switches the fields from NP1 to N. This is more efficient than copying the data from one field the other.

"""
function switch_NP1_to_N() ##TODO check type stability
    active = _get_field("Active")
    for key in keys(data["NP1_to_N"])
        if key == "Bond Damage"
            field_NP1 = get_field(key, "NP1")
            field_N = get_field(key, "N")
            switch_bonds!(field_N, field_NP1)
            continue
        end
        swap_n_np1!(data["NP1_to_N"][key])

        field_NP1 = get_field(key, "NP1")
        fill_field!(field_NP1, active, data["NP1_to_N"][key].value)
    end
end
const NestedArray{T} = Union{AbstractArray{<:AbstractArray{T}},
                             AbstractArray{<:AbstractArray{<:AbstractArray{T}}},
                             Vector{Vector{Vector{T}}}}
function fill_field!(field_NP1::NestedArray{T},
                     active::Vector{Bool}, value::T) where {T<:Union{Int64,Float64,Bool}}
    fill_in_place!(field_NP1, value, active)
end

function fill_field!(field_NP1::AbstractArray{T},
                     active::Vector{Bool}, value::T) where {T<:Union{Int64,Float64,Bool}}
    fill!(field_NP1, value)
end
