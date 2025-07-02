# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

export synch_manager
export switch_NP1_to_N
export set_NP1_to_N

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
function set_NP1_to_N(name::String, type::Type)
    if !has_key(name)
        data["NP1_to_N"][name] = [name * "N", name * "NP1", zero(type)]
    end
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
        data["NP1_to_N"][key][1],
        data["NP1_to_N"][key][2] = data["NP1_to_N"][key][2],
                                   data["NP1_to_N"][key][1]
        field_NP1 = get_field(key, "NP1")
        fill_field!(field_NP1, active, key)
    end
end
const NestedArray{T} = Union{AbstractArray{<:AbstractArray{T}},
                             AbstractArray{<:AbstractArray{<:AbstractArray{T}}},
                             Vector{Vector{Vector{T}}}}
function fill_field!(field_NP1::NestedArray{T},
                     active::Vector{Bool}, key::String) where {T<:Union{Int64,Float64,Bool}}
    fill_in_place!(field_NP1, data["NP1_to_N"][key][3], active)
end

function fill_field!(field_NP1::AbstractArray{T},
                     active::Vector{Bool}, key::String) where {T<:Union{Int64,Float64,Bool}}
    fill!(field_NP1, data["NP1_to_N"][key][3])
end

"""
    synch_manager(synchronise_field, direction::String)

Synchronises the fields.

# Arguments
- `synchronise_field`: The function to synchronise the field.
- `direction::String`: The direction of the synchronisation.
"""
function synch_manager(synchronise_field, direction::String)
    synch_fields = get_synch_fields()
    # @debug synch_fields
    for synch_field in keys(synch_fields)
        synchronise_field(get_comm(),
                          synch_fields,
                          get_overlap_map(),
                          get_field,
                          synch_field,
                          direction)
    end
end
