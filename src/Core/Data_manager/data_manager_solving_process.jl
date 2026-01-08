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
function set_NP1_to_N(name::String, ::Type{T}) where {T<:Union{Float64,Bool,Int64}}
    if !has_key(name)
        data["NP1_to_N"][name] = NP1_to_N(name * "N", name * "NP1", zero(T))
    end
end

function swap_n_np1!(entry::NP1_to_N{T}) where {T<:Union{Float64,Bool,Int64}}
    tmp = entry.N
    entry.N = entry.NP1
    entry.NP1 = tmp
    return nothing
end

@inline function _fill_field_barrier(field_NP1, active::Vector{Bool},
                                     np1_to_n::NP1_to_N{T}) where {T}
    fill_field!(field_NP1, active, np1_to_n.value)
end

"""
    switch_NP1_to_N()

Switches the fields from NP1 to N. This is more efficient than copying the data from one field the other.

"""
function switch_NP1_to_N() ##TODO check type stability
    active::Vector = _get_field("Active")
    if haskey(data["NP1_to_N"], "Bond Damage")
        field_NP1 = get_field("Bond Damage", "NP1")
        field_N = get_field("Bond Damage", "N")
        switch_bonds!(field_N, field_NP1)
    end

    for key::String in keys(data["NP1_to_N"])
        key == "Bond Damage" && continue
        swap_n_np1!(data["NP1_to_N"][key])
        field_NP1 = get_field(key, "NP1")
        _fill_field_barrier(field_NP1, active, data["NP1_to_N"][key])
    end
end

@inline function fill_field!(field_NP1::Vector{Vector{Vector{T}}},
                             active::Vector{Bool},
                             value::T) where {T<:Union{Int64,Float64,Bool}}
    fill_in_place!(field_NP1, value, active)
end

@inline function fill_field!(field_NP1::Vector{Vector{T}},
                             active::Vector{Bool},
                             value::T) where {T<:Union{Int64,Float64,Bool}}
    fill_in_place!(field_NP1, value, active)
end

@inline function fill_field!(field_NP1::Vector{Array{T,3}},
                             active::Vector{Bool},
                             value::T) where {T<:Union{Int64,Float64,Bool}}
    fill_in_place!(field_NP1, value, active)
end

@inline function fill_field!(field_NP1::Array{T,N},
                             active::Vector{Bool},
                             value::T) where {T<:Union{Int64,Float64,Bool},N}
    fill!(field_NP1, value)
end

@inline function fill_in_place!(A::Vector{Vector{T}},
                                value::T,
                                active::Vector{Bool}) where {T<:Union{Int64,Float64,Bool}}
    @inbounds @simd for i in eachindex(active)
        active[i] && fill!(A[i], value)
    end
end

@inline function fill_in_place!(A::Vector{Array{T,3}},
                                value::T,
                                active::Vector{Bool}) where {T<:Union{Int64,Float64,Bool}}
    @inbounds @simd for i in eachindex(active)
        active[i] && fill!(A[i], value)
    end
end
@inline function fill_in_place!(A::Vector{Vector{Vector{T}}},
                                value::T,
                                active::Vector{Bool}) where {T<:Union{Int64,Float64,Bool}}
    @inbounds for i in eachindex(active)
        active[i] || continue
        @simd for j in eachindex(A[i])
            fill!(A[i][j], value)
        end
    end
    return nothing
end
@inline function fill_in_place!(A::Vector{Vector{Array{T,3}}},
                                value::T,
                                active::Vector{Bool}) where {T<:Union{Int64,Float64,Bool}}
    @inbounds for i in eachindex(active)
        active[i] || continue
        @simd for j in eachindex(A[i])
            fill!(A[i][j], value)
        end
    end
    return nothing
end
