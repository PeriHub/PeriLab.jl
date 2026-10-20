# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

export switch_NP1_to_N
export reset_NP1
export set_NP1_to_N
export switch_bonds!
function switch_bonds!(field_N::BondField{E}, field_NP1::BondField{E}) where {E}
    copyto!(flatview(field_N), flatview(field_NP1))
end

"""
    set_NP1_to_N(name::String, type::Type)

Sets the NP1_to_N dataset

# Arguments
- `name`::String: The name of the field.
- `type`::Type The field type
"""
function set_NP1_to_N(name::String, ::Type{T}) where {T<:Union{Float64,Bool,Int64}}
    if !has_key(name * "N")
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
    if isnothing(field_NP1)
        @abort "Field $(np1_to_n.NP1) not found. Cannot fill NP1 field from N."
        return
    end
    fill_field!(field_NP1, active, np1_to_n.value)
end

"""
    switch_NP1_to_N()

Switches the fields from NP1 to N. This is more efficient than copying the data from one field the other.

"""
function switch_NP1_to_N()
    active = _get_field("Active")::Vector{Bool}
    np1_to_n = data["NP1_to_N"]::Dict{String,NP1_to_N}
    if haskey(np1_to_n, "Bond Damage")
        field_NP1 = get_field("Bond Damage", "NP1")::BondScalarState{Float64}
        field_N = get_field("Bond Damage", "N")::BondScalarState{Float64}
        switch_bonds!(field_N, field_NP1)
    end

    for (key, entry) in np1_to_n
        key == "Bond Damage" && continue
        swap_n_np1!(entry)
        _fill_field_barrier(get_field(key, "NP1"), active, entry)
    end
end

"""
    reset_NP1()

Discards the NP1 state of a failed step: every NP1 field is as directly after
`switch_NP1_to_N` (the fill value on active nodes, `Bond Damage` NP1 equal to N).
"""
function reset_NP1()
    active::Vector = _get_field("Active")
    if haskey(data["NP1_to_N"], "Bond Damage")
        switch_bonds!(get_field("Bond Damage", "NP1"), get_field("Bond Damage", "N"))
    end
    for key::String in keys(data["NP1_to_N"])
        key == "Bond Damage" && continue
        _fill_field_barrier(get_field(key, "NP1"), active, data["NP1_to_N"][key])
    end
    return nothing
end

@inline function fill_field!(field_NP1::BondField{E},
                             active::Vector{Bool},
                             value::T) where {E,T<:Union{Int64,Float64,Bool}}
    fill_in_place!(field_NP1, _bond_value(E, value), active)
end

_bond_value(::Type{E}, value) where {E<:Number} = E(value)
_bond_value(::Type{E}, value) where {E<:StaticArray} = E(ntuple(_ -> value, Val(length(E))))

@inline function fill_field!(field_NP1::Array{T,N},
                             active::Vector{Bool},
                             value::T) where {T<:Union{Int64,Float64,Bool},N}
    fill!(field_NP1, value)
end

@inline function fill_in_place!(A::BondField{E},
                                value::E,
                                active::Vector{Bool}) where {E}
    @inbounds for i in eachindex(active)
        active[i] && fill!(A[i], value)
    end
    return nothing
end
