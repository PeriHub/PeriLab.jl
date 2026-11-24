# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

export switch_NP1_to_N
export reset_NP1
export set_NP1_to_N
export switch_bonds!
export track_bond_damage_N
export set_np1_start
function switch_bonds!(field_N::BondField{E}, field_NP1::BondField{E}) where {E}
    copyto!(flatview(field_N), flatview(field_NP1))
end

"""
    track_bond_damage_N()

Keeps `Bond Damage` N equal to NP1 after every step, for the models that read N (damage
models, the matrix based correspondence) and for a fallback solver, which restores NP1
from N. Without it the copy of all bonds is skipped; NP1 keeps the bond damage anyway.
"""
function track_bond_damage_N()
    data["Track Bond Damage N"] = true
    haskey(data["NP1_to_N"], "Bond Damage") &&
        switch_bonds!(get_bond_damage("N"), get_bond_damage("NP1"))
    return nothing
end

"""
    set_np1_start(name::String, mode::Symbol)

How the NP1 state of `name` starts a step (after `switch_NP1_to_N` and `reset_NP1`):
`:fill` with the fill value (default), `:copy` of N (a monotonic field whose model updates
only some nodes, e.g. the node damage), `:keep` as it is (a field every step overwrites on
every node it reads, e.g. the deformed bond geometry; skips the refill).
"""
function set_np1_start(name::String, mode::Symbol)
    mode in (:fill, :copy, :keep) || @abort "Unknown NP1 start mode $mode for $name."
    data["NP1 Start"][name] = mode
    return nothing
end

# starts the NP1 state of `key` as `set_np1_start` defines
function _start_np1!(key::String, entry::NP1_to_N, active::Vector{Bool})
    mode = get(data["NP1 Start"]::Dict{String,Symbol}, key, :fill)
    mode === :keep && return nothing
    if mode === :copy
        _copy_field!(get_field(key, "NP1"), get_field(key, "N"))
        return nothing
    end
    _fill_field_barrier(get_field(key, "NP1"), active, entry)
end

_copy_field!(field_NP1::AbstractArray, field_N::AbstractArray) = copyto!(field_NP1, field_N)
_copy_field!(field_NP1::BondField, field_N::BondField) = switch_bonds!(field_NP1, field_N)

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
    # the dispatch on the (untyped) field gets the entry, not its value: a value would be
    # boxed on every call
    fill_field!(field_NP1, active, np1_to_n)
end

"""
    switch_NP1_to_N()

Switches the fields from NP1 to N. This is more efficient than copying the data from one field the other.

"""
function switch_NP1_to_N()
    active = _get_field("Active")::Vector{Bool}
    np1_to_n = data["NP1_to_N"]::Dict{String,NP1_to_N}
    if haskey(np1_to_n, "Bond Damage") && data["Track Bond Damage N"]::Bool
        field_NP1 = get_field("Bond Damage", "NP1")::BondScalarState{Float64}
        field_N = get_field("Bond Damage", "N")::BondScalarState{Float64}
        switch_bonds!(field_N, field_NP1)
    end

    # every NP1 field is refilled after the swap unless `set_np1_start` says otherwise
    for (key, entry) in np1_to_n
        key == "Bond Damage" && continue
        swap_n_np1!(entry)
        _start_np1!(key, entry, active)
    end
end

"""
    reset_NP1()

Discards the NP1 state of a failed step: every NP1 field is as directly after
`switch_NP1_to_N` (the fill value on active nodes, `Bond Damage` NP1 equal to N if
tracked, see `track_bond_damage_N`).
"""
function reset_NP1()
    active::Vector = _get_field("Active")
    if haskey(data["NP1_to_N"], "Bond Damage") && data["Track Bond Damage N"]::Bool
        switch_bonds!(get_field("Bond Damage", "NP1"), get_field("Bond Damage", "N"))
    end
    for (key, entry) in data["NP1_to_N"]
        key == "Bond Damage" && continue
        _start_np1!(key, entry, active)
    end
    # restored bonds are not flagged by the damage models; the next step must recompute
    # every node (e.g. the shape tensor, see ONLY_FLAGGED_BONDS_CHANGED)
    has_key("Update") && (_get_field("Update") .= true)
    return nothing
end

@inline function fill_field!(field_NP1::BondField{E},
                             active::Vector{Bool},
                             entry::NP1_to_N{T}) where {E,T<:Union{Int64,Float64,Bool}}
    fill_in_place!(field_NP1, _bond_value(E, entry.value), active)
end

_bond_value(::Type{E}, value) where {E<:Number} = E(value)
_bond_value(::Type{E}, value) where {E<:StaticArray} = E(ntuple(_ -> value, Val(length(E))))

@inline function fill_field!(field_NP1::Array{T,N},
                             active::Vector{Bool},
                             entry::NP1_to_N{T}) where {T<:Union{Int64,Float64,Bool},N}
    fill!(field_NP1, entry.value)
end

@inline function fill_in_place!(A::BondField{E},
                                value::E,
                                active::Vector{Bool}) where {E}
    @inbounds for i in eachindex(active)
        active[i] && fill!(A[i], value)
    end
    return nothing
end
