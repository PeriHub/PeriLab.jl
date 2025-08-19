# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# =============================================================================
# ABSTRACT BASE TYPES
# =============================================================================

abstract type DataField{T} end

# =============================================================================
# NodeField Type
# =============================================================================
struct NodeField{T} <: DataField{T}
    name::String
    vartype::Type{T}
    data::Vector{T}
    structured::AbstractArray{T}  # Gespeicherte Matrix-Form

    function NodeField(name::String, vartype::Type{T}, value::T, nnodes::Int, dof::Int,
                       matrix::Bool = false) where {T}
        if matrix
            data = fill(vartype(value), nnodes * dof * dof)
            matrix_form = reshape(data, nnodes, dof, dof)
        else
            data = fill(vartype(value), nnodes * dof)
            if dof != 1
                matrix_form = reshape(data, nnodes, dof)
            else
                matrix_form = @views data
            end
        end
        new{T}(name, vartype, data, matrix_form)
    end
end

# =============================================================================
# BondField Type
# =============================================================================

# Index-Mapping Struktur
abstract type BondField{T} end

struct ScalarBondField{T} <: BondField{T}
    name::String
    data::Vector{T}
    nBonds::Vector{Int}
    cumBonds::Vector{Int}

    function ScalarBondField(name::String, ::Type{T}, value::T,
                             nBonds::Vector{Int}) where {T}
        data = fill(value, sum(nBonds))
        cumBonds = [0; cumsum(nBonds)]
        new{T}(name, data, nBonds, cumBonds)
    end
end

struct BondView{T} <: AbstractVector{T}  # <-- Erbt von AbstractVector!
    parent::ScalarBondField{T}
    element_idx::Int
end

# Erster Index [i] -> BondView (als Vector)
Base.getindex(s::ScalarBondField, i::Int) = BondView(s, i)

# Zweiter Index [j] auf BondView -> Wert
Base.getindex(bv::BondView, j::Int) = begin
    start_idx = bv.parent.cumBonds[bv.element_idx] + 1
    bv.parent.data[start_idx + j - 1]
end

# Setindex für BondView
function Base.setindex!(bv::BondView, val, j::Int)
    start_idx = bv.parent.cumBonds[bv.element_idx] + 1
    bv.parent.data[start_idx + j - 1] = val
end

# Vector Interface für BondView
Base.length(bv::BondView) = bv.parent.nBonds[bv.element_idx]
Base.size(bv::BondView) = (length(bv),)
Base.IndexStyle(::Type{<:BondView}) = IndexLinear()

# =============================================================================
# FreeSizeField Type
# =============================================================================
struct FreeSizeField{T} <: DataField{T}
    name::String
    vartype::Type{T}
    data::Vector{T}
    structured::AbstractArray{T}  # Gespeicherte Matrix-Form

    function FreeSizeField(name::String, vartype::Type{T}, value::T,
                           dof::Tuple{Vararg{Int64}}) where {T}
        if any(dof .< 1)
            @error "DOF smaller than one not allowed in free sized field Input: $dof."
        end
        data = fill(vartype(value), prod(dof))
        matrix_form = reshape(data, dof)
        new{T}(name, vartype, data, matrix_form)
    end
end
