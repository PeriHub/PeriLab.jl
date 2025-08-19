# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
# =============================================================================
# ElementField Type (analog zu BondField)
# =============================================================================

abstract type ElementField{T} <: DataField{T} end

# Forward declare ElementView
struct ElementView{T} <: AbstractVector{T}
    parent::ElementField{T}
    element_idx::Int
end

# Scalar ElementField (dof=1)
struct ScalarElementField{T} <: ElementField{T}
    name::String
    data::Vector{T}
    nElements::Vector{Int}  # Anzahl Neighbors pro Element
    cumElements::Vector{Int}
    element_views::Vector{ElementView{T}}

    function ScalarElementField(name::String, ::Type{T}, value::T,
                                nElements::Vector{Int}) where {T}
        data = fill(value, sum(nElements))
        cumElements = [0; cumsum(nElements)]
        element_views = ElementView{T}[]

        obj = new{T}(name, data, nElements, cumElements, element_views)

        # Precompute alle ElementViews
        for i in 1:length(nElements)
            push!(obj.element_views, ElementView(obj, i))
        end

        return obj
    end
end

struct VectorElementView{T} <: AbstractMatrix{T}
    parent::ElementField{T}
    element_idx::Int
    cached_matrix::AbstractMatrix{T}
end

# Vector ElementField (dof>1)
struct VectorElementField{T} <: ElementField{T}
    name::String
    data::Vector{T}
    nElements::Vector{Int}
    cumElements::Vector{Int}
    dof::Int
    element_views::Vector{AbstractMatrix{T}}
    cached_views::Vector{VectorElementView{T}}

    function VectorElementField(name::String, ::Type{T}, value::T, nElements::Vector{Int},
                                dof::Int) where {T}
        data = fill(value, sum(nElements) * dof)
        cumElements = [0; cumsum(nElements)]

        # Precompute matrix views
        element_views = AbstractMatrix{T}[]
        for i in 1:length(nElements)
            global_start = cumElements[i] * dof + 1
            global_end = global_start + nElements[i] * dof - 1
            data_slice = view(data, global_start:global_end)
            matrix_view = reshape(data_slice, nElements[i], dof)
            push!(element_views, matrix_view)
        end

        cached_views = VectorElementView{T}[]
        obj = new{T}(name, data, nElements, cumElements, dof, element_views, cached_views)

        # Precompute alle VectorElementViews
        for i in 1:length(nElements)
            push!(obj.cached_views, VectorElementView(obj, i, obj.element_views[i]))
        end

        return obj
    end
end

# =============================================================================
# INDEXING IMPLEMENTATION
# =============================================================================

# ScalarElementField indexing
Base.getindex(s::ScalarElementField, i::Int) = s.element_views[i]  # 0 bytes

function Base.getindex(ev::ElementView, j::Int)
    start_idx = ev.parent.cumElements[ev.element_idx] + 1
    @inbounds ev.parent.data[start_idx + j - 1]
end

function Base.setindex!(ev::ElementView, val, j::Int)
    start_idx = ev.parent.cumElements[ev.element_idx] + 1
    @inbounds ev.parent.data[start_idx + j - 1] = val
end

Base.length(ev::ElementView) = ev.parent.nElements[ev.element_idx]
Base.size(ev::ElementView) = (length(ev),)

# VectorElementField indexing
Base.getindex(s::VectorElementField, i::Int) = s.cached_views[i]  # 0 bytes

Base.getindex(vev::VectorElementView, ::Colon, ::Colon) = vev.cached_matrix
function Base.getindex(vev::VectorElementView, neighbor::Int, dof_idx::Int)
    @inbounds vev.cached_matrix[neighbor, dof_idx]
end
function Base.setindex!(vev::VectorElementView, val, neighbor::Int, dof_idx::Int)
    @inbounds vev.cached_matrix[neighbor, dof_idx] = val
end
function Base.getindex(vev::VectorElementView, neighbor::Int, ::Colon)
    @inbounds vev.cached_matrix[neighbor, :]
end
function Base.getindex(vev::VectorElementView, ::Colon, dof_idx::Int)
    @inbounds vev.cached_matrix[:, dof_idx]
end
function Base.setindex!(vev::VectorElementView, val, neighbor::Int, ::Colon)
    @inbounds vev.cached_matrix[neighbor, :] = val
end
function Base.setindex!(vev::VectorElementView, val, ::Colon, dof_idx::Int)
    @inbounds vev.cached_matrix[:, dof_idx] = val
end

Base.size(vev::VectorElementView) = size(vev.cached_matrix)
Base.length(vev::VectorElementView) = length(vev.cached_matrix)

# =============================================================================
# CONSTRUCTOR FUNCTION
# =============================================================================

function ElementField(name::String, ::Type{T}, value::T,
                      nElements::Vector{Int}, dof::Int) where {T<:Union{Float64,Int64,Bool}}
    if dof == 1
        ScalarElementField(name, T, value, nElements)
    else
        VectorElementField(name, T, value, nElements, dof)
    end
end
