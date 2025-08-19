# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# =============================================================================
# BondField Type
# =============================================================================
abstract type BondField{T} <: DataField{T} end
function BondField(name::String, ::Type{T}, value::T,
                   nBonds::Vector{Int}, dof::Int,
                   matrix_style::Bool = false) where {T<:Union{Float64,Int64,Bool}}
    if dof == 1
        ScalarBondField(name, T, value, nBonds)
    elseif !matrix_style
        VectorBondField(name, T, value, nBonds, dof)
    else
        MatrixBondField(name, T, value, nBonds, dof)
    end
end

# Forward declare ScalarBondView
struct ScalarBondView{T} <: AbstractVector{T}
    parent::BondField{T}
    element_idx::Int
end

struct ScalarBondField{T} <: BondField{T}
    name::String
    data::Vector{T}
    nBonds::Vector{Int}
    cumBonds::Vector{Int}
    element_views::Vector{ScalarBondView{T}}

    function ScalarBondField(name::String, ::Type{T}, value::T,
                             nBonds::Vector{Int}) where {T<:Union{Float64,Int64,Bool}}
        data = fill(value, sum(nBonds))
        cumBonds = [0; cumsum(nBonds)]
        element_views = ScalarBondView{T}[]

        obj = new{T}(name, data, nBonds, cumBonds, element_views)

        for i in 1:length(nBonds)
            push!(obj.element_views, ScalarBondView(obj, i))
        end

        return obj
    end
end

Base.getindex(s::ScalarBondField, i::Int) = s.element_views[i]
function Base.getindex(bv::ScalarBondView, j::Int)
    start_idx = bv.parent.cumBonds[bv.element_idx] + 1
    bv.parent.data[start_idx + j - 1]
end
function Base.setindex!(bv::ScalarBondView, val, j::Int)
    start_idx = bv.parent.cumBonds[bv.element_idx] + 1
    bv.parent.data[start_idx + j - 1] = val
end
# Vector Interface
Base.length(bv::ScalarBondView) = bv.parent.nBonds[bv.element_idx]
Base.size(bv::ScalarBondView) = (length(bv),)
Base.IndexStyle(::Type{<:ScalarBondView}) = IndexLinear()

###################
# Vector style
###################

struct VectorBondView{T} <: AbstractMatrix{T}
    parent::BondField{T}
    element_idx::Int
    cached_matrix::AbstractMatrix{T}
end

struct VectorBondField{T} <: BondField{T}
    name::String
    data::Vector{T}
    nBonds::Vector{Int}
    cumBonds::Vector{Int}
    dof::Int
    element_views::Vector{AbstractMatrix{T}}
    cached_views::Vector{VectorBondView{T}}

    function VectorBondField(name::String, ::Type{T}, value::T, nBonds::Vector{Int},
                             dof::Int) where {T<:Union{Float64,Int64,Bool}}
        data = fill(value, sum(nBonds) * dof)
        cumBonds = [0; cumsum(nBonds)]

        element_views = AbstractMatrix{T}[]
        for i in 1:length(nBonds)
            global_start = cumBonds[i] * dof + 1
            global_end = global_start + nBonds[i] * dof - 1
            data_slice = view(data, global_start:global_end)
            matrix_view = reshape(data_slice, nBonds[i], dof)
            push!(element_views, matrix_view)
        end

        cached_views = VectorBondView{T}[]
        obj = new{T}(name, data, nBonds, cumBonds, dof, element_views, cached_views)

        for i in 1:length(nBonds)
            push!(obj.cached_views, VectorBondView(obj, i, obj.element_views[i]))
        end

        return obj
    end
end

Base.getindex(s::VectorBondField, i::Int) = s.cached_views[i]
Base.getindex(vbv::VectorBondView, ::Colon, ::Colon) = vbv.cached_matrix
function Base.getindex(vbv::VectorBondView, bond::Int, dof_idx::Int)
    @inbounds vbv.cached_matrix[bond, dof_idx]
end
function Base.setindex!(vbv::VectorBondView, val, bond::Int, dof_idx::Int)
    @inbounds vbv.cached_matrix[bond, dof_idx] = val
end

function Base.getindex(vbv::VectorBondView, bond::Int, ::Colon)
    @inbounds vbv.cached_matrix[bond, :]
end
function Base.getindex(vbv::VectorBondView, ::Colon, dof_idx::Int)
    @inbounds vbv.cached_matrix[:, dof_idx]
end
function Base.setindex!(vbv::VectorBondView, val, bond::Int, ::Colon)
    @inbounds vbv.cached_matrix[bond, :] = val
end
function Base.setindex!(vbv::VectorBondView, val, ::Colon, dof_idx::Int)
    @inbounds vbv.cached_matrix[:, dof_idx] = val
end

# Size interface
Base.size(vbv::VectorBondView) = size(vbv.cached_matrix)
Base.size(vbv::VectorBondView, dim::Int) = size(vbv.cached_matrix, dim)
Base.length(vbv::VectorBondView) = length(vbv.cached_matrix)
Base.IndexStyle(::Type{<:VectorBondView}) = IndexCartesian()

# Iteration interface
Base.iterate(vbv::VectorBondView, state...) = iterate(vbv.cached_matrix, state...)

###################
# Matrix style
###################
struct MatrixBondView{T} <: AbstractArray{T,3}
    parent::BondField{T}
    element_idx::Int
    cached_tensor::AbstractArray{T,3}  # Precomputed 3D array!
end
struct MatrixBondField{T} <: BondField{T}
    name::String
    data::Vector{T}
    nBonds::Vector{Int}
    cumBonds::Vector{Int}
    dof::Int
    element_views::Vector{AbstractArray{T,3}}
    cached_views::Vector{MatrixBondView{T}}  # ✅ Precomputed Views!

    function MatrixBondField(name::String, ::Type{T}, value::T, nBonds::Vector{Int},
                             dof::Int) where {T<:Union{Float64,Int64,Bool}}
        data = fill(value, sum(nBonds) * dof * dof)
        cumBonds = [0; cumsum(nBonds)]

        # Precompute 3D views
        element_views = AbstractArray{T,3}[]
        for i in 1:length(nBonds)
            global_start = cumBonds[i] * dof * dof + 1
            global_end = global_start + nBonds[i] * dof * dof - 1
            data_slice = view(data, global_start:global_end)
            tensor_view = reshape(data_slice, nBonds[i], dof, dof)
            push!(element_views, tensor_view)
        end

        cached_views = MatrixBondView{T}[]
        obj = new{T}(name, data, nBonds, cumBonds, dof, element_views, cached_views)

        for i in 1:length(nBonds)
            push!(obj.cached_views, MatrixBondView(obj, i, obj.element_views[i]))
        end

        return obj
    end
end

Base.getindex(s::MatrixBondField, i::Int) = s.cached_views[i]

Base.getindex(mbv::MatrixBondView, ::Colon, ::Colon, ::Colon) = mbv.cached_tensor

function Base.getindex(mbv::MatrixBondView, bond::Int, row::Int, col::Int)
    @inbounds mbv.cached_tensor[bond, row, col]
end
function Base.setindex!(mbv::MatrixBondView, val, bond::Int, row::Int, col::Int)
    @inbounds mbv.cached_tensor[bond, row, col] = val
end

function Base.getindex(mbv::MatrixBondView, bond::Int, ::Colon, ::Colon)
    @inbounds mbv.cached_tensor[bond, :, :]
end
function Base.setindex!(mbv::MatrixBondView, val, bond::Int, ::Colon, ::Colon)
    @inbounds mbv.cached_tensor[bond, :, :] = val
end

function Base.getindex(mbv::MatrixBondView, bond::Int, row::Int, ::Colon)
    @inbounds mbv.cached_tensor[bond, row, :]
end
function Base.getindex(mbv::MatrixBondView, bond::Int, ::Colon, col::Int)
    @inbounds mbv.cached_tensor[bond, :, col]
end
function Base.setindex!(mbv::MatrixBondView, val, bond::Int, row::Int, ::Colon)
    @inbounds mbv.cached_tensor[bond, row, :] = val
end
function Base.setindex!(mbv::MatrixBondView, val, bond::Int, ::Colon, col::Int)
    @inbounds mbv.cached_tensor[bond, :, col] = val
end

function Base.getindex(mbv::MatrixBondView, ::Colon, row::Int, col::Int)
    @inbounds mbv.cached_tensor[:, row, col]
end
function Base.getindex(mbv::MatrixBondView, ::Colon, row::Int, ::Colon)
    @inbounds mbv.cached_tensor[:, row, :]
end
function Base.getindex(mbv::MatrixBondView, ::Colon, ::Colon, col::Int)
    @inbounds mbv.cached_tensor[:, :, col]
end
function Base.setindex!(mbv::MatrixBondView, val, ::Colon, row::Int, col::Int)
    @inbounds mbv.cached_tensor[:, row, col] = val
end
function Base.setindex!(mbv::MatrixBondView, val, ::Colon, row::Int, ::Colon)
    @inbounds mbv.cached_tensor[:, row, :] = val
end
function Base.setindex!(mbv::MatrixBondView, val, ::Colon, ::Colon, col::Int)
    @inbounds mbv.cached_tensor[:, :, col] = val
end

# Size interface
Base.size(mbv::MatrixBondView) = size(mbv.cached_tensor)
Base.size(mbv::MatrixBondView, dim::Int) = size(mbv.cached_tensor, dim)
Base.length(mbv::MatrixBondView) = length(mbv.cached_tensor)
Base.IndexStyle(::Type{<:MatrixBondView}) = IndexCartesian()

# Iteration interface
Base.iterate(mbv::MatrixBondView, state...) = iterate(mbv.cached_tensor, state...)
