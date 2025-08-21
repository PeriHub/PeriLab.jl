# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

abstract type NodeField{T} <: DataField{T} end

# =============================================================================
# ScalarNodeField (dof=1) - Single value per node
# =============================================================================

struct ScalarNodeField{T} <: AbstractVector{T}
    name::String
    data::Vector{T}
    function ScalarNodeField(name::String, ::Type{T}, value::T, nnodes::Int) where {T}
        new{T}(name, fill(value, nnodes))
    end
end

# Delegiere alles an den Vector
Base.size(snf::ScalarNodeField) = size(snf.data)
Base.getindex(snf::ScalarNodeField, i...) = getindex(snf.data, i...)
Base.setindex!(snf::ScalarNodeField, v, i...) = setindex!(snf.data, v, i...)
Base.IndexStyle(::Type{<:ScalarNodeField}) = IndexLinear()

# =============================================================================
# VectorNodeField (dof>1, !matrix) - Vector per node
# =============================================================================

struct VectorNodeField{T} <: AbstractMatrix{T}
    name::String
    data::Matrix{T}  # Direkt eine Matrix!

    function VectorNodeField(name::String, ::Type{T}, value::T, nnodes::Int,
                             dof::Int) where {T}
        data = fill(value, nnodes, dof)
        new{T}(name, data)
    end
end

Base.size(vnf::VectorNodeField) = size(vnf.data)
Base.getindex(vnf::VectorNodeField, i...) = getindex(vnf.data, i...)
Base.setindex!(vnf::VectorNodeField, v, i...) = setindex!(vnf.data, v, i...)
Base.IndexStyle(::Type{<:VectorNodeField}) = IndexLinear()

nnodes(vnf::VectorNodeField) = size(vnf.data, 1)
dof(vnf::VectorNodeField) = size(vnf.data, 2)

# =============================================================================
# MatrixNodeField (matrix=true) - Matrix per node
# =============================================================================
struct MatrixNodeField{T} <: AbstractArray{T,3}
    name::String
    data::Array{T,3}  # Direkt ein 3D-Array!

    function MatrixNodeField(name::String, ::Type{T}, value::T, nnodes::Int,
                             dof::Int) where {T}
        data = fill(value, nnodes, dof, dof)
        new{T}(name, data)
    end
end

Base.size(mnf::MatrixNodeField) = size(mnf.data)
Base.getindex(mnf::MatrixNodeField, i...) = getindex(mnf.data, i...)
Base.setindex!(mnf::MatrixNodeField, v, i...) = setindex!(mnf.data, v, i...)
Base.IndexStyle(::Type{<:MatrixNodeField}) = IndexLinear()

nnodes(mnf::MatrixNodeField) = size(mnf.data, 1)
dof(mnf::MatrixNodeField) = size(mnf.data, 2)
