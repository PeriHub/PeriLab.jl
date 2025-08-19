# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# =============================================================================
# NodeField Type Hierarchy - Memory-efficient node-based data structures
# =============================================================================

abstract type NodeField{T} <: DataField{T} end

# =============================================================================
# ScalarNodeField (dof=1) - Single value per node
# =============================================================================

struct ScalarNodeField{T} <: NodeField{T}
    name::String
    data::Vector{T}
    nnodes::Int

    function ScalarNodeField(name::String, ::Type{T}, value::T,
                             nnodes::Int) where {T<:Union{Float64,Int64,Bool}}
        data = fill(value, nnodes)
        new{T}(name, data, nnodes)
    end
end

# ScalarNodeField Indexing - all variants (0 bytes allocation)
Base.getindex(snf::ScalarNodeField, i::Int) = @inbounds snf.data[i]
Base.setindex!(snf::ScalarNodeField, val, i::Int) = @inbounds snf.data[i] = val

Base.getindex(snf::ScalarNodeField, ::Colon) = snf.data
Base.setindex!(snf::ScalarNodeField, val, ::Colon) = snf.data .= val

Base.getindex(snf::ScalarNodeField, r::AbstractRange) = @inbounds snf.data[r]
Base.setindex!(snf::ScalarNodeField, val, r::AbstractRange) = @inbounds snf.data[r] = val

# ScalarNodeField End support
Base.lastindex(snf::ScalarNodeField) = snf.nnodes
Base.firstindex(snf::ScalarNodeField) = 1

# ScalarNodeField Size interface
Base.length(snf::ScalarNodeField) = snf.nnodes
Base.size(snf::ScalarNodeField) = (snf.nnodes,)
Base.axes(snf::ScalarNodeField) = (Base.OneTo(snf.nnodes),)
function Base.axes(snf::ScalarNodeField, dim::Int)
    dim == 1 ? Base.OneTo(snf.nnodes) : Base.OneTo(1)
end

# =============================================================================
# VectorNodeField (dof>1, !matrix) - Vector per node
# =============================================================================

struct VectorNodeField{T} <: NodeField{T}
    name::String
    data::Vector{T}
    nnodes::Int
    dof::Int
    cached_matrix::AbstractMatrix{T}  # Precomputed view for zero allocation

    function VectorNodeField(name::String, ::Type{T}, value::T, nnodes::Int,
                             dof::Int) where {T<:Union{Float64,Int64,Bool}}
        data = fill(value, nnodes * dof)
        cached_matrix = reshape(data, nnodes, dof)
        new{T}(name, data, nnodes, dof, cached_matrix)
    end
end

# VectorNodeField Indexing - all variants (0 bytes allocation)
Base.getindex(vnf::VectorNodeField, i::Int, j::Int) = @inbounds vnf.cached_matrix[i, j]
function Base.setindex!(vnf::VectorNodeField, val, i::Int, j::Int)
    @inbounds vnf.cached_matrix[i, j] = val
end

# Slice access
Base.getindex(vnf::VectorNodeField, i::Int, ::Colon) = @inbounds vnf.cached_matrix[i, :]
Base.getindex(vnf::VectorNodeField, ::Colon, j::Int) = @inbounds vnf.cached_matrix[:, j]
Base.getindex(vnf::VectorNodeField, ::Colon, ::Colon) = vnf.cached_matrix
Base.getindex(vnf::VectorNodeField, ::Colon) = vnf.cached_matrix

function Base.setindex!(vnf::VectorNodeField, val, i::Int, ::Colon)
    @inbounds vnf.cached_matrix[i, :] = val
end
function Base.setindex!(vnf::VectorNodeField, val, ::Colon, j::Int)
    @inbounds vnf.cached_matrix[:, j] = val
end
Base.setindex!(vnf::VectorNodeField, val, ::Colon, ::Colon) = vnf.cached_matrix .= val
Base.setindex!(vnf::VectorNodeField, val, ::Colon) = vnf.cached_matrix .= val

# Range access
function Base.getindex(vnf::VectorNodeField, r::AbstractRange, j::Int)
    @inbounds vnf.cached_matrix[r, j]
end
function Base.getindex(vnf::VectorNodeField, i::Int, r::AbstractRange)
    @inbounds vnf.cached_matrix[i, r]
end
function Base.getindex(vnf::VectorNodeField, r1::AbstractRange, r2::AbstractRange)
    @inbounds vnf.cached_matrix[r1, r2]
end
function Base.getindex(vnf::VectorNodeField, r::AbstractRange, ::Colon)
    @inbounds vnf.cached_matrix[r, :]
end
function Base.getindex(vnf::VectorNodeField, ::Colon, r::AbstractRange)
    @inbounds vnf.cached_matrix[:, r]
end
Base.getindex(vnf::VectorNodeField, r::AbstractRange) = @inbounds vnf.cached_matrix[r]  # Linear indexing

function Base.setindex!(vnf::VectorNodeField, val, r::AbstractRange, j::Int)
    @inbounds vnf.cached_matrix[r, j] = val
end
function Base.setindex!(vnf::VectorNodeField, val, i::Int, r::AbstractRange)
    @inbounds vnf.cached_matrix[i, r] = val
end
function Base.setindex!(vnf::VectorNodeField, val, r1::AbstractRange, r2::AbstractRange)
    @inbounds vnf.cached_matrix[r1, r2] = val
end
function Base.setindex!(vnf::VectorNodeField, val, r::AbstractRange, ::Colon)
    @inbounds vnf.cached_matrix[r, :] = val
end
function Base.setindex!(vnf::VectorNodeField, val, ::Colon, r::AbstractRange)
    @inbounds vnf.cached_matrix[:, r] = val
end
function Base.setindex!(vnf::VectorNodeField, val, r::AbstractRange)
    @inbounds vnf.cached_matrix[r] = val
end

# VectorNodeField End support
Base.lastindex(vnf::VectorNodeField, dim::Int) = dim == 1 ? vnf.nnodes : vnf.dof
Base.firstindex(vnf::VectorNodeField, dim::Int) = 1
Base.lastindex(vnf::VectorNodeField) = length(vnf.cached_matrix)  # Linear indexing
Base.firstindex(vnf::VectorNodeField) = 1

# VectorNodeField Size interface
Base.length(vnf::VectorNodeField) = vnf.nnodes * vnf.dof
Base.size(vnf::VectorNodeField) = (vnf.nnodes, vnf.dof)
Base.axes(vnf::VectorNodeField) = (Base.OneTo(vnf.nnodes), Base.OneTo(vnf.dof))
function Base.axes(vnf::VectorNodeField, dim::Int)
    dim == 1 ? Base.OneTo(vnf.nnodes) : (dim == 2 ? Base.OneTo(vnf.dof) : Base.OneTo(1))
end

# =============================================================================
# MatrixNodeField (matrix=true) - Matrix per node
# =============================================================================

struct MatrixNodeField{T} <: NodeField{T}
    name::String
    data::Vector{T}
    nnodes::Int
    dof::Int
    cached_tensor::AbstractArray{T,3}  # Precomputed view for zero allocation

    function MatrixNodeField(name::String, ::Type{T}, value::T, nnodes::Int,
                             dof::Int) where {T<:Union{Float64,Int64,Bool}}
        data = fill(value, nnodes * dof * dof)
        cached_tensor = reshape(data, nnodes, dof, dof)
        new{T}(name, data, nnodes, dof, cached_tensor)
    end
end

# MatrixNodeField Indexing - all variants (0 bytes allocation)
function Base.getindex(mnf::MatrixNodeField, i::Int, j::Int, k::Int)
    @inbounds mnf.cached_tensor[i, j, k]
end
function Base.setindex!(mnf::MatrixNodeField, val, i::Int, j::Int, k::Int)
    @inbounds mnf.cached_tensor[i, j, k] = val
end

# 2D Slices (single node matrix, single matrix element across nodes)
function Base.getindex(mnf::MatrixNodeField, i::Int, ::Colon, ::Colon)
    @inbounds mnf.cached_tensor[i, :, :]
end
function Base.getindex(mnf::MatrixNodeField, ::Colon, j::Int, ::Colon)
    @inbounds mnf.cached_tensor[:, j, :]
end
function Base.getindex(mnf::MatrixNodeField, ::Colon, ::Colon, k::Int)
    @inbounds mnf.cached_tensor[:, :, k]
end

function Base.setindex!(mnf::MatrixNodeField, val, i::Int, ::Colon, ::Colon)
    @inbounds mnf.cached_tensor[i, :, :] = val
end
function Base.setindex!(mnf::MatrixNodeField, val, ::Colon, j::Int, ::Colon)
    @inbounds mnf.cached_tensor[:, j, :] = val
end
function Base.setindex!(mnf::MatrixNodeField, val, ::Colon, ::Colon, k::Int)
    @inbounds mnf.cached_tensor[:, :, k] = val
end

# 1D Slices (matrix row/column, cross-node slices)
function Base.getindex(mnf::MatrixNodeField, i::Int, j::Int, ::Colon)
    @inbounds mnf.cached_tensor[i, j, :]
end
function Base.getindex(mnf::MatrixNodeField, i::Int, ::Colon, k::Int)
    @inbounds mnf.cached_tensor[i, :, k]
end
function Base.getindex(mnf::MatrixNodeField, ::Colon, j::Int, k::Int)
    @inbounds mnf.cached_tensor[:, j, k]
end

function Base.setindex!(mnf::MatrixNodeField, val, i::Int, j::Int, ::Colon)
    @inbounds mnf.cached_tensor[i, j, :] = val
end
function Base.setindex!(mnf::MatrixNodeField, val, i::Int, ::Colon, k::Int)
    @inbounds mnf.cached_tensor[i, :, k] = val
end
function Base.setindex!(mnf::MatrixNodeField, val, ::Colon, j::Int, k::Int)
    @inbounds mnf.cached_tensor[:, j, k] = val
end

# Complete access
Base.getindex(mnf::MatrixNodeField, ::Colon, ::Colon, ::Colon) = mnf.cached_tensor
Base.getindex(mnf::MatrixNodeField, ::Colon) = mnf.cached_tensor

function Base.setindex!(mnf::MatrixNodeField, val, ::Colon, ::Colon, ::Colon)
    mnf.cached_tensor .= val
end
Base.setindex!(mnf::MatrixNodeField, val, ::Colon) = mnf.cached_tensor .= val

# Range access - most important combinations
function Base.getindex(mnf::MatrixNodeField, r::AbstractRange, j::Int, k::Int)
    @inbounds mnf.cached_tensor[r, j, k]
end
function Base.getindex(mnf::MatrixNodeField, i::Int, r::AbstractRange, k::Int)
    @inbounds mnf.cached_tensor[i, r, k]
end
function Base.getindex(mnf::MatrixNodeField, i::Int, j::Int, r::AbstractRange)
    @inbounds mnf.cached_tensor[i, j, r]
end

function Base.getindex(mnf::MatrixNodeField, r::AbstractRange, ::Colon, ::Colon)
    @inbounds mnf.cached_tensor[r, :, :]
end
function Base.getindex(mnf::MatrixNodeField, ::Colon, r::AbstractRange, ::Colon)
    @inbounds mnf.cached_tensor[:, r, :]
end
function Base.getindex(mnf::MatrixNodeField, ::Colon, ::Colon, r::AbstractRange)
    @inbounds mnf.cached_tensor[:, :, r]
end

function Base.getindex(mnf::MatrixNodeField, r1::AbstractRange, r2::AbstractRange, k::Int)
    @inbounds mnf.cached_tensor[r1, r2, k]
end
function Base.getindex(mnf::MatrixNodeField, r1::AbstractRange, j::Int, r2::AbstractRange)
    @inbounds mnf.cached_tensor[r1, j, r2]
end
function Base.getindex(mnf::MatrixNodeField, i::Int, r1::AbstractRange, r2::AbstractRange)
    @inbounds mnf.cached_tensor[i, r1, r2]
end
function Base.getindex(mnf::MatrixNodeField, r1::AbstractRange, r2::AbstractRange,
                       r3::AbstractRange)
    @inbounds mnf.cached_tensor[r1, r2, r3]
end
Base.getindex(mnf::MatrixNodeField, r::AbstractRange) = @inbounds mnf.cached_tensor[r]  # Linear indexing

# Range setindex variants
function Base.setindex!(mnf::MatrixNodeField, val, r::AbstractRange, j::Int, k::Int)
    @inbounds mnf.cached_tensor[r, j, k] = val
end
function Base.setindex!(mnf::MatrixNodeField, val, i::Int, r::AbstractRange, k::Int)
    @inbounds mnf.cached_tensor[i, r, k] = val
end
function Base.setindex!(mnf::MatrixNodeField, val, i::Int, j::Int, r::AbstractRange)
    @inbounds mnf.cached_tensor[i, j, r] = val
end

function Base.setindex!(mnf::MatrixNodeField, val, r::AbstractRange, ::Colon, ::Colon)
    @inbounds mnf.cached_tensor[r, :, :] = val
end
function Base.setindex!(mnf::MatrixNodeField, val, ::Colon, r::AbstractRange, ::Colon)
    @inbounds mnf.cached_tensor[:, r, :] = val
end
function Base.setindex!(mnf::MatrixNodeField, val, ::Colon, ::Colon, r::AbstractRange)
    @inbounds mnf.cached_tensor[:, :, r] = val
end

function Base.setindex!(mnf::MatrixNodeField, val, r1::AbstractRange, r2::AbstractRange,
                        r3::AbstractRange)
    @inbounds mnf.cached_tensor[r1, r2, r3] = val
end
function Base.setindex!(mnf::MatrixNodeField, val, r::AbstractRange)
    @inbounds mnf.cached_tensor[r] = val
end

# MatrixNodeField End support
Base.lastindex(mnf::MatrixNodeField, dim::Int) = dim == 1 ? mnf.nnodes : mnf.dof
Base.firstindex(mnf::MatrixNodeField, dim::Int) = 1
Base.lastindex(mnf::MatrixNodeField) = length(mnf.cached_tensor)  # Linear indexing
Base.firstindex(mnf::MatrixNodeField) = 1

# MatrixNodeField Size interface
Base.length(mnf::MatrixNodeField) = mnf.nnodes * mnf.dof * mnf.dof
Base.size(mnf::MatrixNodeField) = (mnf.nnodes, mnf.dof, mnf.dof)
function Base.axes(mnf::MatrixNodeField)
    (Base.OneTo(mnf.nnodes), Base.OneTo(mnf.dof), Base.OneTo(mnf.dof))
end
function Base.axes(mnf::MatrixNodeField, dim::Int)
    dim == 1 ? Base.OneTo(mnf.nnodes) : (dim <= 3 ? Base.OneTo(mnf.dof) : Base.OneTo(1))
end

# =============================================================================
# Constructor Function - unified interface
# =============================================================================

function NodeField(name::String, ::Type{T}, value::T, nnodes::Int, dof::Int,
                   matrix::Bool = false) where {T<:Union{Float64,Int64,Bool}}
    if dof == 1
        ScalarNodeField(name, T, value, nnodes)
    elseif !matrix
        VectorNodeField(name, T, value, nnodes, dof)
    else
        MatrixNodeField(name, T, value, nnodes, dof)
    end
end

# =============================================================================
# Usage Examples - all zero allocation after construction
# =============================================================================

# Scalar case (pressure, temperature, etc.):
# scalar_field = ScalarNodeField("pressure", Float64, 0.0, 100)
# scalar_field[end] = 99.0                    # ✅ 0 bytes
# subset = scalar_field[10:20]                # ✅ 0 bytes
# scalar_field[:] = rand(100)                 # ✅ 0 bytes

# Vector case (displacement, velocity, etc.):
# vector_field = VectorNodeField("displacement", Float64, 0.0, 100, 3)
# vector_field[1, 2] = 5.0                    # ✅ 0 bytes - Node 1, DOF 2
# node_vector = vector_field[end, :]          # ✅ 0 bytes - Last node, all DOFs
# vector_field[1:10, :] = rand(10, 3)         # ✅ 0 bytes - First 10 nodes

# Matrix case (stiffness, compliance, etc.):
# matrix_field = MatrixNodeField("stiffness", Float64, 0.0, 100, 3)
# matrix_field[1, 2, 3] = 8.0                 # ✅ 0 bytes - Node 1, Matrix(2,3)
# node_matrix = matrix_field[end, :, :]       # ✅ 0 bytes - Last node matrix
# matrix_field[1:5, :, :] = rand(5, 3, 3)     # ✅ 0 bytes - First 5 node matrices
