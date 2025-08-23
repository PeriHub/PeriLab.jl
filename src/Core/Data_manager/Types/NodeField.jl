# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
abstract type NodeField{T} <: DataField{T} end

function NodeField(name::String, ::Type{T}, value::T,
                   nnodes::Int64,
                   dof::Int64,
                   matrix_style::Bool = false) where {T<:Union{Float64,Int64,Bool}}
    if dof == 1
        ScalarNodeField(name, T, value, nnodes)
    elseif !matrix_style
        VectorNodeField(name, T, value, nnodes, dof)
    else
        MatrixNodeField(name, T, value, nnodes, dof)
    end
end

# =============================================================================
# ScalarNodeField (dof=1) - Single value per node
# =============================================================================

struct ScalarNodeField{T} <: NodeField{T}
    name::String
    data::Vector{T}
    function ScalarNodeField(name::String, ::Type{T}, value::T,
                             nnodes::Int64) where {T<:Union{Float64,Int64,Bool}}
        new{T}(name, fill(value, nnodes))
    end
end

# =============================================================================
# VectorNodeField (dof>1, !matrix) - Vector per node
# =============================================================================

struct VectorNodeField{T} <: NodeField{T}
    name::String
    data::Matrix{T}  # Direkt eine Matrix!

    function VectorNodeField(name::String, ::Type{T}, value::T, nnodes::Int64,
                             dof::Int64) where {T<:Union{Float64,Int64,Bool}}
        data = fill(value, nnodes, dof)
        new{T}(name, data)
    end
end

# =============================================================================
# MatrixNodeField (matrix=true) - Matrix per node
# =============================================================================
struct MatrixNodeField{T} <: NodeField{T}
    name::String
    data::Array{T,3}  # Direkt ein 3D-Array!

    function MatrixNodeField(name::String, ::Type{T}, value::T, nnodes::Int64,
                             dof::Int64) where {T<:Union{Float64,Int64,Bool}}
        data = fill(value, nnodes, dof, dof)
        new{T}(name, data)
    end
end
