# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

"""
    BondField(name, T, value, nBonds, dof, matrix_style=false)

Creates optimized bond data structures based on DOF and matrix_style.

# Arguments
- `name::String`: Field name
- `T::Type`: Data type (Float64, Int64, Bool)
- `value::T`: Initialization value
- `nBonds::Vector{Int}`: Number of bonds per point
- `dof::Int`: Degrees of Freedom
- `matrix_style::Bool`: Whether to use matrix layout

# Returns
- `ScalarBondField` for dof=1
- `VectorBondField` for dof>1 of type [nnodes][nBonds[node],dofs]
- `MatrixBondField` for dof>1 of size [nnodes][nBonds[node],dofs,dofs]
"""

abstract type BondField{T} <: DataField{T} end

function BondField(name::String,
                   ::Type{T},
                   value::T,
                   nBonds::Vector{Int64},
                   dof::Int64,
                   matrix_style::Bool = false) where {T<:Union{Float64,Int64,Bool}}
    if dof == 1
        ScalarBondField(name, T, value, nBonds)
    elseif !matrix_style
        VectorBondField(name, T, value, nBonds, dof)
    else
        MatrixBondField(name, T, value, nBonds, dof)
    end
end

struct ScalarBondField{T} <: BondField{T}
    name::String
    data::Vector{Vector{T}}

    function ScalarBondField(name::String, ::Type{T}, value::T,
                             nBonds::Vector{Int64}) where {T<:Union{Float64,Int64,Bool}}
        data = Vector{Vector{T}}([fill(value, n) for n in nBonds])
        new{T}(name, data)
    end
end

struct VectorBondField{T} <: BondField{T}
    name::String
    data::Vector{T}

    function VectorBondField(name::String, ::Type{T}, value::T,
                             nBonds::Vector{Int64},
                             dof::Int64) where {T<:Union{Float64,Int64,Bool}}
        data = Vector{Matrix{T}}([fill(value, (n, dof)) for n in nBonds])
        new{T}(name, data)
    end
end

struct MatrixBondField{T} <: BondField{T}
    name::String
    data::Vector{Array{T,3}}

    function MatrixBondField(name::String, ::Type{T}, value::T,
                             nBonds::Vector{Int64},
                             dof::Int64) where {T<:Union{Float64,Int64,Bool}}
        data = Vector{Array{T,3}}([fill(value, (n, dof, dof)) for n in nBonds])
        new{T}(name, data)
    end
end
