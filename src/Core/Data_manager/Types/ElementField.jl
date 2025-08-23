# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
# =============================================================================
# ElementField Type (analog zu BondField)
# =============================================================================
abstract type ElementField{T} <: DataField{T} end
function ElementField(name::String,
                      ::Type{T},
                      value::T,
                      nElements::Vector{Int64},
                      dof::Int64) where {T<:Union{Float64,Int64,Bool}}
    if dof == 1
        ScalarElementField(name, T, value, nBonds)
    else
        VectorElementField(name, T, value, nBonds, dof)
    end
end

struct ScalarElementField{T} <: ElementField{T}
    name::String
    data::Vector{Vector{T}}

    function ScalarElementField(name::String, ::Type{T},
                                value::T,
                                nElements::Vector{Int64}) where {T<:Union{Float64,Int64,
                                                                          Bool}}
        data = Vector{Vector{T}}([fill(value, n) for n in nElements])
        new{T}(name, data)
    end
end

struct VectorElementField{T} <: ElementField{T}
    name::String
    data::Vector{T}

    function VectorElementField(name::String, ::Type{T}, value::T,
                                nElements::Vector{Int64},
                                dof::Int64) where {T<:Union{Float64,Int64,Bool}}
        data = Vector{Matrix{T}}([fill(value, (n, dof)) for n in nElements])
        new{T}(name, data)
    end
end
