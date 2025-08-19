# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
struct FreeSizeField{T} <: DataField{T}
    name::String
    vartype::Type{T}
    data::Vector{T}
    structured::AbstractArray{T}

    function FreeSizeField(name::String, vartype::Type{T}, value::T,
                           dof::Tuple{Vararg{Int64}}) where {T<:Union{Float64,Int64,Bool}}
        if any(dof .< 1)
            @error "DOF smaller than one not allowed in free sized field Input: $dof."
        end
        data = fill(vartype(value), prod(dof))
        matrix_form = reshape(data, dof)
        new{T}(name, vartype, data, matrix_form)
    end
end

# =============================================================================
# ALLOCATION FREE INDEXING IMPLEMENTATION
# =============================================================================

# 1D - Scalar indexing
Base.getindex(fsf::FreeSizeField, i::Int) = @inbounds fsf.structured[i]
Base.setindex!(fsf::FreeSizeField, val, i::Int) = @inbounds fsf.structured[i] = val

# 1D - Slice indexing
Base.getindex(fsf::FreeSizeField, ::Colon) = fsf.structured
Base.setindex!(fsf::FreeSizeField, val, ::Colon) = fsf.structured .= val

# 2D - Scalar indexing
Base.getindex(fsf::FreeSizeField, i::Int, j::Int) = @inbounds fsf.structured[i, j]
function Base.setindex!(fsf::FreeSizeField, val, i::Int, j::Int)
    @inbounds fsf.structured[i, j] = val
end

# 2D - Slice indexing
Base.getindex(fsf::FreeSizeField, i::Int, ::Colon) = @inbounds fsf.structured[i, :]
Base.getindex(fsf::FreeSizeField, ::Colon, j::Int) = @inbounds fsf.structured[:, j]
Base.getindex(fsf::FreeSizeField, ::Colon, ::Colon) = fsf.structured
function Base.setindex!(fsf::FreeSizeField, val, i::Int, ::Colon)
    @inbounds fsf.structured[i, :] = val
end
function Base.setindex!(fsf::FreeSizeField, val, ::Colon, j::Int)
    @inbounds fsf.structured[:, j] = val
end
Base.setindex!(fsf::FreeSizeField, val, ::Colon, ::Colon) = fsf.structured .= val

# 3D - Scalar indexing
function Base.getindex(fsf::FreeSizeField, i::Int, j::Int, k::Int)
    @inbounds fsf.structured[i, j, k]
end
function Base.setindex!(fsf::FreeSizeField, val, i::Int, j::Int, k::Int)
    @inbounds fsf.structured[i, j, k] = val
end

# 3D - Slice indexing (alle Kombinationen)
function Base.getindex(fsf::FreeSizeField, i::Int, j::Int, ::Colon)
    @inbounds fsf.structured[i, j, :]
end
function Base.getindex(fsf::FreeSizeField, i::Int, ::Colon, k::Int)
    @inbounds fsf.structured[i, :, k]
end
function Base.getindex(fsf::FreeSizeField, ::Colon, j::Int, k::Int)
    @inbounds fsf.structured[:, j, k]
end
function Base.getindex(fsf::FreeSizeField, i::Int, ::Colon, ::Colon)
    @inbounds fsf.structured[i, :, :]
end
function Base.getindex(fsf::FreeSizeField, ::Colon, j::Int, ::Colon)
    @inbounds fsf.structured[:, j, :]
end
function Base.getindex(fsf::FreeSizeField, ::Colon, ::Colon, k::Int)
    @inbounds fsf.structured[:, :, k]
end
Base.getindex(fsf::FreeSizeField, ::Colon, ::Colon, ::Colon) = fsf.structured

function Base.setindex!(fsf::FreeSizeField, val, i::Int, j::Int, ::Colon)
    @inbounds fsf.structured[i, j, :] = val
end
function Base.setindex!(fsf::FreeSizeField, val, i::Int, ::Colon, k::Int)
    @inbounds fsf.structured[i, :, k] = val
end
function Base.setindex!(fsf::FreeSizeField, val, ::Colon, j::Int, k::Int)
    @inbounds fsf.structured[:, j, k] = val
end
function Base.setindex!(fsf::FreeSizeField, val, i::Int, ::Colon, ::Colon)
    @inbounds fsf.structured[i, :, :] = val
end
function Base.setindex!(fsf::FreeSizeField, val, ::Colon, j::Int, ::Colon)
    @inbounds fsf.structured[:, j, :] = val
end
function Base.setindex!(fsf::FreeSizeField, val, ::Colon, ::Colon, k::Int)
    @inbounds fsf.structured[:, :, k] = val
end
Base.setindex!(fsf::FreeSizeField, val, ::Colon, ::Colon, ::Colon) = fsf.structured .= val

# 4D - Scalar indexing
function Base.getindex(fsf::FreeSizeField, i::Int, j::Int, k::Int, l::Int)
    @inbounds fsf.structured[i, j, k, l]
end
function Base.setindex!(fsf::FreeSizeField, val, i::Int, j::Int, k::Int, l::Int)
    @inbounds fsf.structured[i, j, k, l] = val
end

# 4D - Häufigste Slice-Kombinationen
function Base.getindex(fsf::FreeSizeField, i::Int, ::Colon, ::Colon, ::Colon)
    @inbounds fsf.structured[i, :, :, :]
end
function Base.getindex(fsf::FreeSizeField, ::Colon, ::Colon, ::Colon, l::Int)
    @inbounds fsf.structured[:, :, :, l]
end
Base.getindex(fsf::FreeSizeField, ::Colon, ::Colon, ::Colon, ::Colon) = fsf.structured

function Base.setindex!(fsf::FreeSizeField, val, i::Int, ::Colon, ::Colon, ::Colon)
    @inbounds fsf.structured[i, :, :, :] = val
end
function Base.setindex!(fsf::FreeSizeField, val, ::Colon, ::Colon, ::Colon, l::Int)
    @inbounds fsf.structured[:, :, :, l] = val
end
function Base.setindex!(fsf::FreeSizeField, val, ::Colon, ::Colon, ::Colon, ::Colon)
    fsf.structured .= val
end

# Size delegation
Base.size(fsf::FreeSizeField) = size(fsf.structured)
Base.size(fsf::FreeSizeField, dim::Int) = size(fsf.structured, dim)
Base.length(fsf::FreeSizeField) = length(fsf.structured)

# =============================================================================
# For higher dimensions but not memory free
# =============================================================================

Base.getindex(fsf::FreeSizeField, args...) = @inbounds getindex(fsf.structured, args...)
function Base.setindex!(fsf::FreeSizeField, val, args...)
    @inbounds setindex!(fsf.structured, val, args...)
end
