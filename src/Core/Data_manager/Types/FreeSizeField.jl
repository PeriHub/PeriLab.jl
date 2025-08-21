# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

struct FreeSizeField{T,N} <: DataField{T}
    name::String
    data::Array{T,N}

    function FreeSizeField(name::String, ::Type{T}, value::T, dims::Vararg{Int}) where {T}
        data = fill(value, dims)
        N = length(dims)
        new{T,N}(name, data)
    end
end

Base.size(fsf::FreeSizeField) = size(fsf.data)
Base.getindex(fsf::FreeSizeField, i...) = getindex(fsf.data, i...)
Base.setindex!(fsf::FreeSizeField, v, i...) = setindex!(fsf.data, v, i...)
Base.IndexStyle(::Type{<:FreeSizeField}) = IndexLinear()
