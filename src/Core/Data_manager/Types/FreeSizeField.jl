# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

struct FreeSizeField{T,N} <: DataField{T}
    name::String
    data::Array{T,N}

    function FreeSizeField(name::String, ::Type{T}, value::T,
                           dims::Vararg{Int}) where {T<:Union{Float64,Int64,Bool}}
        data = Array{T,N}(fill(value, dims))
        N = length(dims)
        new{T,N}(name, data)
    end
end
