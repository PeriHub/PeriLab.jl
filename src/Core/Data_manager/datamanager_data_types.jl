# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# =============================================================================
# ABSTRACT BASE TYPES
# =============================================================================

abstract type DataField{T} end

# =============================================================================
# NodeField Type
# =============================================================================

struct NodeField{T} <: DataField{T}
    name::String
    vartype::Type{T}
    data::Vector{T}
    nnodes::Int
    dof::Int
    matrix::Bool
    matrix_data::AbstractArray{T}  # view to matrix Form

    function NodeField(name::String, vartype::Type{T}, value::T, nnodes::Int, dof::Int,
                       matrix::Bool = false) where {T}
        if matrix
            data = fill(vartype(value), nnodes * dof * dof)
            matrix_form = reshape(data, nnodes, dof, dof)
            new{T}(name, vartype, data, nnodes, dof, matrix, matrix_form)
        else
            data = fill(vartype(value), nnodes * dof)
            if dof != 1
                matrix_form = reshape(data, nnodes, dof)
            else
                matrix_form = @views data
            end
            new{T}(name, vartype, data, nnodes, dof, matrix, matrix_form)
        end
    end
end

struct FreeSizeField{T} <: DataField{T}
    name::String
    vartype::Type{T}
    data::Vector{T}
    dof::Tuple{Vararg{Int64}}
    matrix_data::AbstractArray{T}  # Gespeicherte Matrix-Form

    function FreeSizeField(name::String, vartype::Type{T}, value::T,
                           dof::Tuple{Vararg{Int64}}) where {T}
        data = fill(vartype(value), prod(dof))
        matrix_form = reshape(data, dof)
        new{T}(name, vartype, data, dof, matrix_form)
    end
end
