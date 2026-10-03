# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

export req, opt

struct NoDefault end
const NO_DEFAULT = NoDefault()

"""
    FieldDecl

What a module author wrote for one field via `req(...)` or `opt(...)`.
"""
struct FieldDecl
    alias::String
    required::Bool
    default::Any
    min::Union{Nothing,Float64}
    max::Union{Nothing,Float64}
    allowed::Union{Nothing,Vector{Any}}
    quantity::Union{Nothing,Symbol}
    description::String
end

_bound(::Nothing) = nothing
_bound(x::Real) = Float64(x)
_allowed(::Nothing) = nothing
_allowed(values) = Vector{Any}(collect(values))

"""
    req(yaml_key; min, max, allowed, quantity, description)

A required parameter. `quantity` (e.g. `:stress`) is documentation only;
PeriLab has no fixed unit system.
"""
function req(alias::AbstractString; min = nothing, max = nothing, allowed = nothing,
             quantity = nothing, description::AbstractString = "")
    return FieldDecl(String(alias), true, NO_DEFAULT, _bound(min), _bound(max),
                     _allowed(allowed), quantity, String(description))
end

"""
    opt(yaml_key; default, min, max, allowed, quantity, description)

An optional parameter; `default` is required.
"""
function opt(alias::AbstractString; default = NO_DEFAULT, min = nothing, max = nothing,
             allowed = nothing, quantity = nothing, description::AbstractString = "")
    return FieldDecl(String(alias), false, default, _bound(min), _bound(max),
                     _allowed(allowed), quantity, String(description))
end

"""
    FieldSpec

A validated field declaration: struct field name, declared type, and the
`FieldDecl` metadata. `default` is already converted to `type`.
"""
struct FieldSpec
    name::Symbol
    type::Any
    alias::String
    required::Bool
    default::Any
    min::Union{Nothing,Float64}
    max::Union{Nothing,Float64}
    allowed::Union{Nothing,Vector{Any}}
    quantity::Union{Nothing,Symbol}
    description::String
end

function FieldSpec(name::Symbol, type, decl::FieldDecl; default = decl.default)
    return FieldSpec(name, type, decl.alias, decl.required, default, decl.min, decl.max,
                     decl.allowed, decl.quantity, decl.description)
end
