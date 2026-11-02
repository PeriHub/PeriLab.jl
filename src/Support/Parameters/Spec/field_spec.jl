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
    category::Union{Nothing,Symbol}
    evaluated::Bool
    when::Union{Nothing,Dict{String,Vector{String}}}
end

_bound(::Nothing) = nothing
_bound(x::Real) = Float64(x)
_allowed(::Nothing) = nothing
_allowed(values) = Vector{Any}(collect(values))
_when(::Nothing) = nothing
function _when(d::AbstractDict)
    Dict{String,Vector{String}}(String(k) => String[String(v) for v in vs]
                                for (k, vs) in d)
end

"""
    req(yaml_key; min, max, allowed, quantity, description, evaluated, when)

A required parameter. `quantity` (e.g. `:stress`) is documentation only;
PeriLab has no fixed unit system. `evaluated = true` marks a string PeriLab
evaluates as Julia (an expression); the schema flags it with `x-perilab-eval`
so front ends can guard it.

`when = Dict("Symmetry" => ["orthotropic", ...])` documents that the key is read
only if every listed sibling key has one of the listed values (exact, case-
insensitive; an unset key counts with its default, without one the condition
holds); the pseudo key `"\$solver"` is the solver chosen in the same solver
step. The schema emits it as `x-perilab-when`; parsing ignores it.
"""
function req(alias::AbstractString; min = nothing, max = nothing, allowed = nothing,
             quantity = nothing, description::AbstractString = "", evaluated::Bool = false,
             when = nothing)
    return FieldDecl(String(alias), true, NO_DEFAULT, _bound(min), _bound(max),
                     _allowed(allowed), quantity, String(description), nothing, evaluated,
                     _when(when))
end

"""
    opt(yaml_key; default, min, max, allowed, quantity, description, evaluated, when)

An optional parameter; `default` is required.
"""
function opt(alias::AbstractString; default = NO_DEFAULT, min = nothing, max = nothing,
             allowed = nothing, quantity = nothing, description::AbstractString = "",
             evaluated::Bool = false, when = nothing)
    return FieldDecl(String(alias), false, default, _bound(min), _bound(max),
                     _allowed(allowed), quantity, String(description), nothing, evaluated,
                     _when(when))
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
    category::Union{Nothing,Symbol}
    evaluated::Bool
    when::Union{Nothing,Dict{String,Vector{String}}}
end

const DOCS_URL = "https://perihub.github.io/PeriLab.jl/stable/"

"""
    docs_links(description)

Expands markdown links into the PeriLab docs, `[theory](@docs theory/theory_plastic/)`,
to `DOCS_URL` links, so they work in the docs, the schema and `describe` alike.
"""
docs_links(s::AbstractString) = replace(s, "](@docs " => "](" * DOCS_URL)

function FieldSpec(name::Symbol, type, decl::FieldDecl; default = decl.default)
    return FieldSpec(name, type, decl.alias, decl.required, default, decl.min, decl.max,
                     decl.allowed, decl.quantity, docs_links(decl.description),
                     decl.category,
                     decl.evaluated, decl.when)
end
