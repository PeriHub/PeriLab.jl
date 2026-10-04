# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

export @params, derive, check!

const _SELF = @__MODULE__

"`true` for structs declared with `@params`."
is_params(::Type) = false

"""
    parameter_spec(T) -> Vector{FieldSpec}

Field declarations of an `@params` struct, in field order.
"""
function parameter_spec end

"""
    derive(p) -> p′

Hook run once after a parameter struct is built from the input deck. Return a
struct of the same struct type (type parameters may differ) with derived or
normalized values filled in. The default returns `p` unchanged.
"""
derive(p) = p

"""
    check!(p, path, ctx)

Hook for cross-field rules of an `@params` struct (e.g. "exactly one solver",
"Low must not exceed High"). Runs once after `p` was built successfully and
before `derive`; add problems with `add_error!(ctx, path, msg)`. The default
does nothing.
"""
check!(p, path::String, ctx::ParseContext) = nothing

const _SCALAR_UNION_MEMBERS = (Nothing, Int64, Float64, String, Bool)

function _is_scalar_union(T)
    T isa Union || return false
    members = Base.uniontypes(T)
    return all(m -> m in _SCALAR_UNION_MEMBERS, members) &&
           !(Int64 in members && Float64 in members) &&
           any(m -> m !== Nothing, members)
end

_typename(T) = replace(string(T), r"(\w+\.)+" => "")

const _SCALAR_TYPES = (Float64, Int64, Bool, String)
const _VECTOR_TYPES = (Vector{Float64}, Vector{Int64}, Vector{String})

function supported_type(T)
    _is_scalar_union(T) && return true
    if T isa Union
        S = _nonnothing(T)
        return Nothing <: T && !(S isa Union) && supported_type(S)
    end
    (T in _SCALAR_TYPES || T in _VECTOR_TYPES || T === Dependent) && return true
    T isa DataType || return false
    T <: Enum && return true
    if T <: Dict && T.parameters[1] === String
        V = T.parameters[2]
        return V in _SCALAR_TYPES || _is_scalar_union(V) || (V isa DataType && is_params(V))
    end
    return is_params(T)
end

function _numeric_type(T)
    T isa Union && return _numeric_type(_nonnothing(T))
    return T in (Float64, Int64, Dependent, Vector{Float64}, Vector{Int64})
end

"""
    build_spec(T, display, entries)

Turns the `(field name, declared type, FieldDecl)` entries emitted by `@params`
into `FieldSpec`s, raising `ParamsDefinitionError` for invalid declarations.
"""
function build_spec(T, display::String, entries::Vector{Any})
    specs = FieldSpec[]
    used = Dict{String,Symbol}()
    for (fname, ftype, decl) in entries
        location = "$display.$fname"
        if ftype isa UnionAll && is_params(ftype)
            throw(ParamsDefinitionError("$location: nested section type $(_typename(ftype)) contains Dependent fields; this is not supported"))
        end
        supported_type(ftype) ||
            throw(ParamsDefinitionError("$location: unsupported field type $(_typename(ftype)). Supported: Float64, Int64, Bool, String, an @enum, Vector{Float64}, Vector{Int64}, Vector{String}, Dependent, Union{Nothing,T}, a union of Int64/Float64/String/Bool/Nothing (not Int64 and Float64 together), a nested @params struct, or Dict{String,V} of a nested @params struct or a scalar type"))
        if haskey(used, decl.alias)
            throw(ParamsDefinitionError("$location: YAML key \"$(decl.alias)\" is already used by field `$(used[decl.alias])`"))
        end
        used[decl.alias] = fname
        if (decl.min !== nothing || decl.max !== nothing) && !_numeric_type(ftype)
            throw(ParamsDefinitionError("$location: min/max are only allowed on numeric fields"))
        end
        default = NO_DEFAULT
        if !decl.required
            if decl.default === NO_DEFAULT
                throw(ParamsDefinitionError("$location: opt(\"$(decl.alias)\") needs a default, e.g. opt(\"$(decl.alias)\"; default = ...). Use a Union{Nothing,T} field with default = nothing if the value may be absent"))
            end
            ctx = ParseContext()
            v = convert_value(ftype, decl.default, "", ctx; alias = decl.alias)
            v === FAILED || check_constraints!(FieldSpec(fname, ftype, decl), v, "", ctx)
            if has_errors(ctx)
                throw(ParamsDefinitionError("$location: default $(_describe(decl.default)) is invalid: $(ctx.errors[1].message)"))
            end
            default = v
        end
        push!(specs, FieldSpec(fname, ftype, decl; default = default))
    end
    return specs
end

function _is_dependent_type(t)
    return t === :Dependent ||
           (t isa Expr && t.head === :. && t.args[end] == QuoteNode(:Dependent))
end

"`Union{Nothing,Dependent}` in a field declaration (either order)."
function _is_optional_dependent_type(t)
    (t isa Expr && t.head === :curly && t.args[1] === :Union && length(t.args) == 3) ||
        return false
    members = t.args[2:3]
    return any(==(:Nothing), members) && any(_is_dependent_type, members)
end

_definition_error(msg) = :(throw($ParamsDefinitionError($msg)))

function _field_usage(name, fname)
    return "$name: field `$fname` must be written as `$fname::Type = req(\"YAML key\"; ...)` or `$fname::Type = opt(\"YAML key\"; default = ...)`"
end

"""
    @params struct Name
        field::Type = req("YAML key"; min, max, allowed, quantity, description)
        field::Type = opt("YAML key"; default, min, max, allowed, quantity, description)
    end

Declares a parameter struct. Generates the plain immutable struct (type
parameters are added for `Dependent` fields so that every instance is a
concrete type) and its `parameter_spec`. Invalid declarations raise a
`ParamsDefinitionError` naming the struct and field.
"""
macro params(structdef)
    if !(structdef isa Expr && structdef.head === :struct)
        return _definition_error("@params must be applied to a struct definition")
    end
    if structdef.args[1]
        return _definition_error("@params structs must be immutable: use `struct`, not `mutable struct`")
    end
    name = structdef.args[2]
    supertype = nothing
    if name isa Expr && name.head === :<: && name.args[1] isa Symbol
        name, supertype = name.args[1], name.args[2]
    end
    if !(name isa Symbol)
        return _definition_error("@params struct $(name): write a plain name without type parameters")
    end
    fields = Any[]
    typeparams = Any[]
    entries = Any[]
    for line in structdef.args[3].args
        if line isa LineNumberNode
            push!(fields, line)
            continue
        end
        line isa AbstractString && continue
        if !(line isa Expr && line.head === :(=) && line.args[1] isa Expr &&
             line.args[1].head === :(::) && length(line.args[1].args) == 2)
            fname = line isa Expr && line.head === :(::) ? line.args[1] : line
            return _definition_error(_field_usage(name, fname))
        end
        fname, ftype = line.args[1].args
        decl = line.args[2]
        if !(decl isa Expr && decl.head === :call && decl.args[1] in (:req, :opt))
            return _definition_error(_field_usage(name, fname))
        end
        call = Expr(:call, GlobalRef(_SELF, decl.args[1]), decl.args[2:end]...)
        if _is_dependent_type(ftype) || _is_optional_dependent_type(ftype)
            bound = _is_dependent_type(ftype) ? Dependent : Union{Nothing,Dependent}
            typeparam = Symbol("T_", fname)
            push!(typeparams, Expr(:<:, typeparam, bound))
            push!(fields, Expr(:(::), fname, typeparam))
            push!(entries, Expr(:tuple, QuoteNode(fname), bound, call))
        else
            push!(fields, Expr(:(::), fname, ftype))
            push!(entries, Expr(:tuple, QuoteNode(fname), ftype, call))
        end
    end
    head = isempty(typeparams) ? name : Expr(:curly, name, typeparams...)
    supertype === nothing || (head = Expr(:<:, head, supertype))
    specname = Symbol("__params_spec_", name)
    structexpr = Expr(:struct, false, head, Expr(:block, fields...))
    return esc(quote
                   Base.@__doc__ $structexpr
                   const $specname = $(_SELF).build_spec($name, $(string(name)),
                                                         Any[$(entries...)])
                   $(_SELF).parameter_spec(::Type{<:$name}) = $specname
                   $(_SELF).is_params(::Type{<:$name}) = true
                   nothing
               end)
end
