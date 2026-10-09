# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

export NoModel, Composite, WithBase, dependent_tables

"Placeholder for a block that has no model of a category."
struct NoModel end

"Models combined with `+` in the input deck; each part keeps its own struct."
struct Composite{P<:Tuple}
    parts::P
end

"A model together with the base part of its category (see `register_base!`)."
struct WithBase{B,M}
    base::B
    model::M
    extras::Dict{String,Any}   # values of indexed keys (see `key_patterns`), e.g. Property_1
    name::String               # the model name as written, e.g. "A + B"
end
WithBase(base, model) = WithBase(base, model, Dict{String,Any}(), "")
WithBase(base, model, extras) = WithBase(base, model, extras, "")

function _check_key_patterns!(dict::AbstractDict, known::Set{String}, types, path::String,
                              ctx::ParseContext)
    extras = Dict{String,Any}()
    for k in keys(dict)
        key = string(k)
        key in known && continue
        for T in types
            pattern = findfirst(p -> occursin(first(p), key), key_patterns(T))
            pattern === nothing && continue
            v = convert_value(last(key_patterns(T)[pattern]), dict[k], join_path(path, key),
                              ctx)
            v === FAILED || (extras[key] = v)
            push!(known, key)
            break
        end
    end
    return extras
end

function _check_alias_conflicts!(names, types, path::String, ctx::ParseContext)
    seen = Dict{String,Tuple{String,Any}}()
    ok = true
    for (name, T) in zip(names, types), fs in parameter_spec(T)
        if haskey(seen, fs.alias)
            other_name, other_type = seen[fs.alias]
            if other_type !== fs.type
                add_error!(ctx, join_path(path, fs.alias),
                           "declared as $(_typename(other_type)) by \"$other_name\" but as $(_typename(fs.type)) by \"$name\"; models combined with + must agree")
                ok = false
            end
        else
            seen[fs.alias] = (name, fs.type)
        end
    end
    return ok
end

"""
    parse_model(category, dict, path, ctx; name_key)

Builds the model(s) named by `dict[name_key]` (e.g. "Material Model"), which
may combine several registered models with `+`. All parts read their keys from
the same `dict`; a key is unknown only if no part declares it. Returns
`NoModel()` for `dict === nothing`, the single model struct, a `Composite`, or
`nothing` if there were errors.
"""
function parse_model(category::Symbol, dict::Union{Nothing,AbstractDict}, path::String,
                     ctx::ParseContext; name_key::String)
    # invokelatest: models may come from modules loaded at runtime (licensed modules)
    return Base.invokelatest(_parse_model, category, dict, path, ctx, name_key)
end

function _parse_model(category::Symbol, dict::Union{Nothing,AbstractDict}, path::String,
                      ctx::ParseContext, name_key::String)
    dict === nothing && return NoModel()
    name_path = join_path(path, name_key)
    raw = get(dict, name_key, nothing)
    if !(raw isa AbstractString)
        add_error!(ctx, name_path,
                   raw === nothing ? "missing (names the model to use)" :
                   "expected a model name, got $(_describe(raw))")
        return nothing
    end
    names = String.(strip.(split(raw, "+")))
    if any(isempty, names)
        add_error!(ctx, name_path, "empty model name in \"$raw\"")
        return nothing
    end
    types = Any[]
    for name in names
        entry = lookup_model(category, name)
        if entry === nothing
            suggestion = suggest(name, registered_names(category))
            add_error!(ctx, name_path,
                       suggestion === nothing ?
                       "model \"$name\" not found; it may require a licensed module, or its module does not register its parameters" :
                       "model \"$name\" not found — did you mean \"$suggestion\"?")
        elseif entry isa UnavailableModel
            add_error!(ctx, name_path, "model \"$name\" $(entry.reason)")
        else
            push!(types, entry)
        end
    end
    length(types) == length(names) || return nothing
    base = base_model(category)
    all_names = base === nothing ? names : ["base parameters"; names]
    all_types = base === nothing ? types : Any[base; types]
    _check_alias_conflicts!(all_names, all_types, path, ctx) || return nothing
    base_part = nothing
    if base !== nothing
        base_part = build(base, dict, path, ctx; owner = "every $category model")
        base_part === nothing || check!(base_part, path, ctx)
        base_part = base_part === nothing ? nothing : derive(base_part)
    end
    parts = Any[]
    for (name, T) in zip(names, types)
        part = build(T, dict, path, ctx; owner = name)
        part === nothing || check!(part, path, ctx)
        push!(parts, part === nothing ? nothing : derive(part))
    end
    known = Set{String}([name_key])
    for T in all_types
        union!(known, aliases(T))
    end
    extras = _check_key_patterns!(dict, known, all_types, path, ctx)
    check_unknown!(dict, known, path, ctx)
    any(isnothing, parts) && return nothing
    base !== nothing && base_part === nothing && return nothing
    model = length(parts) == 1 ? parts[1] : Composite(Tuple(parts))
    return base === nothing ? model : WithBase(base_part, model, extras, String(strip(raw)))
end

"""
    dependent_tables(x)

The `Table1D` parameters of an `@params` struct or of every part of a `Composite`.
"""
function dependent_tables(x)
    found = Table1D[]
    for part in (x isa Composite ? x.parts : (x,))
        part === nothing && continue
        for fs in parameter_spec(typeof(part))
            v = getfield(part, fs.name)
            v isa Table1D && push!(found, v)
        end
    end
    return found
end
