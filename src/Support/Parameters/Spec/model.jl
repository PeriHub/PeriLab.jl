# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

export NoModel, Composite

"Placeholder for a block that has no model of a category."
struct NoModel end

"Models combined with `+` in the input deck; each part keeps its own struct."
struct Composite{P<:Tuple}
    parts::P
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
                       "model \"$name\" not found; it may require a licensed module" :
                       "model \"$name\" not found — did you mean \"$suggestion\"?")
        elseif entry isa UnavailableModel
            add_error!(ctx, name_path, "model \"$name\" $(entry.reason)")
        else
            push!(types, entry)
        end
    end
    length(types) == length(names) || return nothing
    _check_alias_conflicts!(names, types, path, ctx) || return nothing
    parts = Any[]
    for (name, T) in zip(names, types)
        part = build(T, dict, path, ctx; owner = name)
        part === nothing || check!(part, path, ctx)
        push!(parts, part === nothing ? nothing : derive(part))
    end
    known = Set{String}([name_key])
    for T in types
        union!(known, aliases(T))
    end
    check_unknown!(dict, known, path, ctx)
    any(isnothing, parts) && return nothing
    return length(parts) == 1 ? parts[1] : Composite(Tuple(parts))
end
