# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

"YAML keys declared by an `@params` struct (`one_of` fields have none, see `one_of_keys`)."
aliases(T) = Set{String}(fs.alias for fs in parameter_spec(T) if fs.category === nothing)

"""
    build(T, dict, path, ctx; owner)

Builds an instance of the `@params` struct `T` from a YAML dict. Every field
is checked, so all problems are collected in `ctx`; returns `nothing` if any
field failed. Does not check for unknown keys (see `check_unknown!`), so that
composite models can share one dict. `owner` names the struct or model in
"missing" messages.
"""
function build(T, dict::AbstractDict, path::String, ctx::ParseContext;
               owner::String = string(nameof(T)))
    values = Any[]
    failed = false
    for fs in parameter_spec(T)
        if fs.category !== nothing
            v = _select_one_of(fs, dict, path, ctx)
            v === FAILED ? (failed = true) : push!(values, v)
            continue
        end
        field_path = join_path(path, fs.alias)
        if !haskey(dict, fs.alias)
            if fs.required
                add_error!(ctx, field_path, "missing (required by $owner)")
                failed = true
            else
                push!(values, _fresh(fs.default))
            end
            continue
        end
        v = convert_value(fs.type, dict[fs.alias], field_path, ctx; alias = fs.alias)
        if v === FAILED || !check_constraints!(fs, v, field_path, ctx)
            failed = true
        else
            push!(values, v)
        end
    end
    return failed ? nothing : T(values...)
end

"""
    check_unknown!(dict, known, path, ctx)

Reports keys of `dict` that are not in `known`: errors in strict mode,
warnings otherwise. A key that matches a known key apart from case, spaces and
punctuation is always an error, since its value would silently be replaced by
a default. `Globals` is an explicit escape hatch and never reported.
"""
function check_unknown!(dict::AbstractDict, known, path::String, ctx::ParseContext)
    for k in keys(dict)
        key = string(k)
        (key in known || key == "Globals") && continue
        message = unknown_key_message(key, known)
        key_path = join_path(path, key)
        near_match = any(alias -> _normalize(alias) == _normalize(key), known)
        if ctx.strict || near_match
            add_error!(ctx, key_path, message)
        else
            add_warning!(ctx, key_path, message)
        end
    end
    return nothing
end

"""
    parse_section(T, dict, path, ctx)

Builds `T` from `dict`, reports unknown keys, and runs `derive`. Keys in
`known_extra` are read by the caller itself: they are not reported as unknown
and are offered as suggestions for misspelled keys.
"""
function parse_section(T, dict::AbstractDict, path::String, ctx::ParseContext;
                       known_extra = ())
    # invokelatest: `T` may come from a module loaded at runtime (licensed modules)
    return Base.invokelatest(_parse_section, T, dict, path, ctx, known_extra)
end

function _parse_section(T, dict::AbstractDict, path::String, ctx::ParseContext,
                        known_extra = ())
    p = build(T, dict, path, ctx)
    check_unknown!(dict, union(aliases(T), one_of_keys(T), known_extra), path, ctx)
    p === nothing && return nothing
    check_required_groups!(T, dict, path, ctx)
    check!(p, path, ctx)
    return derive(p)
end

# a key of `dict` with a value (an empty section or `key:` without value counts as missing)
_given(dict::AbstractDict, key::String) = (v = get(dict, key, nothing);
                                           v !== nothing && !(v isa AbstractDict && isempty(v)))

# `when` of a required group holds for `dict` (values match like enums, ignoring case,
# spaces and punctuation)
function _holds(when::Union{Nothing,AbstractDict}, dict::AbstractDict)
    when === nothing && return true
    return all(when) do (key, values)
        _given(dict, key) || return nothing in values
        given = filter(!isnothing, values)
        v = dict[key]
        return isempty(given) ||
               any(x -> x == v || (x isa AbstractString && v isa AbstractString &&
                                   _normalize(x) == _normalize(v)), given)
    end
end

const _COUNT_WORDS = ("one", "two", "three", "four")
_count_word(n::Int) = n <= length(_COUNT_WORDS) ? _COUNT_WORDS[n] : string(n)

_quoted_list(keys) = length(keys) == 1 ? "\"$(only(keys))\"" :
                     join(["\"$k\"" for k in keys[1:(end - 1)]], ", ") * " or \"$(keys[end])\""

"""
    check_required_groups!(T, dict, path, ctx)

Reports the `required_groups` of `T` that `dict` breaks. In a `from_mesh` group a
key also counts as given if the mesh defines a node field of its name with `_`
for spaces (`ctx.node_fields`), e.g. `Bulk_Modulus`; such a group is not checked
if the node fields are not known yet (`nothing`, e.g. Exodus meshes).
"""
function check_required_groups!(T, dict::AbstractDict, path::String, ctx::ParseContext)
    for g in required_groups(T)
        _holds(g.when, dict) || continue
        g.from_mesh && ctx.node_fields === nothing && continue
        n = count(k -> _given(dict, k) ||
                       (g.from_mesh && replace(k, " " => "_") in ctx.node_fields), g.keys)
        if n < g.min
            add_error!(ctx, path,
                       g.min == 1 ? "$(_quoted_list(g.keys)) is required" :
                       "$(_count_word(g.min)) of $(join(["\"$k\"" for k in g.keys], ", ")) are required")
        elseif g.max !== nothing && n > g.max
            add_error!(ctx, path,
                       g.max == 1 ? "give only one of $(_quoted_list(g.keys))" :
                       "at most $(_count_word(g.max)) of $(_quoted_list(g.keys)) may be given")
        end
    end
    return nothing
end
