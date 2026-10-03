# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

"YAML keys declared by an `@params` struct."
aliases(T) = Set{String}(fs.alias for fs in parameter_spec(T))

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

Builds `T` from `dict`, reports unknown keys, and runs `derive`.
"""
function parse_section(T, dict::AbstractDict, path::String, ctx::ParseContext)
    # invokelatest: `T` may come from a module loaded at runtime (licensed modules)
    return Base.invokelatest(_parse_section, T, dict, path, ctx)
end

function _parse_section(T, dict::AbstractDict, path::String, ctx::ParseContext)
    p = build(T, dict, path, ctx)
    check_unknown!(dict, aliases(T), path, ctx)
    return p === nothing ? nothing : derive(p)
end
