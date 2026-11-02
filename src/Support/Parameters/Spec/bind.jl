# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

"""
    bind_dependents!(x, lookup, path, ctx)

Binds every `Table1D` inside `x` (an `@params` struct, a `Composite`, or a
dict of them) to its node field. `lookup(field_name)` returns the field's
`Vector{Float64}` or `nothing`. Call after the fields exist; problems are added
to `ctx`. A table keeps a reference to the array it was bound to, so call it
again whenever the field array is replaced (e.g. after the N/NP1 switch).
"""
function bind_dependents!(x, lookup, path::String, ctx::ParseContext)
    # invokelatest: `x` may come from a module loaded at runtime (licensed modules)
    return Base.invokelatest(_bind!, x, lookup, path, ctx)
end

function _bind!(x, lookup, path::String, ctx::ParseContext)
    is_params(typeof(x)) || return nothing
    for fs in parameter_spec(typeof(x))
        _bind!(getfield(x, fs.name), lookup, join_path(path, fs.alias), ctx)
    end
    return nothing
end

function _bind!(t::Table1D, lookup, path::String, ctx::ParseContext)
    field = lookup(t.field_name)
    if field === nothing
        add_error!(ctx, path, "field \"$(t.field_name)\" required by $(t.source) does not exist")
    elseif !(field isa Vector{Float64})
        add_error!(ctx, path,
                   "field \"$(t.field_name)\" must be a per-node Vector{Float64}, got $(typeof(field))")
    else
        bind_table!(t, field)
    end
    return nothing
end

function _bind!(c::Composite, lookup, path::String, ctx::ParseContext)
    for part in c.parts
        _bind!(part, lookup, path, ctx)
    end
    return nothing
end

function _bind!(d::AbstractDict, lookup, path::String, ctx::ParseContext)
    for (k, v) in d
        _bind!(v, lookup, join_path(path, string(k)), ctx)
    end
    return nothing
end
