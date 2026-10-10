# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

"Returned by `convert_value` when conversion failed; the error is already in the context."
struct Failed end
const FAILED = Failed()

function _fail(ctx::ParseContext, path::AbstractString, msg::AbstractString)
    add_error!(ctx, path, msg)
    return FAILED
end

_describe(::Nothing) = "an empty value"
_describe(raw::AbstractString) = "\"$raw\""
_describe(raw) = repr(raw)

_fresh(x) = (x isa AbstractArray || x isa AbstractDict) ? copy(x) : x

_nonnothing(T::Union) = T.a === Nothing ? T.b : T.a

_fmt(x::Real) = isinteger(x) ? string(Int(x)) : string(x)

"""
    enum_aliases(::Type{E})

Extra YAML spellings for an `@enum`, for values that are not valid Julia
names. Extend it in the module defining the enum, e.g.
`ParameterSpec.enum_aliases(::Type{Symmetry}) = Dict("3D" => Full3D)`.
Without an alias, a YAML string matches an enum instance if both are equal
ignoring case, spaces and punctuation ("plane stress" matches `PlaneStress`).
"""
enum_aliases(::Type{E}) where {E<:Enum} = Dict{String,E}()

"""
    convert_value(T, raw, path, ctx; alias = "")

Converts a raw YAML value to the declared type `T`. On failure the error is
added to `ctx` and `FAILED` is returned. `alias` is the YAML key, needed to
find the column in a `Dependent` data file.
"""
function convert_value(T, raw, path::String, ctx::ParseContext; alias::String = "")
    raw isa T && return _fresh(raw)
    if T === Union{Float64,Vector{Float64}}
        raw isa Real && !(raw isa Bool) && return Float64(raw)
        raw isa AbstractVector && return _convert_vector(Vector{Float64}, raw, path, ctx)
        return _fail(ctx, path, "expected a number or a list of numbers, got $(_describe(raw))")
    end
    if _is_scalar_union(T)
        return _convert_scalar_union(T, raw, path, ctx)
    elseif T isa Union
        # `Section: false` switches an optional section off
        raw === false && is_params(_nonnothing(T)) && return nothing
        return convert_value(_nonnothing(T), raw, path, ctx; alias = alias)
    elseif T === Float64
        raw isa Real && !(raw isa Bool) && return Float64(raw)
        if raw isa AbstractString && endswith(lowercase(raw), ".txt")
            return _fail(ctx, path,
                         "expected a number, got file path \"$raw\" — this parameter does not support dependent values")
        end
        return _fail(ctx, path, "expected a number, got $(_describe(raw))")
    elseif T === Int64
        !(raw isa Bool) && _fits_int64(raw) && return Int64(raw)
        return _fail(ctx, path, "expected an integer, got $(_describe(raw))")
    elseif T === Bool
        return _fail(ctx, path, "expected true or false, got $(_describe(raw))")
    elseif T === String
        raw isa AbstractString && return String(raw)
        return _fail(ctx, path, "expected text, got $(_describe(raw))")
    elseif T === Dependent
        raw isa Real && !(raw isa Bool) && return Constant(Float64(raw))
        if raw isa AbstractString
            table = read_table(joinpath(ctx.directory, raw), alias, path, ctx)
            return table === nothing ? FAILED : table
        end
        return _fail(ctx, path, "expected a number or a data file path, got $(_describe(raw))")
    elseif T isa DataType && T <: Enum
        return _convert_enum(T, raw, path, ctx)
    elseif T isa DataType && T <: Vector
        return _convert_vector(T, raw, path, ctx)
    elseif T isa DataType && T <: Dict && T.parameters[1] === String &&
           supported_type(T)
        return _convert_named_sections(T, raw, path, ctx)
    elseif is_params(T)
        raw isa AbstractDict ||
            return _fail(ctx, path,
                         "expected a section of `key: value` entries, got $(_describe(raw))")
        section = parse_section(T, raw, path, ctx)
        return section === nothing ? FAILED : section
    end
    throw(ArgumentError("convert_value: unsupported type $T"))
end

function _convert_named_sections(::Type{Dict{String,V}}, raw, path::String,
                                 ctx::ParseContext) where {V}
    raw isa AbstractDict ||
        return _fail(ctx, path, "expected named entries, got $(_describe(raw))")
    out = Dict{String,V}()
    ok = true
    for (k, item) in raw
        string(k) == "Globals" && continue
        v = convert_value(V, item, join_path(path, string(k)), ctx)
        if v === FAILED
            ok = false
        else
            out[string(k)] = v
        end
    end
    return ok ? out : FAILED
end

function _convert_enum(::Type{E}, raw, path::String, ctx::ParseContext) where {E<:Enum}
    aliases = enum_aliases(E)
    if raw isa AbstractString
        haskey(aliases, raw) && return aliases[raw]
        key = _normalize(raw)
        for instance in instances(E)
            _normalize(string(instance)) == key && return instance
        end
    end
    names = vcat([string(instance) for instance in instances(E)], collect(keys(aliases)))
    return _fail(ctx, path, "$(_describe(raw)) is not one of: $(join(names, ", "))")
end

function _convert_vector(::Type{Vector{S}}, raw, path::String, ctx::ParseContext) where {S}
    raw isa AbstractVector || return _fail(ctx, path, "expected a list, got $(_describe(raw))")
    out = Vector{S}(undef, length(raw))
    ok = true
    for (i, item) in enumerate(raw)
        v = convert_value(S, item, "$path[$i]", ctx)
        if v === FAILED
            ok = false
        else
            out[i] = v
        end
    end
    return ok ? out : FAILED
end

_numbers(v::Bool) = ()
_numbers(v::Real) = (v,)
_numbers(v::AbstractVector{<:Real}) = v
_numbers(v::Constant) = (v.value,)
_numbers(v::Table1D) = v.y
_numbers(v) = ()

"""
    check_constraints!(fs, v, path, ctx)

Checks `min` / `max` (on every number in `v`, including all table values) and
`allowed`. Adds the first violation to `ctx` and returns `false`.
"""
function check_constraints!(fs::FieldSpec, v, path::String, ctx::ParseContext)
    for x in _numbers(v)
        if fs.min !== nothing && x < fs.min
            add_error!(ctx, path, "$x is below minimum $(_fmt(fs.min))")
            return false
        end
        if fs.max !== nothing && x > fs.max
            add_error!(ctx, path, "$x is above maximum $(_fmt(fs.max))")
            return false
        end
    end
    if fs.allowed !== nothing && v !== nothing && !(v in fs.allowed)
        add_error!(ctx, path,
                   "$(_describe(v)) is not one of: $(join(_describe.(fs.allowed), ", "))")
        return false
    end
    return true
end

const _SCALAR_KIND_NAMES = Dict{Any,String}(Int64 => "an integer", Float64 => "a number",
                                            String => "text", Bool => "true or false")

# an integer value (also a float like 100.0) in the range of Int64
_fits_int64(raw::Integer) = typemin(Int64) <= raw <= typemax(Int64)
_fits_int64(raw::AbstractFloat) = isinteger(raw) && -2.0^63 <= raw < 2.0^63
_fits_int64(raw) = false

function _convert_scalar_union(T, raw, path::String, ctx::ParseContext)
    members = Base.uniontypes(T)
    if raw isa Bool
        Bool in members && return raw
    elseif raw isa Integer
        Int64 in members && _fits_int64(raw) && return Int64(raw)
        Float64 in members && return Float64(raw)
    elseif raw isa AbstractFloat
        Float64 in members && return Float64(raw)
        Int64 in members && _fits_int64(raw) && return Int64(raw)
    elseif raw isa AbstractString
        String in members && return String(raw)
    end
    # fixed order, independent of how Julia orders union members
    expected = join([_SCALAR_KIND_NAMES[m] for m in (Int64, Float64, String, Bool)
                     if m in members], " or ")
    return _fail(ctx, path, "expected $expected, got $(_describe(raw))")
end
