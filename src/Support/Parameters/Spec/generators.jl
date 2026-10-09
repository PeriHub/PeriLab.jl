# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

export json_schema, type_schema, type_label, parameter_rows, nested_params, markdown_table

const _GLOBALS_SCHEMA = Dict{String,Any}("type" => "object")

"""
    type_schema(T)

JSON Schema of a declared field type. `Dependent` is a number or a data file
path; enums are strings (their spellings match loosely, so no JSON `enum`).
"""
function type_schema(T)
    T === Float64 && return Dict{String,Any}("type" => "number")
    T === Int64 && return Dict{String,Any}("type" => "integer")
    T === Bool && return Dict{String,Any}("type" => "boolean")
    T === String && return Dict{String,Any}("type" => "string")
    if T === Dependent
        return Dict{String,Any}("oneOf" => Any[Dict{String,Any}("type" => "number"),
                                               Dict{String,Any}("type" => "string",
                                                                "description" => "data file path")])
    end
    if T === Union{Float64,Vector{Float64}}
        return Dict{String,Any}("oneOf" => Any[Dict{String,Any}("type" => "number"),
                                               Dict{String,Any}("type" => "array",
                                                                "items" => Dict{String,Any}("type" => "number"))])
    end
    if _is_scalar_union(T)
        types = [type_schema(m)["type"] for m in Base.uniontypes(T) if m !== Nothing]
        return Dict{String,Any}("type" => length(types) == 1 ? only(types) : types)
    end
    T isa Union && return type_schema(_nonnothing(T))
    if T isa DataType && T <: Enum
        spellings = vcat(string.(instances(T)), collect(keys(enum_aliases(T))))
        return Dict{String,Any}("type" => "string",
                                "description" => "one of: " * join(spellings, ", ") *
                                                 " (case, spaces and punctuation are ignored)")
    end
    if T isa DataType && T <: Vector
        return Dict{String,Any}("type" => "array", "items" => type_schema(eltype(T)))
    end
    if T isa DataType && T <: Dict
        return Dict{String,Any}("type" => "object", "additionalProperties" => type_schema(valtype(T)))
    end
    is_params(T) && return json_schema(T)
    throw(ArgumentError("type_schema: unsupported type $T"))
end

# min / max on the numbers of a schema (also inside arrays and oneOf branches)
function _bounds!(s::Dict{String,Any}, lo, hi)
    if haskey(s, "oneOf")
        foreach(b -> _bounds!(b, lo, hi), s["oneOf"])
    elseif get(s, "type", "") == "array"
        _bounds!(s["items"], lo, hi)
    elseif get(s, "type", "") in ("number", "integer")
        lo === nothing || (s["minimum"] = lo)
        hi === nothing || (s["maximum"] = hi)
    end
    return s
end

_json_value(x::Union{Real,String}) = x
_json_value(x::Enum) = string(x)
_json_value(x::AbstractVector) = all(v -> _json_value(v) !== nothing, x) ?
                                 Any[_json_value(v) for v in x] : nothing
_json_value(x) = nothing

function _quantity_text(q)
    q === nothing && return ""
    return "quantity: $q (in your consistent unit system)"
end

function field_schema(fs::FieldSpec)
    s = type_schema(fs.type)
    _bounds!(s, fs.min, fs.max)
    fs.allowed === nothing || (s["enum"] = Any[_json_value(v) for v in fs.allowed])
    if !fs.required
        v = _json_value(fs.default)
        v === nothing || (s["default"] = v)
    end
    text = filter(!isempty, [fs.description, get(s, "description", ""), _quantity_text(fs.quantity)])
    isempty(text) || (s["description"] = join(text, "; "))
    return s
end

"""
    json_schema(T)

JSON Schema object of the `@params` struct `T`: its YAML keys with types,
defaults, bounds, allowed values and descriptions; `required` lists the
required keys; indexed keys (`key_patterns`) become `patternProperties`.
"""
function json_schema(T)
    props = Dict{String,Any}("Globals" => copy(_GLOBALS_SCHEMA))
    required = String[]
    for fs in parameter_spec(T)
        props[fs.alias] = field_schema(fs)
        fs.required && push!(required, fs.alias)
    end
    schema = Dict{String,Any}("type" => "object", "properties" => props,
                              "additionalProperties" => false)
    isempty(required) || (schema["required"] = sort!(required))
    patterns = key_patterns(T)
    if !isempty(patterns)
        schema["patternProperties"] = Dict{String,Any}(first(p).pattern => type_schema(last(p))
                                                       for p in patterns)
    end
    return schema
end

"Readable name of a declared field type (docs and `describe`)."
function type_label(T)
    T === Float64 && return "number"
    T === Int64 && return "integer"
    T === Bool && return "true/false"
    T === String && return "text"
    T === Dependent && return "number or data file"
    T === Union{Float64,Vector{Float64}} && return "number or list of numbers"
    if _is_scalar_union(T)
        return join([type_label(m) for m in Base.uniontypes(T) if m !== Nothing], " or ")
    end
    T isa Union && return type_label(_nonnothing(T))
    T isa DataType && T <: Enum && return "one of: " * join(string.(instances(T)), ", ")
    T isa DataType && T <: Vector && return "list of " * type_label(eltype(T))
    T isa DataType && T <: Dict && return "named entries of " * type_label(valtype(T))
    is_params(T) && return "section"
    return string(T)
end

_fmt_value(x::AbstractFloat) = isinteger(x) && abs(x) < 1e15 ? string(Int64(x)) : string(x)
_fmt_value(x) = string(x)

function _range_text(fs::FieldSpec)
    fs.allowed === nothing || return "one of: " * join(_fmt_value.(fs.allowed), ", ")
    fs.min !== nothing && fs.max !== nothing &&
        return _fmt_value(fs.min) * " … " * _fmt_value(fs.max)
    fs.min !== nothing && return "≥ " * _fmt_value(fs.min)
    fs.max !== nothing && return "≤ " * _fmt_value(fs.max)
    return ""
end

function _default_text(fs::FieldSpec)
    fs.required && return ""
    fs.default === nothing && return "—"
    fs.default isa AbstractDict && isempty(fs.default) && return "—"
    return _fmt_value(fs.default)
end

"""
    parameter_rows(T)

One row of strings per YAML key of `T`: key, type, required, default, range,
quantity, description.
"""
function parameter_rows(T)
    return [(key = fs.alias, type = type_label(fs.type),
             required = fs.required ? "required" : "optional",
             default = _default_text(fs), range = _range_text(fs),
             quantity = fs.quantity === nothing ? "" : string(fs.quantity),
             description = fs.description) for fs in parameter_spec(T)]
end

function _params_type(T)
    T isa Union && return _params_type(_nonnothing(T))
    T isa DataType && T <: Dict && return _params_type(valtype(T))
    T isa DataType && T <: Vector && return _params_type(eltype(T))
    return is_params(T) ? T : nothing
end

"The nested `@params` sections of `T` as alias => type, in field order."
function nested_params(T)
    nested = Pair{String,Any}[]
    for fs in parameter_spec(T)
        P = _params_type(fs.type)
        P === nothing || push!(nested, fs.alias => P)
    end
    return nested
end

_md_cell(s::AbstractString) = replace(s, "|" => "\\|", "\n" => " ")

"Markdown table of `parameter_rows`."
function markdown_table(rows)
    lines = ["| YAML key | Type | Required | Default | Range | Quantity | Description |",
             "|---|---|---|---|---|---|---|"]
    for r in rows
        push!(lines,
              "| " * join(_md_cell.([r.key, r.type, r.required, r.default, r.range,
                                     r.quantity, r.description]), " | ") * " |")
    end
    return join(lines, "\n") * "\n"
end
