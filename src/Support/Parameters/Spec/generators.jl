# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

export json_schema, type_schema, type_label, parameter_rows, nested_params, markdown_table

const _GLOBALS_SCHEMA = Dict{String,Any}("type" => "object")

"""
    type_schema(T)

JSON Schema of a declared field type. `Dependent` is a number or a data file
path; enums are strings (their spellings match loosely, so no JSON enum; the
spellings are listed as examples).
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
        spellings = vcat(collect(string.(instances(T))), collect(keys(enum_aliases(T))))
        return Dict{String,Any}("type" => "string", "examples" => spellings,
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

# allowed values of a schema (on the elements of an array, as check_constraints!)
function _allowed!(s::Dict{String,Any}, values::Vector{Any})
    if get(s, "type", "") == "array"
        _allowed!(s["items"], values)
    else
        s["enum"] = values
    end
    return s
end

_json_value(x::Union{Real,String}) = x
_json_value(x::AbstractFloat) = isfinite(x) ? x : nothing   # JSON has no Inf / NaN
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
    fs.allowed === nothing || _allowed!(s, Any[_json_value(v) for v in fs.allowed])
    if !fs.required
        v = _json_value(fs.default)
        v === nothing || (s["default"] = v)
    end
    # a nested section brings its own description; the field's replaces it
    generated = _params_type(fs.type) === nothing ? get(s, "description", "") : ""
    text = filter(!isempty, [fs.description, generated, _quantity_text(fs.quantity)])
    isempty(text) || (s["description"] = join(text, "; "))
    fs.evaluated && (s["x-perilab-eval"] = true)
    fs.when === nothing || (s["x-perilab-when"] = fs.when)
    return s
end

"""
    json_schema(T)

JSON Schema object of the `@params` struct `T`: its description
(`params_description`), its YAML keys with types, defaults, bounds, allowed
values and descriptions, in field order; `required` lists the
required keys; indexed keys (`key_patterns`) become `patternProperties`;
`conditional_required` and `required_groups` become `allOf` entries, the groups
also `x-perilab-required-groups`.
"""
function json_schema(T)
    props = OrderedDict{String,Any}()
    required = String[]
    choices = Any[]
    for fs in parameter_spec(T)
        if fs.category !== nothing
            for (name, S) in registered_models(fs.category)
                props[name] = json_schema(S)
                props[name]["type"] = ["object", "null"]     # `Verlet:` without keys: defaults
                push!(choices, Dict{String,Any}("required" => [name]))
            end
            continue
        end
        props[fs.alias] = field_schema(fs)
        fs.required && push!(required, fs.alias)
    end
    for (alias, values) in schema_choices(T)
        props[alias]["enum"] = values
    end
    props["Globals"] = copy(_GLOBALS_SCHEMA)
    schema = Dict{String,Any}("type" => "object", "properties" => props,
                              "additionalProperties" => false)
    description = params_description(T)
    isempty(description) || (schema["description"] = description)
    isempty(required) || (schema["required"] = sort!(required))
    isempty(choices) || (schema["oneOf"] = choices)
    groups = required_groups(T)
    cases = Any[[_required_case(condition, needed)
                 for (condition, needed) in conditional_required(T)];
                [_group_case(g) for g in groups]]
    isempty(cases) || (schema["allOf"] = cases)
    isempty(groups) ||
        (schema["x-perilab-required-groups"] = Any[_group_json(g) for g in groups])
    patterns = key_patterns(T)
    if !isempty(patterns)
        schema["patternProperties"] = Dict{String,Any}(first(p).pattern => type_schema(last(p))
                                                       for p in patterns)
    end
    return schema
end

# schema of a `conditional_required` / `required_groups` condition: every key has
# one of its values (`nothing`: or is not given; no values: is given with any value)
function _condition_schema(condition::AbstractDict)
    props = OrderedDict{String,Any}()
    required = String[]
    for (k, vs) in condition
        given = filter(!isnothing, vs)
        isempty(given) ||
            (props[k] = length(given) == 1 ? Dict{String,Any}("const" => only(given)) :
                        Dict{String,Any}("enum" => collect(given)))
        nothing in vs || push!(required, k)
    end
    s = Dict{String,Any}()
    isempty(props) || (s["properties"] = props)
    isempty(required) || (s["required"] = required)
    return s
end

# `if` `condition` holds, `then` `required` are required
function _required_case(condition::AbstractDict, required::Vector{String})
    return Dict{String,Any}("if" => _condition_schema(condition),
                            "then" => Dict{String,Any}("required" => required))
end

# the `k`-element subsets of `xs`, in order
function _subsets(xs::Vector{String}, k::Int)
    k == 0 && return [String[]]
    length(xs) < k && return Vector{String}[]
    return [[[xs[1]; s] for s in _subsets(xs[2:end], k - 1)]; _subsets(xs[2:end], k)]
end

_any_of(keys, k) = Any[Dict{String,Any}("required" => s) for s in _subsets(keys, k)]

# a required group: exactly one -> oneOf; else anyOf of the `min`-subsets, and
# `not` anyOf of the `max + 1`-subsets; inside `if`/`then` for a `when`
function _group_case(g)
    s = if g.min == 1 && g.max == 1
        Dict{String,Any}("oneOf" => _any_of(g.keys, 1))
    else
        s = Dict{String,Any}("anyOf" => _any_of(g.keys, g.min))
        g.max === nothing || g.max >= length(g.keys) ||
            (s["not"] = Dict{String,Any}("anyOf" => _any_of(g.keys, g.max + 1)))
        s
    end
    g.when === nothing && return s
    return Dict{String,Any}("if" => _condition_schema(g.when), "then" => s)
end

# "exactly one of A, B", "two of A, B, C" (+ the note)
function _group_text(g)
    text = (g.min == g.max ? "exactly " : "") * _count_word(g.min) * " of " *
           join(g.keys, ", ")
    return g.note === nothing ? text : text * " ($(g.note))"
end

# the `x-perilab-required-groups` entry: the group as data, for consumers that
# do not evaluate the anyOf branches; `when` as declared (null: key not given)
function _group_json(g)
    entry = Dict{String,Any}("keys" => g.keys, "min" => g.min, "max" => g.max,
                             "description" => "give " * _group_text(g))
    g.when === nothing || (entry["when"] = g.when)
    return entry
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
    choices = schema_choices(T)
    groups = Dict{String,Vector{String}}()
    for g in required_groups(T), key in g.keys
        push!(get!(groups, key, String[]), "required: " * _group_text(g))
    end
    rows = [_row(fs) for fs in parameter_spec(T)]
    rows = [haskey(groups, r.key) && r.required == "optional" ?
            merge(r, (required = join(groups[r.key], "; "),)) : r for r in rows]
    return [haskey(choices, r.key) ?
            merge(r, (range = "one of: " * join(choices[r.key], ", "),)) : r
            for r in rows]
end

function _row(fs::FieldSpec)
    if fs.category !== nothing
        return (key = "<" * _category_noun(fs.category) * ">",
                type = "one of: " * join(registered_names(fs.category), ", "),
                required = "required", default = "", range = "", quantity = "",
                description = fs.description)
    end
    return (key = fs.alias, type = type_label(fs.type),
            required = fs.required ? "required" : "optional",
            default = _default_text(fs), range = _range_text(fs),
            quantity = fs.quantity === nothing ? "" : string(fs.quantity),
            description = join(filter(!isempty, [fs.description, _when_text(fs)]), "; "))
end

# "only with Type: Disk" for a field with `when`
function _when_text(fs::FieldSpec)
    fs.when === nothing && return ""
    conditions = [(k == "\$solver" ? "solver" : k) * ": " * join(vs, ", ")
                  for (k, vs) in sort!(collect(fs.when); by = first)]
    return "only with " * join(conditions, " and ")
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
        if fs.category !== nothing
            append!(nested, [name => S for (name, S) in registered_models(fs.category)])
            continue
        end
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
