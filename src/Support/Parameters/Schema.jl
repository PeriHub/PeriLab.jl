# SPDX-License-Identifier: BSD-3-Clause
#
# Core DSL for defining PeriLab's input schema once, and deriving from it:
#   1. runtime validation of loaded YAML/JSON params (validate_params)
#   2. a JSON Schema document for PeriHub / frontend / API (to_json_schema)
#
# This replaces the hand-rolled `expected_structure::Dict{Any,Any}` +
# `validate_structure_recursive` in parameter_handling.jl.

module Schema

using DataStructures: OrderedDict

export SField, SObject, SAny, SOneOf, SArray, SWith,
       requires_together, required_if,
       ValidationError, validate_params, to_json_schema

# ---------------------------------------------------------------------------
# Node types
# ---------------------------------------------------------------------------

abstract type SchemaNode end

"""
    SField(type; required=false, default=nothing, enum=nothing,
           min=nothing, max=nothing, description="")

A leaf value. `type` may be a single Julia type or a `Vector` of types
(mirrors the old `Union{Float64,Int64}` pattern, but stays introspectable
for JSON Schema export, where Julia `Union`s don't translate directly).
"""
struct SField <: SchemaNode
    type::Vector{DataType}
    required::Bool
    default::Any
    enum::Union{Nothing,Vector}
    min::Union{Nothing,Real}
    max::Union{Nothing,Real}
    description::String
end
function SField(type; required::Bool = false, default = nothing,
                enum = nothing, min = nothing, max = nothing, description = "")
    types = type isa AbstractVector ? Vector{DataType}(type) : DataType[type]
    SField(types, required, default, enum, min, max, description)
end

"""
    SObject(fields::AbstractDict; required=false, description="")

A dict with a fixed, known set of keys (e.g. "Discretization").
"""
struct SObject <: SchemaNode
    fields::OrderedDict{String,SchemaNode}
    required::Bool
    description::String
end
function SObject(fields::AbstractDict; required::Bool = false, description = "")
    SObject(OrderedDict{String,SchemaNode}(fields), required, description)
end

"""
    SAny(value::SchemaNode; required=false, min_entries=nothing, description="")

A dict whose *keys are arbitrary, user-chosen names* (block names, output
block names, step names, ...) but whose *values* all share one schema.
Replaces the old `"Any" => [schema, required]` convention.
"""
struct SAny <: SchemaNode
    value::SchemaNode
    required::Bool
    min_entries::Union{Nothing,Int}
    description::String
end
function SAny(value::SchemaNode; required::Bool = false, min_entries = nothing,
             description = "")
    SAny(value, required, min_entries, description)
end

"""
    SOneOf(options::AbstractDict; required=false, description="")

Exactly one of several named alternatives must be present as a key
(e.g. Solver must contain exactly one of Verlet/Static/Newmark/...).
"""
struct SOneOf <: SchemaNode
    options::OrderedDict{String,SchemaNode}
    required::Bool
    description::String
end
function SOneOf(options::AbstractDict; required::Bool = false, description = "")
    SOneOf(OrderedDict{String,SchemaNode}(options), required, description)
end

"""
    SArray(item::SchemaNode; required=false, min_items=nothing, description="")

A list of items sharing one schema (not heavily used yet in PeriLab's
input, but needed for e.g. future list-valued fields).
"""
struct SArray <: SchemaNode
    item::SchemaNode
    required::Bool
    min_items::Union{Nothing,Int}
    description::String
end
function SArray(item::SchemaNode; required::Bool = false, min_items = nothing,
                description = "")
    SArray(item, required, min_items, description)
end

"""
    SWith(base::SObject, constraints::Vector; required=false, description="")

Wraps an `SObject` with cross-field rules that plain key/type checks can't
express — "if C11 is set, C12 must be too", "if Symmetry is anisotropic,
these 21 keys are required together", etc.

Each constraint is a NamedTuple with:
  - `check::Function`     -- (dict, path) -> Vector{ValidationError}
  - `dependent_required::Dict{String,Vector{String}}` (optional, for JSON
    Schema's native `dependentRequired`) OR
  - `if_then::Dict`        (optional, for JSON Schema's native `if`/`then`)

Use `requires_together(...)` / `required_if(...)` below to build these
instead of writing the NamedTuple by hand.
"""
struct SWith <: SchemaNode
    base::SObject
    constraints::Vector{NamedTuple}
    required::Bool
    description::String
end
function SWith(base::SObject, constraints::Vector; required::Bool = base.required,
              description = base.description)
    SWith(base, Vector{NamedTuple}(constraints), required, description)
end

"""
    requires_together(keys::Vector{String}; description="")

If ANY of `keys` is present in the dict, ALL of them must be present.
Typical use: a physical quantity (e.g. a stiffness matrix) that's only
meaningful when fully specified.
"""
function requires_together(keys::Vector{String}; description::String = "")
    check = (dict, path) -> begin
        errs = ValidationError[]
        present = [k for k in keys if haskey(dict, k)]
        if !isempty(present) && length(present) < length(keys)
            missing = setdiff(keys, present)
            push!(errs,
                 ValidationError(path,
                                 "$(join(present, ", ")) given, also requires: $(join(missing, ", "))" *
                                 (isempty(description) ? "" : " ($description)")))
        end
        return errs
    end
    dependent_required = Dict(k => setdiff(keys, [k]) for k in keys)
    return (; check, dependent_required)
end

"""
    required_if(discriminant_key::String, discriminant_value, required_keys::Vector{String})

If `dict[discriminant_key] == discriminant_value`, then all of
`required_keys` must be present. Typical use: "Symmetry" == "anisotropic"
implies a specific subset of Cij entries is required.
"""
function required_if(discriminant_key::String, discriminant_value, required_keys::Vector{String})
    check = (dict, path) -> begin
        errs = ValidationError[]
        if get(dict, discriminant_key, nothing) == discriminant_value
            missing = [k for k in required_keys if !haskey(dict, k)]
            if !isempty(missing)
                push!(errs,
                     ValidationError(path,
                                     "$discriminant_key = $discriminant_value requires: $(join(missing, ", "))"))
            end
        end
        return errs
    end
    if_then = Dict(
        "if" => Dict("properties" => Dict(discriminant_key => Dict("const" => discriminant_value))),
        "then" => Dict("required" => required_keys),
    )
    return (; check, if_then)
end

# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

struct ValidationError
    path::String
    message::String
end
Base.show(io::IO, e::ValidationError) = print(io, "$(e.path): $(e.message)")

_join(path::String, key) = isempty(path) ? string(key) : "$path.$key"

function _type_ok(value, types::Vector{DataType})
    return any(t -> value isa t, types)
end

function _validate!(errors::Vector{ValidationError}, node::SField, value, path::String)
    if !_type_ok(value, node.type)
        push!(errors,
             ValidationError(path,
                             "expected type $(join(node.type, " or ")), got $(typeof(value))"))
        return
    end
    if node.enum !== nothing && !(value in node.enum)
        push!(errors, ValidationError(path, "value $value not in allowed set $(node.enum)"))
    end
    if node.min !== nothing && value isa Real && value < node.min
        push!(errors, ValidationError(path, "value $value below minimum $(node.min)"))
    end
    if node.max !== nothing && value isa Real && value > node.max
        push!(errors, ValidationError(path, "value $value above maximum $(node.max)"))
    end
end

function _validate!(errors::Vector{ValidationError}, node::SObject, value, path::String)
    if !(value isa AbstractDict)
        push!(errors, ValidationError(path, "expected an object/dict, got $(typeof(value))"))
        return
    end
    for (key, child) in node.fields
        child_path = _join(path, key)
        if !haskey(value, key)
            if child_required(child)
                push!(errors, ValidationError(child_path, "missing required key"))
            end
            continue
        end
        _validate!(errors, child, value[key], child_path)
    end
    # Warn (not fail) on truly unknown keys, mirroring old "Key not known" @warn.
    known = Set(keys(node.fields))
    for key in keys(value)
        key == "Globals" && continue   # explicit escape hatch, same as before
        if !(key in known)
            push!(errors, ValidationError(_join(path, key), "unknown key (ignored)"))
        end
    end
end

function _validate!(errors::Vector{ValidationError}, node::SAny, value, path::String)
    if !(value isa AbstractDict)
        push!(errors, ValidationError(path, "expected an object/dict, got $(typeof(value))"))
        return
    end
    if node.min_entries !== nothing && length(value) < node.min_entries
        push!(errors, ValidationError(path, "expected at least $(node.min_entries) entries"))
    end
    for (key, entry) in value
        key == "Globals" && continue
        _validate!(errors, node.value, entry, _join(path, key))
    end
end

function _validate!(errors::Vector{ValidationError}, node::SOneOf, value, path::String)
    if !(value isa AbstractDict)
        push!(errors, ValidationError(path, "expected an object/dict, got $(typeof(value))"))
        return
    end
    present = [k for k in keys(node.options) if haskey(value, k)]
    if length(present) == 0
        push!(errors,
             ValidationError(path,
                             "must contain exactly one of: $(join(keys(node.options), ", "))"))
        return
    end
    if length(present) > 1
        push!(errors, ValidationError(path, "must contain only one of: $(join(present, ", "))"))
        return
    end
    chosen = present[1]
    _validate!(errors, node.options[chosen], value[chosen], _join(path, chosen))
end

function _validate!(errors::Vector{ValidationError}, node::SWith, value, path::String)
    _validate!(errors, node.base, value, path)
    value isa AbstractDict || return  # base check already flagged this
    for c in node.constraints
        append!(errors, c.check(value, path))
    end
end

function _validate!(errors::Vector{ValidationError}, node::SArray, value, path::String)
    if !(value isa AbstractVector)
        push!(errors, ValidationError(path, "expected an array, got $(typeof(value))"))
        return
    end
    if node.min_items !== nothing && length(value) < node.min_items
        push!(errors, ValidationError(path, "expected at least $(node.min_items) items"))
    end
    for (i, item) in enumerate(value)
        _validate!(errors, node.item, item, _join(path, i))
    end
end

child_required(node::SchemaNode) = node.required

"""
    validate_params(schema::SObject, params::AbstractDict) -> Vector{ValidationError}

Validate a fully-loaded params dict (e.g. from YAML.load_file) against
`schema`. Returns an empty vector if valid. Unlike the old validator this
collects *all* errors instead of stopping/warning on the first mismatch.
"""
function validate_params(schema::SObject, params::AbstractDict)
    errors = ValidationError[]
    _validate!(errors, schema, params, "")
    return errors
end

# ---------------------------------------------------------------------------
# JSON Schema export
# ---------------------------------------------------------------------------

_JULIA_TO_JSON_TYPE = Dict(
    Int64 => "integer", Int32 => "integer",
    Float64 => "number", Float32 => "number",
    String => "string",
    Bool => "boolean",
)

function _json_types(types::Vector{DataType})
    jtypes = unique(get(_JULIA_TO_JSON_TYPE, t, "string") for t in types)
    return length(jtypes) == 1 ? jtypes[1] : jtypes
end

function to_json_schema(node::SField)
    d = Dict{String,Any}("type" => _json_types(node.type))
    node.enum !== nothing && (d["enum"] = node.enum)
    node.min !== nothing && (d["minimum"] = node.min)
    node.max !== nothing && (d["maximum"] = node.max)
    node.default !== nothing && (d["default"] = node.default)
    !isempty(node.description) && (d["description"] = node.description)
    return d
end

function to_json_schema(node::SObject)
    props = OrderedDict{String,Any}(k => to_json_schema(v) for (k, v) in node.fields)
    required = [k for (k, v) in node.fields if child_required(v)]
    d = Dict{String,Any}("type" => "object", "properties" => props)
    !isempty(required) && (d["required"] = required)
    !isempty(node.description) && (d["description"] = node.description)
    return d
end

function to_json_schema(node::SAny)
    d = Dict{String,Any}("type" => "object",
                         "patternProperties" => Dict(".+" => to_json_schema(node.value)),
                         "additionalProperties" => false)
    node.min_entries !== nothing && (d["minProperties"] = node.min_entries)
    !isempty(node.description) && (d["description"] = node.description)
    return d
end

function to_json_schema(node::SOneOf)
    alts = [Dict{String,Any}("type" => "object",
                             "properties" => Dict(k => to_json_schema(v)),
                             "required" => [k],
                             "additionalProperties" => true)
            for (k, v) in node.options]
    d = Dict{String,Any}("oneOf" => alts)
    !isempty(node.description) && (d["description"] = node.description)
    return d
end

function to_json_schema(node::SWith)
    d = to_json_schema(node.base)
    dep_req = Dict{String,Any}()
    if_thens = Any[]
    for c in node.constraints
        haskey(c, :dependent_required) && merge!(dep_req, c.dependent_required)
        haskey(c, :if_then) && push!(if_thens, c.if_then)
    end
    !isempty(dep_req) && (d["dependentRequired"] = dep_req)
    !isempty(if_thens) && (d["allOf"] = if_thens)
    return d
end

function to_json_schema(node::SArray)
    d = Dict{String,Any}("type" => "array", "items" => to_json_schema(node.item))
    node.min_items !== nothing && (d["minItems"] = node.min_items)
    !isempty(node.description) && (d["description"] = node.description)
    return d
end

"""
    to_json_schema(root::SObject; title="PeriLab Input", id=nothing) -> Dict

Wraps the object schema with the `schema`/`id`/`title` envelope PeriHub
expects, ready for `JSON3.write(...)` / `JSON.print(...)` to a `.json` file.
"""
function to_json_schema(root::SObject, ::Val{:document}; title = "PeriLab Input",
                        id = nothing)
    d = to_json_schema(root)
    d["\$schema"] = "https://json-schema.org/draft/2020-12/schema"
    id !== nothing && (d["\$id"] = id)
    d["title"] = title
    return d
end

end # module
