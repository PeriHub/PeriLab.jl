# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

export to_json_schema

# a named model entry: the category's base keys, the key naming the model, and
# per registered model (`if` the name is exactly that model) its own keys
function _model_entry_schema(category::Symbol, name_key::String)
    base = ParameterSpec.base_model(category)
    entry = base === nothing ?
            Dict{String,Any}("type" => "object",
                             "properties" => Dict{String,Any}("Globals" => Dict{String,Any}("type" => "object"))) :
            ParameterSpec.json_schema(base)
    delete!(entry, "additionalProperties")      # model keys come from the matching model
    names = ParameterSpec.registered_names(category)
    entry["properties"][name_key] = Dict{String,Any}("type" => "string",
                                                     "description" => "one of: " *
                                                                      join(names, ", ") *
                                                                      "; models can be combined with +")
    entry["required"] = sort!(unique([get(entry, "required", String[]); name_key]))
    cases = Any[]
    for (name, T) in ParameterSpec.registered_models(category)
        model = ParameterSpec.json_schema(T)
        then = Dict{String,Any}("properties" => model["properties"])
        haskey(model, "required") && (then["required"] = model["required"])
        haskey(model, "patternProperties") &&
            (then["patternProperties"] = model["patternProperties"])
        push!(cases,
              Dict{String,Any}("if" => Dict{String,Any}("properties" => Dict{String,Any}(name_key => Dict{String,Any}("const" => name))),
                               "then" => then))
    end
    isempty(cases) || (entry["allOf"] = cases)
    return entry
end

function _pre_calculation_switches_schema()
    props = Dict{String,Any}(name => Dict{String,Any}("type" => "boolean")
                             for name in ParameterSpec.registered_names(:pre_calculation))
    for (old, new) in DEPRECATED_PRE_CALCULATIONS
        props[old] = Dict{String,Any}("type" => "boolean", "const" => false,
                                      "description" => "deprecated, use \"$new\"")
    end
    return Dict{String,Any}("type" => "object", "properties" => props,
                            "additionalProperties" => false)
end

"""
    to_json_schema(PeriLabInput)

JSON Schema (draft 2020-12) of an input deck, generated from the `@params`
declarations and the registered models (installed and licensed modules). Model
entries list the category's shared keys; the keys of a single named model are
added for exactly that name (models combined with `+` are checked only for the
shared keys).
"""
function to_json_schema(::Type{PeriLabInput})
    deck = ParameterSpec.json_schema(PeriLabSections)
    props = deck["properties"]
    models = Dict{String,Any}()
    for (section, category, name_key) in MODEL_SECTIONS
        models[section] = Dict{String,Any}("type" => "object",
                                           "additionalProperties" => _model_entry_schema(category,
                                                                                         name_key))
    end
    switches = _pre_calculation_switches_schema()
    models["Pre Calculation Global"] = switches
    models["Pre Calculation Models"] = Dict{String,Any}("type" => "object",
                                                        "additionalProperties" => switches)
    props["Models"] = Dict{String,Any}("type" => "object", "properties" => models,
                                       "additionalProperties" => false)
    contact = ParameterSpec.json_schema(ContactModelParams)
    props["Contact"] = Dict{String,Any}("type" => "object",
                                        "properties" => Dict{String,Any}("Globals" => ParameterSpec.json_schema(ContactGlobalsParams)),
                                        "additionalProperties" => contact)
    props["Globals"] = Dict{String,Any}("type" => "object")
    deck["required"] = sort!(unique([get(deck, "required", String[]); "Models"]))
    return Dict{String,Any}("\$schema" => "https://json-schema.org/draft/2020-12/schema",
                            "title" => "PeriLab input deck", "type" => "object",
                            "properties" => Dict{String,Any}("PeriLab" => deck),
                            "required" => ["PeriLab"])
end

export describe

# (title, struct, base struct or nothing, (name key, base title) or nothing) of a
# section or registered model name; aborts with a suggestion if unknown
function _described(name::AbstractString)
    for fs in ParameterSpec.parameter_spec(PeriLabSections)
        fs.alias == name || continue
        T = ParameterSpec._params_type(fs.type)
        T === nothing || return ("$name (section)", T, nothing, nothing)
    end
    name == "Contact" && return ("Contact (section, one entry per contact model)",
                                 ContactModelParams, nothing, nothing)
    for (section, category, name_key) in MODEL_SECTIONS
        for (model, T) in ParameterSpec.registered_models(category)
            model == name || continue
            kind = lowercase(replace(section, " Models" => ""))
            return ("$name ($kind model)", T, ParameterSpec.base_model(category),
                    (name_key, "Shared $kind keys"))
        end
    end
    candidates = [[fs.alias for fs in ParameterSpec.parameter_spec(PeriLabSections)];
                  "Contact";
                  [first(m) for (_, c, _) in MODEL_SECTIONS
                   for m in ParameterSpec.registered_models(c)]]
    suggestion = ParameterSpec.suggest(name, candidates)
    @abort "unknown section or model \"$name\"" *
           (suggestion === nothing ? "" : " — did you mean \"$suggestion\"?")
end

function _describe_table(io::IO, T, indent::String)
    rows = ParameterSpec.parameter_rows(T)
    isempty(rows) && return println(io, indent, "(no keys)")
    width = maximum(length(r.key) for r in rows)
    for r in rows
        info = filter(!isempty,
                      [r.type, r.required == "required" ? "required" : "default " * r.default,
                       r.range, r.quantity, r.description])
        println(io, indent, rpad(r.key, width), "  ", join(info, "; "))
    end
    for (alias, P) in ParameterSpec.nested_params(T)
        println(io, indent, alias, ":")
        _describe_table(io, P, indent * "  ")
    end
end

# YAML lines of `T`: required keys with a placeholder, optional keys commented out
function _template(io::IO, T, indent::String)
    for fs in ParameterSpec.parameter_spec(T)
        prefix = fs.required ? "" : "# "
        P = ParameterSpec._params_type(fs.type)
        if P !== nothing && P === ParameterSpec._nonnothing(fs.type)
            println(io, indent, prefix, fs.alias, ":")
            _template(io, P, indent * prefix * "  ")
            continue
        end
        placeholder = "<" * ParameterSpec.type_label(fs.type) * ">"
        default = ParameterSpec._default_text(fs)
        value = fs.required || default == "—" ? placeholder : default
        note = filter(!isempty, [fs.required ? "required" : "optional",
                                 ParameterSpec._range_text(fs),
                                 fs.quantity === nothing ? "" : string(fs.quantity),
                                 fs.description])
        println(io, indent, prefix, fs.alias, ": ", value, "  # ", join(note, "; "))
    end
end

"""
    describe(name; template = false)
    describe(io, name; template = false)

Prints the parameters of a section (e.g. "Solver") or a registered model (e.g.
"Correspondence Plastic"). With `template = true`, prints a YAML block to copy
into an input deck: required keys with a placeholder, optional keys commented
out with their default.
"""
function describe(io::IO, name::AbstractString; template::Bool = false)
    title, T, base, model_key = _described(name)
    if template
        if model_key === nothing
            println(io, name, ":")
        else
            println(io, "My ", lowercase(first(model_key)), ":")
            println(io, "  ", first(model_key), ": \"", name, "\"")
        end
        _template(io, T, "  ")
        base === nothing || _template(io, base, "  ")
        return nothing
    end
    println(io, title)
    _describe_table(io, T, "  ")
    if base !== nothing
        println(io, last(model_key), ":")
        _describe_table(io, base, "  ")
    end
    return nothing
end

describe(name::AbstractString; template::Bool = false) = describe(stdout, name;
                                                                  template = template)
