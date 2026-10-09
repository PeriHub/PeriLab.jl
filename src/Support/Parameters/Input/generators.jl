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
