# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

export to_json_schema

# model sections for the generators: those under `Models` and the contact models
# (named entries of the top-level `Contact` section)
const _GENERATED_MODEL_SECTIONS = (MODEL_SECTIONS..., ("Contact", :contact, "Type"))

# a named model entry: the category's base keys, the key naming the model, and
# per registered model (`if` the name is exactly that model) its own keys, with
# the key set closed for that model (composite names match no `if`: shared keys only)
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
        then = Dict{String,Any}("properties" => merge(entry["properties"], model["properties"]),
                                "additionalProperties" => false)
        haskey(model, "required") && (then["required"] = model["required"])
        patterns = merge(get(entry, "patternProperties", Dict{String,Any}()),
                         get(model, "patternProperties", Dict{String,Any}()))
        isempty(patterns) || (then["patternProperties"] = patterns)
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
    props["Contact"] = Dict{String,Any}("type" => "object",
                                        "properties" => Dict{String,Any}("Globals" => ParameterSpec.json_schema(ContactGlobalsParams)),
                                        "additionalProperties" => _model_entry_schema(:contact,
                                                                                         "Type"))
    props["Globals"] = Dict{String,Any}("type" => "object")
    deck["required"] = sort!(unique([get(deck, "required", String[]); "Models"]))
    return Dict{String,Any}("\$schema" => "https://json-schema.org/draft/2020-12/schema",
                            "title" => "PeriLab input deck", "type" => "object",
                            "properties" => Dict{String,Any}("PeriLab" => deck),
                            "required" => ["PeriLab"])
end

export describe

# the declared type without `Nothing`
_declared(T) = T isa Union && Nothing <: T ? ParameterSpec._nonnothing(T) : T

# `P` if `T` holds named entries of the `@params` struct `P` (`Dict{String,P}`), else nothing
function _named_entries(T)
    D = _declared(T)
    return D isa DataType && D <: Dict && ParameterSpec.is_params(valtype(D)) ? valtype(D) :
           nothing
end

# YAML spelling of a default value
_yaml_value(x::AbstractFloat) = isinf(x) ? (x > 0 ? ".inf" : "-.inf") : ParameterSpec._fmt_value(x)
_yaml_value(x) = ParameterSpec._fmt_value(x)

# placeholder name of one named entry, e.g. "blocks_1"
_entry_name(alias::AbstractString) = replace(lowercase(alias), r"[^a-z0-9]+" => "_") * "_1"

# (title, struct, base struct or nothing, (name key, base title) or nothing, named
# entries?) of a section or registered model name; aborts with a suggestion if unknown
function _described(name::AbstractString)
    for fs in ParameterSpec.parameter_spec(PeriLabSections)
        fs.alias == name || continue
        T = ParameterSpec._params_type(fs.type)
        T === nothing && continue
        named = _named_entries(fs.type) !== nothing
        title = named ? "$name (section, one entry per name)" : "$name (section)"
        return (title, T, nothing, nothing, named)
    end
    name == "Contact" && return ("Contact (section, one entry per contact model; shared keys)",
                                 ContactBaseParams, nothing, nothing, true)
    for (section, category, name_key) in _GENERATED_MODEL_SECTIONS
        for (model, T) in ParameterSpec.registered_models(category)
            model == name || continue
            kind = lowercase(replace(section, " Models" => ""))
            return ("$name ($kind model)", T, ParameterSpec.base_model(category),
                    (name_key, "Shared $kind keys"), false)
        end
    end
    candidates = [[fs.alias for fs in ParameterSpec.parameter_spec(PeriLabSections)
                   if ParameterSpec._params_type(fs.type) !== nothing];
                  "Contact";
                  [first(m) for (_, c, _) in _GENERATED_MODEL_SECTIONS
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

# YAML lines of `T`: required keys with a placeholder, optional keys commented out;
# named entries get one placeholder entry
function _template(io::IO, T, indent::String)
    for fs in ParameterSpec.parameter_spec(T)
        prefix = fs.required ? "" : "# "
        entry = _named_entries(fs.type)
        if entry !== nothing
            println(io, indent, prefix, fs.alias, ":")
            println(io, indent, prefix, "  ", _entry_name(fs.alias), ":")
            _template(io, entry, indent * prefix * "    ")
            continue
        end
        P = _declared(fs.type)
        if ParameterSpec.is_params(P)
            println(io, indent, prefix, fs.alias, ":")
            _template(io, P, indent * prefix * "  ")
            continue
        end
        placeholder = "<" * ParameterSpec.type_label(fs.type) * ">"
        default = fs.required || fs.default === nothing ||
                  (fs.default isa AbstractDict && isempty(fs.default)) ? placeholder :
                  _yaml_value(fs.default)
        note = filter(!isempty, [fs.required ? "required" : "optional",
                                 ParameterSpec._range_text(fs),
                                 fs.quantity === nothing ? "" : string(fs.quantity),
                                 fs.description])
        println(io, indent, prefix, fs.alias, ": ", default, "  # ", join(note, "; "))
    end
end

"""
    describe(name; template = false)
    describe(io, name; template = false)

Prints the parameters of a section (e.g. "Solver") or a registered model (e.g.
"Correspondence Plastic"). With `template = true`, prints a YAML block to copy
into an input deck: required keys with a placeholder, optional keys commented
out with their default; sections of named entries (e.g. "Blocks") get one
placeholder entry.
"""
function describe(io::IO, name::AbstractString; template::Bool = false)
    title, T, base, model_key, named = _described(name)
    if template
        if model_key !== nothing
            println(io, "My ", lowercase(first(model_key)), ":")
            println(io, "  ", first(model_key), ": \"", name, "\"")
            _template(io, T, "  ")
            base === nothing || _template(io, base, "  ")
        elseif name == "Contact"
            # one placeholder entry of the first contact model
            model, M = first(ParameterSpec.registered_models(:contact))
            println(io, name, ":")
            println(io, "  ", _entry_name(name), ":")
            println(io, "    Type: \"", model, "\"")
            _template(io, M, "    ")
            _template(io, T, "    ")
        elseif named
            println(io, name, ":")
            println(io, "  ", _entry_name(name), ":")
            _template(io, T, "    ")
        else
            println(io, name, ":")
            _template(io, T, "  ")
        end
        return nothing
    end
    println(io, title)
    name == "Contact" &&
        println(io, "  Type: one of ", join(ParameterSpec.registered_names(:contact), ", "))
    _describe_table(io, T, "  ")
    if base !== nothing
        println(io, last(model_key), ":")
        _describe_table(io, base, "  ")
    end
    return nothing
end

describe(name::AbstractString; template::Bool = false) = describe(stdout, name;
                                                                  template = template)

export generate_parameter_docs

const _DOCS_HEADER = "<!-- generated by PeriLab.generate_parameter_docs; do not edit -->\n\n"
const _QUANTITY_NOTE = "Quantities name what a value measures; PeriLab has no fixed units, use them in your consistent unit system.\n\n"

function _docs_section(io::IO, title::String, T, level::String)
    println(io, level, " ", title, "\n")
    print(io, ParameterSpec.markdown_table(ParameterSpec.parameter_rows(T)), "\n")
    for (alias, P) in ParameterSpec.nested_params(T)
        _docs_section(io, title * " → " * alias, P, level * "#")
    end
end

"""
    generate_parameter_docs(dir)

Writes the input reference pages `input_sections.md` (the deck's sections) and
`input_models.md` (shared keys and every registered model per category, plus the
pre-calculation switches) into `dir` and returns their paths.
"""
function generate_parameter_docs(dir::AbstractString)
    mkpath(dir)
    sections_file = joinpath(dir, "input_sections.md")
    open(sections_file, "w") do io
        print(io, _DOCS_HEADER, "# Input Sections\n\n", _QUANTITY_NOTE)
        for fs in ParameterSpec.parameter_spec(PeriLabSections)
            T = ParameterSpec._params_type(fs.type)
            T === nothing || _docs_section(io, fs.alias, T, "##")
        end
        _docs_section(io, "Contact → Globals", ContactGlobalsParams, "##")
        println(io, "Every other entry of `Contact` is a contact model, see the input models.\n")
    end
    models_file = joinpath(dir, "input_models.md")
    open(models_file, "w") do io
        print(io, _DOCS_HEADER, "# Input Models\n\n", _QUANTITY_NOTE)
        for (section, category, name_key) in _GENERATED_MODEL_SECTIONS
            println(io, "## ", section, "\n")
            println(io, "Each entry names its model in `", name_key, "`.\n")
            base = ParameterSpec.base_model(category)
            base === nothing || _docs_section(io, "Shared keys", base, "###")
            for (name, T) in ParameterSpec.registered_models(category)
                _docs_section(io, name, T, "###")
            end
        end
        println(io, "## Pre Calculation\n")
        println(io, "`Pre Calculation Global` and each entry of `Pre Calculation Models` switch these on (`true`) or off (`false`):\n")
        for name in ParameterSpec.registered_names(:pre_calculation)
            println(io, "- ", name)
        end
    end
    return [sections_file, models_file]
end
