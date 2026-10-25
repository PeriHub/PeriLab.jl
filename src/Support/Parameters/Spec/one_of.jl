# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

export one_of

"""
    one_of(category; description)

A sub-section whose key is a model name registered in `category` (e.g.
`:solver`); exactly one must be given. The field holds the parsed struct of
that model, so its declared type must be the abstract supertype of the
registered structs (see `register_one_of_type!`). An empty sub-section means
the model's defaults.
"""
function one_of(category::Symbol; description::AbstractString = "")
    return FieldDecl("", true, NO_DEFAULT, nothing, nothing, nothing, nothing,
                     String(description), category, false)
end

"`:ut_drive` → \"ut drive\", used in messages."
_category_noun(category::Symbol) = replace(string(category), "_" => " ")

# all names of `category`, available or not
_one_of_names(category::Symbol) = sort!(collect(keys(get(REGISTRY, category,
                                                         Dict{String,Any}()))))

"YAML keys read by the `one_of` fields of `T`: the names registered in their categories."
function one_of_keys(T)
    found = Set{String}()
    for fs in parameter_spec(T)
        fs.category === nothing || union!(found, _one_of_names(fs.category))
    end
    return found
end

function _select_one_of(fs::FieldSpec, dict::AbstractDict, path::String, ctx::ParseContext)
    noun = _category_noun(fs.category)
    given = sort!([string(k) for k in keys(dict)
                   if lookup_model(fs.category, string(k)) !== nothing])
    if isempty(given)
        return _fail(ctx, path,
                     "one $noun is required: $(join(registered_names(fs.category), ", "))")
    elseif length(given) > 1
        return _fail(ctx, path, "only one $noun may be given, found: $(join(given, ", "))")
    end
    name = only(given)
    key_path = join_path(path, name)
    entry = lookup_model(fs.category, name)
    entry isa UnavailableModel && return _fail(ctx, key_path, "$noun \"$name\" $(entry.reason)")
    raw = first(v for (k, v) in dict if string(k) == name)
    # `Verlet:` without keys is an empty YAML value: the model's defaults
    raw === nothing && (raw = Dict{String,Any}())
    return convert_value(entry, raw, key_path, ctx)
end
