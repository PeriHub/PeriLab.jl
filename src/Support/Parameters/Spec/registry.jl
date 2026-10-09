# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

export register_model!, register_unavailable!, register_material, register_damage,
       register_thermal, register_additive, register_degradation, register_pre_calculation,
       register_contact,
       register_base!, base_model

"A model name that is known (e.g. from a license manifest) but not loaded."
struct UnavailableModel
    reason::String
end

# category => model name => @params struct type or UnavailableModel.
# Only mutated at runtime (module `__init__` functions, runtime-loaded modules),
# never during precompilation.
const REGISTRY = Dict{Symbol,Dict{String,Any}}()

"""
    register_model!(category, name, T)

Makes model `name` with parameter struct `T` available in `category`. Call it
from your module's `__init__()`.
"""
function register_model!(category::Symbol, name::AbstractString, T::Type)
    is_params(T) ||
        throw(ParamsDefinitionError("register_model!: $(_typename(T)) is not an @params struct"))
    models = get!(Dict{String,Any}, REGISTRY, category)
    existing = get(models, name, nothing)
    if existing isa Type && existing !== T
        throw(ParamsDefinitionError("$category model \"$name\" is already registered by $(_typename(existing))"))
    end
    models[String(name)] = T
    return nothing
end

"""
    register_unavailable!(category, name; reason)

Records a model name that exists but cannot be used (e.g. no license), so a
deck referencing it gets a specific message. Never hides a registered model.
"""
function register_unavailable!(category::Symbol, name::AbstractString;
                               reason::AbstractString = "requires a license that is not available")
    models = get!(Dict{String,Any}, REGISTRY, category)
    haskey(models, name) && models[name] isa Type && return nothing
    models[String(name)] = UnavailableModel(String(reason))
    return nothing
end

function lookup_model(category::Symbol, name::AbstractString)
    return get(get(REGISTRY, category, Dict{String,Any}()), name, nothing)
end

function registered_names(category::Symbol)
    models = get(REGISTRY, category, Dict{String,Any}())
    return sort!([name for (name, entry) in models if entry isa Type])
end

"Registered models of `category` as name => parameter type, sorted by name."
function registered_models(category::Symbol)
    models = get(REGISTRY, category, Dict{String,Any}())
    return sort!([name => T for (name, T) in models if T isa Type]; by = first)
end

register_material(name::AbstractString, T::Type) = register_model!(:material, name, T)
register_damage(name::AbstractString, T::Type) = register_model!(:damage, name, T)
register_thermal(name::AbstractString, T::Type) = register_model!(:thermal, name, T)
register_additive(name::AbstractString, T::Type) = register_model!(:additive, name, T)
register_degradation(name::AbstractString, T::Type) = register_model!(:degradation, name, T)
register_contact(name::AbstractString, T::Type) = register_model!(:contact, name, T)
function register_pre_calculation(name::AbstractString, T::Type)
    register_model!(:pre_calculation, name, T)
end

# category => @params struct read by every model of the category
const BASES = Dict{Symbol,Any}()

"""
    register_base!(category, T)

Declares `T` as the base part of `category`: every model of the category reads
`T`'s keys from its YAML block in addition to its own (e.g. the shared material
moduli). Call it from the factory's `__init__()`.
"""
function register_base!(category::Symbol, T::Type)
    is_params(T) ||
        throw(ParamsDefinitionError("register_base!: $(_typename(T)) is not an @params struct"))
    existing = get(BASES, category, nothing)
    if existing !== nothing && existing !== T
        throw(ParamsDefinitionError("$category base parameters are already registered by $(_typename(existing))"))
    end
    BASES[category] = T
    return nothing
end

base_model(category::Symbol) = get(BASES, category, nothing)
