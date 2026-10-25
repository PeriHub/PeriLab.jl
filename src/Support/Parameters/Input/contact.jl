# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# In the Contact section, "Globals" is real configuration (not the generic
# escape hatch); every other entry is a contact model: the shared keys
# (`ContactBaseParams`) and the keys of the model named in `Type`, which the
# contact model module registers with `register_contact`.

@params "Contact search settings shared by all contact models" struct ContactGlobalsParams
    "Steps between two global contact searches; in between only the pairs found last are updated"
    global_search_frequency::Int64 = opt("Global Search Frequency"; default = 1, min = 1)
    "Uses only free surface nodes of the contact blocks as contact nodes, instead of all nodes"
    only_surface_contact_nodes::Bool = opt("Only Surface Contact Nodes"; default = true)
end

@params "A pair of blocks that can come into contact" struct ContactGroupParams
    "Block id of the master side"
    master_block_id::Int64 = req("Master Block ID"; min = 1)
    "Block id of the slave side"
    slave_block_id::Int64 = req("Slave Block ID"; min = 1)
    "Distance within which the global search collects candidate contact pairs; larger than the contact radius"
    search_radius::Float64 = req("Search Radius"; min = 0, quantity = :length)
    "Steps between two global contact searches of this group; overrides the value in Contact Globals"
    global_search_frequency::Union{Nothing,Int64} = opt("Global Search Frequency";
                                                        default = nothing, min = 1)
end

@params "Contact model between pairs of blocks" struct ContactBaseParams
    "Distance below which two nodes are in contact and repel each other"
    contact_radius::Float64 = req("Contact Radius"; min = 0, quantity = :length)
    "Selects the contact stiffness for \"plane stress\", \"plane strain\" or 3D (any other value)"
    symmetry::String = opt("Symmetry"; default = "3D")
    "Named pairs of blocks this contact model acts on"
    contact_groups::Dict{String,ContactGroupParams} = req("Contact Groups")
end

"""
    ContactInput

The typed `Contact` section: `globals` and the contact models by name, each a
`ParameterSpec.WithBase` (shared keys in `base`, the model's own keys in `model`,
its `Type` in `name`).
"""
struct ContactInput
    globals::ContactGlobalsParams
    models::Dict{String,Any}
end

function _contact_section(raw, path::String, ctx::ParseContext)
    raw isa AbstractDict && return true
    add_error!(ctx, path,
               "expected a section of `key: value` entries, got $(ParameterSpec._describe(raw))")
    return false
end

"""
    parse_contact(raw, path, ctx)

Reads the `Contact` section: `Globals` (search settings for all contact
models, defaults if absent) and the named contact models; contact models cannot
be combined with `+`.
"""
function parse_contact(raw, path::String, ctx::ParseContext)
    _contact_section(raw, path, ctx) || return nothing
    globals_path = join_path(path, "Globals")
    raw_globals = get(raw, "Globals", Dict{String,Any}())
    globals = _contact_section(raw_globals, globals_path, ctx) ?
              parse_section(ContactGlobalsParams, raw_globals, globals_path, ctx) : nothing
    models = Dict{String,Any}()
    ok = globals !== nothing
    for (name, entry) in raw
        key = string(name)
        key == "Globals" && continue
        model_path = join_path(path, key)
        if !_contact_section(entry, model_path, ctx)
            ok = false
            continue
        end
        model = ParameterSpec.parse_model(:contact,
                                          Dict{String,Any}(string(k) => v for (k, v) in entry),
                                          model_path, ctx; name_key = "Type")
        if model === nothing
            ok = false
        elseif model.model isa ParameterSpec.Composite
            add_error!(ctx, join_path(model_path, "Type"),
                       "contact models cannot be combined with +")
            ok = false
        else
            models[key] = model
        end
    end
    return ok ? ContactInput(globals, models) : nothing
end

function check!(g::ContactGroupParams, path::String, ctx::ParseContext)
    if g.master_block_id == g.slave_block_id
        add_error!(ctx, path,
                   "Master Block ID and Slave Block ID are equal; self contact is not implemented")
    end
    g.search_radius > 0 ||
        add_error!(ctx, join_path(path, "Search Radius"), "must be greater than zero")
    return nothing
end

"Sorted ids of all blocks that take part in a contact group."
function contact_blocks(contact::ContactInput)
    ids = Int64[]
    for model in values(contact.models), group in values(model.base.contact_groups)
        push!(ids, group.master_block_id, group.slave_block_id)
    end
    return sort!(unique!(ids))
end

"Global search frequency of a contact group; falls back to `Globals`."
contact_search_frequency(group::ContactGroupParams, globals::ContactGlobalsParams) = something(group.global_search_frequency,
                                                                                               globals.global_search_frequency)
