# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# In the Contact section, "Globals" is real configuration (not the generic
# escape hatch); every other entry is a contact model.

@params struct ContactGlobalsParams
    global_search_frequency::Int64 = opt("Global Search Frequency"; default = 1, min = 1)
    only_surface_contact_nodes::Bool = opt("Only Surface Contact Nodes"; default = true)
end

@params struct ContactGroupParams
    master_block_id::Int64 = req("Master Block ID"; min = 1)
    slave_block_id::Int64 = req("Slave Block ID"; min = 1)
    search_radius::Float64 = req("Search Radius"; min = 0, quantity = :length)
    global_search_frequency::Union{Nothing,Int64} = opt("Global Search Frequency";
                                                        default = nothing, min = 1)
end

@params struct ContactModelParams
    type::String = req("Type")
    contact_radius::Float64 = req("Contact Radius"; min = 0, quantity = :length)
    contact_stiffness::Float64 = opt("Contact Stiffness"; default = 1e8, min = 0)
    friction_coefficient::Float64 = opt("Friction Coefficient"; default = 0.0, min = 0)
    symmetry::String = opt("Symmetry"; default = "3D",
                           description = "\"plane stress\", \"plane strain\" or 3D (anything else)")
    contact_groups::Dict{String,ContactGroupParams} = req("Contact Groups")
end

struct ContactInput
    globals::ContactGlobalsParams
    models::Dict{String,ContactModelParams}
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
models, defaults if absent) and the named contact models.
"""
function parse_contact(raw, path::String, ctx::ParseContext)
    _contact_section(raw, path, ctx) || return nothing
    globals_path = join_path(path, "Globals")
    raw_globals = get(raw, "Globals", Dict{String,Any}())
    globals = _contact_section(raw_globals, globals_path, ctx) ?
              parse_section(ContactGlobalsParams, raw_globals, globals_path, ctx) : nothing
    models = Dict{String,ContactModelParams}()
    ok = globals !== nothing
    for (name, entry) in raw
        key = string(name)
        key == "Globals" && continue
        model_path = join_path(path, key)
        if !_contact_section(entry, model_path, ctx)
            ok = false
            continue
        end
        model = parse_section(ContactModelParams, entry, model_path, ctx)
        if model === nothing
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
    for model in values(contact.models), group in values(model.contact_groups)
        push!(ids, group.master_block_id, group.slave_block_id)
    end
    return sort!(unique!(ids))
end

"Global search frequency of a contact group; falls back to `Globals`."
contact_search_frequency(group::ContactGroupParams, globals::ContactGlobalsParams) = something(group.global_search_frequency,
                                                                                               globals.global_search_frequency)
