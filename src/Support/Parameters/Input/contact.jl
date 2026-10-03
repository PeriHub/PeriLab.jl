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
    contact_stiffness::Float64 = req("Contact Stiffness"; min = 0)
    friction_coefficient::Union{Nothing,Float64} = opt("Friction Coefficient";
                                                       default = nothing, min = 0)
    symmetry::Union{Nothing,String} = opt("Symmetry"; default = nothing)
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
