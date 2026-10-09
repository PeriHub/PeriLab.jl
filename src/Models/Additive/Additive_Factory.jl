# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Additive

using .....Data_Manager
using .....PeriLabExceptions: @abort
using .....ModuleLoader: licensed_modules, optional_local_modules

local_user_modules = optional_local_modules(@__DIR__, "additive_name")

for mod in local_user_modules
    include(mod["File"])
end

using .....Helpers: find_inverse_bond_id

export compute_model
export init_model
export init_fields
export fields_for_local_synchronization

"""
    init_fields()

Initialize additive model fields
"""
function init_fields()
    if !Data_Manager.has_key("Activation_Time")
        @abort "'Activation_Time' is missing. Please define an 'Activation_Time' for each point in the mesh file."
        return
    end
    # must be specified, because it might be that no temperature model has been defined
    Data_Manager.create_node_scalar_field("Temperature", Float64)
    Data_Manager.create_node_scalar_field("Heat Flow", Float64)

    bond_damageN = Data_Manager.get_bond_damage("N")
    bond_damageNP1 = Data_Manager.get_bond_damage("NP1")
    nnodes = Data_Manager.get_nnodes()
    if !Data_Manager.has_key("Active")
        active = Data_Manager.create_constant_node_scalar_field("Active", Bool;
                                                                default_value = false)
        for iID in 1:nnodes
            bond_damageN[iID] .= 0
            bond_damageNP1[iID] .= 0
        end
    end
    nlist = Data_Manager.get_nlist()
    inverse_nlist = Data_Manager.set_inverse_nlist(find_inverse_bond_id(nlist))
end

global licensed = Any[]

"""
    load_licensed_models()

Includes the licensed additive modules (if a license is configured) so that they
register their parameter structs before the input deck is read. Runs once.
"""
function load_licensed_models()
    isempty(licensed) || return licensed
    global licensed = licensed_modules(@__MODULE__, "Additive")
    return licensed
end

"""
    load_licensed_models(raw_deck)

Loads the licensed additive modules only if the raw input deck (the YAML dict) has
`Models: Additive Models`, so that runs without additive models never contact the
license server.
"""
function load_licensed_models(raw_deck::AbstractDict)
    models = get(get(raw_deck, "PeriLab", Dict()), "Models", nothing)
    (models isa AbstractDict && haskey(models, "Additive Models")) || return Any[]
    return load_licensed_models()
end

model_module(p) = parentmodule(typeof(p))

"""
    compute_model(nodes::AbstractVector{Int64}, p, block::Int64, time::Float64, dt::Float64)

Computes the additive model of a block.

# Arguments
- `nodes::AbstractVector{Int64}`: The nodes
- `p`: The typed additive model of the block
- `block::Int64`: The block
- `time::Float64`: The current time
- `dt::Float64`: The time step
"""
function compute_model(nodes::AbstractVector{Int64}, p, block::Int64, time::Float64,
                       dt::Float64)
    # invokelatest: licensed models are defined at runtime
    return Base.@invokelatest model_module(p).compute_model(nodes, p, block, time, dt)
end

"""
    init_model(nodes::AbstractVector{Int64}, block::Int64)

Initialize the additive model of a block (`Data_Manager.get_block_models(block).additive`).

# Arguments
- `nodes::AbstractVector{Int64}`: Nodes for the additive model.
- `block::Int64`: Block identifier for the additive model.
"""
function init_model(nodes::AbstractVector{Int64}, block::Int64)
    p = Data_Manager.get_block_models(block).additive
    Base.@invokelatest model_module(p).init_model(nodes, p, block)
end

"""
    fields_for_local_synchronization(model, block)

Defines all synchronization fields for local synchronization

# Arguments
- `model::String`: Model class.
- `block::Int64`: block ID
"""
function fields_for_local_synchronization(model, block)
    p = Data_Manager.get_block_models(block).additive
    return Base.@invokelatest model_module(p).fields_for_local_synchronization(model)
end

end
