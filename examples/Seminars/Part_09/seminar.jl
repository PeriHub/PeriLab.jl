# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Damage_template

using ......Data_Manager

export init_model
export compute_model
export fields_for_local_synchronization

using ......ParameterSpec: @params, register_damage

"""
    DamageTemplateParams

Declare the YAML keys your damage model needs beyond the shared damage keys
(Critical Value, Interblock Damage, Anisotropic Damage, Local Damping are in
`damage.base`). Register it under your model name by uncommenting `__init__` (the
template stays unregistered so that a copy never collides with it).
"""
@params struct DamageTemplateParams
end

# __init__() = register_damage("Damage Template", DamageTemplateParams)

"""
    compute_model(nodes, p, damage, block, time, dt)

Calculates the damage criterion of each bond. This template has to be copied, the file renamed and edited by the user to create a new material. Additional files can be called from here using include and `import .any_module` or `using .any_module`.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `p::DamageTemplateParams`: The model parameters.
- `damage::BlockDamage`: Base and model parameters of the block.
- `block::Int64`: Block number
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
Example:
```julia
```
"""
function compute_model(nodes::AbstractVector{Int64}, p::DamageTemplateParams, damage,
                       block::Int64, time::Float64, dt::Float64)
    @info "Register your model name with register_damage in __init__()."
    @info "You can call your routine within the yaml file."
    @info "Fill the compute_model(nodes, p, damage, block, time, dt) function."
    @info "The Data_Manager, p and damage hold all you need to solve your problem on material level."
    @info "Add own files and refer to them. If a module does not exist. Add it to the project or contact the developer."
end

"""
    fields_for_local_synchronization( model::String)

Returns a user developer defined local synchronization. This happens before each model.



# Arguments

"""
function fields_for_local_synchronization(model::String)
    #download_from_cores = false
    #upload_to_cores = true
    #Data_Manager.set_local_synch(model, "Bond Forces", download_from_cores, upload_to_cores)
end
"""
    init_model(nodes::AbstractVector{Int64}, p::DamageTemplateParams, damage, block::Int64)

Inits the damage model. Should be used to init damage specific fields, etc.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `p::DamageTemplateParams`: The model parameters.
- `damage::BlockDamage`: Base and model parameters of the block.
- `block::Int64`: Block number
Example:
```julia
```
"""
function init_model(nodes::AbstractVector{Int64}, p::DamageTemplateParams, damage,
                    block::Int64)
end
end
