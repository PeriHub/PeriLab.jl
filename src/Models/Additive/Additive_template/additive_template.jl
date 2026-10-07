# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Additive_template

using .....Data_Manager

export compute_model
export additive_name
export init_model
export fields_for_local_synchronization

using ......ParameterSpec: @params, register_additive

"""
    AdditiveTemplateParams

Declare the YAML keys your additive model needs. Register it under your model name
by uncommenting `__init__` (the template stays unregistered so that a copy never
collides with it). Licensed modules register the same way; they are loaded before
the input deck is read (`Additive.load_licensed_models`).
"""
@params struct AdditiveTemplateParams
end

# __init__() = register_additive("Additive Template", AdditiveTemplateParams)

"""
    additive_name()

Gives the additive name. It is needed for comparison with the yaml input deck.

# Arguments

# Returns
- `name::String`: The name of the additive model.

Example:
```julia
println(additive_name())
"additive Template"
```
"""
function additive_name()
    return "Additive Template"
end

"""
    compute_model(nodes::AbstractVector{Int64}, p::AdditiveTemplateParams, block::Int64, time::Float64, dt::Float64)

Calculates the force densities of the additive. This template has to be copied, the file renamed and edited by the user to create a new additive. Additional files can be called from here using include and `import .any_module` or `using .any_module`.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `p::AdditiveTemplateParams`: The model parameters.
- `block::Int64`: The current block.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
Example:
```julia
```
"""
function compute_model(nodes::AbstractVector{Int64}, p::AdditiveTemplateParams, block::Int64,
                       time::Float64, dt::Float64)
    @info "Please write a additive name in additive_name()."
    @info "You can call your routine within the yaml file."
    @info "Fill the compute_model(nodes, p, block, time, dt) function."
    @info "The Data_Manager and p hold all you need to solve your problem on additive level."
    @info "add own files and refer to them. If a module does not exist. Add it to the project or contact the developer."
end

"""
    init_model(nodes, p, block)

Inits the additive model. This template has to be copied, the file renamed and edited by the user to create a new additive. Additional files can be called from here using include and `import .any_module` or `using .any_module`.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `p::AdditiveTemplateParams`: The model parameters.
- `block::Int64`: The current block.

"""
function init_model(nodes::AbstractVector{Int64}, p::AdditiveTemplateParams, block::Int64)
    @info "Please write a additive name in additive_name()."
    @info "You can call your routine within the yaml file."
    @info "Fill the compute_model(nodes, p, block, time, dt) function."
    @info "The Data_Manager and p hold all you need to solve your problem on additive level."
    @info "add own files and refer to them. If a module does not exist. Add it to the project or contact the developer."
end

"""
    fields_for_local_synchronization(model::String)

Returns a user developer defined local synchronization. This happens before each model.



# Arguments

"""
function fields_for_local_synchronization(model::String)
    #download_from_cores = false
    #upload_to_cores = true
    #Data_Manager.set_local_synch(model, "Bond Forces", download_from_cores, upload_to_cores)
end

end
