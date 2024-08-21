# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Additive_template
export compute_additive_model
export additive_name
export init_model

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
    compute_additive_model(datamanager, nodes, additive_parameter, time, dt)

Calculates the additive model. This template has to be copied, the file renamed and edited by the user to create a new additive. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `additive parameter::Dict(String, Any)`: Dictionary with additive parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
# Returns
- `datamanager::Data_manager`: Datamanager.
Example:
```julia
  ```
"""
function compute_additive_model(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    additive_parameter::Dict,
    time::Float64,
    dt::Float64,
)
    @info "Please write a additive name in additive_name()."
    @info "You can call your routine within the yaml file."
    @info "Fill the compute_additive_model(datamanager, nodes, additive_parameter, time, dt) function."
    @info "The datamanager and additive_parameter holds all you need to solve your problem on additive level."
    @info "add own files and refer to them. If a module does not exist. Add it to the project or contact the developer."
    return datamanager
end


"""
    init_model(datamanager, nodes, additive_parameter)

Inits the additive model. This template has to be copied, the file renamed and edited by the user to create a new additive. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `additive parameter::Dict(String, Any)`: Dictionary with additive parameter.
- `block::Int64`: The current block.
# Returns
- `datamanager::Data_manager`: Datamanager.

"""
function init_model(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    additive_parameter::Dict,
    block::Int64,
)
    @info "Please write a additive name in additive_name()."
    @info "You can call your routine within the yaml file."
    @info "Fill the compute_additive_model(datamanager, nodes, additive_parameter, time, dt) function."
    @info "The datamanager and additive_parameter holds all you need to solve your problem on additive level."
    @info "add own files and refer to them. If a module does not exist. Add it to the project or contact the developer."
    return datamanager
end


end
