# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Corrosion_template
export compute_model
export corrosion_name
export init_model

"""
    corrosion_name()

Gives the corrosion name. It is needed for comparison with the yaml input deck.

# Arguments

# Returns
- `name::String`: The name of the corrosion model.

Example:
```julia
println(corrosion_name())
"Corrosion Template"
```
"""
function corrosion_name()
    return "Corrosion Template"
end

"""
    compute_model(datamanager, nodes, corrosion_parameter, time, dt)

Calculates the corrosion model. This template has to be copied, the file renamed and edited by the user to create a new corrosion. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `corrosion parameter::Dict(String, Any)`: Dictionary with corrosion parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
# Returns
- `datamanager::Data_manager`: Datamanager.
Example:
```julia
  ```
"""
function compute_model(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    corrosion_parameter::Dict,
    time::Float64,
    dt::Float64,
)
    @info "Please write a corrosion name in corrosion_name()."
    @info "You can call your routine within the yaml file."
    @info "Fill the compute_model(datamanager, nodes, corrosion_parameter, time, dt) function."
    @info "The datamanager and corrosion_parameter holds all you need to solve your problem on corrosion level."
    @info "add own files and refer to them. If a module does not exist. Add it to the project or contact the developer."
    return datamanager
end


"""
    init_model(datamanager, nodes, corrosion_parameter)

Inits the corrosion model. This template has to be copied, the file renamed and edited by the user to create a new corrosion. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `corrosion parameter::Dict(String, Any)`: Dictionary with corrosion parameter.
- `block::Int64`: The current block.
# Returns
- `datamanager::Data_manager`: Datamanager.

"""
function init_model(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    corrosion_parameter::Dict,
    block::Int64,
)
    @info "Please write a corrosion name in corrosion_name()."
    @info "You can call your routine within the yaml file."
    @info "Fill the compute_model(datamanager, nodes, corrosion_parameter, time, dt) function."
    @info "The datamanager and corrosion_parameter holds all you need to solve your problem on corrosion level."
    @info "add own files and refer to them. If a module does not exist. Add it to the project or contact the developer."
    return datamanager
end


end
