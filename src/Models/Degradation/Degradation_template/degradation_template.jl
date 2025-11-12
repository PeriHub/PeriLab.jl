# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Degradation_template

using .......Data_Manager

export compute_model
export degradation_name
export init_model
export fields_for_local_synchronization
"""
    degradation_name()

Gives the degradation name. It is needed for comparison with the yaml input deck.

# Arguments

# Returns
- `name::String`: The name of the degradation model.

Example:
```julia
println(degradation_name())
"Degradation Template"
```
"""
function degradation_name()
    return "Degradation Template"
end

"""
    compute_model(nodes, degradation_parameter, block::Int64, time, dt)

Calculates the degradation model. This template has to be copied, the file renamed and edited by the user to create a new degradation. Additional files can be called from here using include and `import .any_module` or `using .any_module`.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `degradation parameter::Dict(String, Any)`: Dictionary with degradation parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
Example:
```julia
  ```
"""
function compute_model(nodes::AbstractVector{Int64},
                       degradation_parameter::Dict,
                       block::Int64,
                       time::Float64,
                       dt::Float64)
    @info "Please write a degradation name in degradation_name()."
    @info "You can call your routine within the yaml file."
    @info "Fill the compute_model(nodes, degradation_parameter, time, dt) function."
    @info "The Data_Manager and degradation_parameter holds all you need to solve your problem on degradation level."
    @info "add own files and refer to them. If a module does not exist. Add it to the project or contact the developer."
end

"""
    init_model(nodes, degradation_parameter)

Inits the degradation model. This template has to be copied, the file renamed and edited by the user to create a new degradation. Additional files can be called from here using include and `import .any_module` or `using .any_module`.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `degradation parameter::Dict(String, Any)`: Dictionary with degradation parameter.
- `block::Int64`: The current block.

"""
function init_model(nodes::AbstractVector{Int64},
                    degradation_parameter::Dict,
                    block::Int64)
    @info "Please write a degradation name in degradation_name()."
    @info "You can call your routine within the yaml file."
    @info "Fill the compute_model(nodes, degradation_parameter, time, dt) function."
    @info "The Data_Manager and degradation_parameter holds all you need to solve your problem on degradation level."
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
