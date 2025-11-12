# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Short_Range

using .....Data_Manager

export contact_model_name
export init_contact_model
export compute_model

"""
    contact_model_name()

Gives the contact model name. It is needed for comparison with the yaml input deck.

# Arguments

# Returns
- `name::String`: The name of the contact model.

Example:
```julia
println(contact_model_name())
"Short Range"
```
"""
function contact_model_name()
    return "Short Range"
end

"""
  init_contact_model(
    nodes::AbstractVector{Int64},
    contact_parameter::Dict,
    block::Int64,
)

Inits the contact model. This template has to be copied, the file renamed and edited by the user to create a new contact. Additional files can be called from here using include and `import .any_module` or `using .any_module`.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `contact_parameter::Dict(String, Any)`: Dictionary with contact parameter.
- `block::Int64`: The current block.

"""
function init_contact_model(nodes::AbstractVector{Int64},
                            contact_parameter::Dict,
                            block::Int64)
end

"""
    compute_model(nodes, contact_parameter, block::Int64,, time, dt)

Not yet implemented short range contact model.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `contact_parameter::Dict(String, Any)`: Dictionary with flow parameter.
- `block::Int64`: The current block.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
Example:
```julia
```
"""
function compute_model(nodes::AbstractVector{Int64},
                       contact_parameter::Dict,
                       block::Int64,
                       time::Float64,
                       dt::Float64)
    @info "Please write a contact model name in thermal_name()."
    @info "You can call your routine within the yaml file."
    @info "Fill the compute_model(nodes, contact_parameter, block, time, dt) function."
    @info "The Data_Manager and contact_parameter holds all you need to solve your problem on contact flow level."
    @info "Add own files and refer to them. If a module does not exist. Add it to the project or contact the developer."
end

end
