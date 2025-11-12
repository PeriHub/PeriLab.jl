# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module FEM_template

using ..Data_Manager
#export init_element
#export compute_element
#export element_name
#export shape_function
"""
  element_name()

Gives the element name. It is needed for comparison with the yaml input deck.

# Arguments

# Returns
- `name::String`: The name of the element model.

Example:
```julia
println(element_name())
"element Template"
```
"""
function element_name()
    return "element Template"
end

function init_element(elements::AbstractVector{Int64},
                      element_params::Dict,
                      p::Vector{Int64})
end
"""
  compute_element(nodes, element_parameter, time, dt)

Calculates element model of the material. This template has to be copied, the file renamed and edited by the user to create a new element. Additional files can be called from here using include and `import .any_module` or `using .any_module`.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `element parameter::Dict(String, Any)`: Dictionary with element parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
Example:
```julia
  ```
"""
function compute_element(nodes::AbstractVector{Int64},
                         element_parameter::Dict,
                         time::Float64,
                         dt::Float64)
    @info "Please write a element name in element_name()."
    @info "You can call your routine within the yaml file."
    @info "Fill the compute_element(nodes, element_parameter, time, dt) function."
    @info "The Data_Manager and element_parameter holds all you need to solve your problem on element level."
    @info "add own files and refer to them. If a module does not exist. Add it to the project or contact the developer."
end

end
