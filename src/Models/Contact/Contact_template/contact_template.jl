# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Contact_template
export contact_model_name
export init_contact_model
export compute_model

"""
    contact_model_name()

Gives the contact model name. It is needed for comparison with the yaml input deck.

# Arguments

# Returns
- `name::String`: The name of the contact flow model.

Example:
```julia
println(flow_name())
"Contact Template"
```
"""
function contact_model_name()
    return "Contact Template"
end

"""
   init_contact_model(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    contact_parameter::Dict,
    block::Int64,

Inits the contact model. This template has to be copied, the file renamed and edited by the user to create a new contact. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `contact_parameter::Dict(String, Any)`: Dictionary with contact parameter.
- `block::Int64`: The current block.
# Returns
- `datamanager::Data_manager`: Datamanager.

"""
function init_contact_model(datamanager::Module,
                            nodes::Union{SubArray,Vector{Int64}},
                            contact_parameter::Dict,
                            block::Int64)
    return datamanager
end

"""
    compute_model(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    contact_parameter::Dict,
    block::Int64,
    time::Float64,
    dt::Float64,
)

Calculates the contact behavior of the material. This template has to be copied, the file renamed and edited by the user to create a new flow. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `flow parameter::Dict(String, Any)`: Dictionary with flow parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
# Returns
- `datamanager::Data_manager`: Datamanager.
Example:
```julia
```
"""
function compute_model(datamanager::Module,
                       nodes::Union{SubArray,Vector{Int64}},
                       contact_parameter::Dict,
                       block::Int64,
                       time::Float64,
                       dt::Float64)
    @info "Please write a contact model name in thermal_name()."
    @info "You can call your routine within the yaml file."
    @info "Fill the compute_model(datamanager, nodes, contact_parameter, time, dt) function."
    @info "The datamanager and contact_parameter holds all you need to solve your problem on contact flow level."
    @info "Add own files and refer to them. If a module does not exist. Add it to the project or contact the developer."
    return datamanager
end

end
