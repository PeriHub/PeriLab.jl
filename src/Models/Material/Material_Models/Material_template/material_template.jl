# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Material_template
using TimerOutputs
export fe_support
export init_model
export material_name
export compute_model
export init_model
export fields_for_local_synchronization
"""
  fe_support()

Gives the information if the material supports the FEM part of PeriLab

# Arguments

# Returns
- bool: true - for FEM support; false - for no FEM support

Example:
```julia
println(fe_support())
false
```
"""
function fe_support()
    return false
end

"""
  init_model(datamanager::Module, nodes::AbstractVector{Int64}, material_parameter::Dict)

Initializes the material model.

# Arguments
  - `datamanager::Data_Manager`: Datamanager.
  - `nodes::AbstractVector{Int64}`: List of block nodes.
  - `material_parameter::Dict(String, Any)`: Dictionary with material parameter.

# Returns
  - `datamanager::Data_Manager`: Datamanager.
"""
function init_model(datamanager::Module,
                    nodes::AbstractVector{Int64},
                    material_parameter::Dict)
    return datamanager
end

"""
    material_name()

Gives the material name. It is needed for comparison with the yaml input deck.

# Arguments

# Returns
- `name::String`: The name of the material.

Example:
```julia
println(material_name())
"Material Template"
```
"""
function material_name()
    return "Material Template"
end

"""
    fields_for_local_synchronization(datamanager::Module, model::String)

Returns a user developer defined local synchronization. This happens before each model.



# Arguments

"""
function fields_for_local_synchronization(datamanager::Module, model::String)
    # download_from_cores = false
    # upload_to_cores = true
    # datamanager.set_local_synch(model, "Bond Forces", download_from_cores, upload_to_cores)
    return datamanager
end

"""
    compute_model(datamanager, nodes, material_parameter, time, dt, to::TimerOutput)

Calculates the force densities of the material. This template has to be copied, the file renamed and edited by the user to create a new material. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

# Arguments
- `datamanager::Data_Manager`: Datamanager.
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `material_parameter::Dict(String, Any)`: Dictionary with material parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
# Returns
- `datamanager::Data_Manager`: Datamanager.
Example:
```julia
```
"""
function compute_model(datamanager::Module,
                       nodes::AbstractVector{Int64},
                       material_parameter::Dict,
                       block::Int64,
                       time::Float64,
                       dt::Float64,
                       to::TimerOutput)
    @info "Please write a material name in material_name()."
    @info "You can call your routine within the yaml file."
    @info "Fill the compute_model() and init_model() function."
    @info "The datamanager and material_parameter holds all you need to solve your problem on material level."
    @info "Add own files and refer to them. If a module does not exist. Add it to the project or contact the developer."
    return datamanager
end

end
