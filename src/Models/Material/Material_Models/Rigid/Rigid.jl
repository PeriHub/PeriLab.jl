# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Rigid

using .......Data_Manager

export init_model
export fe_support
export material_name
export compute_model

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
    return true
end

"""
  init_model(nodes::AbstractVector{Int64}, material_parameter::Dict{String, Any})

Initializes the material model.

# Arguments
  - `nodes::AbstractVector{Int64}`: List of block nodes.
  - `material_parameter::Dict(String, Any)`: Dictionary with material parameter.
"""
function init_model(nodes::AbstractVector{Int64},
                    material_parameter::Dict{String,Any})
    @info "Rigid material is applied. No internal forces are calculated. No deformation occurs only rigid body motion."
end

"""
    material_name()

Returns the name of the material model.
"""
function material_name()
    return "Rigid"
end

"""
    compute_model(nodes::AbstractVector{Int64}, material_parameter::Dict{String, Any}, time::Float64, dt::Float64)

Calculate the elastic bond force for each node.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `material_parameter::Dict(String, Any)`: Dictionary with material parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
"""
function compute_model(nodes::AbstractVector{Int64},
                       material_parameter::Dict{String,Any},
                       block::Int64,
                       time::Float64,
                       dt::Float64)
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
