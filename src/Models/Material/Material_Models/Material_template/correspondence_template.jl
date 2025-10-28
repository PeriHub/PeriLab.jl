# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Correspondence_template

using .......Data_Manager
export compute_stresses
export correspondence_name
export fe_support
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
  init_model(nodes::AbstractVector{Int64}, material_parameter::Dict,
    block::Int64)

Initializes the material model.

# Arguments
  - `nodes::AbstractVector{Int64}`: List of block nodes.
  - `material_parameter::Dict(String, Any)`: Dictionary with material parameter.
  - `block::Int64`: Current block
"""
function init_model(nodes::AbstractVector{Int64},
                    material_parameter::Dict,
                    block::Int64)
end

"""
    correspondence_name()

Gives the correspondence material name. It is needed for comparison with the yaml input deck.

# Arguments

# Returns
- `name::String`: The name of the material.

Example:
```julia
println(correspondence_name())
"Material Template"
```
"""
function correspondence_name()
    return "Correspondence Template"
end

"""
    compute_stresses(nodes::AbstractVector{Int64}, dof::Int64, material_parameter::Dict, time::Float64, dt::Float64, strain_increment::SubArray, stress_N::SubArray, stress_NP1::SubArray, iID_jID_nID::Tuple=())

Calculates the stresses of the material. This template has to be copied, the file renamed and edited by the user to create a new material. Additional files can be called from here using include and `import .any_module` or `using .any_module`.

# Arguments
- `iID::Int64`: Node ID.
- `dof::Int64`: Degrees of freedom
- `material_parameter::Dict(String, Any)`: Dictionary with material parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
- `strainInc::Union{Array{Float64,3},Array{Float64,6}}`: Strain increment.
- `stress_N::SubArray`: Stress of step N.
- `stress_NP1::SubArray`: Stress of step N+1.
- `iID_jID_nID::Tuple=(): (optional) are the index and node id information. The tuple is ordered iID as index of the point,  jID the index of the bond of iID and nID the neighborID.
# Returns
- `stress_NP1::SubArray`: updated stresses

Example:
```julia
```
"""
function compute_stresses(iID::Int64,
                          dof::Int64,
                          material_parameter::Dict,
                          time::Float64,
                          dt::Float64,
                          strain_increment::AbstractArray{Float64},
                          stress_N::AbstractArray{Float64},
                          stress_NP1::AbstractArray{Float64},
                          iID_jID_nID::Tuple = ())
    @info "Please write a material name in material_name()."
    @info "You can call your routine within the yaml file."
    @info "Fill the compute_model() and init_model() function."
    @info "The Data_Manager and material_parameter holds all you need to solve your problem on material level."
    @info "Add own files and refer to them. If a module does not exist. Add it to the project or contact the developer."
    return stress_NP1
end

"""
    compute_stresses(dof::Int64, material_parameter::Dict, time::Float64, dt::Float64, strain_increment::SubArray, stress_N::SubArray, stress_NP1::SubArray)

Calculates the stresses of a single node. Needed for FEM. This template has to be copied, the file renamed and edited by the user to create a new material. Additional files can be called from here using include and `import .any_module` or `using .any_module`.

# Arguments
- `dof::Int64`: Degrees of freedom
- `material_parameter::Dict(String, Any)`: Dictionary with material parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
- `strainInc::Union{Array{Float64,3},Array{Float64,6}}`: Strain increment.
- `stress_N::SubArray`: Stress of step N.
- `stress_NP1::SubArray`: Stress of step N+1.
# Returns
- `stress_NP1::SubArray`: updated stresses
Example:
```julia
```
"""
function compute_stresses(dof::Int64,
                          material_parameter::Dict,
                          time::Float64,
                          dt::Float64,
                          strain_increment::Vector{Float64},
                          stress_N::Vector{Float64},
                          stress_NP1::Vector{Float64})
    return stress_NP1
end

function compute_stresses_ba(nodes,
                             nlist,
                             dof::Int64,
                             material_parameter::Dict,
                             time::Float64,
                             dt::Float64,
                             strain_increment::Union{SubArray,Array{Float64,3},
                                                     Vector{Float64}},
                             stress_N::Union{SubArray,Array{Float64,3},Vector{Float64}},
                             stress_NP1::Union{SubArray,Array{Float64,3},
                                               Vector{Float64}})
    @error "$(correspondence_name()) not yet implemented for bond associated."
end

"""
    fields_for_local_synchronization(model::String)

Returns a user developer defined local synchronization. This happens before each model.



# Arguments

"""
function fields_for_local_synchronization(model::String)
    # download_from_cores = false
    # upload_to_cores = true
    # Data_Manager.set_local_synch(model, "Bond Forces", download_from_cores, upload_to_cores)
end

end
