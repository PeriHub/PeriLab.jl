# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Correspondence_template

using .......Data_Manager
using .......PeriLabExceptions: @abort
using .......ParameterSpec: @params, register_material
export compute_stresses
export correspondence_name
export fe_support
export init_model
export fields_for_local_synchronization

"""
    CorrespondenceTemplateParams

Declare the YAML keys your model needs beyond the shared material keys (Symmetry,
moduli, … are in `material.base` / `material.moduli`). Register it under your model
name by uncommenting `__init__`; the name must contain "Correspondence" (that selects
the correspondence formulation). The template stays unregistered so that a copy never
collides with it.
"""
@params struct CorrespondenceTemplateParams
end

# __init__() = register_material("Correspondence Template", CorrespondenceTemplateParams)

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
  init_model(nodes::AbstractVector{Int64}, p::CorrespondenceTemplateParams, material)

Initializes the material model.

# Arguments
  - `nodes::AbstractVector{Int64}`: List of block nodes.
  - `p::CorrespondenceTemplateParams`: The model parameters.
  - `material::BlockMaterial`: The typed block material (base, moduli, symmetry).
"""
function init_model(nodes::AbstractVector{Int64}, p::CorrespondenceTemplateParams, material)
end

"""
    correspondence_name()

Gives the correspondence material name. PeriLab loads the module because it defines this function; the input deck uses the name passed to `register_*` in `__init__()`.

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
    compute_stresses(nodes::AbstractVector{Int64}, dof::Int64, p::CorrespondenceTemplateParams, material, time::Float64, dt::Float64, strain_increment, stress_N, stress_NP1)

Calculates the stresses of the material. This template has to be copied, the file renamed and edited by the user to create a new material. Additional files can be called from here using include and `import .any_module` or `using .any_module`.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `dof::Int64`: Degrees of freedom
- `p::CorrespondenceTemplateParams`: The model parameters.
- `material::BlockMaterial`: The typed block material (base, moduli, symmetry).
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
- `strain_increment`: Strain increment.
- `stress_N`: Stress of step N.
- `stress_NP1`: Stress of step N+1.
# Returns
- `stress_NP1`: updated stresses

Example:
```julia
```
"""
function compute_stresses(nodes::AbstractVector{Int64},
                          dof::Int64,
                          p::CorrespondenceTemplateParams,
                          material,
                          time::Float64,
                          dt::Float64,
                          strain_increment::AbstractArray{Float64},
                          stress_N::AbstractArray{Float64},
                          stress_NP1::AbstractArray{Float64})
    @info "Please write a material name in material_name()."
    @info "You can call your routine within the yaml file."
    @info "Fill the compute_model() and init_model() function."
    @info "The Data_Manager, p and material hold all you need to solve your problem on material level."
    @info "Add own files and refer to them. If a module does not exist. Add it to the project or contact the developer."
    return stress_NP1
end

function compute_stresses_ba(nodes,
                             nlist,
                             dof::Int64,
                             p::CorrespondenceTemplateParams,
                             material,
                             time::Float64,
                             dt::Float64,
                             strain_increment,
                             stress_N,
                             stress_NP1)
    @abort "$(correspondence_name()) not yet implemented for bond associated."
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
