# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Correspondence_Elastic
include("../material_basis.jl")
export compute_stresses
export correspondence_name
export fe_support
export init_material_model
export material_name
export compute_forces

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
  init_material_model(datamanager::Module)

Initializes the material model.

# Arguments
  - `datamanager::Data_manager`: Datamanager.

# Returns
  - `datamanager::Data_manager`: Datamanager.
"""
function init_material_model(datamanager::Module)
   return datamanager
end
"""
   correspondence_name()

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
function correspondence_name()
   return "Correspondence Elastic"
end

"""
   compute_stresses(datamanager, nodes, material_parameter, time, dt)

Calculates the force densities of the material. This template has to be copied, the file renamed and edited by the user to create a new material. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `dof::Int64`: Degrees of freedom
- `material_parameter::Dict(String, Any)`: Dictionary with material parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
- `strainInc::Union{Array{Float64,3},Array{Float64,6}}`: Strain increment.
- `stress_N::SubArray`: Stress of step N.
- `stress_NP1::SubArray`: Stress of step N+1.
# Returns
- `datamanager::Data_manager`: Datamanager.
- `stress_NP1::SubArray`: updated stresses
Example:
```julia
```
"""
function compute_stresses(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, dof::Int64, material_parameter::Dict, time::Float64, dt::Float64, strain_increment::SubArray, stress_N::SubArray, stress_NP1::SubArray)

   hookeMatrix = get_Hooke_matrix(material_parameter, material_parameter["Symmetry"], dof)

   for iID in nodes
      stress_NP1[iID, :, :] = voigt_to_matrix(hookeMatrix * matrix_to_voigt(strain_increment[iID, :, :])) + stress_N[iID, :, :]
   end

   return stress_NP1, datamanager
end
"""
for FEM
"""

function compute_stresses(datamanager::Module, dof::Int64, material_parameter::Dict, time::Float64, dt::Float64, strain_increment::Vector{Float64}, stress_N::Vector{Float64}, stress_NP1::Vector{Float64})

   hookeMatrix = get_Hooke_matrix(material_parameter, material_parameter["Symmetry"], dof)

   return hookeMatrix * strain_increment + stress_N, datamanager
end


end