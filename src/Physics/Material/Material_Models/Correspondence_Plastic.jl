# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Correspondence_Plastic
include("../material_basis.jl")
using LinearAlgebra
export fe_support
export init_material_model
export init_material_model
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
  init_material_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, material_parameter::Dict)

Initializes the material model.

# Arguments
  - `datamanager::Data_manager`: Datamanager.
  - `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
  - `material_parameter::Dict(String, Any)`: Dictionary with material parameter.

# Returns
  - `datamanager::Data_manager`: Datamanager.
"""
function init_material_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, material_parameter::Dict)
  if !haskey(material_parameter, "Shear Modulus")
    @error "Shear Modulus must be defined to be able to run this plastic material"
    return nothing
  end
  if !haskey(material_parameter, "Yield Stress")
    @error "No ''Yield Stress'' is defined."
    return nothing
  end

  # datamanager.create_node_field("von Mises Stress", Float64, 1)
  datamanager.create_node_field("Plastic Strain", Float64, 1)

  return datamanager
end

"""
   correspondence_name()

   Gives the material name. It is needed for comparison with the yaml input deck.

   Parameters:

   Returns:
   - `name::String`: The name of the material.

   Example:
   ```julia
   println(material_name())
   "Material Template"
   ```
   """
function correspondence_name()
  return "Correspondence Plastic"
end
"""
   compute_force(datamanager, nodes, material_parameter, time, dt)

   Calculates the force densities of the material. This template has to be copied, the file renamed and edited by the user to create a new material. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

   Parameters:
        - `datamanager::Data_manager`: Datamanager.
        - `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
        - `material_parameter::Dict(String, Any)`: Dictionary with material parameter.
        - `time::Float64`: The current time.
        - `dt::Float64`: The current time step.
   Returns:
        - - `datamanager::Data_manager`: Datamanager.
   Example:
   ```julia
     ```
   """
function compute_stresses(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, dof::Int64, material_parameter::Dict, time::Float64, dt::Float64, strain_increment::SubArray, stress_N::SubArray, stress_NP1::SubArray)

  von_Mises_stress = datamanager.get_field("von Mises Stress", "NP1")
  plastic_strain_N = datamanager.get_field("Plastic Strain", "N")
  plastic_strain_NP1 = datamanager.get_field("Plastic Strain", "NP1")
  coordinates = datamanager.get_field("Coordinates")
  yield_stress::Float64 = material_parameter["Yield Stress"]
  spherical_stress_N::Float64 = 0
  deviatoric_stress_N::Matrix{Float64} = zeros(dof, dof)

  spherical_stress_NP1::Float64 = 0
  deviatoric_stress_NP1::Matrix{Float64} = zeros(dof, dof)

  temp_A::Matrix{Float64} = zeros(dof, dof)
  temp_B::Matrix{Float64} = zeros(dof, dof)

  sqrt23::Float64 = sqrt(2 / 3)
  for iID in nodes
    von_Mises_stress[iID] = get_von_mises_stress(stress_NP1[iID, :, :], dof)
    reduced_yield_stress = yield_stress

    reduced_yield_stress = flaw_function(material_parameter, coordinates[iID, :], yield_stress)
    if von_Mises_stress[iID] < reduced_yield_stress
      # material is elastic and nothing happens
      plastic_strain_NP1[iID] = plastic_strain_N[iID]
      continue
    end
    deviatoric_stress_magnitude_NP1 = maximum([1.0e-20, von_Mises_stress[iID] / sqrt23])
    deviatoric_stress_NP1 .*= sqrt23 * reduced_yield_stress / deviatoric_stress_magnitude_NP1
    stress_NP1[iID, :, :] = deviatoric_stress_NP1 + spherical_stress_NP1 .* I(dof)
    von_Mises_stress[iID] = sqrt(3.0 / 2.0 * sum(deviatoric_stress_NP1[:, :]))
    #############################
    # comment taken from Peridigm elastic_plastic_correspondence.cxx
    #############################
    # Update the equivalent plastic strain
    #
    # The algorithm below is generic and should not need to be modified for
    # any J2 plasticity yield surface.  It uses the difference in the yield
    # surface location at the NP1 and N steps to increment eqps regardless
    # of how the plastic multiplier was found in the yield surface
    # evaluation.
    #
    # First go back to step N and compute deviatoric stress and its
    # magnitude.  We didn't do this earlier because it wouldn't be necassary
    # if the step is elastic.

    spherical_stress_N = sum(stress_N[iID, i, i] for i in 1:dof) / 3
    deviatoric_stress_N = stress_N[iID, :, :] - spherical_stress_N .* I(dof)

    deviatoric_stress_magnitude_N = maximum([1.0e-20, sqrt(sum(deviatoric_stress_N .* deviatoric_stress_N))])
    # Contract the two tensors. This represents a projection of the plastic
    # strain increment tensor onto the "direction" of deviatoric stress
    # increment
    temp_A = (deviatoric_stress_NP1 - deviatoric_stress_N) ./ 2 / material_parameter["Shear Modulus"]
    temp_B = (deviatoric_stress_NP1 ./ deviatoric_stress_magnitude_NP1 + deviatoric_stress_N ./ deviatoric_stress_magnitude_N) ./ 2
    temp_scalar = sum(temp_A .* temp_B)
    plastic_strain_NP1[iID] = plastic_strain_N[iID] + maximum([0, sqrt23 * temp_scalar])
  end

  return stress_NP1, datamanager
end



end