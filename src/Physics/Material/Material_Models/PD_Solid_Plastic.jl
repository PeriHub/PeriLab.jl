# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module PD_Solid_Elastic_Plastic
using TimerOutputs
using StaticArrays
import .Ordinary

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

  if !haskey(material_parameter, "Yield Stress")
    @error "Yield Stress is not defined in input deck"
    return nothing
  end
  yield_stress = material_parameter["Yield Stress"]
  yield = datamanager.create_constant_node_field("Yield Value", Float64, 1)
  datamanager.create_constant_bond_field("Deviatoric Plastic Extension State", Float64, 1)
  if get_symmmetry(material_parameter) == "3D"
    yield = 25 * yield_stress * yield_stress ./ (8 * pi .* horizon .^ 5)
  else
    yield = 225 * yield_stress * yield_stress ./ (24 * thickness * pi ./ horizon .^ 4)
  end

  bond_force_deviatoric_part = datamanager.create_constant_bond_field("Bond Forces Deviatoric", Float64, 1)
  bond_force_isotropic_part = datamanager.create_constant_bond_field("Bond Forces Isotropic", Float64, 1)

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
  return "PD Solid Elastic Plastic"
end

"""
    compute_forces(datamanager, nodes, material_parameter, time, dt, to::TimerOutput)

Calculates the force densities of the material. This template has to be copied, the file renamed and edited by the user to create a new material. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `material_parameter::Dict(String, Any)`: Dictionary with material parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
# Returns
- `datamanager::Data_manager`: Datamanager.
Example:
```julia
```
"""
function compute_forces(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, material_parameter::Dict, time::Float64, dt::Float64, to::TimerOutput)
  horizon = datamanager.get_field("Horizon")
  symmetry::String = get_symmmetry(material_parameter)
  undeformed_bond = datamanager.get_field("Bond Geometry")
  deformed_bond = datamanager.get_field("Deformed Bond Geometry", "NP1")
  bond_force_elastic = datamanager.get_field
  bond_damage = datamanager.get_bond_damage("NP1")
  bond_force = datamanager.get_field("Bond Forces")
  shear_modulus = material_parameter["Shear Modulus"]
  bulk_modulus = material_parameter["Bulk Modulus"]
  bond_force_deviatoric_part = datamanager.get_field("Bond Forces Deviatoric")
  bond_force_isotropic_part = datamanager.get_field("Bond Forces Isotropic")

  alpha, gamma, kappa = Ordinary.calculate_symmetry_params(symmetry, shear_modulus, bulk_modulus)
  td_norm = compute_deviatoric_force_state_norm(nodes, nlist, alpha, bond_force_deviatoric_part, bond_damage, omega, volume, deviatoric_plastic_extension_state)
  return datamanager
end




function compute_deviatoric_force_state_norm(nodes::nodes::Union{SubArray,Vector{Int64}}, nlist, alpha::Float64, bond_force_deviatoric::SubArray, bond_damage::SubArray, omega::SubArray, volume::SubArray, deviatoric_plastic_extension_state)
  # not optimal allocation of memory, but not check of indices is needed
  norm = @MMatrix zeros(Float64, maximum(nodes))
  for iID in nodes
    td_trial = bond_force_deviatoric[iID] - alpha .* bond_damage[iID][:] * omega[iID][:] * deviatoric_plastic_extension_state[iID][:]
    norm[iID] = sum(td_trial * td_trial * volume[nlist[iID]])
  end

  return sqrt.(norm)
end


function plastic()


  for iID in nodes


    if td_norm[iID] * td_norm[iID] / 2 - yield_value[iID] < 0
      lambdaNP1[iID] = lambdaN[iID]

    end
    deltaLambda = (td_norm / sqrt(2.0 * yield_value[iID]) - 1.0) / alpha
    lambdaNP1[iID] = lambdaN[iID] + delta_lambda

    for jID in nneighbors[jID]
      td_trial = bond_force_elastic[iID] - alpha .* bond_damage[iID][:] * omega[iID][:] * deviatoric_plastic_extension_state[iID][:]
      bond_force_deviatoric[iID][jID] = sqrt(2.0 * yield_value) * td_trial / td_norm


      # Update deviatoric plastic deformation state

      deviatoric_plastic_extension_state += bond_force_deviatoric_part[iID][jID] * deltaLambda


      # Compute isotropic part of force state

      #ti = c * zeta
      bond_force[iID][jID] = bond_damage[iID][jID] * (bond_force_isotropic_part[iID][jID] + bond_force_deviatoric_part[iID][jID])

    end
    return bond_force, deviatoric_plastic_extension_state
  end
end
end