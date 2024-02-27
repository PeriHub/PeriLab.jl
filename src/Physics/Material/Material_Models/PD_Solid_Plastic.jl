# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module PD_Solid_Elastic_Plastic
using TimerOutputs
using StaticArrays
include("./Ordinary/Ordinary.jl")
using .Ordinary: calculate_symmetry_params, get_bond_forces

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
  datamanger.create_node_field("Lambda Plastic", Float64, 1)
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
  volume = datamanager.get_field("Volume")
  symmetry::String = get_symmmetry(material_parameter)
  deformed_bond = datamanager.get_field("Deformed Bond Geometry", "NP1")
  omega = datamanager.get_field("Influence Function")
  bond_damage = datamanager.get_bond_damage("NP1")
  bond_force = datamanager.get_field("Bond Forces")
  shear_modulus = material_parameter["Shear Modulus"]
  bulk_modulus = material_parameter["Bulk Modulus"]
  bond_force_deviatoric_part = datamanager.get_field("Bond Forces Deviatoric")
  bond_force_isotropic_part = datamanager.get_field("Bond Forces Isotropic")
  deviatoric_plastic_extension_state = datamanager.get_field("Deviatoric Plastic Extension State")
  yield_value = datamanager.get_field("Yield Value")
  lambdaN = datamanager.get_field("Lambda Plastic", "N")
  lambdaNP1 = datamanager.get_field("Lambda Plastic", "NP1")

  alpha, gamma, kappa = Ordinary.calculate_symmetry_params(symmetry, shear_modulus, bulk_modulus)
  td_norm = compute_deviatoric_force_state_norm(nodes, nlist, alpha, bond_force_deviatoric_part, bond_damage, omega, volume, deviatoric_plastic_extension_state)
  lambdaNP1 = copy(lambdaN)

  bond_force_deviatoric_part, deviatoric_plastic_extension_state = plastic(nodes, td_norm, yield_value, lambdaNP1, alpha, omega, bond_damage, deviatoric_plastic_extension_state, bond_force_deviatoric_part)
  bond_force = get_bond_forces(nodes, bond_force_deviatoric_part + bond_force_isotropic_part, deformed_bond, bond_force)
  return datamanager
end

"""
    compute_deviatoric_force_state_norm(nodes::Union{SubArray,Vector{Int64}}, 
                                        nlist::SubArray, 
                                        alpha::Float64, 
                                        bond_force_deviatoric::SubArray, 
                                        bond_damage::SubArray, 
                                        omega::SubArray, 
                                        volume::SubArray, 
                                        deviatoric_plastic_extension_state)

Compute the norm of the deviatoric force state for each node.

# Arguments
- `nodes::Union{SubArray,Vector{Int64}}`: Vector of node indices or a subarray representing the indices of the nodes.
- `nlist::SubArray`: Subarray representing the neighbor list for each node.
- `alpha::Float64`: Alpha parameter.
- `bond_force_deviatoric::SubArray`: Subarray representing the deviatoric bond forces.
- `bond_damage::SubArray`: Subarray representing the bond damage.
- `omega::SubArray`: Subarray representing the weights for each bond.
- `volume::SubArray`: Subarray representing the volume for each node.
- `deviatoric_plastic_extension_state`: Deviatoric plastic extension state.

# Returns
- `td_norm::Vector{Float64}`: Vector containing the norm of the deviatoric force state for each node.
"""

function compute_deviatoric_force_state_norm(nodes::Union{SubArray,Vector{Int64}}, nlist::SubArray, alpha::Float64, bond_force_deviatoric::SubArray, bond_damage::SubArray, omega::SubArray, volume::SubArray, deviatoric_plastic_extension_state)
  # not optimal allocation of memory, but not check of indices is needed
  td_norm = @MMatrix zeros(Float64, maximum(nodes))
  for iID in nodes
    td_trial = bond_force_deviatoric[iID] - alpha .* bond_damage[iID][:] * omega[iID][:] * deviatoric_plastic_extension_state[iID][:]
    td_norm[iID] = sum(td_trial * td_trial * volume[nlist[iID]])
  end

  return sqrt.(td_norm)
end

"""
    plastic(nodes::Union{SubArray,Vector{Int64}}, 
            td_norm::Vector{Float64}, 
            yield_value::Vector{Float64}, 
            lambdaNP1::Vector{Float64}, 
            alpha::Float64, 
            omega::SubArray, 
            bond_damage::SubArray, 
            deviatoric_plastic_extension_state::SubArray, 
            bond_force_deviatoric::SubArray)

Update the plastic state based on the deviatoric force norm.

# Arguments
- `nodes::Union{SubArray,Vector{Int64}}`: Vector of node indices or a subarray representing the indices of the nodes.
- `td_norm::Vector{Float64}`: Vector containing the norm of the deviatoric force state for each node.
- `yield_value::Vector{Float64}`: Vector containing the yield values for each node.
- `lambdaNP1::Vector{Float64}`: Vector containing the plastic multipliers for each node.
- `alpha::Float64`: Alpha parameter.
- `omega::SubArray`: Subarray representing the weights for each bond.
- `bond_damage::SubArray`: Subarray representing the bond damage.
- `deviatoric_plastic_extension_state::SubArray`: Subarray representing the deviatoric plastic extension state.
- `bond_force_deviatoric::SubArray`: Subarray representing the deviatoric bond forces.

# Returns
- `bond_force_deviatoric::SubArray`: Updated deviatoric bond forces.
- `deviatoric_plastic_extension_state::SubArray`: Updated deviatoric plastic extension state.
"""

function plastic(nodes::Union{SubArray,Vector{Int64}}, td_norm, yield_value::SubArray, lambdaNP1::SubArray, alpha::Float64, omega::SubArray, bond_damage::SubArray, deviatoric_plastic_extension_state::SubArray, bond_force_deviatoric::SubArray)


  for iID in nodes
    if td_norm[iID] * td_norm[iID] / 2 - yield_value[iID] < 0
      continue
    end
    delta_lambda = (td_norm / sqrt(2.0 * yield_value[iID]) - 1.0) / alpha
    lambdaNP1[iID] .+= delta_lambda
    td_trial = bond_force_deviatoric[iID] - alpha .* bond_damage[iID][:] * omega[iID][:] * deviatoric_plastic_extension_state[iID][:]
    bond_force_deviatoric[iID] = sqrt(2.0 * yield_value) * td_trial ./ td_norm[iID]
    deviatoric_plastic_extension_state += bond_force_deviatoric_part[iID][:] .* delta_lambda

    return bond_force_deviatoric, deviatoric_plastic_extension_state
  end
end
end