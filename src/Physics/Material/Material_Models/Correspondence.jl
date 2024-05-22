# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Correspondence
using LinearAlgebra
using TensorOperations
using TimerOutputs
include("./Zero_Energy_Control/global_control.jl")
include("./Bond_Associated_Correspondence.jl")
using .Bond_Associated_Correspondence
include("../material_basis.jl")
include("../../../Support/geometry.jl")
using .Global_zero_energy_control

include("../../../Core/Module_inclusion/set_Modules.jl")
using .Set_modules
global module_list = Set_modules.find_module_files(@__DIR__, "correspondence_name")
Set_modules.include_files(module_list)
using .Geometry
export init_material_model
export material_name
export compute_forces
export init_material_model

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
  # global dof
  # global rotation
  # global angles
  if !haskey(material_parameter, "Symmetry")
    @error "Symmetry for correspondence material is missing; options are 'isotropic plane strain', 'isotropic plane stress', 'anisotropic plane stress', 'anisotropic plane stress','isotropic' and 'anisotropic'. For 3D the plane stress or plane strain option is ignored."
    return nothing
  end
  dof = datamanager.get_dof()
  datamanager.create_node_field("Strain", Float64, "Matrix", dof)
  datamanager.create_constant_node_field("Strain Increment", Float64, "Matrix", dof)
  datamanager.create_node_field("Cauchy Stress", Float64, "Matrix", dof)
  datamanager.create_node_field("von Mises Stress", Float64, 1)

  material_models = split(material_parameter["Material Model"], "+")
  material_models = map(r -> strip(r), material_models)
  #occursin("Correspondence", material_name)
  for material_model in material_models

    mod = Set_modules.create_module_specifics(material_model, module_list, "correspondence_name")
    if isnothing(mod)
      @error "No correspondence material of name " * material_model * " exists."
      return nothing
    end
    datamanager.set_model_module(material_model, mod)
    datamanager = mod.init_material_model(datamanager, nodes, material_parameter)

  end
  if haskey(material_parameter, "Bond Associated") && material_parameter["Bond Associated"]
    return Bond_Associated_Correspondence.init_material_model(datamanager, nodes, material_parameter)
  end
  material_parameter["Bond Associated"] = false
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
  return "Correspondence"
end

"""
    compute_forces(datamanager, nodes, material_parameter, time, dt, to::TimerOutput)

Calculates the force densities of the material. This template has to be copied, the file renamed and edited by the user to create a new material. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray, Vector{Int64}}`: List of block nodes.
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

  if material_parameter["Bond Associated"]
    return Bond_Associated_Correspondence.compute_forces(datamanager, nodes, material_parameter, time, dt, to)
  end

  rotation::Bool, angles = datamanager.rotation_data()
  dof = datamanager.get_dof()
  deformation_gradient = datamanager.get_field("Deformation Gradient")
  bond_force = datamanager.get_field("Bond Forces")
  bond_damage = datamanager.get_bond_damage("NP1")

  undeformed_bond = datamanager.get_field("Bond Geometry")
  inverse_shape_tensor = datamanager.get_field("Inverse Shape Tensor")

  strain_N = datamanager.get_field("Strain", "N")
  strain_NP1 = datamanager.get_field("Strain", "NP1")
  stress_N = datamanager.get_field("Cauchy Stress", "N")
  stress_NP1 = datamanager.get_field("Cauchy Stress", "NP1")
  strain_increment = datamanager.get_field("Strain Increment")
  strain_NP1 = Geometry.strain(nodes, deformation_gradient, strain_NP1)
  strain_increment[nodes, :, :] = strain_NP1[nodes, :, :] - strain_N[nodes, :, :]

  if rotation
    stress_N = rotate(nodes, dof, stress_N, angles, false)
    strain_increment = rotate(nodes, dof, strain_increment, angles, false)
  end

  material_models = split(material_parameter["Material Model"], "+")
  material_models = map(r -> strip(r), material_models)
  for material_model in material_models
    mod = datamanager.get_model_module(material_model)

    for iID in nodes
      stress_NP1, datamanager = mod.compute_stresses(datamanager, iID, dof, material_parameter, time, dt, strain_increment, stress_N, stress_NP1)
    end
  end
  if rotation
    stress_NP1 = rotate(nodes, dof, stress_NP1, angles, true)
  end
  bond_force = calculate_bond_force(nodes, deformation_gradient, undeformed_bond, bond_damage, inverse_shape_tensor, stress_NP1, bond_force)
  # general interface, because it might be a flexbile Set_modules interface in future
  datamanager = zero_energy_mode_compensation(datamanager, nodes, material_parameter, time, dt)
  return datamanager
end

"""
    zero_energy_mode_compensation(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, material_parameter::Dict, time::Float64, dt::Float64)

Global - J. Wan et al., "Improved method for zero-energy mode suppression in peridynamic correspondence model in Acta Mechanica Sinica https://doi.org/10.1007/s10409-019-00873-y

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `material_parameter::Dict(String, Any)`: Dictionary with material parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function zero_energy_mode_compensation(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, material_parameter::Dict, time::Float64, dt::Float64)
  if !haskey(material_parameter, "Zero Energy Control")
    if time == 0
      @warn "No zero energy control activated for corresponcence."
    end
    return datamanager
  end
  if material_parameter["Zero Energy Control"] == Global_zero_energy_control.control_name()
    datamanager = Global_zero_energy_control.compute_control(datamanager, nodes, material_parameter, time, dt)
  end
  return datamanager
end

"""
    calculate_bond_force(nodes::Union{SubArray,Vector{Int64}}, deformation_gradient::SubArray, undeformed_bond::SubArray, bond_damage::SubArray, inverse_shape_tensor::SubArray, stress_NP1::SubArray, bond_force::SubArray)

Calculate bond forces for specified nodes based on deformation gradients.

# Arguments
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `deformation_gradient::SubArray`: Deformation gradient.
- `undeformed_bond::SubArray`: Undeformed bond geometry.
- `bond_damage::SubArray`: Bond damage.
- `inverse_shape_tensor::SubArray`: Inverse shape tensor.
- `stress_NP1::SubArray`: Stress at time step n+1.
- `bond_force::SubArray`: Bond force.
# Returns
- `bond_force::SubArray`: Bond force.
"""
function calculate_bond_force(nodes::Union{SubArray,Vector{Int64}}, deformation_gradient::SubArray, undeformed_bond::SubArray, bond_damage::SubArray, inverse_shape_tensor::SubArray, stress_NP1::SubArray, bond_force::SubArray)
  for iID in nodes
    jacobian = det(deformation_gradient[iID, :, :])
    if jacobian <= 1e-8
      @error "Deformation Gradient is singular and cannot be inverted.\n - Check if your mesh is 3D, but has only one layer of nodes\n - Check number of damaged bonds."
    end
    # taken from corresponcence.cxx -> computeForcesAndStresses
    invdeformation_gradient::Matrix{Float64} = inv(deformation_gradient[iID, :, :])
    piolaStress::Matrix{Float64} = jacobian .* invdeformation_gradient * stress_NP1[iID, :, :]
    temp::Matrix{Float64} = piolaStress * inverse_shape_tensor[iID, :, :]

    bond_force[iID][:, :] = (bond_damage[iID] .* undeformed_bond[iID]) * temp

  end
  return bond_force
end

"""
    rotate(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, matrix::Union{SubArray,Array{Float64,3}}, angles::SubArray, back::Bool)

Rotates the matrix.

# Arguments
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `dof::Int64`: Degree of freedom.
- `matrix::Union{SubArray,Array{Float64,3}}`: Matrix.
- `angles::SubArray`: Angles.
- `back::Bool`: Back.
# Returns
- `matrix::SubArray`: Matrix.
"""
function rotate(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, matrix::Union{SubArray,Array{Float64,3}}, angles::SubArray, back::Bool)
  for iID in nodes
    matrix[iID, :, :] = rotate_second_order_tensor(angles[iID, :], matrix[iID, :, :], dof, back)
  end
  return matrix
end

"""
    rotate_second_order_tensor(angles::Union{Vector{Float64},Vector{Int64}}, tensor::Matrix{Float64}, dof::Int64, back::Bool)

Rotates the second order tensor.

# Arguments
- `angles::Union{Vector{Float64},Vector{Int64}}`: Angles.
- `tensor::Matrix{Float64}`: Second order tensor.
- `dof::Int64`: Degree of freedom.
- `back::Bool`: Back.
# Returns
- `tensor::Matrix{Float64}`: Second order tensor.
"""
function rotate_second_order_tensor(angles::Union{Vector{Float64},Vector{Int64}}, tensor::Matrix{Float64}, dof::Int64, back::Bool)
  rot = Geometry.rotation_tensor(angles)

  R = rot[1:dof, 1:dof]

  if back
    @tensor begin
      tensor[m, n] = tensor[i, j] * R[m, i] * R[n, j]
    end
  else
    @tensor begin
      tensor[m, n] = tensor[i, j] * R[i, m] * R[j, n]
    end
  end
end


end