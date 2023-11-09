# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Correspondence
using LinearAlgebra
using TensorOperations
include("Correspondence_Elastic.jl")
include("./Zero_Energy_Control/global_control.jl")
include("../material_basis.jl")
include("../../../Support/geometry.jl")
using .global_zero_energy_control
using .Correspondence_Elastic
using .Geometry
export compute_force
export material_name
"""
   material_name()

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
function material_name()
  return Correspondence_Elastic.correspondence_name()
end
"""
   compute_force(datamanager, nodes, material_parameter, time, dt)

   Calculates the force densities of the material. This template has to be copied, the file renamed and edited by the user to create a new material. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

   Parameters:
        - `datamanager::Data_manager`: Datamanager.
        - `nodes::Union{SubArray, Vector{Int64}}`: List of block nodes.
        - `material_parameter::Dict(String, Any)`: Dictionary with material parameter.
        - `time::Float64`: The current time.
        - `dt::Float64`: The current time step.
   Returns:
        - - `datamanager::Data_manager`: Datamanager.
   Example:
   ```julia
     ```
   """
function compute_forces(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, material_parameter::Dict, time::Float64, dt::Float64)
  rotation::Bool = false
  if "Angles" in datamanager.get_all_field_keys()
    rotation = true
    angles = datamanager.get_field("Angles")
  end
  dof = datamanager.get_dof()
  defGrad = datamanager.get_field("Deformation Gradient")
  bond_force = datamanager.get_field("Bond Forces")

  bond_damage = datamanager.get_field("Bond Damage", "NP1")

  bondGeom = datamanager.get_field("Bond Geometry")
  invShapeTensor = datamanager.get_field("Inverse Shape Tensor")

  strainN, strainNP1 = datamanager.create_node_field("Strain", Float64, "Matrix", dof)
  stressN, stressNP1 = datamanager.create_node_field("Cauchy Stress", Float64, "Matrix", dof)
  strainNP1 = Geometry.strain(nodes, defGrad, strainNP1)
  strainInc = strainNP1 - strainN

  if rotation
    stressN = rotate(nodes, dof, stressN, angles, false)
    strainInc = rotate(nodes, dof, strainInc, angles, false)
  end

  # in future this part must be changed -> using set Modules
  stressNP1, datamanager = Correspondence_Elastic.compute_stresses(datamanager, nodes, dof, material_parameter, time, dt, strainInc, stressN, stressNP1)

  if rotation
    stressNP1 = rotate(nodes, dof, stressNP1, angles, true)
  end
  bond_force[:] = calculate_bond_force(nodes, defGrad, bondGeom, bond_damage, invShapeTensor, stressNP1, bond_force)
  # general interface, because it might be a flexbile Set_modules interface in future
  datamanager = zero_energy_mode_compensation(datamanager, nodes, material_parameter, time, dt)

  return datamanager
end

"""
Global - J. Wan et al., "Improved method for zero-energy mode suppression in peridynamic correspondence model in Acta Mechanica Sinica https://doi.org/10.1007/s10409-019-00873-y
"""
function zero_energy_mode_compensation(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, material_parameter, time::Float64, dt::Float64)
  if !haskey(material_parameter, "Zero Energy Mode Control")
    return datamanager
  end
  if material_parameter["Zero Energy Mode Control"] == global_zero_energy_control.control_name()
    datamanager = global_zero_energy_control.compute_control(datamanager, nodes, material_parameter, time, dt)
  end
  return datamanager
end



function calculate_bond_force(nodes::Union{SubArray,Vector{Int64}}, defGrad, bondGeom, bond_damage, invShapeTensor, stressNP1, bond_force)
  for iID in nodes
    jacobian = det(defGrad[iID, :, :])
    if jacobian <= 1e-8
      @error "Deformation Gradient is singular and cannot be inverted.\n - Check if your mesh is 3D, but has only one layer of nodes\n - Check number of damaged bonds."
    end
    # taken from corresponcence.cxx -> computeForcesAndStresses
    invDefGrad = inv(defGrad[iID, :, :])
    piolaStress = jacobian .* invDefGrad * stressNP1[iID, :, :]
    temp = piolaStress * invShapeTensor[iID, :, :]
    for jID in eachindex(bond_force[iID][:, 1])
      bond_force[iID][jID, :] = temp * (bond_damage[iID][jID] .* bondGeom[iID][jID, 1:end-1])
    end
  end

  return bond_force
end


function rotate(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, matrix, angles, back::Bool)
  for iID in nodes
    matrix[iID, :, :] = rotate_second_order_tensor(angles[iID, :], matrix[iID, :, :], dof, back)
  end
  return matrix

end
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