# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Correspondence
using LinearAlgebra
include("Correspondence_Elastic.jl")
include("../material_basis.jl")
include("../../../Support/geometry.jl")
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
        - `nodes::Vector{Int64}`: List of block nodes.
        - `material_parameter::Dict(String, Any)`: Dictionary with material parameter.
        - `time::Float32`: The current time.
        - `dt::Float32`: The current time step.
   Returns:
        - - `datamanager::Data_manager`: Datamanager.
   Example:
   ```julia
     ```
   """
function compute_forces(datamanager::Module, nodes::Vector{Int64}, material_parameter, time::Float32, dt::Float32)
  rotation::Bool = false
  if "Angles" in datamanager.get_all_field_keys()
    rotation = true
    angles = datamanager.get_field("Angles")
  end
  dof = datamanager.get_dof()
  nlist = datamanager.get_nlist()
  volume = datamanager.get_field("Volume")
  strainInc = datamanager.create_constant_node_field("Strain Increment", Float32, "Matrix", dof)
  defGradN = datamanager.get_field("Deformation Gradient", "N")
  defGradNP1 = datamanager.get_field("Deformation Gradient", "NP1")
  stressN, stressNP1 = datamanager.create_node_field("Cauchy Stress", Float32, "Matrix", dof)
  bond_force = datamanager.create_constant_bond_field("Bond Forces", Float32, dof)
  force_densities = datamanager.get_field("Force Densities", "NP1")
  bondGeom = datamanager.get_field("Bond Geometry")
  invShapeTensor = datamanager.get_field("Inverse Shape Tensor")
  strainInc = Geometry.strain_increment(nodes, defGradNP1, defGradN, strainInc)

  if rotation
    @info "not implemented"
    # rotation of strain increments and stresses N
  end

  # in future this part can be changed
  stressNP1, datamanager = Correspondence_Elastic.compute_stresses(datamanager, nodes, dof, material_parameter, time, dt, strainInc, stressN, stressNP1)
  bond_force[:] = calculate_bond_force(nodes, defGradNP1, bondGeom, invShapeTensor, stressNP1, bond_force)

  if rotation
    # rotation of stresses NP1
  end

  force_densities[:] = distribute_forces(nodes, nlist, bond_force, volume, force_densities)

  return datamanager
end

function calculate_bond_force(nodes::Vector{Int64}, defGrad, bondGeom, invShapeTensor, stressNP1, bond_force)
  for iID in nodes
    jacobian = det(defGrad[iID, :, :])
    if jacobian <= 1e-8
      @error "Deformation Gradient is singular and cannot be inverted.\n - Check if your mesh is 3D, but has only one layer of nodes\n - Check number of damaged bonds."
    end
    invDefGrad = inv(defGrad[iID, :, :])

    bond_force[iID][:, :] = transpose(jacobian .* invDefGrad * stressNP1[iID, :, :] * invShapeTensor[iID, :, :] * transpose(bondGeom[iID][:, 1:end-1]))
  end

  return bond_force
end
end