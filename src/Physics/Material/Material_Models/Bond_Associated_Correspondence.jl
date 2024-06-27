# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Bond_Associated_Correspondence
using LinearAlgebra
include("../../../Support/geometry.jl")
using .Geometry: strain, compute_left_stretch_tensor, compute_weighted_deformation_gradient
include("../../../Support/helpers.jl")
using .Helpers: qdim


using TimerOutputs
using NearestNeighbors: BallTree, inrange
function init_material_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, material_parameter::Dict)
  if !haskey(material_parameter, "Symmetry")
    @error "Symmetry for correspondence material is missing; options are 'isotropic plane strain', 'isotropic plane stress', 'anisotropic plane stress', 'anisotropic plane stress','isotropic' and 'anisotropic'. For 3D the plane stress or plane strain option is ignored."
    return nothing
  end
  if !haskey(material_parameter, "Accuracy Order")
    @warn "No Accuracy Order has been defined it is set to 1."
    merge!(material_parameter, Dict("Accuracy Order" => 1))
  end
  accuracy_order = material_parameter["Accuracy Order"]
  if typeof(accuracy_order) != Int
    @error "Accuracy Order must be an integer."
    return nothing
  elseif accuracy_order < 1
    @error "Accuracy Order must be an greater than zero."
    return nothing
  end
  dof = datamanager.get_dof()
  datamanager.create_bond_field("Bond Strain", Float64, "Matrix", dof)
  datamanager.create_bond_field("Bond Cauchy Stress", Float64, "Matrix", dof)
  datamanager.create_constant_bond_field("Bond Strain Increment", Float64, "Matrix", dof)
  weighted_volume = datamanager.create_constant_node_field("Weighted Volume", Float64, 1)
  gradient_weights = datamanager.create_constant_bond_field("Lagrangian Gradient Weights", Float64, dof)
  integral_stress = datamanager.create_constant_bond_field("Integral Nodal Stress", Float64, "Matrix", dof)

  deformation_gradient = datamanager.create_constant_bond_field("Bond Deformation Gradient", Float64, "Matrix", dof)
  deformation_gradient_dot = datamanager.create_constant_bond_field("Bond Rated Deformation Gradient", Float64, "Matrix", dof)

  return datamanager
end

function compute_weighted_volume(nodes::Union{SubArray,Vector{Int64}}, nlist::Union{Vector{Vector{Int64}},SubArray}, volume::SubArray, bond_damage::SubArray, omega::SubArray, weighted_volume::SubArray)
  for iID in nodes
    weighted_volume = sum(bond_damage[iID][:] .* omega[iID][:] .* volume[nlist[iID]])
  end
  return weighted_volume
end


function compute_stress_integral(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, nlist::Union{Vector{Vector{Int64}},SubArray}, omega::SubArray, bond_damage::SubArray, volume::SubArray, weighted_volume::SubArray, bond_geometry::SubArray, bond_length::SubArray, bond_stresses::SubArray, stress_integral::SubArray)
  temp::Matrix{Float64} = zeros(dof, dof)
  for iID in nodes
    stress_integral[iID, :, :] .= 0.0
    for (jID, nID) in enumerate(nlist[iID])
      for i in 1:dof
        temp[i, i] = 1 - bond_geometry[iID][jID, i] * bond_geometry[iID][jID, i]
        for j in i:dof
          temp[i, j] = bond_geometry[iID][jID, i] * bond_geometry[iID][jID, j]
          temp[j, i] = temp[i, j]
        end
      end
      temp ./= bond_length[iID][jID]
      stress_integral[iID, :, :] += (volume[nID] * omega[iID] * bond_damage[iID][jID] * (0.5 / weighted_volume[iID] + 0.5 / weighted_volume[nID])) .* bond_stresses[iID][jID, :, :] * temp
    end
  end
  return stress_integral
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
  return "Correspondence Bond-Associated"
end

"""
https://link.springer.com/article/10.1007/s10409-021-01055-5


- Bond associated neighborhood is the overlap between nlist[iID] and nlist[nlist[iID][jID]]
- Filter equal nodes and create a new neighborhoodlist for bond -> bond_nlist
- calculate K, Kinv and defGrad -> already there if the neighborhood loop is in a function
- weighted volume (sum(volume(bond_nlist))/sum(volume[nlist[iID]]))
- global local IDs to be checked 
  -> all neighbors search for neighbors at each core
  -> numbers are correct and it allows a change in size -> local ID is correct
"""



function find_local_neighbors(nID::Int64, coordinates::Union{SubArray,Matrix{Float64},Matrix{Int64}}, nlist::Union{Vector{Int64},SubArray{Int64}}, bond_horizon::Union{Float64,Int64})
  # excludes right now iID node in the coordinates list. Because it is an abritrary sublist it should be fine.
  # saving faster than recalculation?
  nlist_without_neighbor = view(nlist[nlist.!=nID], :)
  balltree = BallTree(transpose(coordinates[nlist_without_neighbor, :]))
  return nlist_without_neighbor[inrange(balltree, coordinates[nID, :], bond_horizon, true)]

end


function compute_forces(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, material_parameter::Dict, time::Float64, dt::Float64, to::TimerOutput)

  rotation::Bool, angles = datamanager.rotation_data()

  dof = datamanager.get_dof()
  nlist = datamanager.get_field("Neighborhoodlist")

  bond_damage = datamanager.get_bond_damage("NP1")

  omega = datamanager.get_field("Influence Function")
  volume = datamanager.get_field("Volume")

  bond_length = datamanager.get_field("Bond Length")

  undeformed_bond = datamanager.get_field("Bond Geometry")

  strain_N = datamanager.get_field("Bond Strain", "N")
  strain_NP1 = datamanager.get_field("Bond Strain", "NP1")

  integral_stress = datamanager.get_field("Integral Nodal Stress")

  stress_N = datamanager.get_field("Bond Cauchy Stress", "N")
  stress_NP1 = datamanager.get_field("Bond Cauchy Stress", "NP1")

  strain_increment = datamanager.get_field("Bond Strain Increment")
  bond_force = datamanager.get_field("Bond Forces")

  gradient_weight = datamanager.get_field("Lagrangian Gradient Weights")
  weighted_volume = datamanager.get_field("Weighted Volume")

  displacments = datamanager.get_field("Displacement", "NP1")
  velocity = datamanager.get_field("Velocity", "NP1")
  accuracy_order = material_parameter["Accuracy Order"]

  deformation_gradient = datamanager.get_field("Bond Deformation Gradient")
  deformation_gradient_dot = datamanager.get_field("Bond Rated Deformation Gradient")


  gradient_weights = compute_Lagrangian_gradient_weights(nodes, dof, accuracy_order, volume, nlist, horizon, bond_damage, omega, undeformed_bond, gradient_weights)
  weighted_volume = compute_weighted_volume(nodes, nlist, volume, bond_damage, omega, weighted_volume)

  deformation_gradient, deformation_gradient_dot = compute_weighted_deformation_gradient(nodes, dof, nlist, volume, gradient_weight, displacement, velocity, deformation_gradient, deformation_gradient_dot)


  strain_increment = update_Green_Langrange_bond_strain_increment(nodes, nlist, dt, deformation_gradient, deformation_gradient_dot, strain_increment)
  strain_NP1[:][:, :, :] = strain_N[:][:, :, :] .+ strain_increment[:][:, :, :]
  if rotation
    stress_N = rotate(nodes, dof, stress_N, angles, false)
    strain_increment = rotate(nodes, dof, strain_increment, angles, false)
  end

  material_models = split(material_parameter["Material Model"], "+")
  material_models = map(r -> strip(r), material_models)

  for material_model in material_models
    mod = datamanager.get_model_module(material_model)
    for iID in nodes
      for (jID, nID) in enumerate(nlist)


        strain_increment = compute_stra

        stress_NP1[iID][jID, :, :], datamanager = mod.compute_stresses(datamanager, nID, dof, material_parameter, time, dt, strain_increment[iID][:, :, :], stress_N[iID][:, :, :], stress_NP1[iID][:, :, :], (iID, jID))

      end
    end
  end
  if rotation
    stress_NP1 = rotate(nodes, dof, stress_NP1, angles, true)
  end
  bond_force = compute_bond_forces(nodes, nlist, bond_geometry, bond_length, bond_stress, integratal_nodal_stress, weighted_volume, gradient_weights, omega, bond_damage, bond_force)
  return datamanager

end

"""
function compute_bond_associated_weighted_volume(nodes::Union{SubArray, Vector{Int64}}, nlist::SubArray, coordinates::Union{SubArray,Vector{Float64}}, bond_damage::Union{SubArray,Vector{Float64}}, omega::Union{SubArray,Vector{Float64}}, volume::Union{SubArray,Vector{Float64}}, bond_horizon::Float64,
bond_horizon::Float64, weighted_volume::SubArray)

Compute the bond-associated weighted volume for given nodes. This function computes the bond-associated weighted volume for a given set of nodes. It iterates over each node in `nodes` and calculates the weighted volume associated with each bond connected to the node.

The function first computes the neighborhood volume of the node, which is the sum of the product of volume, bond damage, and omega for each neighboring node. Then, for each neighboring node, it calculates the local neighborhood list excluding the current neighbor, and computes the weighted volume for the bond between the current node and the neighbor. The weighted volume is calculated as the sum of the product of volume, bond damage, and omega for all bonds in the local neighborhood, divided by the neighborhood volume.

# Arguments
- `nodes::Union{SubArray, Vector{Int64}}`: A vector of integers representing node IDs.
- `nlist::SubArray`: Neighborhood list.
- `coordinates::Union{SubArray, Vector{Float64}}`: Node coordinates.
- `bond_damage::Union{SubArray, Vector{Float64}}`: Damage values for bonds.
- `omega::Union{SubArray, Vector{Float64}}`: Influence function values for bonds.
- `volume::Union{SubArray, Vector{Float64}}`: Volumes of nodes.
- `bond_horizon::Float64`: Local horizon around the neighbor
- `weighted_volume::SubArray`: Array to store the computed weighted volumes.

# Output
- `weighted_volume::SubArray`: Updated array containing the computed weighted volumes.

"""

function compute_bond_associated_weighted_volume(nodes::Union{SubArray,Vector{Int64}}, nlist::Union{Vector{Vector{Int64}},SubArray}, coordinates::Union{SubArray,Matrix{Float64}}, bond_damage::Union{SubArray,Vector{Vector{Float64}}}, omega::Union{SubArray,Vector{Vector{Float64}}}, volume::Union{SubArray,Vector{Float64}}, bond_horizon::Float64, weighted_volume::Union{SubArray,Vector{Vector{Float64}}})
  neighborhood_volume::Float64 = 0
  for iID in nodes

    neighborhood_volume = sum(volume[nlist[iID]] .* bond_damage[iID] .* omega[iID])

    for (jID, nID) in enumerate(nlist[iID])
      local_nlist = find_local_neighbors(nID, coordinates, nlist[iID], bond_horizon)
      sub_bond_list = [findfirst(x -> x == elem, nlist[iID]) for elem in local_nlist]
      weighted_volume[iID][jID] = sum(volume[local_nlist] .* bond_damage[iID][sub_bond_list] .* omega[iID][sub_bond_list] / neighborhood_volume)
    end
  end
  return weighted_volume
end
"""
accuracy_order::Int64 - needs a number of bonds which are linear independent

"""

function calculate_Q(accuracy_order::Int64, dof::Int64, undeformed_bond::Vector{Float64}, horizon::Union{Int64,Float64})

  Q = ones(Float64, qdim(accuracy_order, dof))  # Initialize Q with ones
  counter = 1
  p = zeros(Int64, dof)
  for this_order in 1:accuracy_order
    for p[1] in this_order:-1:0
      if dof == 3
        for p[2] in this_order-p[1]:-1:0
          p[3] = this_order - p[1] - p[2]
          # Calculate the product for Q[counter]
          Q[counter] = prod((undeformed_bond ./ horizon) .^ p)
          counter += 1
        end
      else
        p[2] = this_order - p[1]
        Q[counter] = prod((undeformed_bond ./ horizon) .^ p)
        counter += 1
      end
    end
  end
  return Q
end
function compute_Lagrangian_gradient_weights(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, accuracy_order::Int64, volume::Union{SubArray,Vector{Float64}}, nlist::Union{Vector{Vector{Int64}},SubArray}, horizon::Union{SubArray,Vector{Float64}}, bond_damage::Union{SubArray,Vector{Vector{Float64}}}, omega::Union{SubArray,Vector{Vector{Float64}}}, undeformed_bond, gradient_weights)
  #https://arxiv.org/pdf/2004.11477
  # maybe as static array
  dim = qdim(accuracy_order, dof)
  Minv = zeros(Float64, dim, dim)
  for iID in nodes
    M = zeros(Float64, dim, dim)
    for (jID, nID) in enumerate(nlist[iID])
      Q = calculate_Q(accuracy_order, dof, undeformed_bond[iID][jID, :], horizon[iID])
      M += omega[iID][jID] * bond_damage[iID][jID] * volume[nID] .* Q * Q'
    end
    try
      Minv = inv(M)
    catch
      @error "In compute_Lagrangian_gradient_weights the matrix M is singular and cannot be inverted. To many bond damages or a to small horizon might cause this."
      return nothing
    end
    for (jID, nID) in enumerate(nlist[iID])
      Q = calculate_Q(accuracy_order, dof, undeformed_bond[iID][jID, :], horizon[iID])
      # this comes from Eq(19) in 10.1007/s40571-019-00266-9
      # or example 1 in https://arxiv.org/pdf/2004.11477
      for idof in 1:dof
        gradient_weights[iID][jID, idof] = omega[iID][jID] / horizon[iID] .* (Minv[idof, :]' * Q)
      end
    end
  end
  return gradient_weights
end

function update_Green_Langrange_nodal_strain_increment(nodes::Union{SubArray,Vector{Int64}}, dt::Float64, deformation_gradient::SubArray, deformation_gradient_dot::SubArray, strain_increment::SubArray)

  for iID in nodes
    strain_increment[iID, :, :] = update_Green_Langrange_strain(dt, deformation_gradient[iID, :, :], deformation_gradient_dot[iID, :, :])
  end

end

function update_Green_Langrange_bond_strain_increment(nodes::Union{SubArray,Vector{Int64}}, nlist::Union{Vector{Vector{Int64}},SubArray}, dt::Float64, deformation_gradient::SubArray, deformation_gradient_dot::SubArray, strain_increment::SubArray)

  for iID in nodes
    for jID in nlist[iID]
      strain_increment[iID][jID, :, :] = update_Green_Langrange_strain(dt, deformation_gradient[iID][jID, :, :], deformation_gradient_dot[iID][jID, :, :])
    end
  end

end

function update_Green_Langrange_strain(dt::Float64, deformation_gradient::Matrix{Float64}, deformation_gradient_dot::Matrix{Float64})
  # later in Geometry.jl
  A = dt * 0.5 .* (deformation_gradient * deformation_gradient_dot)
  return A + A'
end

function compute_bond_forces(nodes, nlist, bond_geometry, bond_length, bond_stress, integratal_nodal_stress, weighted_volume, gradient_weights, omega, bond_damage, bond_forces)

  for iID in nodes
    for (jID, nID) in enumerate(nlist[iID])
      if bond_damage[iID][jID] == 0
        continue
      end
      bond_forces[iID][jID] = bond_damage[iID][jID] * omega[iID][jID] / (weighted_volume[iID] * bond_length[iID][jID]) .* bond_stress[iID][jID, :, :] * bond_geometry[iID][jID, :]
      bond_forces[iID][jID] += integratal_nodal_stress[iID, :, :] * gradient_weights[iID][jID, :]
    end
  end
  return bond_forces
end


end