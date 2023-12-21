# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Thermal_Flow
using LinearAlgebra
export compute_thermal_model
export thermal_model_name

"""
    thermal_model_name()

Gives the model name. It is needed for comparison with the yaml input deck.

# Arguments

# Returns
- `name::String`: "Thermal Flow"
"""
function thermal_model_name()
  return "Thermal Flow"
end

"""
    compute_thermal_model(datamanager, nodes, thermal_parameter, time, dt)

Calculates the thermal behavior of the material. This template has to be copied, the file renamed and edited by the user to create a new flow. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `flow parameter::Dict(String, Any)`: Dictionary with flow parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
# Returns
- `datamanager::Data_manager`: Datamanager.
Example:
```julia
```
"""
function compute_thermal_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, thermal_parameter::Dict, time::Float64, dt::Float64)

  if !haskey(thermal_parameter, "Type")
    @error "No model type has beed defined; ''Type'': ''Bond based'' or Type: ''Correspondence''"
    return datamanager
  end
  dof = datamanager.get_dof()
  nlist = datamanager.get_nlist()
  bond_damage = datamanager.get_field("Bond Damage", "NP1")
  bond_heat_flow = datamanager.get_field("Bond Heat Flow")
  undeformed_bond = datamanager.get_field("Bond Geometry")
  volume = datamanager.get_field("Volume")
  temperature = datamanager.get_field("Temperature", "NP1")
  lambda = thermal_parameter["Heat Transfer Coefficient"]

  if thermal_parameter["Type"] == "Bond based"
    horizon = datamanager.get_field("Horizon")
    if length(lambda) > 1
      lambda = lambda[1]
    end
    bond_heat_flow = compute_heat_flow_state_bond_based(nodes, dof, nlist, lambda, bond_damage, undeformed_bond, horizon, temperature, bond_heat_flow)
    return datamanager

  elseif thermal_parameter["Type"] == "Correspondence"

    lambda_matrix = zeros(Float64, dof, dof)
    Kinv = datamanager.get_field("Inverse Shape Tensor")
    if length(lambda) == 1
      for i in 1:dof
        lambda_matrix[i, i] = lambda
      end
    else
      for i in 1:dof
        lambda_matrix[i, i] = lambda[i]
      end
    end

    bond_heat_flow = compute_heat_flow_state_correspondence(nodes, dof, nlist, lambda_matrix, bond_damage, undeformed_bond, Kinv, temperature, volume, bond_heat_flow)
  else
    @error "No model valid type has beed defined; ''Bond based'' or ''Correspondence''"
  end
  return datamanager
end


"""
[BrighentiR2021](@cite)
is a prototype with some errors
"""
function compute_heat_flow_state_correspondence(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, nlist::SubArray, lambda::Matrix{Float64}, bond_damage::SubArray, undeformed_bond::SubArray, Kinv::SubArray, temperature::SubArray, volume::SubArray, bond_heat_flow::SubArray)

  nablaT = zeros(Float64, dof)
  H = zeros(Float64, dof)
  for iID in nodes
    H .= 0
    for (jID, neighborID) in enumerate(nlist[iID])
      tempState = (temperature[neighborID] - temperature[iID]) * volume[neighborID] * bond_damage[iID][jID]
      H += tempState .* undeformed_bond[iID][jID, 1:dof]
    end
    nablaT = Kinv[iID, :, :] * H
    """
if (MATRICES::vectorNorm(angles, 3)!=0){  
  MATRICES::tensorRotation(angles,lambda,true,rotatedLambda);
  for (int i=0 ; i<3 ; ++i) {
    q[i] = 0.0;
    for (int j=0 ; j<3 ; ++j) {
      q[i] += rotatedLambda[3*i + j] * nablaT[j];
    }
  }
}"""
    q = lambda * nablaT
    for (jID, neighborID) in enumerate(nlist[iID])
      temp = Kinv[iID, :, :] * undeformed_bond[iID][jID, 1:dof]
      bond_heat_flow[iID][jID] = dot(temp, q)
    end
  end
  return bond_heat_flow

end

"""
    compute_heat_flow_state_bond_based(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, nlist::SubArray,
      lambda::Union{Float64, Int64}, bond_damage::SubArray, undeformed_bond::SubArray, horizon::SubArray,
      temperature::SubArray, volume::SubArray, bond_heat_flow::SubArray)

Calculate heat flow based on a bond-based model for thermal analysis.

# Arguments
- `nodes::Union{SubArray,Vector{Int64}}`: An array of node indices for which heat flow should be computed.
- `dof::Int64`: The degree of freedom, either 2 or 3, indicating whether the analysis is 2D or 3D.
- `nlist::SubArray`: A SubArray representing the neighbor list for each node.
- `lambda::Union{Float64, Int64}`: The thermal conductivity.
- `bond_damage::SubArray`: A SubArray representing the damage state of bonds between nodes.
- `undeformed_bond::SubArray`: A SubArray representing the geometry of the bonds.
- `horizon::SubArray`: A SubArray representing the horizon for each node.
- `temperature::SubArray`: A SubArray representing the temperature at each node.
- `volume::SubArray`: A SubArray representing the volume at each node.
- `bond_heat_flow::SubArray`: A SubArray where the computed bond heat flow values will be stored.

## Returns
- `bond_heat_flow`: updated bond heat flow values will be stored.

## Description
This function calculates the heat flow between neighboring nodes based on a bond-based model for thermal analysis [OterkusS2014b](@cite). It considers various parameters, including thermal conductivity, damage state of bonds, geometry of bonds, horizons, temperature, and volume. The calculated bond heat flow values are stored in the `bond_heat_flow` array.

"""
function compute_heat_flow_state_bond_based(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, nlist::SubArray, lambda::Union{Float64,Int64}, bond_damage::SubArray, undeformed_bond::SubArray, horizon::SubArray, temperature::SubArray, bond_heat_flow::SubArray)
  kernel::Float64 = 0.0
  for iID in nodes
    #nlist
    if dof == 2
      kernel = 6.0 / (pi * horizon[iID] * horizon[iID] * horizon[iID])
    else
      kernel = 6.0 / (pi * horizon[iID] * horizon[iID] * horizon[iID] * horizon[iID])
    end

    for (jID, neighborID) in enumerate(nlist[iID])
      if bond_damage[iID][jID] == 0
        continue
      end
      tempState = bond_damage[iID][jID] * (temperature[neighborID] - temperature[iID])
      bond_heat_flow[iID][jID] = lambda * kernel * tempState / undeformed_bond[iID][jID, end]
    end
  end
  return bond_heat_flow
end
end