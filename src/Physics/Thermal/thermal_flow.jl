# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Thermal_Flow
export compute_thermal_model
export thermal_model_name
"""
   thermal_model_name()

   Gives the model name. It is needed for comparison with the yaml input deck.

   Parameters:

   Returns:
   - `name::String`: "Thermal Flow"

    """
function thermal_model_name()
  return "Thermal Flow"
end
"""
   compute_thermal_model(datamanager, nodes, thermal_parameter, time, dt)

   Calculates the thermal behavior of the material. This template has to be copied, the file renamed and edited by the user to create a new flow. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

   Parameters:
        - `datamanager::Data_manager`: Datamanager.
        - `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
        - `flow parameter::Dict(String, Any)`: Dictionary with flow parameter.
        - `time::Float64`: The current time.
        - `dt::Float64`: The current time step.
   Returns:
        - - `datamanager::Data_manager`: Datamanager.
   Example:
   ```julia
     ```
   """
function compute_thermal_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, thermal_parameter::Dict, time::Float64, dt::Float64)

  if !haskey(thermal_parameter, "Type")
    @error "No model type has beed defined; Type: ''Bond based'' or Type: ''Correspondence''"
    return datamanager
  end
  dof = datamanager.get_dof()
  nlist = datamanager.get_nlist()
  bond_damage = datamanager.get_field("Bond Damage", "NP1")
  bond_geometry = datamanager.get_field("Bond Geometry")
  horizon = datamanager.get_field("Horizon")
  volume = datamanager.get_field("Volume")

  lamda = thermal_parameter["Lambda"]
  if thermal_parameter["Type"] == "Bond based"
    if length(lamda) > 1
      lambda = lambda[1]
    end
    return compute_heat_flow_state_bondbased(nodes, dof, nlist, lambda, bond_damage, bond_geometry, horizon, temperature, volume, bond_heat_flow)

  elseif thermal_parameter["Type"] == "Correspondence"
    lambda_matrix = zeros(dof, dof)
    if length(lamda) == 1
      for i in 1:dof
        lambda_matrix[i, i] = lambda
      end
    else
      for i in 1:dof
        lambda_matrix[i, i] = lambda[i]
      end
    end

    return compute_heat_flow_state_correspondence()
  else
    @error "No model valid type has beed defined; ''Bond based'' or ''Correspondence''"
  end
  return datamanager
end




"""
[BrighentiR2021](@cite)
"""

function compute_heat_flow_state_correspondence()


  for iID in nodes
    H = zeros(Float64, dof)

    tempState = (temperature[nlist[iID]] .- temperature[iID]) * volume[nlist[iID]] * bond_damage[iID]


    H = sum(tempState .* undeformed_bond[iID][:, 1:dof])
    nablaT = Kinv[iID, :, :] * H
    """    for(iNID=0 ; iNID<numNeighbors ; ++iNID, bondDamage++){
      for (int i=0 ; i<3 ; ++i) 
        X[i] = modelCoord[3*neighborID+i] - Xp[i];
      tempState = (temperature[neighborID] - temperature[iID]) * volume[neighborID] * (1 - *bondDamage);
      // sum_j (Tj-Ti)*rij*Vj -> EQ. (8)
      for (int i=0 ; i<3 ; ++i) H[i] += tempState * X[i] ;
      // std::cout<<*bondDamage<<std::endl;
    }
    // Ki * H -> EQ. (8)
    for (int i=0 ; i<3 ; ++i) {
      nablaT[i] = 0.0;
      for (int j=0 ; j<3 ; ++j) {
        nablaT[i] += KInv[3*i + j] * H[j];
        //std::cout<< nablaT[i]<<std::endl;
      }
    }
    if (MATRICES::vectorNorm(angles, 3)!=0){  
      MATRICES::tensorRotation(angles,lambda,true,rotatedLambda);
      for (int i=0 ; i<3 ; ++i) {
        q[i] = 0.0;
        for (int j=0 ; j<3 ; ++j) {
          q[i] += rotatedLambda[3*i + j] * nablaT[j];
        }
      }
    }
    else{
      for (int i=0 ; i<3 ; ++i) {
        q[i] = 0.0;
        for (int j=0 ; j<3 ; ++j) {
          q[i] += lambda[3*i + j] * nablaT[j];
        }
      }
    }
    for(iNID=0 ; iNID<numNeighbors ; ++iNID){
      neighborID = neighborhoodList[secondNeighborhoodListIndex++];
      for (int i=0 ; i<3 ; ++i)X[i] = modelCoord[3*neighborID+i] - Xp[i];
      for (int i=0 ; i<3 ; ++i){
        temp[i] = 0.0;
        for (int j=0 ; j<3 ; ++j) {
          temp[i] += KInv[3*i + j] * X[j]; // K * rij -> Eq. (7)
        }
      }
      // bond damage muss hier noch rein
       // qj * temp -> Eq (6) both must be negative, because rij is equal
      //for (int i=0 ; i<3 ; ++i)  heatFlowState[iID] += temp[i] * q[i] * volume[neighborID];
      //for (int i=0 ; i<3 ; ++i)  heatFlowState[neighborID] += temp[i] * q[i] * volume[iID];
      for (int i=0 ; i<3 ; ++i)  heatFlowState[iID] -= temp[i] * q[i] * volume[neighborID];
      for (int i=0 ; i<3 ; ++i)  heatFlowState[neighborID] += temp[i] * q[i] * volume[iID];
    }
    
  }
  
}"""
  end
end
end

function compute_heat_flow_state_bondbased(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, nlist::SubArray, lambda::Float64, bond_damage::SubArray, bond_geometry::SubArray, horizon::SubArray, temperature::SubArray, volume::SubArray, bond_heat_flow::SubArray)
  kernel::Float64 = 0.0
  for iID in nodes
    #nlist
    if dof == 2
      kernel = 6.0 / (pi * horizon[iID] * horizon[iID] * horizon[iID] * bond_geometry[iID][:, end])
    else
      kernel = 6.0 / (pi * horizon[iID] * horizon[iID] * horizon[iID] * horizon[iID])
    end

    for (jID, neighborID) in enumerate(nlist[iID])
      if bond_damage[iID][jID] == 0
        continue
      end
      tempState = bond_damage[iID][jID] * (temperature[neighborID] - temperature[iID])
      bond_heat_flow[iID][jID] = lambda * kernel * tempState * volume[neighborID] / bond_geometry[iID][jID, end]
    end
  end
end