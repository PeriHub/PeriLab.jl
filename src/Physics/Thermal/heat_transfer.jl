# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Heat_transfer
export compute_thermal_model
export thermal_model_name

"""
    thermal_model_name()

Gives the flow name. It is needed for comparison with the yaml input deck.

# Arguments

# Returns
- `name::String`: The name of the thermal flow model.

Example:
```julia
println(flow_name())
"Thermal Template"
```
"""
function thermal_model_name()
  return "Heat Transfer"
end

"""
    compute_thermal_model(datamanager, nodes, thermal_parameter, time, dt)

Calculates the heat transfer to the environment. [BrighentiR2021](@cite)

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
  dof = datamanager.get_dof()
  volume = datamanager.get_field("Volume")
  kappa = thermal_parameter["Thermal Conductivity"]
  Tenv = thermal_parameter["Environmental Temperature"]
  heat_flow = datamanager.get_field("Heat Flow", "NP1")
  temperature = datamanager.get_field("Temperature", "NP1")
  surface_nodes = datamanager.get_field("Surface_Nodes")
  specific_volume = datamanager.get_field("Specific Volume", "NP1")
  bond_damage = datamanager.get_field("Bond Damage", "NP1")
  horizon = datamanager.get_field("Horizon")
  nlist = datamanager.get_nlist()
  dx = 1.0
  area = 1.0

  specific_volume = calculate_specific_volume(nodes, nlist, volume, bond_damage, specific_volume, dof, horizon)

  for iID in nodes

    if dof == 2
      dx = sqrt(volume[iID])
    elseif dof == 3
      dx = volume[iID]^(1 / 3)
    end

    if surface_nodes[iID] && specific_volume[iID] > 1.0 # could be more as sometimes specific_volume is sometimes weird
      heat_flow[iID] += (kappa * (temperature[iID] - Tenv)) / dx * specific_volume[iID]
    end
  end

  return datamanager
end

"""
  calculate_specific_volume(iID::Int64, nlist::SubArray, coordinates::Union{SubArray,Vector{Float64}}, volume::SubArray, surface_nodes::Union{SubArray,Vector{Bool}})

Calculates the specific volume.

# Arguments
- `iID::Int64`: The index of the node.
- `nlist::SubArray`: The neighbor list.
- `coordinates::Union{SubArray,Vector{Float64}}`: The coordinates of the nodes.
- `volume::SubArray`: The volume of the nodes.
- `surface_nodes::Union{SubArray,Vector{Bool}}`: The surface nodes.
# Returns
- `specific_volume::Union{SubArray,Vector{Bool}}`: The surface nodes.
"""
function calculate_specific_volume(nodes, nlist, volume, bond_damage, specific_volume, dof, horizon)

  for iID in nodes
    neighbor_volume = 0.0
    for (jID, neighborID) in enumerate(nlist[iID])
      if bond_damage[iID][jID] == 0.0
        continue
      end
      neighbor_volume += volume[nlist[iID][jID]]
    end
    if dof == 2
      horizon_volume = pi * horizon[iID]^2
    elseif dof == 3
      horizon_volume = 4 / 3 * pi * horizon[iID]^3
    end
    specific_volume[iID] = horizon_volume / neighbor_volume
  end
  return specific_volume
end
end