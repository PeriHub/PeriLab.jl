# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Heat_transfer
export compute_thermal_model
export thermal_model_name
"""
   thermal_flow_name()

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
   compute_force(datamanager, nodes, thermal_parameter, time, dt)

   Calculates the heat transfer to the environment. [BrighentiR2021](@cite)

   # Arguments
   - `datamanager::Data_manager`: Datamanager.
   - `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
   - `flow parameter::Dict(String, Any)`: Dictionary with flow parameter.
   - `time::Float64`: The current time.
   - `dt::Float64`: The current time step.
   # Returns
        - - `datamanager::Data_manager`: Datamanager.
   Example:
   ```julia
     ```
   """
function compute_thermal_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, thermal_parameter::Dict, time::Float64, dt::Float64)
  dof = datamanager.get_dof()
  volume = datamanager.get_field("Volume")
  alpha = thermal_parameter["Alpha"]
  Tenv = thermal_parameter["Environmental Temperature"]
  heat_flow = datamanager.get_field("Heat Flow", "NP1")
  temperature = datamanager.get_field("Temperature", "NP1")
  surface_nodes = datamanager.get_field("Surface_Nodes")
  dx = 1.0
  area = 1.0


  for iID in nodes

    if dof == 2
      dx = sqrt(volume[iID])
    elseif dof == 3
      dx = volume[iID]^(1 / 3)
    end

    if surface_nodes[iID]
      heat_flow[iID] += (alpha * (temperature[iID] - Tenv)) / dx * area
    end
  end

  return datamanager
end

"""
   get_surface_nodes(iID::Int64, nlist::SubArray, coordinates::Union{SubArray,Vector{Float64}}, volume::SubArray, surface_nodes::Union{SubArray,Vector{Bool}})

   Returns the surface nodes.

   # Arguments
   - `iID::Int64`: The index of the node.
   - `nlist::SubArray`: The neighbor list.
   - `coordinates::Union{SubArray,Vector{Float64}}`: The coordinates of the nodes.
   - `volume::SubArray`: The volume of the nodes.
   - `surface_nodes::Union{SubArray,Vector{Bool}}`: The surface nodes.
   # Returns
   - `surface_nodes::Union{SubArray,Vector{Bool}}`: The surface nodes.
"""
function get_surface_nodes(iID::Int64, nlist::SubArray, coordinates::Union{SubArray,Vector{Float64}}, volume::SubArray, surface_nodes::Union{SubArray,Vector{Bool}})
  compare_neighbor = 0
  neighbor_volume = 0.0
  right, left, front, back, above, below = false, false, false, false, false, false, false, false

  for jID in 1:nneighbors[iID]
    if bond_damage[iID][jID] == 0.0
      continue
    end

    if equal_discretized
      if !right && coordinates[iID, 1] > coordinates[nlist[iID][jID], 1] && coordinates[iID, 2] == coordinates[nlist[iID][jID], 2] && coordinates[iID, 3] == coordinates[nlist[iID][jID], 3]
        right = true
        compare_neighbor += 1
      end
      if !left && coordinates[iID, 1] < coordinates[nlist[iID][jID], 1] && coordinates[iID, 2] == coordinates[nlist[iID][jID], 2] && coordinates[iID, 3] == coordinates[nlist[iID][jID], 3]
        left = true
        compare_neighbor += 1
      end
      if !front && coordinates[iID, 2] > coordinates[nlist[iID][jID], 2] && coordinates[iID, 1] == coordinates[nlist[iID][jID], 1] && coordinates[iID, 3] == coordinates[nlist[iID][jID], 3]
        front = true
        compare_neighbor += 1
      end
      if !back && coordinates[iID, 2] < coordinates[nlist[iID][jID], 2] && coordinates[iID, 1] == coordinates[nlist[iID][jID], 1] && coordinates[iID, 3] == coordinates[nlist[iID][jID], 3]
        back = true
        compare_neighbor += 1
      end
      if !above && coordinates[iID, 3] > coordinates[nlist[iID][jID], 3] && coordinates[iID, 1] == coordinates[nlist[iID][jID], 1] && coordinates[iID, 2] == coordinates[nlist[iID][jID], 2]
        above = true
        compare_neighbor += 1
      end
      if !below && coordinates[iID, 3] < coordinates[nlist[iID][jID], 3] && coordinates[iID, 1] == coordinates[nlist[iID][jID], 1] && coordinates[iID, 2] == coordinates[nlist[iID][jID], 2]
        below = true
        compare_neighbor += 1
      end

    else
      neighbor_volume += volume[nlist[iID][jID]]
    end

  end

  if equal_discretized

    if dof == 2 && compare_neighbor != 4
      area = 4 - compare_neighbor
    elseif dof == 3 && compare_neighbor != 6
      area = 6 - compare_neighbor
    end

    specific_volume = compare_neighbor

  else

    specific_volume = neighbor_volume / dx

    if dof == 2
      area = specific_volume / (1 / 4)
    elseif dof == 3
      area = specific_volume / (1 / 6)
    end

  end
  return area
end
end