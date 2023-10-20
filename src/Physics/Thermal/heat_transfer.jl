# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Heat_transfer
export compute_thermal_model
export thermal_model_name
"""
   thermal_flow_name()

   Gives the flow name. It is needed for comparison with the yaml input deck.

   Parameters:

   Returns:
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
  dof = datamanager.get_dof()
  volume = datamanager.get_field("Volume")

  numNeighbors = datamanager.get_field("Number of Neighbors")
  alpha = thermal_parameter["Alpha"]
  Tenv = thermal_parameter["Environment Temperature"]
  equal_discretized = thermal_parameter["Equal Discretized"]
  heat_flow = datamanager.get_field("Heat Flow", "NP1")
  temperature = datamanager.get_field("Temperature", "NP1")
  specific_volume = datamanager.get_field("Specific Volume")
  bond_damage = datamanager.get_field("Bond Damage")
  coordinates = datamanager.get_field("Coordinates")

  dx = 1.0
  area = 1.0

  if dof == 2
    dx = sqrt(volume[iID])
  elseif dof == 3
    dx = volume[iID]^(1 / 3)
  end

  for iID in nodes
    compareNeighbor = 0
    neighbor_volume = 0.0
    right, left, front, back, above, below = false, false, false, false, false, false, false, false
    for jID in numNeighbors[iID]
      if bond_damage[iID, jID] == 0.0
        continue
      end

      if equal_discretized
        if !right && coordinates[iID, 1] > coordinates[jID, 1] && coordinates[iID, 2] == coordinates[jID, 2] && coordinates[iID, 3] == coordinates[jID, 3]
          right = true
          compare_neighbor += 1
        end
        if !left && coordinates[iID, 1] < coordinates[jID, 1] && coordinates[iID, 2] == coordinates[jID, 2] && coordinates[iID, 3] == coordinates[jID, 3]
          left = true
          compare_neighbor += 1
        end
        if !front && coordinates[iID, 2] > coordinates[jID, 2] && coordinates[iID, 1] == coordinates[jID, 1] && coordinates[iID, 3] == coordinates[jID, 3]
          front = true
          compare_neighbor += 1
        end
        if !back && coordinates[iID, 2] < coordinates[jID, 2] && coordinates[iID, 1] == coordinates[jID, 1] && coordinates[iID, 3] == coordinates[jID, 3]
          back = true
          compare_neighbor += 1
        end
        if !above && coordinates[iID, 3] > coordinates[jID, 3] && coordinates[iID, 1] == coordinates[jID, 1] && coordinates[iID, 2] == coordinates[jID, 2]
          above = true
          compare_neighbor += 1
        end
        if !below && coordinates[iID, 3] < coordinates[jID, 3] && coordinates[iID, 1] == coordinates[jID, 1] && coordinates[iID, 2] == coordinates[jID, 2]
          below = true
          compare_neighbor += 1
        end

      else
        neighbor_volume += volume[jID]
      end

    end

    if equal_discretized

      if dof == 2 && compare_neighbor != 4
        area = 4 - compareNeighbor
      elseif dof == 3 && compare_neighbor != 6
        area = 6 - compareNeighbor
      end

      specific_volume[iID] = compareNeighbor

    else

      specific_volume[iID] = neighbor_volume / dx

      if dof == 2
        area = specific_volume[iID] / (1 / 4)
      elseif dof == 3
        area = specific_volume[iID] / (1 / 6)
      end

    end

    heat_flow[iID] += (alpha * (temperature[iID] - Tenv)) / dx * area

  end

  return datamanager
end

end