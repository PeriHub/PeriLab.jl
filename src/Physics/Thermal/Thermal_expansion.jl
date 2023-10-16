# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause


# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Thermal_expansion
export compute_thermal_model
export thermal_model_name
"""
   thermal_model_name()

   Gives the expansion model name. It is needed for comparison with the yaml input deck.

   Parameters:

   Returns:
   - `name::String`: The name of the thermal expansion model.

   Example:
   ```julia
   println(flow_name())
   "Thermal Expansion"
   ```
   """
function thermal_model_name()
    return "Thermal Expansion"
end
"""
   compute_thermal_model(datamanager, nodes, flow_parameter, time, dt)

   Calculates the thermal expansion of the material. 

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
function compute_thermal_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, flow_parameter::Dict, time::Float64, dt::Float64)
    @info "Please write a thermal flow model name in thermal_flow_name()."
    @info "You can call your routine within the yaml file."
    @info "Fill the compute_force(datamanager, nodes, flow_parameter, time, dt) function."
    @info "The datamanger and flow_parameter holds all you need to solve your problem on thermal flow level."
    @info "add own files and refer to them. If a module does not exist. Add it to the project or contact the developer."
    return datamanager
end

function thermal_stretch(nodes::Union{Subarray,Vector{Int64}}, alpha::Float64, temperature::SubArray, deformed_bond::SubArray, undeformed_bond::SubArray, thermal_stretch::SubArray)

    for iID in nnodes
        thermal_stretch[iID][:] -= alpha .* temperature[iID] .* (deformed_bond[iID][:, end] - undeformed_bond[iID][:, end])
    end
end
return thermal_stretch
end

function thermal_strain(nodes::Union{Subarray,Vector{Int64}}, alpha::Matrix{Float64}, temperature::SubArray, thermal_strain::SubArray)
    for iID in nodes
        thermal_strain[iID] = alpha .* temperature[iID]
    end
    return thermal_strain
end
