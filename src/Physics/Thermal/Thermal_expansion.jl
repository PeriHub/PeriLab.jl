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
   compute_thermal_model(datamanager, nodes, thermal_parameter, time, dt)

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
function compute_thermal_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, thermal_parameter::Dict, time::Float64, dt::Float64)

    temperature = datamanager.get_field("Temperature", "NP1")
    alpha = thermal_parameter["Heat expansion"]
    if haskey(thermal_parameter, "Thermal Stretch") && thermal_parameter["Thermal Stretch"]
        undeformed_bond = datamanager.get_field("Bond Geoemtry")
        deformed_bond = datamanager.get_field("Deformed Bond Geometry", "NP1")
        thermal_bond_stretch = datamanager.create_constant_bond_field(Float64, "Thermal Stretch", 1)
        if length(alpha) > 1
            @warn "Matrix is defined and first entry is used in ''Thermal Stretch''."
            alpha = alpha[1]
        end
        thermal_bond_stretch = thermal_stretch(nodes, alpha, temperature, deformed_bond, undeformed_bond, thermal_bond_stretch)
    elseif haskey(thermal_parameter, "Thermal Strain") && thermal_parameter["Thermal Strain"]
        thermal_strain_tensor = datamanager.get_field("Thermal Strain")
        dof = datamanager.get_dof()
        alpha_mat = zeros(Float64, dof, dof)
        if length(alpha) == 1
            for i in 1:dof
                alpha_mat[i, i] = alpha
            end
        elseif length(alpha) == dof || length(alpha) == 3
            for i in 1:dof
                alpha_mat[i, i] = alpha[i]
            end
        elseif length(alpha) == dof * dof || length(alpha) == 9
            @error "not yet implemented for full heat expansion matrix"
        end
        thermal_strain_tensor = thermal_strain(nodes, alpha, temperature, thermal_strain_tensor)
    else
        @error "No valid thermal expansion measure defined ''Thermal Stretch'' or ''Thermal Strain''"
    end

    return datamanager
end

function thermal_stretch(nodes::Union{SubArray,Vector{Int64}}, alpha::Float64, temperature::SubArray, deformed_bond::SubArray, undeformed_bond::SubArray, thermal_stretch::SubArray)

    for iID in nnodes
        thermal_stretch[iID][:] -= alpha .* temperature[iID] .* (deformed_bond[iID][:, end] - undeformed_bond[iID][:, end])
    end
end
return thermal_stretch
end

function thermal_strain(nodes::Union{SubArray,Vector{Int64}}, alpha::Matrix{Float64}, temperature::SubArray, thermal_strain::SubArray)
    for iID in nodes
        thermal_strain[iID] = alpha .* temperature[iID]
    end
    return thermal_strain
end
