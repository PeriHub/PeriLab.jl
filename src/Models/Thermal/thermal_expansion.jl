# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Thermal_Expansion
using LinearAlgebra
using StaticArrays
using ...Pre_Calculation.Deformation_Gradient
export fields_for_local_synchronization
export compute_model
export thermal_model_name
export init_model
"""
	thermal_model_name()

Gives the expansion model name. It is needed for comparison with the yaml input deck.

# Arguments

# Returns
- `name::String`: The name of the thermal expansion model.

Example:
```julia
"Thermal Expansion"
```
"""
function thermal_model_name()
    return "Thermal Expansion"
end

function thermal_expansion_matrix(alpha::T,
                                  ::Val{2}) where {T<:Union{Float64,Matrix{Float64}}}
    if length(alpha) == 1
        return alpha_mat = SMatrix{2,2,Float64}(alpha, 0, 0, alpha)
    elseif length(alpha) == 2
        return alpha_mat = SMatrix{2,2,Float64}(alpha[1], 0, 0, alpha[2])
    elseif length(alpha) == 4
        @error "Full heat expansion matrix is not implemented yet."
    end
end
function thermal_expansion_matrix(alpha::T,
                                  ::Val{3}) where {T<:Union{Float64,Matrix{Float64}}}
    if length(alpha) == 1
        return alpha_mat = SMatrix{3,3,Float64}(alpha, 0, 0, 0, alpha, 0, 0, 0, alpha)
    elseif length(alpha) == 3
        return alpha_mat = SMatrix{3,3,Float64}(alpha[1], 0, 0, 0, alpha[2], 0, 0, 0,
                                                alpha[3])
    elseif length(alpha) == 9
        @error "Full heat expansion matrix is not implemented yet."
    end
end

"""
	init_model(datamanager, nodes, thermal_parameter)

Inits the thermal model. This template has to be copied, the file renamed and edited by the user to create a new thermal. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

# Arguments
- `datamanager::Data_Manager`: Datamanager.
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `thermal parameter::Dict(String, Any)`: Dictionary with thermal parameter.
- `block::Int64`: The current block.
# Returns
- `datamanager::Data_Manager`: Datamanager.

"""
function init_model(datamanager::Module,
                    nodes::AbstractVector{Int64},
                    thermal_parameter::Dict)
    if !haskey(thermal_parameter, "Reference Temperature")
        @warn "No reference temperature defined. Assuming 0"
    end
    return datamanager
end

"""
	compute_model(datamanager, nodes, thermal_parameter, block::Int64, time, dt)

Calculates the thermal expansion of the material.

# Arguments
- `datamanager::Data_Manager`: Datamanager.
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `flow parameter::Dict(String, Any)`: Dictionary with flow parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
# Returns
- `datamanager::Data_Manager`: Datamanager.
Example:
```julia
```
"""
function compute_model(datamanager::Module,
                       nodes::AbstractVector{Int64},
                       thermal_parameter::Dict,
                       block::Int64,
                       time::Float64,
                       dt::Float64)
    thermal_parameter["Matrix based"]
    if get(thermal_parameter, "Matrix based", false)
        return compute_model_matrix_based(datamanager, nodes, thermal_parameter, block,
                                          time, dt)
    end

    temperature_NP1 = datamanager.get_field("Temperature", "NP1")
    dof = datamanager.get_dof()

    alpha_mat=thermal_expansion_matrix(thermal_parameter["Thermal Expansion Coefficient"],
                                       Val(dof))

    ref_temp = 0.0
    ref_temp = get(thermal_parameter, "Reference Temperature", 0.0)

    undeformed_bond = datamanager.get_field("Bond Geometry")
    undeformed_bond_length = datamanager.get_field("Bond Length")
    deformed_bond = datamanager.get_field("Deformed Bond Geometry", "NP1")
    deformed_bond_length = datamanager.get_field("Deformed Bond Length", "NP1")

    temp_diff = 0.0

    for iID in nodes
        temp_diff = temperature_NP1[iID] - ref_temp
        @inbounds @fastmath @views for jID in eachindex(undeformed_bond[iID])
            for j in 1:dof
                deformed_bond[iID][jID][j] -= temp_diff * alpha_mat[j, j] *
                                              undeformed_bond[iID][jID][j]
            end
        end
        deformed_bond_length[iID] .-= sum(alpha_mat) / dof * temp_diff .*
                                      undeformed_bond_length[iID]
    end

    if datamanager.has_key("Deformation Gradient")
        datamanager = Deformation_Gradient.compute(datamanager, nodes, Dict(), block)
    end

    return datamanager
end

"""
	fields_for_local_synchronization(datamanager::Module, model::String)

Returns a user developer defined local synchronization. This happens before each model.



# Arguments

"""
function fields_for_local_synchronization(datamanager::Module, model::String)
    #download_from_cores = false
    #upload_to_cores = true
    #datamanager.set_local_synch(model, "Bond Forces", download_from_cores, upload_to_cores)
    return datamanager
end

function compute_model_matrix_based(datamanager::Module,
                                    nodes::AbstractVector{Int64},
                                    thermal_parameter::Dict,
                                    block::Int64,
                                    time::Float64,
                                    dt::Float64)
    thermal_factor = datamanager.get_field("Thermal Expension Factor")
    temperatureN = datamanager.get_field("Temperature", "N")
    temperatureNP1 = datamanager.get_field("Temperature", "NP1")
    dof = datamanager.get_dof()

    ref_temp = get(thermal_parameter, "Reference Temperature", 0.0)

    alpha_mat=thermal_expansion_matrix(thermal_parameter["Thermal Expansion Coefficient"],
                                       Val(dof))
    for iID in nodes
        @views thermal_factor[iID, :]=(temperatureNP1[iID]-ref_temp) .* diag(alpha_mat) .+ 1
    end
    return datamanager
end

end # Module end
