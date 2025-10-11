# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Thermal_expansion
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
println(flow_name())
"Thermal Expansion"
```
"""
function thermal_model_name()
    return "Thermal Expansion"
end

"""
    init_model(datamanager, nodes, thermal_parameter)

Inits the thermal model. This template has to be copied, the file renamed and edited by the user to create a new thermal. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `thermal parameter::Dict(String, Any)`: Dictionary with thermal parameter.
- `block::Int64`: The current block.
# Returns
- `datamanager::Data_manager`: Datamanager.

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
- `datamanager::Data_manager`: Datamanager.
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `flow parameter::Dict(String, Any)`: Dictionary with flow parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
# Returns
- `datamanager::Data_manager`: Datamanager.
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
    temperature_NP1 = datamanager.get_field("Temperature", "NP1")
    dof = datamanager.get_dof()

    alpha = thermal_parameter["Thermal Expansion Coefficient"]
    ref_temp = 0.0
    if haskey(thermal_parameter, "Reference Temperature")
        ref_temp = thermal_parameter["Reference Temperature"]
    end

    alpha_mat::Matrix{Float64} = @MMatrix zeros(Float64, dof, dof)
    if length(alpha) == 1
        alpha_mat = alpha .* I(dof)
    elseif length(alpha) == dof || length(alpha) == 3
        for i in 1:dof
            alpha_mat[i, i] = alpha[i]
        end
    elseif length(alpha) == dof * dof || length(alpha) == 9
        @error "Full heat expansion matrix is not implemented yet."
    end
    undeformed_bond = datamanager.get_field("Bond Geometry")
    undeformed_bond_length = datamanager.get_field("Bond Length")
    deformed_bond = datamanager.get_field("Deformed Bond Geometry", "NP1")
    deformed_bond_length = datamanager.get_field("Deformed Bond Length", "NP1")

    temp_diff = 0.0

    for iID in nodes
        temp_diff = temperature_NP1[iID] - ref_temp
        for j in 1:dof
            deformed_bond[iID][:][j] .-= temp_diff * alpha_mat[j, j] .*
                                         undeformed_bond[iID][:][j]
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

end # Module end
