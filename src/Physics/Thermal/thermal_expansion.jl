# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause


# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Thermal_expansion
using LinearAlgebra
using StaticArrays
include("../Pre_calculation/deformation_gradient.jl")
include("../Pre_calculation/bond_deformation_gradient.jl")
using .Deformation_Gradient
using .Bond_Deformation_Gradient

export compute_thermal_model
export thermal_model_name
export init_thermal_model
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
    init_thermal_model(datamanager, nodes, thermal_parameter)

Inits the thermal model. This template has to be copied, the file renamed and edited by the user to create a new thermal. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `thermal parameter::Dict(String, Any)`: Dictionary with thermal parameter.
- `block::Int64`: The current block.
# Returns
- `datamanager::Data_manager`: Datamanager.

"""
function init_thermal_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, thermal_parameter::Dict)

    return datamanager
end

"""
    compute_thermal_model(datamanager, nodes, thermal_parameter, time, dt)

Calculates the thermal expansion of the material. 

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
    temperature_NP1 = datamanager.get_field("Temperature", "NP1")
    dof = datamanager.get_dof()

    alpha = thermal_parameter["Thermal Expansion Coefficient"]

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

    for iID in nodes
        for j in 1:dof
            deformed_bond[iID][:, j] .-= temperature_NP1[iID] * alpha_mat[j, j] .* undeformed_bond[iID][:, j]
        end
        deformed_bond_length[iID] .-= sum(alpha_mat) / dof * temperature_NP1[iID] .* undeformed_bond_length[iID]
    end

    if datamanager.has_key("Deformation Gradient")
        block_id = datamanger.get_field("Block_Id")
        # is a work around to get the block ID information. Because all points in this routine are in the same block, the first node give the information
        datamanager = Deformation_Gradient.compute(datamanager, nodes, block_id[nodes[1]])
    end

    return datamanager
end

"""
    thermal_deformation(nodes, alpha, temperature, undeformed_bond, thermal_deformation)

Calculate thermal deformation for a set of nodes.

This function calculates thermal deformation for a specified set of nodes based on the given inputs, including the nodes, thermal expansion coefficient (alpha), nodal temperatures, undeformed bond information, and the resulting thermal deformation.

## Arguments

- `nodes::Union{SubArray, Vector{Int64}}`: A collection of node indices where thermal deformation is to be calculated.

- `alpha::Union{Matrix{Float64},Matrix{Int64}}`: The thermal expansion matrix, representing the material's response to temperature changes.

- `temperature::SubArray`: A SubArray containing nodal temperatures for the specified nodes.

- `undeformed_bond::SubArray`: A SubArray containing information about the undeformed bond geometry.

- `thermal_deformation::SubArray`: A SubArray to store the calculated thermal deformation for each node.

## Returns

- `thermal_deformation::SubArray`: A SubArray containing the calculated thermal deformation for the specified nodes.

## Notes

- The thermal deformation is calculated based on the formula: `thermal_deformation[iID] = alpha * temperature[iID] * undeformed_bond[iID]`.

## Example

```julia
nodes = [1, 2, 3]
alpha = [1.3 0.0; 0.0 1.3] # Example thermal expansion coefficient
temperature = SubArray(...) # Provide temperature data
undeformed_bond = SubArray(...) # Provide undeformed bond geometry data
thermal_deformation = SubArray(zeros(3, 3)) # Initialize thermal_deformation with zeros

result = thermal_deformation(nodes, alpha, temperature, undeformed_bond, thermal_deformation)
"""
function thermal_deformation(nodes::Union{SubArray,Vector{Int64}}, alpha::Union{Matrix{Float64},Matrix{Int64}}, temperature::Union{Vector{Float64},SubArray}, undeformed_bond::SubArray, thermal_deformation::SubArray)
    for iID in nodes
        for jID in eachindex(undeformed_bond[iID][:, 1])
            thermal_deformation[iID][jID, :] = -thermal_strain(alpha, temperature[iID]) * undeformed_bond[iID][jID, :]
        end
    end
    return thermal_deformation
end

"""
    thermal_strain(alpha, temperature)

Calculate thermal strain using thermal expansion coefficients and temperature.

This function calculates thermal strain based on the given thermal expansion coefficients (alpha) and temperature values. The thermal strain is computed as the element-wise product of alpha and temperature.

## Arguments

- `alpha::Union{Matrix{Float64}, Matrix{Int64}}`: A matrix of thermal expansion coefficients, representing the material's response to temperature changes.

- `temperature::Union{Float64, Int64}`: The temperature value or a scalar representing the temperature change.

## Returns

- `thermal_strain::Union{Matrix{Float64}, Matrix{Int64}}`: A matrix representing the calculated thermal strain based on alpha and temperature.

## Notes

- The thermal strain is calculated as the element-wise product of alpha and temperature.

## Example

```julia
alpha = [0.001 0.002 0.003; 0.004 0.005 0.006; 0.007 0.008 0.009] # Example thermal expansion coefficients
temperature = 100.0 # Example temperature value

result = thermal_strain(alpha, temperature)
"""
function thermal_strain(alpha::Union{Matrix{Float64},Matrix{Int64}}, temperature::Union{Float64,Int64})
    return alpha .* temperature
end

end # Module end