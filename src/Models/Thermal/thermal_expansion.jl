# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Thermal_Expansion
using LinearAlgebra
using StaticArrays
using TimerOutputs: @timeit
using .....Data_Manager
using ...Pre_Calculation.Deformation_Gradient: compute
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
        return alpha_mat = SMatrix{2,2,Float64}(alpha, 0.0, 0.0, alpha)
    elseif length(alpha) == 2
        return alpha_mat = SMatrix{2,2,Float64}(alpha[1], 0.0, 0.0, alpha[2])
    elseif length(alpha) == 4
        @error "Full heat expansion matrix is not implemented yet."
        return alpha_mat = SMatrix{2,2,Float64}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    end
end
function thermal_expansion_matrix(alpha::T,
                                  ::Val{3}) where {T<:Union{Float64,Matrix{Float64}}}
    if length(alpha) == 1
        return alpha_mat = SMatrix{3,3,Float64}(alpha, 0.0, 0.0, 0.0, alpha, 0.0, 0.0, 0.0,
                                                alpha)
    elseif length(alpha) == 3
        return alpha_mat = SMatrix{3,3,Float64}(alpha[1], 0.0, 0.0, 0.0, alpha[2], 0.0, 0.0,
                                                0.0,
                                                alpha[3])
    elseif length(alpha) == 9
        @error "Full heat expansion matrix is not implemented yet."
        return alpha_mat = SMatrix{3,3,Float64}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    end
end

"""
    init_model(nodes, thermal_parameter)

Inits the thermal model. This template has to be copied, the file renamed and edited by the user to create a new thermal. Additional files can be called from here using include and `import .any_module` or `using .any_module`.
# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `thermal parameter::Dict(String, Any)`: Dictionary with thermal parameter.
- `block::Int64`: The current block.

"""
function init_model(nodes::AbstractVector{Int64},
                    thermal_parameter::Dict)
    if !haskey(thermal_parameter, "Reference Temperature")
        @warn "No reference temperature defined. Assuming 0"
    end
end

"""
    compute_model(nodes, thermal_parameter, block::Int64, time, dt)

Calculates the thermal expansion of the material.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `flow parameter::Dict(String, Any)`: Dictionary with flow parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
Example:
```julia
```
"""
function compute_model(nodes::AbstractVector{Int64},
                       thermal_parameter::Dict,
                       block::Int64,
                       time::Float64,
                       dt::Float64)
    temperature_NP1::NodeScalarField{Float64} = Data_Manager.get_field("Temperature",
                                                                       "NP1")
    dof::Int64 = Data_Manager.get_dof()

    @timeit "thermal_expansion_matrix" alpha_mat::AbstractMatrix{Float64}=thermal_expansion_matrix(thermal_parameter["Thermal Expansion Coefficient"],
                                                                                                   Val(dof))
    ref_temp::Float64 = get(thermal_parameter, "Reference Temperature", 0.0)

    undeformed_bond::BondVectorState{Float64} = Data_Manager.get_field("Bond Geometry")
    undeformed_bond_length::BondScalarState{Float64} = Data_Manager.get_field("Bond Length")
    deformed_bond::BondVectorState{Float64} = Data_Manager.get_field("Deformed Bond Geometry",
                                                                     "NP1")
    deformed_bond_length::BondScalarState{Float64} = Data_Manager.get_field("Deformed Bond Length",
                                                                            "NP1")

    @timeit "apply_thermal_expansion!" apply_thermal_expansion!(deformed_bond,
                                                                deformed_bond_length,
                                                                nodes,
                                                                temperature_NP1,
                                                                ref_temp,
                                                                alpha_mat,
                                                                undeformed_bond,
                                                                undeformed_bond_length,
                                                                dof)

    if Data_Manager.has_key("Deformation Gradient")
        #TODO all forces computed are from the original configuration
        @timeit "Deformation_Gradient" compute(nodes, Dict(), block)
    end
end

"""
    apply_thermal_expansion!(deformed_bond, deformed_bond_length, nodes,
                            temperature_NP1, ref_temp, alpha_mat,
                            undeformed_bond, undeformed_bond_length, dof)

Apply thermal expansion corrections to deformed bonds and bond lengths.

# Arguments
- `deformed_bond`: Array of deformed bond vectors (modified in-place)
- `deformed_bond_length`: Array of deformed bond lengths (modified in-place)
- `nodes`: Node indices to process
- `temperature_NP1`: Current temperature field
- `ref_temp`: Reference temperature
- `alpha_mat`: Thermal expansion coefficient matrix
- `undeformed_bond`: Array of undeformed bond vectors
- `undeformed_bond_length`: Array of undeformed bond lengths
- `dof`: Degrees of freedom
"""
function apply_thermal_expansion!(deformed_bond::BondVectorState{Float64},
                                  deformed_bond_length::BondScalarState{Float64},
                                  nodes::AbstractVector{Int64},
                                  temperature_NP1::NodeScalarField{Float64},
                                  ref_temp::Float64,
                                  alpha_mat::AbstractMatrix{Float64},
                                  undeformed_bond::BondVectorState{Float64},
                                  undeformed_bond_length::BondScalarState{Float64},
                                  dof::Int64)
    # Precompute average thermal expansion coefficient
    avg_alpha = sum(alpha_mat) / dof

    for iID in nodes
        temp_diff::Float64 = temperature_NP1[iID] - ref_temp

        # Skip if no temperature difference
        if iszero(temp_diff)
            continue
        end

        # Apply thermal expansion to bond vectors
        @inbounds @fastmath @views for jID in eachindex(undeformed_bond[iID])
            for j in 1:dof
                deformed_bond[iID][jID][j] -= temp_diff * alpha_mat[j, j] *
                                              undeformed_bond[iID][jID][j]
            end
        end

        # Apply thermal expansion to bond lengths
        deformed_bond_length[iID] .-= avg_alpha * temp_diff .*
                                      undeformed_bond_length[iID]
    end

    return nothing
end

"""
    fields_for_local_synchronization(model::String)

Returns a user developer defined local synchronization. This happens before each model.



# Arguments

"""
function fields_for_local_synchronization(model::String)
    #download_from_cores = false
    #upload_to_cores = true
    #Data_Manager.set_local_synch(model, "Bond Forces", download_from_cores, upload_to_cores)
end

end # Module end
