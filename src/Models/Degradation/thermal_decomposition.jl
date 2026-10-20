# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Thermal_Decomposition

using .......Data_Manager

export compute_model
export init_model
export fields_for_local_synchronization

using ......ParameterSpec: @params, register_degradation
@params struct ThermalDecompositionParams
    decomposition_temperature::Float64 = req("Decomposition Temperature";
                                             quantity = :temperature)
end
__init__() = register_degradation("Thermal Decomposition", ThermalDecompositionParams)

"""
    compute_model(nodes, p, degradation, block, time, dt)

Calculates the bond-based degradation model. This template has to be copied, the file renamed and edited by the user to create a new degradation. Additional files can be called from here using include and `import .any_module` or `using .any_module`.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `p::ThermalDecompositionParams`: The model parameters.
- `degradation`: `nothing` (degradation models have no shared parameters).
- `block::Int64`: The current block.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
Example:
```julia
  ```
"""
function compute_model(nodes::AbstractVector{Int64}, p::ThermalDecompositionParams,
                       degradation, block::Int64, time::Float64, dt::Float64)
    nlist = Data_Manager.get_nlist()::BondScalarState{Int64}
    inverse_nlist = Data_Manager.get_inverse_nlist()::Vector{Dict{Int64,Int64}}
    bond_damage = Data_Manager.get_bond_damage("NP1")::BondScalarState{Float64}
    temperature = Data_Manager.get_field("Temperature", "NP1")::NodeScalarField{Float64}
    active = Data_Manager.get_field("Active")::NodeScalarField{Bool}
    number_of_neighbors = Data_Manager.get_field("Number of Neighbors")::NodeScalarField{Int64}
    decomp_temp = p.decomposition_temperature

    for iID in nodes
        if !active[iID]
            continue
        end
        if temperature[iID] >= decomp_temp
            active[iID] = false
            neighbors = nlist[iID]
            for jID in 1:number_of_neighbors[iID]
                neighborID = neighbors[jID]
                bond_damage[iID][jID] = 0.0
                if haskey(inverse_nlist[neighborID], iID)
                    bond_damage[neighborID][inverse_nlist[neighborID][iID]] = 0.0
                end
            end
        end
    end
end

"""
    init_model(nodes, p, degradation, block)

Inits the bond-based degradation model. This template has to be copied, the file renamed and edited by the user to create a new degradation. Additional files can be called from here using include and `import .any_module` or `using .any_module`.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `p::ThermalDecompositionParams`: The model parameters.
- `degradation`: `nothing` (degradation models have no shared parameters).
- `block::Int64`: The current block.

"""
function init_model(nodes::AbstractVector{Int64}, p::ThermalDecompositionParams,
                    degradation, block::Int64)
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

end
