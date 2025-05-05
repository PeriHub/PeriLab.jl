# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Thermal_Decomposition
using TimerOutputs
export compute_model
export degradation_name
export init_model
export fields_for_local_synchronization

"""
    degradation_name()

Gives the degradation name. It is needed for comparison with the yaml input deck.

# Arguments

# Returns
- `name::String`: The name of the bond-based degradation model.

Example:
```julia
println(degradation_name())
"Bond-based Corrosion"
```
"""
function degradation_name()
    return "Thermal Decomposition"
end

"""
    compute_model(datamanager, nodes, degradation_parameter, block::Int64, time, dt)

Calculates the bond-based degradation model. This template has to be copied, the file renamed and edited by the user to create a new degradation. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `degradation parameter::Dict(String, Any)`: Dictionary with degradation parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
# Returns
- `datamanager::Data_manager`: Datamanager.
Example:
```julia
  ```
"""
function compute_model(datamanager::Module,
                       nodes::Union{SubArray,Vector{Int64}},
                       degradation_parameter::Dict,
                       block::Int64,
                       time::Float64,
                       dt::Float64)
    nlist = datamanager.get_nlist()
    inverse_nlist = datamanager.get_inverse_nlist()
    bond_damage = datamanager.get_bond_damage("NP1")
    temperature = datamanager.get_field("Temperature", "NP1")
    active = datamanager.get_field("Active")
    number_of_neighbors = datamanager.get_field("Number of Neighbors")

    decomp_temp = degradation_parameter["Decomposition Temperature"]
    nn::Int64 = 1

    for iID in nodes
        if !active[iID]
            continue
        end
        if temperature[iID] >= decomp_temp
            active[iID] = false
            nlist_temp = nlist[iID]
            nn = number_of_neighbors[iID]
            for jID in 1:nn
                @views neighborID = nlist_temp[jID]
                bond_damage[iID][jID] = 0.0
                if haskey(inverse_nlist[neighborID], iID)
                    bond_damage[neighborID][inverse_nlist[neighborID][iID]] = 0.0
                end
            end
        end
    end

    return datamanager
end

"""
    init_model(datamanager, nodes, block::Int64, degradation_parameter)

Inits the bond-based degradation model. This template has to be copied, the file renamed and edited by the user to create a new degradation. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `degradation parameter::Dict(String, Any)`: Dictionary with degradation parameter.
- `block::Int64`: The current block.
# Returns
- `datamanager::Data_manager`: Datamanager.

"""
function init_model(datamanager::Module,
                    nodes::Union{SubArray,Vector{Int64}},
                    degradation_parameter::Dict,
                    block::Int64)
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

end
