# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Damage_Based
export compute_model
export additive_name
export init_model
export fields_for_local_synchronization
"""
    additive_name()

Gives the additive name. It is needed for comparison with the yaml input deck.

# Arguments

# Returns
- `name::String`: The name of the additive model.

Example:
```julia
println(additive_name())
"additive Template"
```
"""
function additive_name()
    return "Simple"
end

"""
    compute_model(
    datamanager::Module,
    nodes::AbstractVector{Int64},
    additive_parameter::Dict,
    block::Int64,
    time::Float64,
    dt::Float64,
)

Calculates the force densities of the additive. This template has to be copied, the file renamed and edited by the user to create a new additive. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

# Arguments
- `datamanager::Data_Manager`: Datamanager.
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `additive parameter::Dict(String, Any)`: Dictionary with additive parameter.
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
                       additive_parameter::Dict,
                       block::Int64,
                       time::Float64,
                       dt::Float64)
    nlist = datamanager.get_nlist()
    inverse_nlist = datamanager.get_inverse_nlist()
    activation_time = datamanager.get_field("Activation_Time")
    bond_damage = datamanager.get_bond_damage("NP1")
    active = datamanager.get_field("Active")
    number_of_neighbors = datamanager.get_field("Number of Neighbors")
    # must be specified, because it might be that no temperature model has been defined
    flux = datamanager.get_field("Heat Flow", "NP1")
    add_flux = datamanager.get_field("Additive Heat Flux")
    nn::Int64 = 1
    ###########
    for iID in nodes
        if active[iID]
            continue
        end
        if time - dt <= activation_time[iID] < time
            active[iID] = true
            flux[iID] = add_flux[iID] / dt
            nlist_temp = nlist[iID]
            nn = number_of_neighbors[iID]
            for jID in 1:nn
                @views neighborID = nlist_temp[jID]
                if activation_time[neighborID] > time
                    continue
                end
                bond_damage[iID][jID] = 1.0
                if haskey(inverse_nlist[neighborID], iID)
                    bond_damage[neighborID][inverse_nlist[neighborID][iID]] = 1.0
                end
            end
        end
    end

    return datamanager
end

"""
    init_model(datamanager, nodes, additive_parameter)

Inits the simple additive model.

# Arguments
- `datamanager::Data_Manager`: Datamanager.
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `additive parameter::Dict(String, Any)`: Dictionary with additive parameter.
- `block::Int64`: The current block.
# Returns
- `datamanager::Data_Manager`: Datamanager.

"""
function init_model(datamanager::Module,
                    nodes::AbstractVector{Int64},
                    additive_parameter::Dict,
                    block::Int64)
    add_flux = datamanager.create_constant_node_field("Additive Heat Flux", Float64, 1)
    heat_capacity = datamanager.get_field("Specific Heat Capacity")
    density = datamanager.get_field("Density")

    printTemperature = additive_parameter["Print Temperature"]

    for iID in nodes
        add_flux[iID] = -printTemperature * heat_capacity[iID] * density[iID]
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

end
