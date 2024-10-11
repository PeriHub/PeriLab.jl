# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Heat_transfer
export compute_model
export thermal_model_name
export init_model
"""
    thermal_model_name()

Gives the flow name. It is needed for comparison with the yaml input deck.

# Arguments

# Returns
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
    init_model(datamanager, nodes, thermal_parameter)

Inits the thermal model. This template has to be copied, the file renamed and edited by the user to create a new thermal. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `thermal parameter::Dict(String, Any)`: Dictionary with thermal parameter.
- `block::Int64`: The current block.
# Returns
- `datamanager::Data_manager`: Datamanager.

"""
function init_model(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    thermal_parameter::Dict,
)

    return datamanager
end

"""
    compute_model(datamanager, nodes, thermal_parameter, time, dt)

Calculates the heat transfer to the environment. [BrighentiR2021](@cite)

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
function compute_model(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    thermal_parameter::Dict,
    block::Int64,
    time::Float64,
    dt::Float64,
)
    dof = datamanager.get_dof()
    volume = datamanager.get_field("Volume")
    kappa = thermal_parameter["Heat Transfer Coefficient"]
    Tenv = thermal_parameter["Environmental Temperature"]
    req_specific_volume = get(thermal_parameter, "Required Specific Volume", 1.1)
    heat_flow = datamanager.get_field("Heat Flow", "NP1")
    temperature = datamanager.get_field("Temperature", "NP1")
    surface_nodes = datamanager.get_field("Surface_Nodes")
    specific_volume = datamanager.get_field("Specific Volume")
    active = datamanager.get_field("Active")
    horizon = datamanager.get_field("Horizon")
    nlist = datamanager.get_nlist()
    dx = 1.0

    specific_volume = calculate_specific_volume(
        nodes,
        nlist,
        volume,
        active,
        specific_volume,
        dof,
        horizon,
    )
    for iID in nodes

        if surface_nodes[iID] && specific_volume[iID] > req_specific_volume

            if dof == 2
                dx = sqrt(volume[iID])
            elseif dof == 3
                dx = volume[iID]^(1 / 3)
            end

            heat_flow[iID] +=
                (kappa * (temperature[iID] - Tenv)) / dx * floor(specific_volume[iID])
        end
    end

    return datamanager
end

#TODO @Jan-Timo update documentation
"""
  calculate_specific_volume(nodes::Int64, nlist::SubArray, coordinates::Union{SubArray,Vector{Float64}}, volume::SubArray, surface_nodes::Union{SubArray,Vector{Bool}})

Calculates the specific volume.

# Arguments
- `iID::Int64`: The index of the node.
- `nlist::SubArray`: The neighbor list.
- `coordinates::Union{SubArray,Vector{Float64}}`: The coordinates of the nodes.
- `volume::SubArray`: The volume of the nodes.
- `surface_nodes::Union{SubArray,Vector{Bool}}`: The surface nodes.
# Returns
- `specific_volume::Union{SubArray,Vector{Bool}}`: The surface nodes.
"""
function calculate_specific_volume(
    nodes::Union{SubArray,Vector{Int64}},
    nlist::Union{SubArray,Vector{Vector{Int64}}},
    volume::Vector{Float64},
    active::Vector{Bool},
    specific_volume::Vector{Float64},
    dof::Int64,
    horizon::Vector{Float64},
)

    for iID in nodes
        neighbor_volume = 0.0
        for (jID, neighborID) in enumerate(nlist[iID])
            if active[neighborID]
                neighbor_volume += volume[neighborID]
            end
        end
        if dof == 2
            horizon_volume = pi * horizon[iID]^2
        elseif dof == 3
            horizon_volume = 4 / 3 * pi * horizon[iID]^3
        end
        if neighbor_volume != 0.0
            specific_volume[iID] = horizon_volume / neighbor_volume
        end
    end

    if dof == 2
        scaling_factor = 2.0 / maximum(specific_volume)
    elseif dof == 3
        scaling_factor = 4.0 / maximum(specific_volume)
    end
    if !isinf(scaling_factor)
        specific_volume .*= scaling_factor
    end

    return specific_volume
end

"""
    fields_for_local_synchronization()

Returns a user developer defined local synchronization. This happens before each model.

The structure of the Dict must because

    synchfield = Dict(
        "Field name" =>
            Dict("upload_to_cores" => true, "dof" => datamanager.get_dof()),
    )

or

    synchfield = Dict(
        "Field name" =>
            Dict("download_from_cores" => true, "dof" => datamanager.get_dof()),
    )

# Arguments

"""
function fields_for_local_synchronization()
    return Dict()
end

end
