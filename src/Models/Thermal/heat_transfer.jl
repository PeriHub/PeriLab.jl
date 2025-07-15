# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Heat_transfer
using LinearAlgebra: dot
export compute_model
export thermal_model_name
export init_model
include("../../Support/Helpers.jl")
using .Helpers: normalize_in_place!
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
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `thermal parameter::Dict(String, Any)`: Dictionary with thermal parameter.
- `block::Int64`: The current block.
# Returns
- `datamanager::Data_manager`: Datamanager.

"""
function init_model(datamanager::Module,
                    nodes::AbstractVector{Int64},
                    thermal_parameter::Dict)
    nlist = datamanager.get_nlist()
    dof = datamanager.get_dof()

    datamanager.create_constant_node_field("Specific Volume Check", Bool, 1, true)

    undeformed_bond = datamanager.get_field("Bond Geometry")
    bond_norm_field = datamanager.create_constant_bond_field("Bond Norm", Float64, dof, 1)
    for iID in nodes
        for (jID, neighborID) in enumerate(nlist[iID])
            normalize_in_place!(bond_norm_field[iID][jID], undeformed_bond[iID][jID])
        end
    end
    return datamanager
end

"""
    compute_model(datamanager, nodes, thermal_parameter, time, dt)

Calculates the heat transfer to the environment. [BrighentiR2021](@cite)

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
    dof = datamanager.get_dof()
    volume = datamanager.get_field("Volume")
    kappa = thermal_parameter["Heat Transfer Coefficient"]
    Tenv = thermal_parameter["Environmental Temperature"]
    allow_surface_change = get(thermal_parameter, "Allow Surface Change", true)
    additive_enabled = haskey(datamanager.get_active_models(), "Additive Model")
    heat_flow = datamanager.get_field("Heat Flow", "NP1")
    temperature = datamanager.get_field("Temperature", "NP1")
    surface_nodes = datamanager.get_field("Surface_Nodes")
    specific_volume = datamanager.get_field("Specific Volume")
    active = datamanager.get_field("Active")
    bond_norm = datamanager.get_field("Bond Norm")
    rotation_tensor = datamanager.get_rotation() ?
                      datamanager.get_field("Rotation Tensor") : nothing
    specific_volume_check = datamanager.get_field("Specific Volume Check")
    nlist = datamanager.get_nlist()
    dx = 1.0

    kappa = thermal_parameter["Heat Transfer Coefficient"]

    if allow_surface_change
        calculate_specific_volume!(specific_volume,
                                   nodes,
                                   nlist,
                                   active,
                                   bond_norm,
                                   rotation_tensor,
                                   specific_volume_check,
                                   dof)
    end
    for iID in nodes
        if !surface_nodes[iID] && (additive_enabled || !allow_surface_change)
            continue
        end
        if specific_volume[iID] > 0 || !allow_surface_change
            if !additive_enabled
                surface_nodes[iID] = true
            end
            if dof == 2
                dx = sqrt(volume[iID])
            elseif dof == 3
                dx = volume[iID]^(1 / 3)
            end
            heat_flow[iID] += (kappa * (temperature[iID] - Tenv)) / dx *
                              specific_volume[iID]
        else
            surface_nodes[iID] = false
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
function calculate_specific_volume!(specific_volume::Vector{Int64},
                                    nodes::AbstractVector{Int64},
                                    nlist::Union{SubArray,Vector{Vector{Int64}}},
                                    active::Vector{Bool},
                                    bond_norm::Vector{Vector{Vector{Float64}}},
                                    rotation_tensor::Union{Array{Float64,3},Nothing},
                                    specific_volume_check::Vector{Bool},
                                    dof::Int64)
    directions = [[1, 0], [-1, 0], [0, 1], [0, -1]]
    if dof == 3
        directions = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]]
    end
    directions_free = fill(true, dof*2)
    for iID in nodes
        if !specific_volume_check[iID]
            continue
        end
        if !isnothing(rotation_tensor)
            directions = [rotation_tensor[iID, :, :] * direction
                          for direction in directions]
        end
        directions_free .= true
        for (jID, neighborID) in enumerate(nlist[iID])
            if !active[neighborID]
                continue
            end
            for (kID, direction) in enumerate(directions)
                if !directions_free[kID]
                    continue
                end
                if dot(bond_norm[iID][jID], direction) > 0.5
                    directions_free[kID] = false
                    break
                end
            end
        end
        # If the up direction is not free, the specific_volume should not change
        if dof == 3 && !directions_free[5]
            specific_volume_check[iID] = false
        end
        specific_volume[iID] = sum(directions_free)
    end
end

"""
    fields_for_local_synchronization(datamanager::Module, model::String)

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
function fields_for_local_synchronization(datamanager::Module, model::String)
    return Dict()
end

end
