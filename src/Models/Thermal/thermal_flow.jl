# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Thermal_Flow
using LinearAlgebra
using StaticArrays
export compute_model
export thermal_model_name
export init_thermal_model
"""
    thermal_model_name()

Gives the model name. It is needed for comparison with the yaml input deck.

# Arguments

# Returns
- `name::String`: "Thermal Flow"
"""
function thermal_model_name()
    return "Thermal Flow"
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
function init_thermal_model(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    thermal_parameter::Dict,
)
    dof = datamanager.get_dof()
    if !haskey(thermal_parameter, "Type") || (
        thermal_parameter["Type"] != "Bond based" &&
        thermal_parameter["Type"] != "Correspondence"
    )
        @error "No model type has beed defined; ''Type'': ''Bond based'' or Type: ''Correspondence''"
        return nothing
    end

    if haskey(thermal_parameter, "Print Bed Temperature")
        if dof < 3
            @warn "Print bed temperature can only be defined for 3D problems. Its deactivated."
            delete!(thermal_parameter, "Print Bed Temperature")
        end
    end
    if !haskey(thermal_parameter, "Thermal Conductivity")
        @error "Thermal Conductivity not defined."
        return nothing
    end

    return datamanager
end

"""
    compute_model(datamanager, nodes, thermal_parameter, time, dt)

Calculates the thermal behavior of the material. This template has to be copied, the file renamed and edited by the user to create a new flow. Additional files can be called from here using include and `import .any_module` or `using .any_module`. Make sure that you return the datamanager.

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
    time::Float64,
    dt::Float64,
)

    dof = datamanager.get_dof()
    nlist = datamanager.get_nlist()
    coordinates = datamanager.get_field("Coordinates")
    bond_damage = datamanager.get_bond_damage("NP1")
    heat_flow = datamanager.get_field("Heat Flow", "NP1")
    undeformed_bond = datamanager.get_field("Bond Geometry")
    undeformed_bond_length = datamanager.get_field("Bond Length")
    volume = datamanager.get_field("Volume")
    temperature = datamanager.get_field("Temperature", "NP1")
    lambda = thermal_parameter["Thermal Conductivity"]
    apply_print_bed = false

    t_bed = 0.0
    lambda_bed = 0.0
    print_bed = nothing

    if haskey(thermal_parameter, "Print Bed Temperature")
        apply_print_bed = true
        t_bed = thermal_parameter["Print Bed Temperature"]
        lambda_bed = thermal_parameter["Thermal Conductivity Print Bed"]
        # print_bed = datamanager.get_field("Print_bed")
    end

    lambda = thermal_parameter["Thermal Conductivity"]

    if thermal_parameter["Type"] == "Bond based"
        horizon = datamanager.get_field("Horizon")
        if length(lambda) > 1
            lambda = lambda[1]
        end
        heat_flow = compute_heat_flow_state_bond_based(
            nodes,
            dof,
            nlist,
            lambda,
            apply_print_bed,
            t_bed,
            lambda_bed,
            print_bed,
            coordinates,
            bond_damage,
            undeformed_bond,
            undeformed_bond_length,
            horizon,
            temperature,
            volume,
            heat_flow,
        )
        return datamanager

    elseif thermal_parameter["Type"] == "Correspondence"

        lambda_matrix = @MMatrix zeros(Float64, dof, dof)
        Kinv = datamanager.get_field("Inverse Shape Tensor")
        if length(lambda) == 1
            for i = 1:dof
                lambda_matrix[i, i] = lambda
            end
        else
            for i = 1:dof
                lambda_matrix[i, i] = lambda[i]
            end
        end
        heat_flow = compute_heat_flow_state_correspondence(
            nodes,
            dof,
            nlist,
            lambda_matrix,
            bond_damage,
            undeformed_bond,
            Kinv,
            temperature,
            volume,
            heat_flow,
        )
    end
    return datamanager
end


"""
[BrighentiR2021](@cite)
is a prototype with some errors
"""
function compute_heat_flow_state_correspondence(
    nodes::Union{SubArray,Vector{Int64}},
    dof::Int64,
    nlist::SubArray,
    lambda::Union{Matrix{Float64},StaticArraysCore.MMatrix},
    bond_damage::SubArray,
    undeformed_bond::SubArray,
    Kinv::SubArray,
    temperature::SubArray,
    volume::SubArray,
    heat_flow::SubArray,
)

    nablaT = @MVector zeros(Float64, dof)
    H = @MVector zeros(Float64, dof)
    for iID in nodes
        H .= 0
        for (jID, neighborID) in enumerate(nlist[iID])
            temp_state =
                (temperature[neighborID] - temperature[iID]) *
                volume[neighborID] *
                bond_damage[iID][jID]
            H .+= temp_state .* undeformed_bond[iID][jID, 1:dof]
        end
        nablaT = Kinv[iID, :, :] * H
        # -> rotation must be included sometime #144
        q = lambda * nablaT
        for (jID, neighborID) in enumerate(nlist[iID])
            temp = Kinv[iID, :, :] * undeformed_bond[iID][jID, 1:dof]
            heat_flow[iID] -= dot(temp, q) * volume[neighborID]
            heat_flow[neighborID] += dot(temp, q) * volume[iID]
        end
    end
    return heat_flow

end

"""
    compute_heat_flow_state_bond_based(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, nlist::SubArray,
      lambda::Union{Float64, Int64}, bond_damage::SubArray, undeformed_bond::SubArray, horizon::SubArray,
      temperature::SubArray, heat_flow::SubArray)

Calculate heat flow based on a bond-based model for thermal analysis.

# Arguments
- `nodes::Union{SubArray,Vector{Int64}}`: An array of node indices for which heat flow should be computed.
- `dof::Int64`: The degree of freedom, either 2 or 3, indicating whether the analysis is 2D or 3D.
- `nlist::SubArray`: A SubArray representing the neighbor list for each node.
- `lambda::Union{Float64, Int64}`: The thermal conductivity.
- `apply_print_bed::Bool`: A boolean indicating whether the print bed should be applied to the thermal conductivity.
- `t_bed::Float64`: The thickness of the print bed.
- `lambda_bed::Float64`: The thermal conductivity of the print bed.
- `bond_damage::SubArray`: A SubArray representing the damage state of bonds between nodes.
- `undeformed_bond::SubArray`: A SubArray representing the geometry of the bonds.
- `undeformed_bond_length::SubArray`: A SubArray representing the undeformed bond length for each bond.
- `horizon::SubArray`: A SubArray representing the horizon for each node.
- `temperature::SubArray`: A SubArray representing the temperature at each node.
- `heat_flow::SubArray`: A SubArray where the computed heat flow values will be stored.

## Returns
- `heat_flow`: updated bond heat flow values will be stored.

## Description
This function calculates the heat flow between neighboring nodes based on a bond-based model for thermal analysis [OterkusS2014b](@cite). It considers various parameters, including thermal conductivity, damage state of bonds, geometry of bonds, horizons, temperature, and volume. The calculated bond heat flow values are stored in the `heat_flow` array.

"""
function compute_heat_flow_state_bond_based(
    nodes::Union{SubArray,Vector{Int64}},
    dof::Int64,
    nlist::SubArray,
    lambda::Union{Float64,Int64},
    apply_print_bed::Bool,
    t_bed::Float64,
    lambda_bed::Float64,
    print_bed,
    coordinates::SubArray,
    bond_damage::SubArray,
    undeformed_bond::SubArray,
    undeformed_bond_length::SubArray,
    horizon::SubArray,
    temperature::SubArray,
    volume::SubArray,
    heat_flow::SubArray,
)
    kernel::Float64 = 0.0
    for iID in nodes
        if dof == 2
            kernel = 6.0 / (pi * horizon[iID]^3)
        else
            kernel = 6.0 / (pi * horizon[iID]^4)
        end
        # if apply_print_bed && print_bed[iID] != 0
        #   temp_state = t_bed - temperature[iID]
        #   # heat_flow[iID] -= lambda_bed * temp_state * 3 * print_bed[iID]
        #   heat_flow[iID] -= lambda_bed * temp_state * 3 / print_bed[iID]
        # end
        for (jID, neighborID) in enumerate(nlist[iID])
            if bond_damage[iID][jID] == 0
                continue
            end
            if apply_print_bed
                # check if node is near print bed and if the mirror neighbor is in a higher layer
                if coordinates[iID, 3] < horizon[iID] && undeformed_bond[iID][jID, 3] > 0
                    # check if mirrored neighbor would be lower z=0
                    if coordinates[iID, 3] - undeformed_bond[iID][jID, 3] <= 0
                        temp_state = bond_damage[iID][jID] * (t_bed - temperature[iID])
                        if coordinates[iID, 3] == 0.0
                            heat_flow[iID] -=
                                lambda_bed * kernel * temp_state /
                                undeformed_bond_length[iID][jID] * volume[neighborID]
                        else
                            heat_flow[iID] -=
                                lambda_bed * kernel * temp_state / coordinates[iID, 3] *
                                volume[neighborID]
                        end
                    end
                end
            end
            temp_state =
                bond_damage[iID][jID] * (temperature[neighborID] - temperature[iID])
            heat_flow[iID] -=
                lambda * kernel * temp_state / undeformed_bond_length[iID][jID] *
                volume[neighborID]
        end
    end
    return heat_flow
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
