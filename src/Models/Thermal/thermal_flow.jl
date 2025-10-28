# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Thermal_Flow
using LinearAlgebra
using StaticArrays

using .....Data_Manager
using .....Helpers: rotate_second_order_tensor
export compute_model
export thermal_model_name
export init_model
export fields_for_local_synchronization
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
    init_model(nodes, thermal_parameter, block)

Inits the thermal model. This template has to be copied, the file renamed and edited by the user to create a new thermal. Additional files can be called from here using include and `import .any_module` or `using .any_module`.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `thermal parameter::Dict(String, Any)`: Dictionary with thermal parameter.
- `block::Int64`: The current block.

"""
function init_model(nodes::AbstractVector{Int64},
                    thermal_parameter::Dict)
    dof = Data_Manager.get_dof()
    if !haskey(thermal_parameter, "Type") || (thermal_parameter["Type"] != "Bond based" &&
        thermal_parameter["Type"] != "Correspondence")
        @error "No model type has beed defined; ''Type'': ''Bond based'' or Type: ''Correspondence''"
        return nothing
    end

    if haskey(thermal_parameter, "Print Bed Temperature")
        if dof < 3
            @warn "Print bed temperature can only be defined for 3D problems. Its deactivated."
            delete!(thermal_parameter, "Print Bed Temperature")
        else
            coordinates = Data_Manager.get_field("Coordinates")
            print_bed_z_coord = get(thermal_parameter, "Print Bed Z Coordinate", 0.0)
            if print_bed_z_coord >= minimum(coordinates[:, 3])
                @error "The Print Bed Z Coordinate needs to be smaller than the minimum Z coordinate."
            end
        end
    end
    if !haskey(thermal_parameter, "Thermal Conductivity")
        @error "Thermal Conductivity not defined."
        return nothing
    end
end

"""
    compute_model(nodes, thermal_parameter, time, dt)

Calculates the thermal behavior of the material. This template has to be copied, the file renamed and edited by the user to create a new flow. Additional files can be called from here using include and `import .any_module` or `using .any_module`.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `thermal_parameter::Dict(String, Any)`: Dictionary with flow parameter.
- `block::Int64`: Current block
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
    dof = Data_Manager.get_dof()
    nlist = Data_Manager.get_nlist()
    coordinates = Data_Manager.get_field("Coordinates")
    bond_damage = Data_Manager.get_bond_damage("NP1")
    heat_flow = Data_Manager.get_field("Heat Flow", "NP1")
    undeformed_bond = Data_Manager.get_field("Bond Geometry")
    undeformed_bond_length = Data_Manager.get_field("Bond Length")
    volume = Data_Manager.get_field("Volume")
    temperature = Data_Manager.get_field("Temperature", "NP1")
    active = Data_Manager.get_field("Active")
    lambda = thermal_parameter["Thermal Conductivity"]
    rotation::Bool = Data_Manager.get_element_rotation()
    rotation_tensor = nothing
    if rotation
        rotation_tensor = Data_Manager.get_field("Rotation Tensor")
    end
    apply_print_bed = false

    t_bed = 0.0
    lambda_bed = 0.0

    if haskey(thermal_parameter, "Print Bed Temperature")
        apply_print_bed = true
        t_bed = thermal_parameter["Print Bed Temperature"]
        lambda_bed = thermal_parameter["Thermal Conductivity Print Bed"]
    end

    print_bed_z_coord = get(thermal_parameter, "Print Bed Z Coordinate", 0.0)
    lambda = thermal_parameter["Thermal Conductivity"]

    if thermal_parameter["Type"] == "Bond based"
        horizon = Data_Manager.get_field("Horizon")
        if length(lambda) > 1
            lambda = lambda[1]
        end
        heat_flow = compute_heat_flow_state_bond_based(nodes,
                                                       dof,
                                                       nlist,
                                                       lambda,
                                                       apply_print_bed,
                                                       t_bed,
                                                       lambda_bed,
                                                       print_bed_z_coord,
                                                       coordinates,
                                                       bond_damage,
                                                       active,
                                                       undeformed_bond,
                                                       undeformed_bond_length,
                                                       horizon,
                                                       temperature,
                                                       volume,
                                                       heat_flow)
        return

    elseif thermal_parameter["Type"] == "Correspondence"
        lambda_matrix = @MMatrix zeros(Float64, dof, dof)
        Kinv = Data_Manager.get_field("Inverse Shape Tensor")
        if length(lambda) == 1
            for i in 1:dof
                lambda_matrix[i, i] = lambda
            end
        else
            for i in 1:dof
                lambda_matrix[i, i] = lambda[i]
            end
        end
        heat_flow = compute_heat_flow_state_correspondence(nodes,
                                                           dof,
                                                           nlist,
                                                           lambda_matrix,
                                                           rotation_tensor,
                                                           bond_damage,
                                                           undeformed_bond,
                                                           Kinv,
                                                           temperature,
                                                           volume,
                                                           heat_flow)
    end
end

"""
[BrighentiR2021](@cite)
is a prototype with some errors
"""
function compute_heat_flow_state_correspondence(nodes::AbstractVector{Int64},
                                                dof::Int64,
                                                nlist::Vector{Vector{Int64}},
                                                lambda::Union{Matrix{Float64},MMatrix},
                                                rotation_tensor,
                                                bond_damage::Vector{Vector{Float64}},
                                                undeformed_bond::Vector{Vector{Vector{Float64}}},
                                                Kinv::Array{Float64,3},
                                                temperature::Vector{Float64},
                                                volume::Vector{Float64},
                                                heat_flow::Vector{Float64})
    nablaT = @MVector zeros(Float64, dof)
    H = @MVector zeros(Float64, dof)
    for iID in nodes
        H .= 0
        for (jID, neighborID) in enumerate(nlist[iID])
            temp_state = (temperature[neighborID] - temperature[iID]) *
                         volume[neighborID] *
                         bond_damage[iID][jID]
            H .+= temp_state .* undeformed_bond[iID][jID]
        end
        nablaT = Kinv[iID, :, :] * H

        if isnothing(rotation_tensor)
            q = lambda * nablaT
        else
            q = rotate_second_order_tensor(rotation_tensor[iID, :, :], lamba[:, :], false)
        end
        for (jID, neighborID) in enumerate(nlist[iID])
            temp = Kinv[iID, :, :] * undeformed_bond[iID][jID]
            heat_flow[iID] -= dot(temp, q) * volume[neighborID]
            heat_flow[neighborID] += dot(temp, q) * volume[iID]
        end
    end
    return heat_flow
end

"""
    compute_heat_flow_state_bond_based(nodes::AbstractVector{Int64}, dof::Int64, nlist::Vector{Vector{Int64},
      lambda::Union{Float64, Int64}, bond_damage::Vector{Vector{Float64}}, undeformed_bond::Vector{Matrix{Float64}}, horizon::Vector{Float64},
      temperature::Vector{Float64}, heat_flow::Vector{Float64})

Calculate Heat Flow based on a bond-based model for thermal analysis.

# Arguments
- `nodes::AbstractVector{Int64}`: An array of node indices for which Heat Flow should be computed.
- `dof::Int64`: The degree of freedom, either 2 or 3, indicating whether the analysis is 2D or 3D.
- `nlist::Vector{Vector{Int64}`: A Vector representing the neighbor list for each node.
- `lambda::Union{Float64, Int64}`: The thermal conductivity.
- `apply_print_bed::Bool`: A boolean indicating whether the print bed should be applied to the thermal conductivity.
- `t_bed::Float64`: The thickness of the print bed.
- `lambda_bed::Float64`: The thermal conductivity of the print bed.
- `bond_damage::Vector{Vector{Float64}}`: A Vector representing the damage state of bonds between nodes.
- `undeformed_bond::Vector{Matrix{Float64}}`: A Vector representing the geometry of the bonds.
- `undeformed_bond_length::Vector{Vector{Float64}}`: A Vector representing the undeformed bond length for each bond.
- `horizon::Vector{Float64}`: A Vector representing the horizon for each node.
- `temperature::Vector{Float64}`: A Vector representing the temperature at each node.
- `heat_flow::Vector{Float64}`: A Vector where the computed Heat Flow values will be stored.

## Returns
- `heat_flow`: updated bond Heat Flow values will be stored.

## Description
This function calculates the Heat Flow between neighboring nodes based on a bond-based model for thermal analysis [OterkusS2014b](@cite). It considers various parameters, including thermal conductivity, damage state of bonds, geometry of bonds, horizons, temperature, and volume. The calculated bond Heat Flow values are stored in the `heat_flow` array.

"""
function compute_heat_flow_state_bond_based(nodes::AbstractVector{Int64},
                                            dof::Int64,
                                            nlist::Vector{Vector{Int64}},
                                            lambda::Union{Float64,Int64},
                                            apply_print_bed::Bool,
                                            t_bed::Float64,
                                            lambda_bed::Float64,
                                            print_bed_z_coord::Float64,
                                            coordinates::Matrix{Float64},
                                            bond_damage::Vector{Vector{Float64}},
                                            active::Vector{Bool},
                                            undeformed_bond::Vector{Vector{Vector{Float64}}},
                                            undeformed_bond_length::Vector{Vector{Float64}},
                                            horizon::Vector{Float64},
                                            temperature::Vector{Float64},
                                            volume::Vector{Float64},
                                            heat_flow::Vector{Float64})
    kernel::Float64 = 0.0
    for iID in nodes
        if !active[iID]
            continue
        end
        if dof == 2
            kernel = 6.0 / (pi * horizon[iID]^3)
        else
            kernel = 6.0 / (pi * horizon[iID]^4)
        end
        if apply_print_bed
            print_bed_distance = coordinates[iID, 3] - print_bed_z_coord
            if print_bed_distance < horizon[iID]
                temp_state = t_bed - temperature[iID]
                print_bed_volume = (pi*print_bed_distance^2/3)*(3*horizon[iID]-print_bed_distance) #Partial spherical volume below print bed
                heat_flow[iID] -= lambda_bed * kernel * temp_state * print_bed_volume /
                                  print_bed_distance
            end
        end
        for (jID, neighborID) in enumerate(nlist[iID])
            if bond_damage[iID][jID] == 0
                continue
            end
            temp_state = bond_damage[iID][jID] *
                         (temperature[neighborID] - temperature[iID])
            heat_flow[iID] -= lambda * kernel * temp_state *
                              volume[neighborID]/undeformed_bond_length[iID][jID]
        end
    end
    return heat_flow
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
