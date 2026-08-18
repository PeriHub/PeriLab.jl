# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Heat_Transfer
using LinearAlgebra: dot
export compute_model
export thermal_model_name
export init_model
using TimerOutputs: @timeit
using .....Data_Manager
using .....Helpers: normalize_in_place!

const directions2D = [[1, 0], [-1, 0], [0, 1], [0, -1]]
const directions3D = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]]
const directions_free2D = fill(true, 4)
const directions_free3D = fill(true, 6)
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
    init_model(nodes, thermal_parameter)

Inits the thermal model. This template has to be copied, the file renamed and edited by the user to create a new thermal. Additional files can be called from here using include and `import .any_module` or `using .any_module`.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `thermal parameter::Dict(String, Any)`: Dictionary with thermal parameter.
- `block::Int64`: The current block.

"""
function init_model(nodes::AbstractVector{Int64},
                    thermal_parameter::Dict)
    nlist = Data_Manager.get_nlist()
    dof = Data_Manager.get_dof()

    Data_Manager.create_constant_node_scalar_field("Specific Volume Check", Bool;
                                                   default_value = true)

    undeformed_bond = Data_Manager.get_field("Bond Geometry")
    bond_norm_field = Data_Manager.create_constant_bond_vector_state("Bond Norm", Float64,
                                                                     dof; default_value = 1)
    for iID in nodes
        for (jID, neighborID) in enumerate(nlist[iID])
            normalize_in_place!(bond_norm_field[iID][jID], undeformed_bond[iID][jID])
        end
    end
end

"""
    compute_model(nodes, thermal_parameter, time, dt)

Calculates the heat transfer to the environment. [BrighentiR2021](@cite)

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
    dof = Data_Manager.get_dof()
    volume::NodeScalarField{Float64} = Data_Manager.get_field("Volume")
    kappa::Float64 = thermal_parameter["Heat Transfer Coefficient"]
    Tenv::Float64 = 0.0
    if thermal_parameter["Environmental Temperature"] isa String
        global t = time
        Tenv = eval(Meta.parse(thermal_parameter["Environmental Temperature"]))
    else
        Tenv = thermal_parameter["Environmental Temperature"]
    end
    allow_surface_change = get(thermal_parameter, "Allow Surface Change", true)
    additive_enabled = haskey(Data_Manager.get_active_models(), "Additive Model")
    heat_flow::NodeScalarField{Float64} = Data_Manager.get_field("Heat Flow", "NP1")
    temperature::NodeScalarField{Float64} = Data_Manager.get_field("Temperature", "NP1")
    surface_nodes::NodeScalarField{Bool} = Data_Manager.get_field("Surface_Nodes")
    specific_volume::NodeScalarField{Float64} = Data_Manager.get_field("Specific Volume")
    active::NodeScalarField{Bool} = Data_Manager.get_field("Active")
    bond_norm::BondVectorState{Float64} = Data_Manager.get_field("Bond Norm")
    rotation_tensor = Data_Manager.get_rotation() ?
                      Data_Manager.get_field("Rotation Tensor") : nothing
    specific_volume_check::NodeScalarField{Bool} = Data_Manager.get_field("Specific Volume Check")
    nlist::BondScalarState{Int64} = Data_Manager.get_nlist()
    dx = 1.0

    if allow_surface_change
        @timeit "calculate_specific_volume!" calculate_specific_volume!(specific_volume,
                                                                        nodes,
                                                                        nlist,
                                                                        active,
                                                                        bond_norm,
                                                                        rotation_tensor,
                                                                        specific_volume_check,
                                                                        dof)
    end
    @timeit "apply_surface_heat_flow!" begin
        apply_surface_heat_flow!(heat_flow, temperature, volume, specific_volume,
                                 surface_nodes, nodes, kappa, Tenv,
                                 additive_enabled, allow_surface_change, Val(dof))
    end
end

"""
    apply_surface_heat_flow!(..., ::Val{DOF})

Surface heat flow kernel, spezialisiert auf dof == 2 bzw. dof == 3.
Die `dx`-Berechnung wird über den Val-Parameter zur Compile-Zeit aufgelöst.
"""
function apply_surface_heat_flow!(heat_flow::Vector{Float64},
                                  temperature::Vector{Float64},
                                  volume::Vector{Float64},
                                  specific_volume::Vector{Float64},
                                  surface_nodes::Vector{Bool},
                                  nodes::AbstractVector{Int64},
                                  kappa::Float64,
                                  Tenv::T,
                                  additive_enabled::Bool,
                                  allow_surface_change::Bool,
                                  ::Val{2}) where {T<:Union{Int64,Float64}}
    @inbounds for iID in nodes
        if !surface_nodes[iID] && (additive_enabled || !allow_surface_change)
            continue
        end
        if specific_volume[iID] > 0 || !allow_surface_change
            if !additive_enabled
                surface_nodes[iID] = true
            end
            dx = sqrt(volume[iID])
            heat_flow[iID] += (kappa * (temperature[iID] - Tenv)) / dx *
                              specific_volume[iID]
        else
            surface_nodes[iID] = false
        end
    end
    return nothing
end

function apply_surface_heat_flow!(heat_flow::Vector{Float64},
                                  temperature::Vector{Float64},
                                  volume::Vector{Float64},
                                  specific_volume::Vector{Float64},
                                  surface_nodes::Vector{Bool},
                                  nodes::AbstractVector{Int64},
                                  kappa::Float64,
                                  Tenv::T,
                                  additive_enabled::Bool,
                                  allow_surface_change::Bool,
                                  ::Val{3}) where {T<:Union{Int64,Float64}}
    @inbounds for iID in nodes
        if !surface_nodes[iID] && (additive_enabled || !allow_surface_change)
            continue
        end
        if specific_volume[iID] > 0 || !allow_surface_change
            if !additive_enabled
                surface_nodes[iID] = true
            end
            dx = cbrt(volume[iID])
            heat_flow[iID] += (kappa * (temperature[iID] - Tenv)) / dx *
                              specific_volume[iID]
        else
            surface_nodes[iID] = false
        end
    end
    return nothing
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
function calculate_specific_volume!(specific_volume::NodeScalarField{Float64},
                                    nodes::AbstractVector{Int64},
                                    nlist::Union{SubArray,BondScalarState{Int64}},
                                    active::NodeScalarField{Bool},
                                    bond_norm::BondVectorState{Float64},
                                    rotation_tensor::Array{Float64,3},
                                    specific_volume_check::NodeScalarField{Bool},
                                    dof::Int64)
    if dof == 2
        directions = directions2D
        directions_free = directions_free2D
    else
        directions = directions3D
        directions_free = directions_free3D
    end
    rot_directions = copy(directions)

    for iID in nodes
        !specific_volume_check[iID] && continue

        rot_directions = [rotation_tensor[iID, :, :] * direction
                          for direction in directions]

        directions_free .= true
        for (jID, neighborID) in enumerate(nlist[iID])
            !active[neighborID] && continue

            for (kID, direction) in enumerate(rot_directions)
                !directions_free[kID] && continue

                if dot(bond_norm[iID][jID], direction) >= 0.6
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

function calculate_specific_volume!(specific_volume::NodeScalarField{Float64},
                                    nodes::AbstractVector{Int64},
                                    nlist::Union{SubArray,BondScalarState{Int64}},
                                    active::NodeScalarField{Bool},
                                    bond_norm::BondVectorState{Float64},
                                    rotation_tensor::Nothing,
                                    specific_volume_check::NodeScalarField{Bool},
                                    dof::Int64)
    if dof == 2
        directions = directions2D
        directions_free = directions_free2D
    else
        directions = directions3D
        directions_free = directions_free3D
    end

    for iID in nodes
        !specific_volume_check[iID] && continue
        directions_free .= true
        for (jID, neighborID) in enumerate(nlist[iID])
            !active[neighborID] && continue

            for (kID, direction) in enumerate(directions)
                !directions_free[kID] && continue

                if dot(bond_norm[iID][jID], direction) >= 0.6
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
    fields_for_local_synchronization(model::String)

Returns a user developer defined local synchronization. This happens before each model.

The structure of the Dict must because

    synchfield = Dict(
        "Field name" =>
            Dict("upload_to_cores" => true, "dof" => Data_Manager.get_dof()),
    )

or

    synchfield = Dict(
        "Field name" =>
            Dict("download_from_cores" => true, "dof" => Data_Manager.get_dof()),
    )

# Arguments

"""
function fields_for_local_synchronization(model::String)
end

end
