# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Heat_Transfer
using LinearAlgebra: dot, normalize
using StaticArrays: SVector, MVector
export compute_model
export init_model

using ......ParameterSpec: @params, register_thermal
@params struct HeatTransferParams
    heat_transfer_coefficient::Float64 = req("Heat Transfer Coefficient"; min = 0)
    environmental_temperature::Union{Float64,String} = req("Environmental Temperature";
                                                           description = "a number, or an expression in the time t",
                                                           evaluated = true)
    allow_surface_change::Bool = opt("Allow Surface Change"; default = true)
end
__init__() = register_thermal("Heat Transfer", HeatTransferParams)

"Environmental temperature: a number, or an expression in the time `t` (evaluated at `time`)."
environmental_temperature(T::Float64, time::Float64) = T
function environmental_temperature(expression::String, time::Float64)
    global t = time
    return Float64(eval(Meta.parse(expression)))
end

using TimerOutputs: @timeit
using .....Data_Manager
using .....Helpers: node_tensor

# axis directions; a direction is free if no active neighbor lies (roughly) in it
const DIRECTIONS_2D = (SVector(1.0, 0.0), SVector(-1.0, 0.0), SVector(0.0, 1.0),
                       SVector(0.0, -1.0))
const DIRECTIONS_3D = (SVector(1.0, 0.0, 0.0), SVector(-1.0, 0.0, 0.0),
                       SVector(0.0, 1.0, 0.0), SVector(0.0, -1.0, 0.0),
                       SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -1.0))
"""
    init_model(nodes, p, thermal, block)

Inits the thermal model. This template has to be copied, the file renamed and edited by the user to create a new thermal. Additional files can be called from here using include and `import .any_module` or `using .any_module`.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `p`: The model parameters.
- `thermal::WithBase`: The typed thermal model of the block.
- `block::Int64`: The current block.

"""
function init_model(nodes::AbstractVector{Int64}, p::HeatTransferParams, thermal,
                    block::Int64)
    nlist = Data_Manager.get_nlist()
    dof = Data_Manager.get_dof()

    Data_Manager.create_constant_node_scalar_field("Specific Volume Check", Bool;
                                                   default_value = true)

    undeformed_bond = Data_Manager.get_field("Bond Geometry")
    bond_norm_field = Data_Manager.create_constant_bond_vector_state("Bond Norm", Float64,
                                                                     dof; default_value = 1)
    for iID in nodes
        for (jID, neighborID) in enumerate(nlist[iID])
            bond_norm_field[iID][jID] = normalize(undeformed_bond[iID][jID])
        end
    end
end

"""
    compute_model(nodes, p, thermal, block, time, dt)

Calculates the heat transfer to the environment. [BrighentiR2021](@cite)

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `p::HeatTransferParams`: The model parameters.
- `thermal::WithBase`: The typed thermal model of the block.
- `block::Int64`: The current block.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
Example:
```julia
```
"""
function compute_model(nodes::AbstractVector{Int64}, p::HeatTransferParams, thermal,
                       block::Int64, time::Float64, dt::Float64)
    dof = Data_Manager.get_dof()
    volume::NodeScalarField{Float64} = Data_Manager.get_field("Volume")
    kappa::Float64 = p.heat_transfer_coefficient
    Tenv::Float64 = environmental_temperature(p.environmental_temperature, time)
    allow_surface_change = p.allow_surface_change
    additive_enabled = haskey(Data_Manager.get_active_models(), "Additive Model")
    heat_flow::NodeScalarField{Float64} = Data_Manager.get_field("Heat Flow", "NP1")
    temperature::NodeScalarField{Float64} = Data_Manager.get_field("Temperature", "NP1")
    surface_nodes::NodeScalarField{Bool} = Data_Manager.get_field("Surface_Nodes")
    specific_volume::NodeScalarField{Float64} = Data_Manager.get_field("Specific Volume")
    active::NodeScalarField{Bool} = Data_Manager.get_field("Active")
    bond_norm::BondVectorField{Float64} = Data_Manager.get_field("Bond Norm")
    rotation_tensor = Data_Manager.get_rotation() ?
                      Data_Manager.get_field("Rotation Tensor")::Array{Float64,3} : nothing
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

"""
    calculate_specific_volume!(specific_volume::NodeScalarField{Int64}, nodes::AbstractVector{Int64}, nlist::BondScalarState{Int64}, active::NodeScalarField{Bool}, bond_norm::BondVectorState{Float64}, rotation_tensor::Union{Array{Float64,3},Nothing}, specific_volume_check::NodeScalarField{Bool}, dof::Int64)

Counts the number of free surface directions for each node to determine its specific volume.
A direction is considered "free" if no neighbor exists along that axis (within a cosine threshold).

# Arguments
- `specific_volume`: Output field storing the count of free directions per node.
- `nodes`: Vector of active node indices.
- `nlist`: Neighbor list for each node.
- `active`: Boolean mask indicating which nodes are active.
- `bond_norm`: Normalized bond vectors to neighbors.
- `rotation_tensor`: Optional rotation tensor (for rotated meshes). If nothing, directions are not rotated.
- `specific_volume_check`: Boolean field indicating which nodes need volume recalculation.
- `dof`: Spatial dimensionality (2 or 3).

# Returns
- Updates `specific_volume` in-place with the count of free surface directions.
"""
function calculate_specific_volume!(specific_volume::NodeScalarField{Float64},
                                    nodes::AbstractVector{Int64},
                                    nlist::Union{SubArray,BondScalarState{Int64}},
                                    active::NodeScalarField{Bool},
                                    bond_norm::BondVectorState{Float64,D},
                                    rotation_tensor::Union{Nothing,Array{Float64,3}},
                                    specific_volume_check::NodeScalarField{Bool},
                                    dof::Int64) where {D}
    directions = D == 2 ? DIRECTIONS_2D : DIRECTIONS_3D
    for iID in nodes
        !specific_volume_check[iID] && continue
        rotated = rotate_directions(directions, rotation_tensor, iID)
        free = MVector(map(_ -> true, directions))
        norms = bond_norm[iID]
        for (jID, neighborID) in enumerate(nlist[iID])
            !active[neighborID] && continue
            for kID in eachindex(rotated)
                !free[kID] && continue
                if dot(norms[jID], rotated[kID]) >= 0.6
                    free[kID] = false
                    break
                end
            end
        end
        # If the up direction is not free, the specific_volume should not change
        if D == 3 && !free[5]
            specific_volume_check[iID] = false
        end
        specific_volume[iID] = count(free)
    end
end

rotate_directions(directions, ::Nothing, iID::Int64) = directions
function rotate_directions(directions::Tuple{SVector{D,Float64},Vararg{SVector{D,Float64}}},
                           rotation_tensor::Array{Float64,3}, iID::Int64) where {D}
    R = node_tensor(rotation_tensor, iID, Val(D))
    return map(direction -> R * direction, directions)
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
