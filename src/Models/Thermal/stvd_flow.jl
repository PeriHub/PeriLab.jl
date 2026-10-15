# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module STVD_Flow
using NearestNeighbors: KDTree, nn
using .....Data_Manager
using .....PeriLabExceptions: @abort
export compute_model
export thermal_model_name
export init_model
export fields_for_local_synchronization

"""
	thermal_model_name()

Gives the model name. It is needed for comparison with the yaml input deck.

# Returns
- `name::String`: "STVD Flow"
"""
function thermal_model_name()
    return "STVD Flow"
end

"""
    init_model(nodes, thermal_parameter)

Inits the spatio-temporal variable diffusivity (STVD) model [AtiehM2021](@cite). The
temperature field is interpreted as concentration of the diffusant. Nodes of the
diffusant blocks form an inert core with zero diffusivity. The core is consumed by a
reaction front which advances from the diffusant surface with the speed `Reaction Rate`.
For each diffusant node the distance to the diffusant surface is stored in "STVD Depth".

The front position is defined by the `Reaction Law`:
- `"Linear"` (default): front depth = `Reaction Rate` * t [AtiehM2022](@cite)
- `"Exponential"`: the core shrinks with exp(-`Reaction Rate` * t), i.e.
  front depth = `Core Size` * (1 - exp(-`Reaction Rate` * t)) [AtiehM2021](@cite)

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `thermal_parameter::AbstractDict(String, Any)`: Dictionary with thermal parameter.
"""
function init_model(nodes::AbstractVector{Int64},
                    thermal_parameter::AbstractDict)
    if !haskey(thermal_parameter, "Thermal Conductivity")
        @abort "Thermal Conductivity not defined."
    end
    if !haskey(thermal_parameter, "Reaction Rate")
        @abort "Reaction Rate not defined."
    end
    if !haskey(thermal_parameter, "Diffusant Block ID")
        @abort "Diffusant Block ID not defined."
    end
    reaction_law = get!(thermal_parameter, "Reaction Law", "Linear")
    if !(reaction_law in ["Linear", "Exponential"])
        @abort "Reaction Law $reaction_law not supported. Use ''Linear'' or ''Exponential''."
    end
    if reaction_law == "Exponential" && !haskey(thermal_parameter, "Core Size")
        @abort "Core Size has to be defined for the exponential reaction law."
    end
    diffusant_blocks = parse_block_ids(thermal_parameter["Diffusant Block ID"])
    depth = Data_Manager.create_constant_node_scalar_field("STVD Depth", Float64)
    compute_depth!(depth,
                   Data_Manager.get_field("Coordinates"),
                   Data_Manager.get_field("Block_Id"),
                   diffusant_blocks,
                   Data_Manager.get_field("Bond Length"))
end

parse_block_ids(ids::Int64) = [ids]
parse_block_ids(ids::AbstractString) = parse.(Int64, split(ids))

"""
    compute_depth!(depth, coordinates, block_ids, diffusant_blocks, bond_length)

Computes the distance of each diffusant node to the diffusant surface. The surface is
assumed half a point spacing outside the outermost diffusant nodes. Substrate nodes get a
negative depth and are therefore active from the beginning.
"""
function compute_depth!(depth::AbstractVector{Float64},
                        coordinates::AbstractMatrix{Float64},
                        block_ids::AbstractVector{Int64},
                        diffusant_blocks::AbstractVector{Int64},
                        bond_length)
    is_diffusant = [block_id in diffusant_blocks for block_id in block_ids]
    depth .= -1.0
    substrate = findall(.!is_diffusant)
    if isempty(substrate)
        @abort "STVD Flow needs at least one substrate node outside of the diffusant blocks."
    end
    tree = KDTree(Matrix(transpose(coordinates[substrate, :])))
    for iID in findall(is_diffusant)
        _, dist = nn(tree, coordinates[iID, :])
        depth[iID] = dist - 0.5 * minimum(bond_length[iID])
    end
end

"""
    compute_model(nodes, thermal_parameter, block, time, dt)

Calculates the concentration flow of the STVD model. A bond transports concentration only
if both nodes are reached by the reaction front, i.e. `depth <= Reaction Rate * time`.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `thermal_parameter::AbstractDict(String, Any)`: Dictionary with flow parameter.
- `block::Int64`: Current block
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
"""
function compute_model(nodes::AbstractVector{Int64},
                       thermal_parameter::AbstractDict,
                       block::Int64,
                       time::Float64,
                       dt::Float64)
    front_depth = compute_front_depth(thermal_parameter, time)
    compute_stvd_flow!(nodes,
                       Data_Manager.get_dof(),
                       Data_Manager.get_nlist(),
                       thermal_parameter["Thermal Conductivity"],
                       front_depth,
                       Data_Manager.get_field("STVD Depth"),
                       Data_Manager.get_bond_damage("NP1"),
                       Data_Manager.get_field("Bond Length"),
                       Data_Manager.get_field("Horizon"),
                       Data_Manager.get_field("Temperature", "NP1"),
                       Data_Manager.get_field("Volume"),
                       Data_Manager.get_field("Heat Flow", "NP1"))
end

"""
    compute_front_depth(thermal_parameter, time)

Depth of the reaction front below the initial diffusant surface.
"""
function compute_front_depth(thermal_parameter::AbstractDict, time::Float64)
    k = thermal_parameter["Reaction Rate"]
    if get(thermal_parameter, "Reaction Law", "Linear") == "Exponential"
        return thermal_parameter["Core Size"] * (1 - exp(-k * time))
    end
    return k * time
end

"""
    compute_stvd_flow!(nodes, dof, nlist, diffusivity, front_depth, depth, bond_damage,
                       bond_length, horizon, concentration, volume, flow)

Bond-based diffusion [OterkusS2014b](@cite) with the diffusivity of a bond set to zero as
long as one of its nodes lies in the not yet reacted core.
"""
function compute_stvd_flow!(nodes::AbstractVector{Int64},
                            dof::Int64,
                            nlist,
                            diffusivity::Float64,
                            front_depth::Float64,
                            depth::AbstractVector{Float64},
                            bond_damage,
                            bond_length,
                            horizon::AbstractVector{Float64},
                            concentration::AbstractVector{Float64},
                            volume::AbstractVector{Float64},
                            flow::AbstractVector{Float64})
    for iID in nodes
        if depth[iID] > front_depth
            continue
        end
        kernel = dof == 2 ? 6.0 / (pi * horizon[iID]^3) : 6.0 / (pi * horizon[iID]^4)
        for (jID, neighborID) in enumerate(nlist[iID])
            if depth[neighborID] > front_depth || bond_damage[iID][jID] == 0
                continue
            end
            flow[iID] -= diffusivity * kernel * bond_damage[iID][jID] *
                         (concentration[neighborID] - concentration[iID]) *
                         volume[neighborID] / bond_length[iID][jID]
        end
    end
end

function fields_for_local_synchronization(model::String)
end

end
