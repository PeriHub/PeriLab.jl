# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module STVD_Flow
using NearestNeighbors: KDTree, nn
using .....Data_Manager
using .....PeriLabExceptions: @abort
export compute_model
export init_model
export fields_for_local_synchronization

using ......ParameterSpec: @params, register_thermal, ParseContext, add_error!, join_path
import ......ParameterSpec: check!

@enum ReactionLaw Linear Exponential
@enum Release Node Continuous

@params "Spatio-temporal variable diffusivity: diffusion of a diffusant released by a reaction front" struct STVDFlowParams
    "Speed of the reaction front (Linear) or decay rate of the core (Exponential)"
    reaction_rate::Float64 = req("Reaction Rate"; min = 0)
    "Block ID(s) of the diffusant, several separated by spaces"
    diffusant_block_id::Union{Int64,String} = req("Diffusant Block ID")
    "Front position over time: Linear or Exponential"
    reaction_law::ReactionLaw = opt("Reaction Law"; default = Linear)
    "Release of the core: per Node or Continuous"
    release::Release = opt("Release"; default = Node)
    "Initial core size; required for the exponential reaction law"
    core_size::Union{Nothing,Float64} = opt("Core Size"; default = nothing, min = 0)
end

function check!(p::STVDFlowParams, path::String, ctx::ParseContext)
    if p.reaction_law == Exponential && p.core_size === nothing
        add_error!(ctx, join_path(path, "Core Size"),
                   "has to be defined for the exponential reaction law")
    end
    return nothing
end

__init__() = register_thermal("STVD Flow", STVDFlowParams)

"""
    init_model(nodes, p, thermal, block)

Inits the spatio-temporal variable diffusivity (STVD) model [AtiehM2021](@cite). The
temperature field is interpreted as concentration of the diffusant. Nodes of the
diffusant blocks form an inert core with zero diffusivity. The core is consumed by a
reaction front which advances from the diffusant surface with the speed `Reaction Rate`.
For each diffusant node the distance to the diffusant surface is stored in "STVD_Depth".
It is estimated from the nearest substrate node, unless the mesh file provides an
`STVD_Depth` column (negative for substrate nodes). An exact depth from the mesh is
recommended for touching particles, where no substrate node lies at the contact points.

The release of the core is defined by `Release`:
- `"Node"` (default): a node is released completely when the front passes it.
- `"Continuous"`: a node is released in proportion to the crossed part of its spacing w,
  with the inert fraction f = clamp((depth + w/2 - s) / w, 0, 1). The bonds transport the
  mobile concentration u - f. This avoids release pulses of one node volume, which disturb
  threshold-based quantities when the diffusion is fast compared to the front.

The front position is defined by the `Reaction Law`:
- `"Linear"` (default): front depth = `Reaction Rate` * t [AtiehM2022](@cite)
- `"Exponential"`: the core shrinks with exp(-`Reaction Rate` * t), i.e.
  front depth = `Core Size` * (1 - exp(-`Reaction Rate` * t)) [AtiehM2021](@cite)

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `p::STVDFlowParams`: The model parameters.
- `thermal::WithBase`: The typed thermal model of the block (`thermal.base`: Thermal Conductivity).
- `block::Int64`: The current block.
"""
function init_model(nodes::AbstractVector{Int64}, p::STVDFlowParams, thermal,
                    block::Int64)
    if Data_Manager.get_mpi_active()
        @abort "STVD Flow does not support MPI yet, the depth below the diffusant surface needs all substrate nodes."
    end
    if thermal.base.thermal_conductivity === nothing
        @abort "Thermal Conductivity not defined."
    end
    if !Data_Manager.has_key("STVD_Spacing")
        spacing = Data_Manager.create_constant_node_scalar_field("STVD_Spacing", Float64)
        if p.release == Continuous
            bond_length = Data_Manager.get_field("Bond Length")
            for iID in eachindex(spacing)
                spacing[iID] = isempty(bond_length[iID]) ? 0.0 : minimum(bond_length[iID])
            end
        end
    end
    if Data_Manager.has_key("STVD_Depth")
        return
    end
    diffusant_blocks = parse_block_ids(p.diffusant_block_id)
    depth = Data_Manager.create_constant_node_scalar_field("STVD_Depth", Float64)
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
    compute_model(nodes, p, thermal, block, time, dt)

Calculates the concentration flow of the STVD model. A bond transports the mobile
concentration only if neither node is completely inert.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `p::STVDFlowParams`: The model parameters.
- `thermal::WithBase`: The typed thermal model of the block.
- `block::Int64`: Current block
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
"""
function compute_model(nodes::AbstractVector{Int64}, p::STVDFlowParams, thermal,
                       block::Int64, time::Float64, dt::Float64)
    front_depth = compute_front_depth(p, time)
    compute_stvd_flow!(nodes,
                       Data_Manager.get_dof(),
                       Data_Manager.get_nlist(),
                       thermal.base.thermal_conductivity,
                       front_depth,
                       p.release == Continuous,
                       Data_Manager.get_field("STVD_Depth"),
                       Data_Manager.get_field("STVD_Spacing"),
                       Data_Manager.get_bond_damage("NP1"),
                       Data_Manager.get_field("Bond Length"),
                       Data_Manager.get_field("Horizon"),
                       Data_Manager.get_field("Temperature", "NP1"),
                       Data_Manager.get_field("Volume"),
                       Data_Manager.get_field("Heat Flow", "NP1"))
end

"""
    compute_front_depth(p, time)

Depth of the reaction front below the initial diffusant surface.
"""
function compute_front_depth(p::STVDFlowParams, time::Float64)
    k = p.reaction_rate
    if p.reaction_law == Exponential
        return p.core_size * (1 - exp(-k * time))
    end
    return k * time
end

"""
    inert_fraction(depth, spacing, front_depth, continuous)

Part of the node concentration that is still bound in the core. Node-wise release: 1 while
the front has not passed the node, 0 afterwards. Continuous release: the crossed part of the
node spacing is released.
"""
function inert_fraction(depth::Float64, spacing::Float64, front_depth::Float64,
                        continuous::Bool)
    if !continuous || spacing == 0
        return depth > front_depth ? 1.0 : 0.0
    end
    return clamp((depth + 0.5 * spacing - front_depth) / spacing, 0.0, 1.0)
end

"""
    compute_stvd_flow!(nodes, dof, nlist, diffusivity, front_depth, continuous, depth,
                       spacing, bond_damage, bond_length, horizon, concentration, volume,
                       flow)

Bond-based diffusion [OterkusS2014b](@cite) of the mobile concentration u - f, where f is
the inert fraction. A bond is switched off as long as one of its nodes is completely inert.
"""
function compute_stvd_flow!(nodes::AbstractVector{Int64},
                            dof::Int64,
                            nlist,
                            diffusivity::Float64,
                            front_depth::Float64,
                            continuous::Bool,
                            depth::AbstractVector{Float64},
                            spacing::AbstractVector{Float64},
                            bond_damage,
                            bond_length,
                            horizon::AbstractVector{Float64},
                            concentration::AbstractVector{Float64},
                            volume::AbstractVector{Float64},
                            flow::AbstractVector{Float64})
    for iID in nodes
        f_i = inert_fraction(depth[iID], spacing[iID], front_depth, continuous)
        if f_i == 1
            continue
        end
        mobile_i = concentration[iID] - f_i
        kernel = dof == 2 ? 6.0 / (pi * horizon[iID]^3) : 6.0 / (pi * horizon[iID]^4)
        for (jID, neighborID) in enumerate(nlist[iID])
            if bond_damage[iID][jID] == 0
                continue
            end
            f_j = inert_fraction(depth[neighborID], spacing[neighborID], front_depth,
                                 continuous)
            if f_j == 1
                continue
            end
            flow[iID] -= diffusivity * kernel * bond_damage[iID][jID] *
                         (concentration[neighborID] - f_j - mobile_i) *
                         volume[neighborID] / bond_length[iID][jID]
        end
    end
end

function fields_for_local_synchronization(model::String)
end

end
