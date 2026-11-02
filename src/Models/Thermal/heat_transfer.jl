# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Heat_Transfer
using StaticArrays: SVector, MVector, SMatrix
export compute_model
export init_model

using ......ParameterSpec: @params, register_thermal, ParseContext, add_error!, join_path
import ......ParameterSpec: check!
@params "Heat exchange of surface nodes with the environment (convection) and with the build plate" struct HeatTransferParams
    "Heat flow per surface area and temperature difference between surface and environment"
    heat_transfer_coefficient::Float64 = req("Heat Transfer Coefficient"; min = 0)
    "Temperature of the surrounding medium: a number, or an expression in the time t"
    environmental_temperature::Union{Float64,String} = req("Environmental Temperature";
                                                           evaluated = true)
    "Recomputes the surface in every step, so cracks and added material exchange heat too"
    allow_surface_change::Bool = opt("Allow Surface Change"; default = true)
    "Temperature of the build plate: a number, or an expression in the time t; without it, faces on the plate exchange heat with the environment"
    bed_temperature::Union{Nothing,Float64,String} = opt("Bed Temperature";
                                                         default = nothing,
                                                         evaluated = true)
    "Heat flow per contact area and temperature difference between node and build plate; Heat Transfer Coefficient if not given"
    bed_heat_transfer_coefficient::Union{Nothing,
                                         Float64} = opt("Bed Heat Transfer Coefficient";
                                                        default = nothing, min = 0)
    "Thermal conductivity of the contact to the build plate; the heat is conducted across half the node height, i.e. Bed Heat Transfer Coefficient = 2 Bed Thermal Conductivity / height"
    bed_thermal_conductivity::Union{Nothing,Float64} = opt("Bed Thermal Conductivity";
                                                           default = nothing, min = 0)
    "Position of the build plate along the last coordinate axis (z in 3D, y in 2D)"
    bed_position::Float64 = opt("Bed Position"; default = 0.0)
end
__init__() = register_thermal("Heat Transfer", HeatTransferParams)

function check!(p::HeatTransferParams, path::String, ctx::ParseContext)
    if p.bed_heat_transfer_coefficient !== nothing && p.bed_thermal_conductivity !== nothing
        add_error!(ctx, join_path(path, "Bed Thermal Conductivity"),
                   "give either Bed Heat Transfer Coefficient or Bed Thermal Conductivity")
    end
    if p.bed_temperature === nothing && (p.bed_heat_transfer_coefficient !== nothing ||
        p.bed_thermal_conductivity !== nothing)
        add_error!(ctx, join_path(path, "Bed Temperature"),
                   "required for the heat exchange with the build plate")
    end
    return nothing
end

"Environmental temperature: a number, or an expression in the time `t` (evaluated at `time`)."
environmental_temperature(T::Float64, time::Float64) = T
function environmental_temperature(expression::String, time::Float64)
    global t = time
    return Float64(eval(Meta.parse(expression)))
end

using TimerOutputs: @timeit
using .....Data_Manager
using .....Helpers: node_tensor

# Each node is a box with the edge lengths "Extent" (from the mesh, e.g. sampling x width x
# height of a G-code strand in its own frame; else a cube of the node volume), oriented by
# the rotation tensor. Face f = 2k-1 is the +k face, f = 2k the -k face (local axes).
# A neighbor covers part of the face it lies against; the overlap of its box with the face
# is computed once from the reference geometry ("Bond Face", "Bond Face Area"). In every
# step the faces are covered by the active, intact bonds only, so added material and
# cracks change the surface. "Specific Volume" is the exposed area / node volume.

"""
    init_model(nodes, p, thermal, block)

Inits the heat transfer model: finds for every bond the face of the node the neighbor
covers and the covered area, and the face of the node lying on the build plate.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `p`: The model parameters.
- `thermal::WithBase`: The typed thermal model of the block.
- `block::Int64`: The current block.

"""
function init_model(nodes::AbstractVector{Int64}, p::HeatTransferParams, thermal,
                    block::Int64)
    dof = Data_Manager.get_dof()
    volume::NodeScalarField{Float64} = Data_Manager.get_field("Volume")
    has_extent = Data_Manager.has_key("Extent")
    extent::NodeVectorField{Float64} = Data_Manager.create_constant_node_vector_field("Extent",
                                                                                      Float64,
                                                                                      dof)
    if !has_extent
        for iID in eachindex(volume)
            extent[iID, :] .= volume[iID]^(1 / dof)
        end
    end
    rotation_tensor = Data_Manager.get_rotation() ?
                      Data_Manager.get_field("Rotation Tensor")::Array{Float64,3} : nothing
    bond_face = Data_Manager.create_constant_bond_scalar_state("Bond Face", Int64)
    bond_face_area = Data_Manager.create_constant_bond_scalar_state("Bond Face Area",
                                                                    Float64)
    bed_face = Data_Manager.create_constant_node_scalar_field("Bed Face", Int64)
    Data_Manager.create_constant_node_scalar_field("Bed Specific Volume", Float64)
    Data_Manager.create_constant_node_scalar_field("Surface Active", Bool)
    Data_Manager.create_constant_node_scalar_field("Surface Update", Bool)
    Data_Manager.create_constant_node_scalar_field("Surface Update Nodes", Int64)
    # nodes with neighbors on other cores; their surface changes with activated ghost nodes
    border = Data_Manager.create_constant_node_scalar_field("Surface Border", Bool)
    nlist = Data_Manager.get_nlist()
    nnodes = Data_Manager.get_nnodes()
    for iID in nodes
        border[iID] = any(>(nnodes), nlist[iID])
    end
    init_bond_faces!(bond_face, bond_face_area, nodes, Data_Manager.get_nlist(),
                     Data_Manager.get_field("Bond Geometry"), extent, rotation_tensor,
                     Val(dof))
    if p.bed_temperature !== nothing
        init_bed_faces!(bed_face, nodes, Data_Manager.get_field("Coordinates"), extent,
                        rotation_tensor, p.bed_position, Val(dof))
    end
    if !p.allow_surface_change
        compute_specific_volume!(nodes, false)
    end
end

"""
    compute_model(nodes, p, thermal, block, time, dt)

Calculates the heat transfer to the environment and to the build plate. [BrighentiR2021](@cite)

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `p::HeatTransferParams`: The model parameters.
- `thermal::WithBase`: The typed thermal model of the block.
- `block::Int64`: The current block.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
"""
function compute_model(nodes::AbstractVector{Int64}, p::HeatTransferParams, thermal,
                       block::Int64, time::Float64, dt::Float64)
    Tenv::Float64 = environmental_temperature(p.environmental_temperature, time)
    kappa_bed = something(p.bed_heat_transfer_coefficient, p.heat_transfer_coefficient)
    Tbed::Float64 = p.bed_temperature === nothing ? Tenv :
                    environmental_temperature(p.bed_temperature, time)
    if p.allow_surface_change
        # without damage a node stays covered, so it is not checked again
        freeze = !haskey(Data_Manager.get_active_models(), "Damage Model")
        @timeit "compute_specific_volume!" compute_specific_volume!(nodes, freeze)
    end
    @timeit "apply_surface_heat_flow!" begin
        apply_surface_heat_flow!(Data_Manager.get_field("Heat Flow", "NP1"),
                                 Data_Manager.get_field("Temperature", "NP1"),
                                 Data_Manager.get_field("Specific Volume"),
                                 Data_Manager.get_field("Bed Specific Volume"),
                                 Data_Manager.get_field("Surface_Nodes"), nodes,
                                 p.heat_transfer_coefficient, Tenv, kappa_bed, Tbed,
                                 Data_Manager.get_field("Density"),
                                 Data_Manager.get_field("Specific Heat Capacity"), dt,
                                 p.bed_thermal_conductivity,
                                 Data_Manager.get_field("Extent"),
                                 Data_Manager.get_field("Bed Face"))
    end
end

function compute_specific_volume!(nodes::AbstractVector{Int64}, freeze::Bool)
    extent::NodeVectorField{Float64} = Data_Manager.get_field("Extent")
    # "Active" is created after the model init if no additive model defines it; all nodes are active then
    active = Data_Manager.has_key("Active") ? Data_Manager.get_field("Active") :
             fill(true, size(extent, 1))
    if freeze
        # without damage only activated nodes change the surface of themselves and their neighbors
        nodes = surface_update_nodes!(Data_Manager.get_field("Surface Update Nodes"),
                                      Data_Manager.get_field("Surface Update"),
                                      Data_Manager.get_field("Surface Active"),
                                      Data_Manager.get_field("Surface Border"), active,
                                      nodes, Data_Manager.get_nlist(),
                                      Data_Manager.get_nnodes())
    end
    compute_specific_volume!(Data_Manager.get_field("Specific Volume"),
                             Data_Manager.get_field("Bed Specific Volume"),
                             Data_Manager.get_field("Surface_Nodes"), nodes,
                             Data_Manager.get_nlist(), active,
                             Data_Manager.get_bond_damage("NP1"),
                             Data_Manager.get_field("Bond Face"),
                             Data_Manager.get_field("Bond Face Area"), extent,
                             Data_Manager.get_field("Bed Face"), freeze,
                             Val(size(extent, 2)))
end

"""
    surface_update_nodes!(buffer, update, surface_active, border, active, nodes, nlist, nnodes)

The nodes of `nodes` whose surface has to be computed: the nodes activated since the last
call and their neighbors. `surface_active` remembers the active nodes, `update` keeps the
flags of nodes not in `nodes` until they are. An activated ghost node (index > `nnodes`)
updates all `border` nodes, as the neighbor lists of ghost nodes are not complete.
"""
function surface_update_nodes!(buffer::AbstractVector{Int64},
                               update::AbstractVector{Bool},
                               surface_active::AbstractVector{Bool},
                               border::AbstractVector{Bool},
                               active::AbstractVector{Bool},
                               nodes::AbstractVector{Int64}, nlist, nnodes::Int64)
    ghost_activated = false
    @inbounds for iID in eachindex(active)
        (active[iID] && !surface_active[iID]) || continue
        surface_active[iID] = true
        update[iID] = true
        if iID > nnodes
            ghost_activated = true
            continue
        end
        for neighborID in nlist[iID]
            update[neighborID] = true
        end
    end
    if ghost_activated
        @inbounds for iID in nodes
            border[iID] && (update[iID] = true)
        end
    end
    count = 0
    @inbounds for iID in nodes
        update[iID] || continue
        update[iID] = false
        count += 1
        buffer[count] = iID
    end
    return view(buffer, 1:count)
end

"""
    apply_surface_heat_flow!(heat_flow, temperature, specific_volume, bed_specific_volume, surface_nodes, nodes, kappa, Tenv, kappa_bed, Tbed, density, heat_capacity, dt, bed_conductivity, extent, bed_face)

Adds the convection `kappa (T - Tenv) A/V` and the build plate exchange `kappa_bed (T - Tbed) A_bed/V`.
With a `bed_conductivity`, `kappa_bed = 2 bed_conductivity / h` with the node height `h`
normal to the plate.

Both are a linear relaxation `q = a (T - Tw)` of the node towards a mean wall temperature
`Tw`. The explicit update `dT = -q dt / (rho c)` is unstable for `a dt / (rho c) > 2`, so the
flow is scaled by `(1 - exp(-x)) / x`, `x = a dt / (rho c)`: the explicit update then gives
the exact solution of the relaxation over the step, for any time step.
"""
function apply_surface_heat_flow!(heat_flow::AbstractVector{Float64},
                                  temperature::AbstractVector{Float64},
                                  specific_volume::AbstractVector{Float64},
                                  bed_specific_volume::AbstractVector{Float64},
                                  surface_nodes::AbstractVector{Bool},
                                  nodes::AbstractVector{Int64},
                                  kappa::Float64, Tenv::Float64, kappa_bed::Float64,
                                  Tbed::Float64,
                                  density::AbstractVector{Float64},
                                  heat_capacity::AbstractVector{Float64}, dt::Float64,
                                  bed_conductivity::Union{Nothing,Float64} = nothing,
                                  extent::AbstractMatrix{Float64} = zeros(0, 0),
                                  bed_face::AbstractVector{Int64} = Int64[])
    @inbounds for iID in nodes
        surface_nodes[iID] || continue
        T = temperature[iID]
        a_env = kappa * specific_volume[iID]
        a_bed = bed_specific_volume[iID] > 0 ?
                bed_coefficient(kappa_bed, bed_conductivity, extent, bed_face, iID) *
                bed_specific_volume[iID] : 0.0
        x = (a_env + a_bed) * dt / (density[iID] * heat_capacity[iID])
        heat_flow[iID] += relaxation_factor(x) *
                          (a_env * (T - Tenv) + a_bed * (T - Tbed))
    end
    return nothing
end

"`(1 - exp(-x)) / x`, 1 for `x = 0`."
relaxation_factor(x::Float64) = x > 0 ? -expm1(-x) / x : 1.0

bed_coefficient(kappa_bed::Float64, ::Nothing, extent, bed_face, iID::Int64) = kappa_bed
function bed_coefficient(::Float64, bed_conductivity::Float64, extent, bed_face,
                         iID::Int64)
    return 2 * bed_conductivity / extent[iID, (bed_face[iID] + 1) >> 1]
end

"""
    compute_specific_volume!(specific_volume, bed_specific_volume, surface_nodes, nodes, nlist, active, bond_damage, bond_face, bond_face_area, extent, bed_face, freeze, ::Val{D})

Exposed area per node volume: sum over the faces of the uncovered face area / volume,
split into the faces exchanging heat with the environment and the face on the build
plate. A face is covered by the active neighbors with intact bonds lying against it.
With `freeze`, a node without exposed area is marked as no surface node and skipped.
"""
function compute_specific_volume!(specific_volume::AbstractVector{Float64},
                                  bed_specific_volume::AbstractVector{Float64},
                                  surface_nodes::AbstractVector{Bool},
                                  nodes::AbstractVector{Int64},
                                  nlist,
                                  active::AbstractVector{Bool},
                                  bond_damage,
                                  bond_face,
                                  bond_face_area,
                                  extent::AbstractMatrix{Float64},
                                  bed_face::AbstractVector{Int64},
                                  freeze::Bool,
                                  ::Val{D}) where {D}
    @inbounds for iID in nodes
        surface_nodes[iID] || continue
        covered = zero(MVector{2D,Float64})
        faces = bond_face[iID]
        areas = bond_face_area[iID]
        damage = bond_damage[iID]
        for (jID, neighborID) in enumerate(nlist[iID])
            f = faces[jID]
            (f == 0 || !active[neighborID] || damage[jID] == 0) && continue
            covered[f] += areas[jID]
        end
        h = SVector{D,Float64}(ntuple(k -> extent[iID, k], Val(D)))
        box = prod(h)
        free = 0.0
        bed = 0.0
        for f in 1:(2D)
            k = (f + 1) >> 1
            face_area = box / h[k]
            exposed = max(face_area - covered[f], 0.0) / box
            if f == bed_face[iID]
                bed += exposed
            else
                free += exposed
            end
        end
        specific_volume[iID] = free
        bed_specific_volume[iID] = bed
        if freeze && free + bed == 0
            surface_nodes[iID] = false
        end
    end
    return nothing
end

"""
    init_bond_faces!(bond_face, bond_face_area, nodes, nlist, bond_geometry, extent, rotation_tensor, ::Val{D})

For every bond: the face of node i the neighbor box lies against (0 if none) and the area
of the face it covers. Both boxes are compared in the frame of node i, the neighbor box by
its bounding box in this frame. The face is the axis along which the neighbor is farthest
relative to the touching distance; it lies against the face if this ratio is in
[0.5, 1.5], so gaps and overlaps of up to half the box sizes count as contact.
"""
function init_bond_faces!(bond_face, bond_face_area, nodes::AbstractVector{Int64}, nlist,
                          bond_geometry, extent::AbstractMatrix{Float64},
                          rotation_tensor::Union{Nothing,Array{Float64,3}},
                          ::Val{D}) where {D}
    for iID in nodes
        Ri = node_rotation(rotation_tensor, iID, Val(D))
        hi = SVector{D,Float64}(ntuple(k -> extent[iID, k], Val(D))) / 2
        for (jID, neighborID) in enumerate(nlist[iID])
            Rj = node_rotation(rotation_tensor, neighborID, Val(D))
            hj = abs.(Ri' * Rj) *
                 SVector{D,Float64}(ntuple(k -> extent[neighborID, k], Val(D))) / 2
            d = Ri' * SVector{D,Float64}(bond_geometry[iID][jID])
            ratio = abs.(d) ./ (hi + hj)
            k = argmax(ratio)
            area = 1.0
            for m in 1:D
                m == k && continue
                area *= max(min(hi[m], d[m] + hj[m]) - max(-hi[m], d[m] - hj[m]), 0.0)
            end
            if 0.5 <= ratio[k] <= 1.5 && area > 0
                bond_face[iID][jID] = d[k] > 0 ? 2k - 1 : 2k
                bond_face_area[iID][jID] = area
            else
                bond_face[iID][jID] = 0
                bond_face_area[iID][jID] = 0.0
            end
        end
    end
end

"""
    init_bed_faces!(bed_face, nodes, coordinates, extent, rotation_tensor, bed_position, ::Val{D})

Marks the downward face (normal along -e_D) of the nodes whose box reaches the build plate
at `bed_position` within half the box height, so nodes at the strand center and at the
nozzle height (G-code: top of the strand) both lie on the plate.
"""
function init_bed_faces!(bed_face, nodes::AbstractVector{Int64},
                         coordinates::AbstractMatrix{Float64},
                         extent::AbstractMatrix{Float64},
                         rotation_tensor::Union{Nothing,Array{Float64,3}},
                         bed_position::Float64, ::Val{D}) where {D}
    for iID in nodes
        R = node_rotation(rotation_tensor, iID, Val(D))
        h = SVector{D,Float64}(ntuple(k -> extent[iID, k], Val(D)))
        height = (abs.(R) * h)[D]
        bed_face[iID] = 0
        coordinates[iID, D] - height / 2 <= bed_position + height / 2 || continue
        # local axis k, whose -k face points down the most
        k = argmax(abs.(R[D, :]))
        abs(R[D, k]) < 0.9 && continue
        bed_face[iID] = R[D, k] > 0 ? 2k : 2k - 1
    end
end

node_rotation(::Nothing, iID::Int64, ::Val{D}) where {D} = one(SMatrix{D,D,Float64})
function node_rotation(rotation_tensor::Array{Float64,3}, iID::Int64,
                       ::Val{D}) where {D}
    return node_tensor(rotation_tensor, iID, Val(D))
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
