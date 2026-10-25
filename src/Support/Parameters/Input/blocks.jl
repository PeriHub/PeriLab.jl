# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

@params "A block of nodes with its density, horizon and the models it uses" struct BlockParams
    "Block id of the nodes in the mesh file"
    block_id::Int64 = req("Block ID"; min = 1)
    "Mass density of the block"
    density::Float64 = req("Density"; min = 0, quantity = :density)
    "Radius of the neighborhood in which a node interacts with other nodes"
    horizon::Float64 = req("Horizon"; min = 0, quantity = :length)
    "Heat capacity per unit mass; required with thermal models"
    specific_heat_capacity::Union{Nothing,Float64} = opt("Specific Heat Capacity";
                                                         default = nothing, min = 0)
    "Name of an entry in Material Models"
    material_model::Union{Nothing,String} = opt("Material Model"; default = nothing)
    "Name of an entry in Damage Models; no damage if not given"
    damage_model::Union{Nothing,String} = opt("Damage Model"; default = nothing)
    "Name of an entry in Thermal Models"
    thermal_model::Union{Nothing,String} = opt("Thermal Model"; default = nothing)
    "Name of an entry in Additive Models"
    additive_model::Union{Nothing,String} = opt("Additive Model"; default = nothing)
    "Name of an entry in Pre Calculation Models"
    pre_calculation_model::Union{Nothing,String} = opt("Pre Calculation Model";
                                                       default = nothing)
    "Name of an entry in Degradation Models"
    degradation_model::Union{Nothing,String} = opt("Degradation Model"; default = nothing)
    "Rotation of the block's material axes about x, in degrees; the only angle in 2D (rotation in the x-y plane); overrides angles from the mesh"
    angle_x::Union{Nothing,Float64} = opt("Angle X"; default = nothing)
    "Rotation of the block's material axes about y, in degrees (3D, required there with Angle X)"
    angle_y::Union{Nothing,Float64} = opt("Angle Y"; default = nothing)
    "Rotation of the block's material axes about z, in degrees (3D, required there with Angle X)"
    angle_z::Union{Nothing,Float64} = opt("Angle Z"; default = nothing)
    "Lists the block as FEM in the block summary of the log; the FEM nodes themselves come from the external topology"
    fem::Union{Nothing,Bool} = opt("FEM"; default = nothing)
end

@params "Coupling between the finite element and the peridynamic region" struct FEMCouplingParams
    "Coupling method, e.g. \"Arlequin\""
    coupling_type::String = req("Coupling Type")
    "Arlequin weight in the coupling zone: forces of the coupled finite element nodes are scaled by this value, those of the coupled peridynamic nodes by one minus it"
    pd_weight::Float64 = opt("PD Weight"; default = 0.5)
    "Stiffness of the coupling between peridynamic nodes and finite elements"
    kappa::Float64 = opt("Kappa"; default = 1.0)
    "Block id of the coupling zone; without it every peridynamic node inside an element is coupled"
    coupling_block::Union{Nothing,Int64} = opt("Coupling Block"; default = nothing)
end

@params "Finite elements for the nodes of the external topology, optionally coupled to peridynamics" struct FEMParams
    "Element formulation, e.g. \"Lagrange\""
    element_type::String = req("Element Type")
    "Polynomial degree of the elements: one value, or one per direction separated by spaces"
    degree::Union{Int64,String} = req("Degree")
    "Name of the entry in Material Models used by the finite elements"
    material_model::String = req("Material Model")
    coupling::Union{Nothing,FEMCouplingParams} = opt("Coupling"; default = nothing)
end

function check!(f::FEMParams, path::String, ctx::ParseContext)
    degrees = f.degree isa Int64 ? [f.degree] :
              [tryparse(Int64, part) for part in split(f.degree)]
    if isempty(degrees) || any(d -> d === nothing || d < 1, degrees)
        add_error!(ctx, join_path(path, "Degree"),
                   "expected a positive integer or positive integers separated by spaces, got \"$(f.degree)\"")
    end
    return nothing
end

"Polynomial degrees of the FEM elements as written: one value, or one per direction."
fem_degree(f::FEMParams) = f.degree isa Int64 ? [f.degree] : parse.(Int64, split(f.degree))

@params "Corrects the missing neighbors of nodes near free surfaces" struct SurfaceCorrectionParams
    "Correction method"
    type::String = req("Type"; allowed = ["Volume Correction"])
    "Recomputes the correction in every step, so cracks and added material change the surface"
    update::Bool = opt("Update"; default = false)
end

"""
    block_by_id(blocks, block_id) -> (name, block)

The block with `Block ID` `block_id`.
"""
function block_by_id(blocks::Dict{String,BlockParams}, block_id::Integer)
    for (name, block) in blocks
        block.block_id == block_id && return name, block
    end
    @abort "Block with ID $block_id is not defined"
end

"""
    block_angles(name, block, dof)

Block rotation: `nothing` if `Angle X` is not given, `Angle X` in 2D, and
`[Angle X, Angle Y, Angle Z]` in 3D (all three are required then).
"""
function block_angles(name::AbstractString, block::BlockParams, dof::Int64)
    block.angle_x === nothing && return nothing
    dof == 2 && return block.angle_x
    dof == 3 || return nothing
    for (key, value) in (("Angle Y", block.angle_y), ("Angle Z", block.angle_z))
        value === nothing && @abort "$key of $name is not defined"
    end
    return [block.angle_x, block.angle_y, block.angle_z]
end

"""
    block_names_and_ids(blocks, mesh_block_ids, mpi) -> (names, ids)

Names and IDs of the blocks present in the mesh, ordered by ID. Blocks of the
input deck that are missing in the mesh are skipped (with a warning unless
running under MPI).
"""
function block_names_and_ids(blocks::Dict{String,BlockParams},
                             mesh_block_ids::AbstractVector{<:Integer}, mpi::Bool)
    names = String[]
    ids = Int64[]
    for id in 1:maximum(block.block_id for block in values(blocks))
        if !(id in mesh_block_ids)
            mpi || @warn "Block with ID $id is not defined in the provided mesh"
            continue
        end
        for (name, block) in blocks
            if block.block_id == id
                push!(names, name)
                push!(ids, id)
            end
        end
    end
    return names, ids
end
