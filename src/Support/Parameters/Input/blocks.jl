# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

@params struct BlockParams
    block_id::Int64 = req("Block ID"; min = 1)
    density::Float64 = req("Density"; min = 0, quantity = :density)
    horizon::Float64 = req("Horizon"; min = 0, quantity = :length)
    specific_heat_capacity::Union{Nothing,Float64} = opt("Specific Heat Capacity";
                                                         default = nothing, min = 0)
    material_model::Union{Nothing,String} = opt("Material Model"; default = nothing)
    damage_model::Union{Nothing,String} = opt("Damage Model"; default = nothing)
    thermal_model::Union{Nothing,String} = opt("Thermal Model"; default = nothing)
    additive_model::Union{Nothing,String} = opt("Additive Model"; default = nothing)
    pre_calculation_model::Union{Nothing,String} = opt("Pre Calculation Model";
                                                       default = nothing)
    degradation_model::Union{Nothing,String} = opt("Degradation Model"; default = nothing)
    angle_x::Union{Nothing,Float64} = opt("Angle X"; default = nothing)
    angle_y::Union{Nothing,Float64} = opt("Angle Y"; default = nothing)
    angle_z::Union{Nothing,Float64} = opt("Angle Z"; default = nothing)
    step_id::Union{Nothing,Int64,String} = opt("Step ID"; default = nothing)
    fem::Union{Nothing,Bool} = opt("FEM"; default = nothing)
end

@params struct FEMCouplingParams
    coupling_type::String = req("Coupling Type")
    pd_weight::Union{Nothing,Float64} = opt("PD Weight"; default = nothing)
    kappa::Union{Nothing,Float64} = opt("Kappa"; default = nothing)
    coupling_block::Union{Nothing,Int64} = opt("Coupling Block"; default = nothing)
end

@params struct FEMParams
    element_type::String = req("Element Type")
    degree::Union{Int64,String} = req("Degree")
    material_model::String = req("Material Model")
    coupling::Union{Nothing,FEMCouplingParams} = opt("Coupling"; default = nothing)
end

@params struct SurfaceCorrectionParams
    type::String = req("Type")
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
