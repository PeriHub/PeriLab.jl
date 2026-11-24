# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Damage

using TimerOutputs: @timeit
using LinearAlgebra: dot
using Polyester: @batch
using ....Data_Manager
using ....PeriLabExceptions: @abort
using ....ModuleLoader: find_registered_modules
using .....ParameterSpec: @params, Dependent, Constant, register_base!, ParseContext,
                          add_error!, join_path, WithBase, Table1D, dependent_tables, value,
                          model_module, model_function
import .....ParameterSpec: check!

@params "Direction-dependent critical values for Critical Energy Anisotropic, along the material axes" struct AnisotropicDamageParams
    "Critical bond energy for the bond share along the material x axis"
    critical_value_x::Float64 = req("Critical Value X")
    "Critical bond energy for the bond share along the material y axis"
    critical_value_y::Float64 = req("Critical Value Y")
    "Critical bond energy for the bond share along the material z axis (3D); the y value if not given"
    critical_value_z::Union{Nothing,Float64} = opt("Critical Value Z"; default = nothing)
    "How a bond fails: Axis (each material axis fails on its own, the bond keeps the share of the intact axes), Axis Force (Axis, and the bond fails completely if a failed axis has a weight >= 0.5), Directional (the bond fails as a whole at the critical value interpolated for its direction), Axis Stretch and Directional Stretch (the same with critical stretches as critical values)"
    bond_criterion::String = opt("Bond Criterion"; default = "Axis",
                                 allowed = ["Axis", "Axis Force", "Directional",
                                     "Axis Stretch", "Directional Stretch"])
    "Axis Force only: a failed axis removes the bond force component along it (correspondence models)"
    force_filter::Bool = opt("Force Filter"; default = true)
    "Exponent n of the directional interpolation of the critical values (Directional, Directional Stretch)"
    exponent::Float64 = opt("Exponent"; default = 1.0, min = 0)
end

@params "Damps the waves released when bonds break (local damping after Silling)" struct LocalDampingParams
    "Young's modulus that sets the wave speed of the damping"
    representative_youngs_modulus::Float64 = req("Representative Young's modulus"; min = 0,
                                                 quantity = :stress)
    "Dimensionless strength of the damping; the stable range shrinks with the time step (checked in the first step); 0 switches it off"
    damping_coefficient::Float64 = req("Damping coefficient"; min = 0)
end

"""
    DamageBaseParams

Keys every damage model may use (read from the same YAML block as the model's
own keys).
"""
@params "Damage model of the blocks that use this entry: the failure limit of a bond and options shared by all damage models" struct DamageBaseParams
    "Failure limit of a bond: the critical stretch or the critical energy release rate, depending on the model; a number or a table over a node field"
    critical_value::Dependent = req("Critical Value"; min = 0)
    "Critical values of bonds between two blocks, as \"Interblock Critical Value i_j\" entries for bonds from block i to block j"
    interblock_damage::Union{Nothing,Dict{String,Float64}} = opt("Interblock Damage";
                                                                 default = nothing)
    anisotropic_damage::Union{Nothing,AnisotropicDamageParams} = opt("Anisotropic Damage";
                                                                     default = nothing)
    local_damping::Union{Nothing,LocalDampingParams} = opt("Local Damping";
                                                           default = nothing)
    "Scales the critical value of every bond by its Bonding Degree (Critical Stretch: B, Energy Release: B²); needs the Interlayer Bonding thermal model in an earlier step"
    interlayer_bonding::Bool = opt("Interlayer Bonding"; default = false)
end

const INTERBLOCK_KEY = r"^Interblock Critical Value \d+_\d+$"

function check!(p::DamageBaseParams, path::String, ctx::ParseContext)
    p.interblock_damage === nothing && return nothing
    for name in keys(p.interblock_damage)
        occursin(INTERBLOCK_KEY, name) ||
            add_error!(ctx, join_path(join_path(path, "Interblock Damage"), name),
                       "unknown key — expected \"Interblock Critical Value <block>_<block>\"")
    end
    p.critical_value isa Constant ||
        add_error!(ctx, join_path(path, "Critical Value"),
                   "must be a number when Interblock Damage is used")
    return nothing
end

# registration runs at load time, never during precompilation
__init__() = register_base!(:damage, DamageBaseParams)

"""
    BlockDamage

The typed damage model of a block: the shared base part, the model's own part,
and the field-dependent tables to bind before each compute.
"""
# mutable only for heap identity (never mutated), see BlockMaterial
mutable struct BlockDamage{B,M}
    base::B
    model::M
    tables::Vector{Table1D}
end

function block_damage(wb::WithBase)
    BlockDamage(wb.base, wb.model,
                vcat(dependent_tables(wb.base),
                     dependent_tables(wb.model)))
end

function init_interface_crit_values(damage::BlockDamage, block_id::Int64)
    interblock = damage.base.interblock_damage
    interblock === nothing && return
    max_block_id = maximum(Data_Manager.get_block_id_list())
    inter_critical_value = Data_Manager.get_crit_values_matrix()
    if inter_critical_value == fill(-1, (1, 1, 1))
        inter_critical_value = fill(value(damage.base.critical_value, 1),
                                    (max_block_id, max_block_id, max_block_id))
    end
    for block_iId in 1:max_block_id, block_jId in 1:max_block_id
        name = "Interblock Critical Value $(block_iId)_$block_jId"
        haskey(interblock, name) &&
            (inter_critical_value[block_iId, block_jId, block_id] = interblock[name])
    end
    Data_Manager.set_crit_values_matrix(inter_critical_value)
end

function init_aniso_crit_values(aniso::AnisotropicDamageParams, block_id::Int64,
                                dof::Int64)
    aniso_crit::Dict{Int64,Any} = Data_Manager.get_aniso_crit_values()
    aniso_crit[block_id] = dof == 2 ?
                           [aniso.critical_value_x, aniso.critical_value_y] :
                           [aniso.critical_value_x, aniso.critical_value_y,
        something(aniso.critical_value_z, aniso.critical_value_y)]
    Data_Manager.set_aniso_crit_values(aniso_crit)
end

"""
    init_inverse_bond_id()

Creates (once) the bond field "Inverse Bond ID", the bond-indexed inverse neighborhood
list (see `find_inverse_bond_index!`), and returns it.
"""
function init_inverse_bond_id()
    Data_Manager.has_key("Inverse Bond ID") &&
        return Data_Manager.get_field("Inverse Bond ID")
    inverse_bond_id = Data_Manager.create_constant_bond_scalar_state("Inverse Bond ID",
                                                                     Int64)
    return find_inverse_bond_index!(inverse_bond_id, Data_Manager.get_nlist(),
                                    Data_Manager.get_inverse_nlist())
end

"""
    critical_value_lookup(critical_value, inter_block_damage::Bool, block::Int64;
                          bonding_exponent::Int = 0)

Returns `(iID, jID, neighborID) -> critical value` of bond `jID` of node `iID`: the
`Critical_Value` node field if it exists, else the interblock matrix, else the
`Critical Value` parameter. With `bonding_exponent > 0` the value is multiplied by
`Bonding Degree[iID][jID]^bonding_exponent` (the bond field of the Interlayer Bonding
thermal model). Each case is its own closure type, so a bond loop taking it as an
argument is type stable.
"""
function critical_value_lookup(critical_value, inter_block_damage::Bool, block::Int64;
                               bonding_exponent::Int = 0)
    return with_critical_value(identity, critical_value, inter_block_damage, block;
                               bonding_exponent)
end

"""
    with_critical_value(f, critical_value, inter_block_damage::Bool, block::Int64;
                        bonding_exponent::Int = 0)

`f(crit_of)` with the lookup of `critical_value_lookup`. Every case calls `f` with its
own closure type, so the lookup is not returned (a returned closure of one of several
types is boxed in every call); the damage models run their bond loop in `f`.
"""
function with_critical_value(f::F, critical_value, inter_block_damage::Bool, block::Int64;
                             bonding_exponent::Int = 0) where {F}
    if Data_Manager.has_key("Critical_Value")
        crit_field::NodeScalarField{Float64} = Data_Manager.get_field("Critical_Value")
        return with_bonding(f, (iID, jID, neighborID) -> crit_field[iID],
                            bonding_exponent)
    elseif inter_block_damage
        block_ids::NodeScalarField{Int64} = Data_Manager.get_field("Block_Id")
        crit_matrix::Array{Float64,3} = Data_Manager.get_crit_values_matrix()
        return with_bonding(f,
                            (iID, jID, neighborID) -> crit_matrix[block_ids[iID],
                                                                  block_ids[neighborID],
                                                                  block],
                            bonding_exponent)
    end
    return with_bonding(f, (iID, jID, neighborID) -> value(critical_value, iID),
                        bonding_exponent)
end

# the critical value scaled by `Bonding Degree^bonding_exponent` (Interlayer Bonding)
function with_bonding(f::F, base::B, bonding_exponent::Int) where {F,B}
    bonding_exponent == 0 && return f(base)
    if !Data_Manager.has_key("Bonding Degree")
        @abort "Interlayer Bonding needs the bond field Bonding Degree: add the Interlayer Bonding thermal model to an earlier solver step"
    end
    bonding::BondScalarState{Float64} = Data_Manager.get_field("Bonding Degree")
    return f((iID, jID, neighborID) -> base(iID, jID, neighborID) *
                                       bonding[iID][jID]^bonding_exponent)
end

"Values of the `Compression` option: how a shortened bond counts (full, skipped, only its transverse part)."
const COMPRESSION_MODES = ["full", "none", "shear"]

"Part of `v` normal to the unit vector `n`."
transverse(v, n) = v - dot(v, n) * n

for file in find_registered_modules(@__DIR__, "register_damage")
    include(file)
end

using .....Helpers: find_inverse_bond_id, find_inverse_bond_index!, MINBATCH
export fields_for_local_synchronization
export compute_model
export init_interface_crit_values
export init_model
export init_fields

"""
    init_fields()

Creates the damage field and the inverse neighbor list.
"""
function init_fields()
    dof = Data_Manager.get_dof()
    Data_Manager.create_node_scalar_field("Damage", Float64)
    # the damage index only grows and is recomputed for the flagged nodes only
    Data_Manager.set_np1_start("Damage", :copy)
    # energy released by bond failure (critical energy models); created here, because the
    # damage models of later solver steps are initialized after the output
    Data_Manager.create_constant_node_scalar_field("Dissipated Energy", Float64)

    anisotropic_damage = false

    nlist = Data_Manager.get_nlist()
    inverse_nlist = Data_Manager.set_inverse_nlist(find_inverse_bond_id(nlist))
    Data_Manager.track_bond_damage_N()
end

"""
    compute_model(nodes, damage, block, time, dt)

Binds the damage's dependent values, computes the block's damage model and the
damage index.

# Arguments
- `nodes::AbstractVector{Int64}`: The nodes
- `damage::BlockDamage`: The typed damage model of the block
- `block::Int64`: The block
- `time::Float64`: The current time
- `dt::Float64`: The time step
"""
function compute_model(nodes::AbstractVector{Int64}, damage::BlockDamage, block::Int64,
                       time::Float64, dt::Float64)
    Data_Manager.bind_dependent_tables!(damage.tables)
    model_function(damage.model, Val(:compute_model))(nodes, damage.model, damage, block,
                                                      time, dt)
    # the first call of a block covers all nodes (bonds broken at the start, e.g. by a bond
    # filter), later ones only the nodes whose bonds the damage model changed
    all_nodes = !(block in INDEXED_BLOCKS)
    push!(INDEXED_BLOCKS, block)
    if isnothing(Data_Manager.get_filtered_nlist())
        @timeit "compute index" return damage_index(nodes, all_nodes)
    end
    @timeit "compute index" return damage_index(nodes, all_nodes,
                                                Data_Manager.get_filtered_nlist())
end

# blocks whose damage index was computed for all nodes once; reset by `init_model`
const INDEXED_BLOCKS = Set{Int64}()

"""
    fields_for_local_synchronization(model, block)

Defines all synchronization fields for local synchronization

# Arguments
- `model::String`: Model class.
- `block::Int64`: block ID
"""
function fields_for_local_synchronization(model, block)
    damage = Data_Manager.get_block_models(block).damage
    model_module(damage.model).fields_for_local_synchronization(model)
end
"""
    damage_index(::Union{SubArray, Vector{Int64})

Function calculates the damage index related to the neighborhood volume for a set of corresponding nodes.
The damage index is defined as damaged volume in relation the neighborhood volume.
damageIndex = sum_i (brokenBonds_i * volume_i) / volumeNeighborhood

# Arguments
- `nodes::AbstractVector{Int64}`: corresponding nodes to this model
"""
function damage_index(nodes::AbstractVector{Int64}, all_nodes::Bool,
                      nlist_filtered_ids::NestedNeighborhoodList)
    # contact bonds of a bond filter stay broken; a damage model must not heal them
    bond_damageNP1 = Data_Manager.get_bond_damage("NP1")::BondScalarState{Float64}
    update = all_nodes ? Bool[] : Data_Manager.get_field("Update")::NodeScalarField{Bool}
    _keep_filtered_broken!(bond_damageNP1, nodes, all_nodes, update, nlist_filtered_ids)
    return damage_index(nodes, all_nodes)
end

function _keep_filtered_broken!(bond_damage::BondScalarState{Float64},
                                nodes::AbstractVector{Int64}, all_nodes::Bool,
                                update::NodeScalarField{Bool}, nlist_filtered_ids)
    for iID in nodes
        (all_nodes || update[iID]) || continue
        # element-wise: a view with an index vector of the bond view copies the indices
        damage = bond_damage[iID]
        for jID in nlist_filtered_ids[iID]
            damage[jID] = 0.0
        end
    end
end

function damage_index(nodes::AbstractVector{Int64}, all_nodes::Bool = true)
    nlist = Data_Manager.get_nlist()::BondScalarState{Int64}
    volume = Data_Manager.get_field("Volume")::NodeScalarField{Float64}
    bond_damageNP1 = Data_Manager.get_bond_damage("NP1")::BondScalarState{Float64}
    damage = Data_Manager.get_damage("NP1")::NodeScalarField{Float64}
    update = all_nodes ? Bool[] : Data_Manager.get_field("Update")::NodeScalarField{Bool}
    compute_index(damage, nodes, volume, nlist, bond_damageNP1, update, all_nodes)
end

function compute_index(damage::NodeScalarField{Float64},
                       nodes::AbstractVector{Int64},
                       volume::NodeScalarField{Float64},
                       nlist::BondScalarState{Int64},
                       bond_damage::BondScalarState{Float64},
                       update::AbstractVector{Bool} = Bool[],
                       all_nodes::Bool = true)::Nothing
    @batch minbatch=MINBATCH for iID in nodes
        # unflagged nodes keep their bonds, their damage NP1 is the copy of N
        (all_nodes || update[iID]) || continue
        undamaged_volume = 0.0  # More explicit than zero(Float64)
        totalDamage = 0.0

        # Cache the vectors to help type inference
        neighbors = nlist[iID]
        bonds = bond_damage[iID]

        @inbounds @fastmath for j in eachindex(neighbors)
            jID = neighbors[j]
            vol_j = volume[jID]
            undamaged_volume += vol_j
            totalDamage += (1.0 - bonds[j]) * vol_j
        end

        # More explicit type handling
        current_damage = damage[iID]
        threshold = @fastmath current_damage * undamaged_volume
        if threshold < totalDamage
            damage[iID] = @fastmath totalDamage / undamaged_volume
        end
    end
    return nothing
end

"""
    init_model(nodes::AbstractVector{Int64}, block::Int64)

Initializes the damage model of a block (`Data_Manager.get_block_models(block).damage`),
its interface and anisotropic critical values.

# Arguments
- `nodes::AbstractVector{Int64}`: Nodes of the block.
- `block::Int64`: Block identifier.
"""
function init_model(nodes::AbstractVector{Int64}, block::Int64)
    delete!(INDEXED_BLOCKS, block)
    damage = Data_Manager.get_block_models(block).damage
    model_module(damage.model).init_model(nodes, damage.model, damage, block)
    model_module(damage.model).fields_for_local_synchronization("Damage Model")
    init_interface_crit_values(damage, block)
    damage.base.anisotropic_damage === nothing ||
        init_aniso_crit_values(damage.base.anisotropic_damage, block,
                               Data_Manager.get_dof())
end
end
