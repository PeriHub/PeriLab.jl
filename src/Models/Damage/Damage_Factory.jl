# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Damage

using TimerOutputs: @timeit
using ....Data_Manager
using ....PeriLabExceptions: @abort
using ....ModuleLoader: find_registered_modules
using .....ParameterSpec: @params, Dependent, Constant, register_base!, ParseContext,
                          add_error!, join_path, WithBase, Table1D, dependent_tables, value,
                          model_module
import .....ParameterSpec: check!

@params struct AnisotropicDamageParams
    critical_value_x::Float64 = req("Critical Value X")
    critical_value_y::Float64 = req("Critical Value Y")
    critical_value_z::Union{Nothing,Float64} = opt("Critical Value Z"; default = nothing,
                                                   description = "defaults to Critical Value Y")
end

@params struct LocalDampingParams
    representative_youngs_modulus::Float64 = req("Representative Young's modulus"; min = 0,
                                                 quantity = :stress)
    damping_coefficient::Float64 = req("Damping coefficient"; min = 0)
end

"""
    DamageBaseParams

Keys every damage model may use (read from the same YAML block as the model's
own keys).
"""
@params struct DamageBaseParams
    critical_value::Dependent = req("Critical Value"; min = 0)
    interblock_damage::Union{Nothing,Dict{String,Float64}} = opt("Interblock Damage";
                                                                 default = nothing,
                                                                 description = "Interblock Critical Value <block>_<block> entries")
    anisotropic_damage::Union{Nothing,AnisotropicDamageParams} = opt("Anisotropic Damage";
                                                                     default = nothing)
    local_damping::Union{Nothing,LocalDampingParams} = opt("Local Damping"; default = nothing)
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
struct BlockDamage{B,M}
    base::B
    model::M
    tables::Vector{Table1D}
end

block_damage(wb::WithBase) = BlockDamage(wb.base, wb.model,
                                         vcat(dependent_tables(wb.base),
                                              dependent_tables(wb.model)))

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
for file in find_registered_modules(@__DIR__, "register_damage")
    include(file)
end

using LoopVectorization
using .....Helpers: find_inverse_bond_id
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

    anisotropic_damage = false

    nlist = Data_Manager.get_nlist()
    inverse_nlist = Data_Manager.set_inverse_nlist(find_inverse_bond_id(nlist))
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
    model_module(damage.model).compute_model(nodes, damage.model, damage, block, time, dt)
    if isnothing(Data_Manager.get_filtered_nlist())
        @timeit "compute index" return damage_index(nodes)
    end
    @timeit "compute index" return damage_index(nodes, Data_Manager.get_filtered_nlist())
end

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
function damage_index(nodes::AbstractVector{Int64},
                      nlist_filtered_ids::BondScalarState{Int64})
    bond_damageNP1 = Data_Manager.get_bond_damage("NP1")
    for iID in nodes
        bond_damageNP1[iID][nlist_filtered_ids[iID]] .= 1
    end
    return damage_index(nodes)
end

function damage_index(nodes::AbstractVector{Int64})
    nlist = Data_Manager.get_nlist()::BondScalarState{Int64}
    volume = Data_Manager.get_field("Volume")::NodeScalarField{Float64}
    bond_damageNP1 = Data_Manager.get_bond_damage("NP1")::BondScalarState{Float64}
    damage = Data_Manager.get_damage("NP1")::NodeScalarField{Float64}
    compute_index(damage, nodes, volume, nlist, bond_damageNP1)
end

function compute_index(damage::NodeScalarField{Float64},
                       nodes::AbstractVector{Int64},
                       volume::NodeScalarField{Float64},
                       nlist::BondScalarState{Int64},
                       bond_damage::BondScalarState{Float64})::Nothing
    @inbounds @fastmath for iID in nodes
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
        threshold = current_damage * undamaged_volume
        if threshold < totalDamage
            damage[iID] = totalDamage / undamaged_volume
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
    damage = Data_Manager.get_block_models(block).damage
    model_module(damage.model).init_model(nodes, damage.model, damage, block)
    model_module(damage.model).fields_for_local_synchronization("Damage Model")
    init_interface_crit_values(damage, block)
    damage.base.anisotropic_damage === nothing ||
        init_aniso_crit_values(damage.base.anisotropic_damage, block, Data_Manager.get_dof())
end
end
