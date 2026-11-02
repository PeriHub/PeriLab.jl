# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Pre_Calculation

using TimerOutputs: @timeit
using ....Data_Manager
using ....Helpers: compute_label
using ....ModuleLoader: find_registered_modules
using .....ParameterSpec
# true during the first model pass of `compute_models`: there only the nodes flagged "Update"
# can have changed bonds since their last computation, the damage models flag them and the
# update pass recomputes them. Not true in the other solver paths.
const ONLY_FLAGGED_BONDS_CHANGED = Ref(false)
for file in find_registered_modules(@__DIR__, "register_pre_calculation")
    include(file)
end

export init_fields
export init_model
export fields_for_local_synchronization
export compute_model
export check_dependencies

"""
    init_fields()

Initializes the fields.
"""
function init_fields()
    dof = Data_Manager.get_dof()
    Data_Manager.create_bond_vector_state("Deformed Bond Geometry", Float64, dof)
    Data_Manager.create_bond_scalar_state("Deformed Bond Length", Float64)
    Data_Manager.create_node_vector_field("Displacements", Float64, dof)
end

"The parameters of the registered pre-calculation `name` (a switch: they have no keys)."
pre_calculation(name::String) = ParameterSpec.lookup_model(:pre_calculation, name)()

"Active names of `switches` (on/off per pre-calculation), in run order."
function active_pre_calculations(switches::AbstractDict{String,Bool})
    order_pre_calculations([name
                            for (name, on) in switches
                            if on])
end

"`names` in run order: the fixed order first, then the others sorted."
function order_pre_calculations(names)
    order = Data_Manager.get_pre_calculation_order()
    return vcat([n for n in order if n in names],
                sort!([n for n in names if !(n in order)]))
end

"""
    init_model(nodes::AbstractVector{Int64}, block::Int64)

Initializes the block's active pre-calculations
(`Data_Manager.get_block_models(block).pre_calculation`).

# Arguments
- `nodes::AbstractVector{Int64}`: The nodes.
- `block::Int64`: Block.
"""
function init_model(nodes::AbstractVector{Int64}, block::Int64)
    for name in Data_Manager.get_block_models(block).pre_calculation
        p = pre_calculation(name)
        model_module(p).init_model(nodes, p, nothing, block)
    end
end

"""
    compute_model(nodes::AbstractVector{Int64}, names::Vector{String}, block::Int64, time::Float64, dt::Float64)

Computes the block's active pre-calculations in run order.

# Arguments
- `nodes::AbstractVector{Int64}`: The nodes
- `names::Vector{String}`: The active pre-calculations of the block
- `block::Int64`: The block
- `time::Float64`: The current time
- `dt::Float64`: The time step
"""
function compute_model(nodes::AbstractVector{Int64}, names::Vector{String}, block::Int64,
                       time::Float64, dt::Float64)
    return compute_each(nodes, typed_pre_calculations(names), block, time, dt)
end

# `name => parameters` of `names`, built once per list: the type of the tuple lets
# `compute_each` dispatch every pre-calculation statically.
# ponytail: keyed by the identity of the block's (never mutated) name list
const TYPED_PRE_CALCULATIONS = IdDict{Vector{String},Any}()
function typed_pre_calculations(names::Vector{String})
    return get!(() -> Tuple(name => pre_calculation(name) for name in names),
                TYPED_PRE_CALCULATIONS, names)
end

compute_each(nodes, ::Tuple{}, block, time, dt) = nothing
function compute_each(nodes, pre::Tuple, block, time, dt)
    name, p = first(pre)
    @timeit compute_label(name) model_function(p, Val(:compute_model))(nodes, p, nothing,
                                                                       block, time, dt)
    return compute_each(nodes, Base.tail(pre), block, time, dt)
end

"""
    fields_for_local_synchronization(model, block)

Defines all synchronization fields for local synchronization

# Arguments
- `model::String`: Model class.
- `block::Int64`: block ID
"""
function fields_for_local_synchronization(model, block)
    for name in Data_Manager.get_block_models(block).pre_calculation
        model_module(pre_calculation(name)).fields_for_local_synchronization(model)
    end
end

"""
    check_dependencies(block_nodes::Dict{Int64,Vector{Int64}})

Adds the pre-calculations the block's material needs, and the ones the active
pre-calculations depend on, and stores the block's list in run order.

# Arguments
- `block_nodes::Dict{Int64,Vector{Int64}}`: block nodes.
"""
function check_dependencies(block_nodes::Dict{Int64,Vector{Int64}})
    for block_id in eachindex(block_nodes)
        models = Data_Manager.get_block_models(block_id)
        material = models.material
        material === nothing && continue
        names = Set{String}(models.pre_calculation)
        push!(names, "Deformed Bond Geometry")
        if material.correspondence
            if material.base.bond_associated
                push!(names, "Bond Associated Correspondence")
            else
                push!(names, "Shape Tensor", "Deformation Gradient")
            end
        end
        "Deformation Gradient" in names && push!(names, "Shape Tensor")
        Data_Manager.set_block_models(block_id,
                                      BlockModels(models.material, models.damage,
                                                  models.thermal, models.additive,
                                                  models.degradation,
                                                  order_pre_calculations(collect(names))))
    end
end

end
