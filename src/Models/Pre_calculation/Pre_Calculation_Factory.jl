# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Pre_Calculation

using TimerOutputs: @timeit
using ....Data_Manager
using ....ModuleLoader: find_module_files
using .....ParameterSpec
global module_list = find_module_files(@__DIR__, "pre_calculation_name")
for mod in module_list
    include(mod["File"])
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

"The module of the registered pre-calculation `name`."
pre_calculation_module(name::String) = parentmodule(ParameterSpec.lookup_model(:pre_calculation,
                                                                               name))

"Active names of `switches` (on/off per pre-calculation), in run order."
active_pre_calculations(switches::AbstractDict{String,Bool}) = order_pre_calculations([name
                                                                                       for (name, on) in switches
                                                                                       if on])

"`names` in run order: the fixed order first, then the others sorted."
function order_pre_calculations(names)
    order = Data_Manager.get_pre_calculation_order()
    return vcat([n for n in order if n in names], sort!([n for n in names if !(n in order)]))
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
        pre_calculation_module(name).init_model(nodes, block)
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
    for name in names
        @timeit "compute $name" pre_calculation_module(name).compute(nodes, block)
    end
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
        pre_calculation_module(name).fields_for_local_synchronization(model)
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
