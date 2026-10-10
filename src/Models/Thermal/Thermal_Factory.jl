# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Thermal

using ....Data_Manager
using ....PeriLabExceptions: @abort
using TimerOutputs: @timeit
using ....ModuleLoader: find_registered_modules

using .....ParameterSpec: @params, register_base!, WithBase, model_parts, model_module

"""
    ThermalBaseParams

Keys every thermal model may use: `Thermal Conductivity` (read by the critical
time step for any thermal model, and by Thermal Flow).
"""
@params struct ThermalBaseParams
    thermal_conductivity::Union{Nothing,Float64} = opt("Thermal Conductivity";
                                                       default = nothing, min = 0)
end

# registration runs at load time, never during precompilation
__init__() = register_base!(:thermal, ThermalBaseParams)
for file in find_registered_modules(@__DIR__, "register_thermal")
    include(file)
end

export init_model
export compute_model
export init_fields
export fields_for_local_synchronization

"""
    init_fields()

Initialize thermal model fields
"""
function init_fields()
    Data_Manager.create_node_scalar_field("Temperature", Float64)
    Data_Manager.create_constant_node_scalar_field("Delta Temperature", Float64)
    Data_Manager.create_node_scalar_field("Heat Flow", Float64)
    Data_Manager.create_constant_node_scalar_field("Specific Volume", Float64)
    # if it is already initialized via mesh file no new field is created here
    Data_Manager.create_constant_node_scalar_field("Surface_Nodes", Bool;
                                                   default_value = true)
end

"""
    compute_model(nodes::AbstractVector{Int64}, thermal::WithBase, block::Int64, time::Float64, dt::Float64)

Computes every part of the block's thermal model, in deck order.

# Arguments
- `nodes::AbstractVector{Int64}`: The nodes
- `thermal::WithBase`: The typed thermal model of the block (`thermal.base`: Thermal Conductivity)
- `block::Int64`: The block
- `time::Float64`: The current time
- `dt::Float64`: The time step
"""
function compute_model(nodes::AbstractVector{Int64}, thermal::WithBase, block::Int64,
                       time::Float64, dt::Float64)
    for part in model_parts(thermal.model)
        @timeit "$(nameof(model_module(part)))" model_module(part).compute_model(nodes, part,
                                                                                  thermal, block,
                                                                                  time, dt)
    end
end

"""
    init_model(nodes::AbstractVector{Int64}, block::Int64)

Initializes every part of the block's thermal model (`Data_Manager.get_block_models(block).thermal`).

# Arguments
- `nodes::AbstractVector{Int64}`: The nodes.
- `block::Int64`: Block.
"""
function init_model(nodes::AbstractVector{Int64}, block::Int64)
    thermal = Data_Manager.get_block_models(block).thermal
    for part in model_parts(thermal.model)
        model_module(part).init_model(nodes, part, thermal, block)
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
    thermal = Data_Manager.get_block_models(block).thermal
    for part in model_parts(thermal.model)
        model_module(part).fields_for_local_synchronization(model)
    end
end

end
