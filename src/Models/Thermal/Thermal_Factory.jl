# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Thermal
using TimerOutputs
include("../../Core/Module_inclusion/set_Modules.jl")
using .Set_modules
global module_list = Set_modules.find_module_files(@__DIR__, "thermal_model_name")
Set_modules.include_files(module_list)
using TimerOutputs
export init_model
export compute_model
export init_fields
export fields_for_local_synchronization

"""
    init_fields(datamanager::Module)

Initialize thermal model fields

# Arguments
- `datamanager::Data_manager`: Datamanager.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function init_fields(datamanager::Module)
    datamanager.create_node_field("Temperature", Float64, 1)
    datamanager.create_constant_node_field("Delta Temperature", Float64, 1)
    datamanager.create_node_field("Heat Flow", Float64, 1)
    datamanager.create_constant_node_field("Specific Volume", Float64, 1)
    # if it is already initialized via mesh file no new field is created here
    datamanager.create_constant_node_field("Surface_Nodes", Bool, 1, true)

    return datamanager
end


"""
    compute_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, model_param::Dict, block::Int64, time::Float64, dt::Float64,to::TimerOutput,)

Computes the thermal models

# Arguments
- `datamanager::Module`: The datamanager
- `nodes::Union{SubArray,Vector{Int64}}`: The nodes
- `model_param::Dict`: The model parameters
- `block::Int64`: The block
- `time::Float64`: The current time
- `dt::Float64`: The time step
# Returns
- `datamanager::Module`: The datamanager
"""
function compute_model(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    model_param::Dict,
    block::Int64,
    time::Float64,
    dt::Float64,
    to::TimerOutput,
)

    thermal_models = split(model_param["Thermal Model"], "+")
    thermal_models = map(r -> strip(r), thermal_models)
    for thermal_model in thermal_models
        mod = datamanager.get_model_module(thermal_model)
        datamanager = mod.compute_model(datamanager, nodes, model_param, block, time, dt)
    end

    return datamanager
end

"""
    init_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}, block::Int64)

Initializes the thermal model.

# Arguments
- `datamanager::Data_manager`: Datamanager
- `nodes::Union{SubArray,Vector{Int64}}`: The nodes.
- `block::Int64`: Block.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function init_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, block::Int64)
    model_param = datamanager.get_properties(block, "Thermal Model")
    thermal_models = split(model_param["Thermal Model"], "+")
    thermal_models = map(r -> strip(r), thermal_models)
    for thermal_model in thermal_models
        mod = Set_modules.create_module_specifics(
            thermal_model,
            module_list,
            "thermal_model_name",
        )
        if isnothing(mod)
            @error "No thermal model of name " * thermal_model * " exists."
            return nothing
        end
        datamanager.set_model_module(thermal_model, mod)
        datamanager = mod.init_model(datamanager, nodes, model_param)
    end
    return datamanager
end

"""
    fields_for_local_synchronization(datamanager, model, block)

Defines all synchronization fields for local synchronization

# Arguments
- `datamanager::Module`: datamanager.
- `model::String`: Model class.
- `block::Int64`: block ID
# Returns
- `datamanager::Module`: Datamanager.
"""

function fields_for_local_synchronization(datamanager, model, block)
    model_param = datamanager.get_properties(block, "Thermal Model")
    thermal_models = split(model_param["Thermal Model"], "+")
    thermal_models = map(r -> strip(r), thermal_models)
    for thermal_model in thermal_models
        mod = datamanager.get_model_module(thermal_model)
        mod.fields_for_local_synchronization(datamanager, model)
    end
    return datamanager
end

end
