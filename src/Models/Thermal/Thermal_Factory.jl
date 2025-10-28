# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Thermal

using ....Data_Manager
using ...Solver_Manager: find_module_files, create_module_specifics
global module_list = find_module_files(@__DIR__, "thermal_model_name")
for mod in module_list
    include(mod["File"])
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
    Data_Manager.create_node_field("Temperature", Float64, 1)
    Data_Manager.create_constant_node_field("Delta Temperature", Float64, 1)
    Data_Manager.create_node_field("Heat Flow", Float64, 1)
    Data_Manager.create_constant_node_field("Specific Volume", Int64, 1)
    # if it is already initialized via mesh file no new field is created here
    Data_Manager.create_constant_node_field("Surface_Nodes", Bool, 1, true)
end

"""
    compute_model(nodes::AbstractVector{Int64}, model_param::Dict, block::Int64, time::Float64, dt::Float64)

Computes the thermal models

# Arguments
- `nodes::AbstractVector{Int64}`: The nodes
- `model_param::Dict`: The model parameters
- `block::Int64`: The block
- `time::Float64`: The current time
- `dt::Float64`: The time step
"""
function compute_model(nodes::AbstractVector{Int64},
                       model_param::Dict,
                       block::Int64,
                       time::Float64,
                       dt::Float64)
    thermal_models = split(model_param["Thermal Model"], "+")
    thermal_models = map(r -> strip(r), thermal_models)
    for thermal_model in thermal_models
        mod = get_model_module(thermal_model)
        mod.compute_model(nodes, model_param, block, time, dt)
    end
end

"""
    init_model(nodes::Union{SubArray,Vector{Int64}, block::Int64)

Initializes the thermal model.

# Arguments
- `nodes::AbstractVector{Int64}`: The nodes.
- `block::Int64`: Block.
"""
function init_model(nodes::AbstractVector{Int64},
                    block::Int64)
    model_param = Data_Manager.get_properties(block, "Thermal Model")
    thermal_models = split(model_param["Thermal Model"], "+")
    thermal_models = map(r -> strip(r), thermal_models)
    for thermal_model in thermal_models
        mod = create_module_specifics(thermal_model,
                                      module_list,
                                      @__MODULE__,
                                      "thermal_model_name")
        if isnothing(mod)
            @error "No thermal model of name " * thermal_model * " exists."
            return nothing
        end
        Data_Manager.set_model_module(thermal_model, mod)
        mod.init_model(nodes, model_param)
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
    model_param = Data_Manager.get_properties(block, "Thermal Model")
    thermal_models = split(model_param["Thermal Model"], "+")
    thermal_models = map(r -> strip(r), thermal_models)
    for thermal_model in thermal_models
        mod = Data_Manager.get_model_module(thermal_model)
        mod.fields_for_local_synchronization(model)
    end
end

end
