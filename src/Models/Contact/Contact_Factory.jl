# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Contact

include("../../Core/Module_inclusion/set_Modules.jl")
using .Set_modules
global module_list = Set_modules.find_module_files(@__DIR__, "contact_model_name")
Set_modules.include_files(module_list)
export init_contact_model
export compute_model

"""
    init_contact_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}, block::Int64)

Initializes the contact model.

# Arguments
- `datamanager::Data_manager`: Datamanager
- `nodes::Union{SubArray,Vector{Int64}}`: The nodes.
- `block::Int64`: Block.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function init_contact_model(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    block::Int64,
)
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
            @error "No contact of name " * material_model * " exists."
        end
        datamanager.set_model_module(thermal_model, mod)
        datamanager = mod.init_thermal_model(datamanager, nodes, model_param)
    end
    return datamanager
end

"""
    compute_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, model_param::Dict, time::Float64, dt::Float64, to::TimerOutput)

Compute the forces.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: The nodes.
- `model_param::Dict`: The material parameter.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function compute_model(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    model_param::Dict,
    time::Float64,
    dt::Float64,
    to::TimerOutput,
)
    material_models = split(model_param["Material Model"], "+")
    material_models = map(r -> strip(r), material_models)
    if occursin("Correspondence", model_param["Material Model"])
        mod = datamanager.get_model_module("Correspondence")
        datamanager = mod.compute_model(datamanager, nodes, model_param, time, dt, to)
        return datamanager
    end
    for material_model in material_models
        mod = datamanager.get_model_module(material_model)
        datamanager = mod.compute_model(datamanager, nodes, model_param, time, dt, to)
    end
    return datamanager
end
end
