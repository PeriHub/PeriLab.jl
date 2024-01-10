# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Thermal

include("../../Core/Module_inclusion/set_Modules.jl")
using .Set_modules
global module_list = Set_modules.find_module_files(@__DIR__, "thermal_model_name")
Set_modules.include_files(module_list)
export init_thermal_model
export compute_thermal_model

"""
    compute_thermal_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, model_param::Dict, time::Float64, dt::Float64)

Compute the thermal model

# Arguments
- `datamanager::Module`: The datamanager
- `nodes::Union{SubArray,Vector{Int64}}`: The nodes
- `model_param::Dict`: The model parameters
- `time::Float64`: The current time
- `dt::Float64`: The time step
# Returns
- `datamanager::Module`: The datamanager
"""
function compute_thermal_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, model_param::Dict, time::Float64, dt::Float64)

    specifics = Dict{String,String}("Call Function" => "compute_thermal_model", "Name" => "thermal_model_name")

    thermal_models = split(model_param["Thermal Model"], "+")
    thermal_models = map(r -> strip(r), thermal_models)
    for thermal_model in thermal_models
        datamanager = Set_modules.create_module_specifics(thermal_model, module_list, specifics, (datamanager, nodes, model_param, time, dt))
        if isnothing(datamanager)
            @error "No thermal model of name " * model_name * " exists."
        end
    end
    datamanager = distribute_heat_flows(datamanager, nodes)
    return datamanager
end

"""
    compute_thermal_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, model_param::Dict, time::Float64, dt::Float64)

Compute the thermal model

# Arguments
- `datamanager::Module`: The datamanager
- `nodes::Union{SubArray,Vector{Int64}}`: The nodes
- `model_param::Dict`: The model parameters
- `time::Float64`: The current time
- `dt::Float64`: The time step
# Returns
- `datamanager::Module`: The datamanager
"""
function init_thermal_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, block::Int64)
    model_param = datamanager.get_properties(block, "Thermal Model")
    specifics = Dict{String,String}("Call Function" => "init_thermal_model", "Name" => "thermal_model_name")

    thermal_models = split(model_param["Thermal Model"], "+")
    thermal_models = map(r -> strip(r), thermal_models)
    for thermal_model in thermal_models
        datamanager = Set_modules.create_module_specifics(thermal_model, module_list, specifics, (datamanager, nodes, model_param, time, dt))
        if isnothing(datamanager)
            @error "No thermal model of name " * model_name * " exists."
        end
    end
    return datamanager
end

"""
    distribute_heat_flows(datamanager::Module, nodes::Union{SubArray,Vector{Int64}})

Distribute the heat flow
Note: is included, because also additional heat flow influences can be included as well here and it might be important for bond-associated formulations.

# Arguments
- `datamanager::Module`: The datamanager
- `nodes::Union{SubArray,Vector{Int64}}`: The nodes
# Returns
- `datamanager::Module`: The datamanager
"""
function distribute_heat_flows(datamanager::Module, nodes::Union{SubArray,Vector{Int64}})
    bond_heat_flow = datamanager.get_field("Bond Heat Flow")
    heat_flow = datamanager.get_field("Heat Flow", "NP1")
    volume = datamanager.get_field("Volume")
    nlist = datamanager.get_nlist()
    for iID in nodes
        for (jID, neighborID) in enumerate(nlist[iID])
            heat_flow[iID] -= bond_heat_flow[iID][jID] * volume[neighborID]
        end
    end
    return datamanager
end
end