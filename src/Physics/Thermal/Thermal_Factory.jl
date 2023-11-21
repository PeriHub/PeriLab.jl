# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Thermal

include("../../Core/Module_inclusion/set_Modules.jl")
using .Set_modules
global module_list = Set_modules.find_module_files(@__DIR__, "thermal_model_name")
Set_modules.include_files(module_list)

export compute_thermal_model

function compute_thermal_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, model_param::Dict, time::Float64, dt::Float64)

    specifics = Dict{String,String}("Call Function" => "compute_thermal_model", "Name" => "thermal_model_name")

    thermal_models = split(model_param["Thermal Model"], "+")
    for thermal_model in thermal_models
        datamanager = Set_modules.create_module_specifics(thermal_model, module_list, specifics, (datamanager, nodes, model_param, time, dt))
        if isnothing(datamanager)
            @error "No thermal model of name " * model_name * " exists."
        end
        datamanager = distribute_heat_flows(datamanager, nodes)
    end
    return datamanager
end
"""

Note: is included, because also additional heat flow influences can be included as well here and it might be important for bond-associated formulations.

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