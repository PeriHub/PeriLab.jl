# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Thermal

include("../../Core/Module_inclusion/set_Modules.jl")
using .Set_module
global module_list = Set_modules.find_module_files(@__DIR__, "thermal_model_name")
Set_modules.include_files(module_list)

export compute_thermal_model

function compute_thermal_model(datamanager, nodes, model_param, time, dt)

    specifics = Dict{String,String}("Call Function" => "compute_thermal_model", "Name" => "thermal_model_name")

    datamanager = Set_modules.create_module_specifics(model_param["Thermal Model"], module_list, specifics, (datamanager, nodes, model_param, time, dt))
    if isnothing(datamanager)
        @error "No thermal model of name " * model_param["Thermal Model"] * " exists."
    end

    return datamanager
end

end