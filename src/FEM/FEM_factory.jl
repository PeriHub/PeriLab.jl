module FEM
include("../../Core/Module_inclusion/set_Modules.jl")
# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause


using .Set_modules

global module_list = Set_modules.find_module_files(@__DIR__, "element_name")
Set_modules.include_files(module_list)

datamanager = Set_modules.create_module_specifics(fem_model, module_list, specifics, (datamanager, nodes, model_param, time, dt))
if isnothing(datamanager)
    @error "No shape function of name " * fem_model * " exists."
end



end