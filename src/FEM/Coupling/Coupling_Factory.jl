# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Coupling_PF_FEM
include("../../Core/Module_inclusion/set_Modules.jl")
using .Set_modules
global module_list = Set_modules.find_module_files(@__DIR__, "coupling_name")
Set_modules.include_files(module_list)

export init_coupling
export compute_coupling

function init_coupling(datamanager::Module, complete_params::Dict)
    datamanager.create_constant_node_field("PD Nodes", Int64, 1)
    if !haskey(complete_params["FEM"], "Coupling")
        return datamanager
    end
    if !haskey(complete_params["FEM"]["Coupling"], "Coupling Type")
        @error "''Coupling Type'' is missing and must be defined, options are Alrequin. "
        return nothing
    end
    coupling_model = complete_params["FEM"]["Coupling"]["Coupling Type"]

    mod = Set_modules.create_module_specifics(coupling_model, module_list, "coupling_name")
    if isnothing(mod)
        @error "No material of name " * material_model * " exists."
    end
    datamanager.set_model_module(coupling_model, mod)

    ###TODO nodes and blcoks
    datamanager =
        mod.init_coupling_model(datamanager, nodes, complete_params["FEM"]["Coupling"])
    return datamanager
end

function compute_coupling(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    fem_params::Dict,
)

    if !haskey(fem_params, "Coupling")
        return datamanager
    end
    coupling_model = fem_params["Coupling"]["Coupling Type"]

    mod = datamanager.get_model_module(coupling_model)
    return mod.compute_coupling(datamanager, nodes, fem_params)

end



end
