# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Coupling

using ...Solver_Manager: find_module_files, create_module_specifics
global module_list = find_module_files(@__DIR__, "coupling_name")
for mod in module_list
    include(mod["File"])
end

export init_coupling
export compute_coupling

function init_coupling(datamanager::Module, nodes, complete_params::Dict)
    datamanager.create_constant_node_field("PD Nodes", Int64, 1)
    if !haskey(complete_params["FEM"], "Coupling")
        return datamanager
    end
    if !haskey(complete_params["FEM"]["Coupling"], "Coupling Type")
        @error "''Coupling Type'' is missing and must be defined, options are Alrequin. "
        return nothing
    end
    coupling_model = complete_params["FEM"]["Coupling"]["Coupling Type"]

    mod = create_module_specifics(coupling_model, module_list, "coupling_name")
    if isnothing(mod)
        @error "No material of name " * material_model * " exists."
    end
    datamanager.set_model_module(coupling_model, mod)

    ###TODO nodes and blocks
    datamanager = mod.init_coupling_model(datamanager, nodes,
                                          convert(Dict{String,Any}, complete_params["FEM"]))
    return datamanager
end

function compute_coupling(datamanager::Module, fem_params::Dict)
    if !haskey(fem_params, "Coupling")
        return datamanager
    end
    coupling_model = fem_params["Coupling"]["Coupling Type"]
    mod = datamanager.get_model_module(coupling_model)
    return mod.compute_coupling(datamanager, fem_params)
end

end
