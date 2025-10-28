# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Coupling

using ....Data_Manager
using ...Solver_Manager: find_module_files, create_module_specifics

global module_list = find_module_files(@__DIR__, "coupling_name")
for mod in module_list
    include(mod["File"])
end

export init_coupling
export compute_coupling

function init_coupling(nodes, complete_params::Dict)
    Data_Manager.create_constant_node_field("PD Nodes", Int64, 1)
    if !haskey(complete_params["FEM"], "Coupling")
        return
    end
    if !haskey(complete_params["FEM"]["Coupling"], "Coupling Type")
        @error "''Coupling Type'' is missing and must be defined, options are Alrequin. "
        return nothing
    end
    coupling_model = complete_params["FEM"]["Coupling"]["Coupling Type"]

    mod = create_module_specifics(coupling_model, module_list,
                                  @__MODULE__, "coupling_name")
    if isnothing(mod)
        @error "No material of name " * material_model * " exists."
    end
    Data_Manager.set_model_module(coupling_model, mod)

    ###TODO nodes and blocks
    mod.init_coupling_model(nodes,
                            convert(Dict{String,Any}, complete_params["FEM"]))
end

function compute_coupling(fem_params::Dict)
    if !haskey(fem_params, "Coupling")
        return
    end
    coupling_model = fem_params["Coupling"]["Coupling Type"]
    mod = Data_Manager.get_model_module(coupling_model)
    return mod.compute_coupling(fem_params)
end

end
