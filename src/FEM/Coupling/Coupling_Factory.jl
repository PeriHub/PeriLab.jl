# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Coupling

using ....Data_Manager
using ....InputDeck: FEMParams
using ....PeriLabExceptions: @abort
using ....ModuleLoader: find_module_files, create_module_specifics

global module_list = find_module_files(@__DIR__, "coupling_name")
for mod in module_list
    include(mod["File"])
end

export init_coupling
export compute_coupling

function init_coupling(nodes, fem::FEMParams)
    Data_Manager.create_constant_node_scalar_field("PD Nodes", Int64)
    fem.coupling === nothing && return
    coupling_model = fem.coupling.coupling_type

    mod = create_module_specifics(coupling_model, module_list,
                                  @__MODULE__, "coupling_name")
    if isnothing(mod)
        @abort "No coupling model of name " * coupling_model * " exists."
        return
    end
    Data_Manager.set_model_module(coupling_model, mod)

    ###TODO nodes and blocks
    mod.init_coupling_model(nodes, fem)
end

function compute_coupling(fem::FEMParams)
    fem.coupling === nothing && return
    mod = Data_Manager.get_model_module(fem.coupling.coupling_type)
    return mod.compute_coupling(fem.coupling)
end

end
