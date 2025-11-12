# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Influence_Function
using ...Data_Manager

function init_influence_function(nodes::AbstractVector{Int64},
                                 params::Dict)
    if !haskey(params, "Influence Function")
        return
    end
    bond_length = Data_Manager.get_field("Bond Length")
    omega = Data_Manager.get_field("Influence Function")
    if params["Influence Function"] == "1/xi^2"
        omega = influence_function_1_div_xi_squared(omega, nodes, bond_length)
    end
end

function influence_function_1_div_xi_squared(omega, nodes, bond_length)
    for iID in nodes
        omega[iID][:] = 1 ./ (bond_length[iID][:] .* bond_length[iID][:])
    end
    return omega
end
end
