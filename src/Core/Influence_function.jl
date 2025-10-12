# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Influence_Function

function init_influence_function(nodes::AbstractVector{Int64},
                                 datamanager::Module,
                                 params::Dict)
    if !haskey(params, "Influence Function")
        return datamanager
    end
    bond_length = datamanager.get_field("Bond Length")
    omega = datamanager.get_field("Influence Function")
    if params["Influence Function"] == "1/xi^2"
        omega = influence_function_1_div_xi_squared(omega, nodes, bond_length)
    end
    return datamanager
end

function influence_function_1_div_xi_squared(omega, nodes, bond_length)
    for iID in nodes
        omega[iID][:] = 1 ./ (bond_length[iID][:] .* bond_length[iID][:])
    end
    return omega
end
end
