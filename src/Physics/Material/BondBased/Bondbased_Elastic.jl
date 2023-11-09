# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Bondbased_Elastic
export compute_force
export material_name


function material_name()
    return "Bond-based Elastic"
end
function compute_forces(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, material_parameter::Dict, time::Float64, dt::Float64)

    bond_geometry = datamanager.get_field("Bond Geometry")
    deformed_bond = datamanager.get_field("Deformed Bond Geometry", "NP1")
    bond_damage = datamanager.get_field("Bond Damage", "NP1")
    bond_force = datamanager.get_field("Bond Forces")
    horizon = datamanager.get_field("Horizon")
    dof = datamanager.get_dof()

    K = material_parameter["Bulk Modulus"]

    for iID in nodes
        constant = 18.0 * K / (pi * horizon[iID]^4)
        if deformed_bond[iID][:, end] == 0
            @error "Length of bond is zero due to its deformation."
        end
        # Calculate the bond force
        bond_force[iID] = (0.5 .* constant .* bond_damage[iID][:] .* (deformed_bond[iID][:, end] .- bond_geometry[iID][:, end]) ./ bond_geometry[iID][:, end]) .* deformed_bond[iID][:, 1:dof] ./ deformed_bond[iID][:, end]

    end

    return datamanager
end

end