# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Bondbased_Elastic
include("../material_basis.jl")
export init_material_model
export material_name
export compute_force

# global dof::Int64
# global horizon::Vector{Float64}
"""
    init_material_model(datamanager::Module)

    Initializes the material model.

    Parameters:
    - `datamanager::Data_manager`: Datamanager.

    Returns:
    - `datamanager::Data_manager`: Datamanager.
"""
function init_material_model(datamanager::Module)
    # global dof
    # global horizon

    # dof = datamanager.get_dof()
    # horizon = datamanager.get_field("Horizon")

    return datamanager
end


function material_name()
    return "Bond-based Elastic"
end
function compute_forces(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, material_parameter::Dict, time::Float64, dt::Float64)
    # global dof
    # global horizon
    dof = datamanager.get_dof()
    horizon = datamanager.get_field("Horizon")
    symmetry::String = get_symmmetry(material_parameter)
    undeformed_bond = datamanager.get_field("Bond Geometry")
    deformed_bond = datamanager.get_field("Deformed Bond Geometry", "NP1")
    bond_damage = datamanager.get_field("Bond Damage", "NP1")
    bond_force = datamanager.get_field("Bond Forces")

    K = material_parameter["Bulk Modulus"]
    E = material_parameter["Young's Modulus"]
    for iID in nodes
        if symmetry == "plane stress"
            constant = 12.0 * E / (pi * (1 + 1.0 / 3) * horizon[iID]^3)
        elseif symmetry == "plane strain"
            constant = 12.0 * E / (pi * (1.25) * horizon[iID]^3)
        else
            constant = 18.0 * K / (pi * horizon[iID]^4)
        end
        if deformed_bond[iID][:, end] == 0
            @error "Length of bond is zero due to its deformation."
        end
        # Calculate the bond force
        bond_force[iID] = (0.5 .* constant .* bond_damage[iID][:] .* (deformed_bond[iID][:, end] .- undeformed_bond[iID][:, end]) ./ undeformed_bond[iID][:, end]) .* deformed_bond[iID][:, 1:dof] ./ deformed_bond[iID][:, end]

    end

    return datamanager
end

end