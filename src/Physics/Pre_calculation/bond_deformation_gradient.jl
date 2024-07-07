# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Bond_Deformation_Gradient

include("../Material/Material_Models/Bond_Associated_Correspondence.jl")
using .Bond_Associated_Correspondence: find_local_neighbors
include("../../Support/geometry.jl")
using .Geometry: calculate_bond_length, bond_associated_deformation_gradient
export compute

"""
    compute(datamanager, nodes)

Compute the bond deformation gradient.

# Arguments
- `datamanager`: Datamanager.
- `nodes`: List of nodes.
# Returns
- `datamanager`: Datamanager.
"""


function compute(datamanager::Module, nodes::Union{SubArray,Vector{Int64}})

    dof = datamanager.get_dof()
    nlist = datamanager.get_nlist()
    volume = datamanager.get_field("Volume")
    omega = datamanager.get_field("Influence Function")
    bond_damage = datamanager.get_bond_damage("NP1")
    undeformed_bond = datamanager.get_field("Bond Geometry")
    deformed_bond = datamanager.get_field("Deformed Bond Geometry", "NP1")
    deformation_gradient = datamanager.get_field("Deformation Gradient")
    gradient_weights = datamanager.get_field("Lagrangian Gradient Weights")




    return datamanager
end





end