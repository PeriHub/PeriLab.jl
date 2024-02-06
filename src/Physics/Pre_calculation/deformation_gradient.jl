# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Deformation_Gradient
include("../../Support/geometry.jl")
using .Geometry
export compute
"""
    compute(datamanager::Module, nodes::Union{SubArray,Vector{Int64}})

Compute the deformation gradient.

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
    inverse_shape_tensor = datamanager.get_field("Inverse Shape Tensor")

    deformation_gradient = Geometry.deformation_gradient(nodes, dof, nlist, volume, omega, bond_damage, deformed_bond, undeformed_bond, inverse_shape_tensor, deformation_gradient)

    return datamanager
end

end