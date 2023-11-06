# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Deformation_Gradient
include("../../Support/geometry.jl")
using .Geometry
export compute

function compute(datamanager::Module, nodes::Union{SubArray,Vector{Int64}})
    dof = datamanager.get_dof()
    nlist = datamanager.get_nlist()
    volume = datamanager.get_field("Volume")
    omega = datamanager.get_field("Influence Function")
    bondDamage = datamanager.get_field("Bond Damage", "NP1")
    bondGeometry = datamanager.get_field("Bond Geometry")
    deformed_bondN = datamanager.get_field("Deformed Bond Geometry", "N")
    deformed_bondNP1 = datamanager.get_field("Deformed Bond Geometry", "NP1")
    defGrad = datamanager.get_field("Deformation Gradient")
    invShapeTensor = datamanager.get_field("Inverse Shape Tensor")

    defGrad = Geometry.deformation_gradient(nodes, dof, nlist, volume, omega, bondDamage, deformed_bondNP1, bondGeometry, invShapeTensor, defGrad)

    return datamanager
end

end