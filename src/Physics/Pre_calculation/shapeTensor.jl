# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Shape_Tensor
include("../../Support/geometry.jl")
using .Geometry
export compute

function compute(datamanager::Module, nodes::Vector{Int64})
    dof = datamanager.get_dof()
    nlist = datamanager.get_nlist()
    volume = datamanager.get_field("Volume")
    omega = datamanager.get_field("Influence Function")
    bondDamage = datamanager.get_field("Bond Damage", "NP1")
    bondGeometry = datamanager.get_field("Bond Geometry")
    shapeTensor = datamanager.get_field("Shape Tensor")
    invShapeTensor = datamanager.get_field("Inverse Shape Tensor")
    shapeTensor, invShapeTensor = Geometry.shape_tensor(nodes, dof, nlist, volume, omega, bondDamage, bondGeometry, shapeTensor, invShapeTensor)
    return datamanager
end


end