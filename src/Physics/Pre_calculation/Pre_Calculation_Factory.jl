# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Pre_calculation
include("bond_deformation.jl")
include("bond_deformation_gradient.jl")
include("bond_shapeTensor.jl")
include("deformation_gradient.jl")
include("shapeTensor.jl")

using .Bond_Deformation
using .Bond_Deformation_Gradient
using .Bond_Shape_Tensor
using .Deformation_Gradient
using .Shape_Tensor
export compute
export init_pre_calculation

function compute(datamanager, nodes, options, time, dt)

    if options["Deformed Bond Geometry"]
        datamanager = Bond_Deformation.compute(datamanager, nodes)
    end
    if options["Shape Tensor"]
        datamanager = Shape_Tensor.compute(datamanager, nodes)
    end
    if options["Deformation Gradient"]
        datamanager = Deformation_Gradient.compute(datamanager, nodes)
    end
    if options["Bond Associated Shape Tensor"]
        datamanager = Bond_Shape_Tensor.compute(datamanager, nodes)
    end
    if options["Bond Associated Deformation Gradient"]
        datamanager = Bond_Deformation_Gradient.compute(datamanager, nodes)
    end
    return datamanager
end

function init_pre_calculation(datamanager, options)
    dof = datamanager.get_dof()
    if options["Deformed Bond Geometry"]
        bond_defN, bond_defNP1 = datamanager.create_bond_field("Deformed Bond Geometry", Float32, dof + 1)
    end
    if options["Shape Tensor"]
        shapeTensor = datamanager.create_constant_node_field("Shape Tensor", Float32, "Matrix", dof)
        invShapeTensor = datamanager.create_constant_node_field("Inverse Shape Tensor", Float32, "Matrix", dof)
    end
    if options["Deformation Gradient"]
        defGradN, defGradNP1 = datamanager.create_node_field("Deformation Gradient", Float32, "Matrix", dof)
    end
    if options["Bond Associated Shape Tensor"]
        bondShapeTensor = datamanager.create_constant_bond_field("Bond Associated Shape Tensor", Float32, "Matrix", dof)
        invBondShapeTensor = datamanager.create_constant_bond_field("Inverse Bond Associated Shape Tensor", Float32, "Matrix", dof)
    end
    if options["Bond Associated Deformation Gradient"]
        bondDefGradN, bondDefGradNP1 = datamanager.create_bond_field("Bond Associated Deformation Gradient", Float32, "Matrix", dof)
    end
    return datamanager
end

end