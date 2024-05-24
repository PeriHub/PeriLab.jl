# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Pre_calculation
include("bond_deformation.jl")
include("bond_deformation_gradient.jl")
include("bond_shape_tensor.jl")
include("deformation_gradient.jl")
include("shape_tensor.jl")

using TimerOutputs
using .Bond_Deformation
using .Bond_Deformation_Gradient
using .Bond_Shape_Tensor
using .Deformation_Gradient
using .Shape_Tensor
export compute
export init_pre_calculation

"""
    compute(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, options, time::Float64, dt::Float64)

Compute the pre-calculation.

# Arguments
- `datamanager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `options`: Options.
- `time::Float64`: Time.
- `dt::Float64`: Time step.
# Returns
- `datamanager`: Datamanager.
"""
function compute(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, options, time::Float64, dt::Float64, to::TimerOutput)

    if options["Deformed Bond Geometry"]
        @timeit to "Deformed Bond Geometry" datamanager = Bond_Deformation.compute(datamanager, nodes, time)
    end
    if options["Shape Tensor"]
        @timeit to "Shape Tensor" datamanager = Shape_Tensor.compute(datamanager, nodes)
    end
    if options["Deformation Gradient"]
        @timeit to "Deformation Gradient" datamanager = Deformation_Gradient.compute(datamanager, nodes)
    end
    if options["Bond Associated Shape Tensor"]
        @timeit to "Bond Associated Shape Tensor" datamanager = Bond_Shape_Tensor.compute(datamanager, nodes)
    end
    if options["Bond Associated Deformation Gradient"]
        @timeit to "Bond Associated Deformation Gradient" datamanager = Bond_Deformation_Gradient.compute(datamanager, nodes)
    end
    return datamanager
end

"""
    init_pre_calculation(datamanager::Module, options::Dict)

Initialize the pre-calculation.

# Arguments
- `datamanager`: Datamanager.
- `options::Dict`: Options.
# Returns
- `datamanager`: Datamanager.
"""
function init_pre_calculation(datamanager::Module, options::Dict)
    dof = datamanager.get_dof()
    if options["Deformed Bond Geometry"]
        datamanager.create_bond_field("Deformed Bond Geometry", Float64, dof)
        datamanager.create_bond_field("Deformed Bond Length", Float64, 1)
    end

    if options["Deformation Gradient"]
        datamanager.create_constant_node_field("Deformation Gradient", Float64, "Matrix", dof)
        options["Shape Tensor"] = true
    end

    if options["Bond Associated Deformation Gradient"]
        datamanager.create_constant_bond_field("Bond Associated Deformation Gradient", Float64, "Matrix", dof)
        options["Bond Associated Shape Tensor"] = true
    end
    if options["Shape Tensor"]
        datamanager.create_constant_node_field("Shape Tensor", Float64, "Matrix", dof)
        datamanager.create_constant_node_field("Inverse Shape Tensor", Float64, "Matrix", dof)
    end
    if options["Bond Associated Shape Tensor"]
        datamanager.create_constant_bond_field("Bond Associated Shape Tensor", Float64, "Matrix", dof)
        datamanager.create_constant_bond_field("Inverse Bond Associated Shape Tensor", Float64, "Matrix", dof)

        for block_id in data_managerget_block_list()
            if isnothing(datamanager.get_property(block_id, "Material Model", "Bond Horizon"))
                horizon = datamanager.get_field("Horizon")
                datamanager.set_property(block_id, "Material Model", "Bond Horizon", maximum(horizon))
            end
        end
    end
    return datamanager
end

end