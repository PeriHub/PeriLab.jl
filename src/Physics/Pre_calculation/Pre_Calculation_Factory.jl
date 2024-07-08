# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Pre_calculation
include("bond_deformation.jl")
include("bond_deformation_gradient.jl")

include("deformation_gradient.jl")
include("shape_tensor.jl")

using TimerOutputs
using .Bond_Deformation
using .Bond_Deformation_Gradient
using .Deformation_Gradient
using .Shape_Tensor
export compute
export init_pre_calculation
export synchronize

function synchronize(datamanager::Module, options::Dict, synchronise_field)
    if options["Bond Associated Deformation Gradient"]
        synchronise_field(datamanager.get_comm(), datamanager.get_synch_fields(), datamanager.get_overlap_map(), datamanager.get_field, "Deformation Gradient", "upload_to_cores")
        synchronise_field(datamanager.get_comm(), datamanager.get_synch_fields(), datamanager.get_overlap_map(), datamanager.get_field, "Weighted Volume", "upload_to_cores")
    end
end

"""
    compute(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, options::Dict, time::Float64, dt::Float64, block_id::Int64, to::TimerOutput)

Compute the pre-calculation.

# Arguments
- `datamanager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `options::Dict`: Options.
- `time::Float64`: Time.
- `dt::Float64`: Time step.
- `block_id::Int64`: Block ID
- `to::TimerOutput`::TimerOutput 
# Returns
- `datamanager`: Datamanager.
"""
function compute(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, options::Dict, time::Float64, dt::Float64, block_id::Int64, to::TimerOutput)

    if options["Deformed Bond Geometry"]
        @timeit to "Deformed Bond Geometry" datamanager = Bond_Deformation.compute(datamanager, nodes, block_id)
    end
    if options["Shape Tensor"]
        @timeit to "Shape Tensor" datamanager = Shape_Tensor.compute(datamanager, nodes, block_id)
    end
    if options["Deformation Gradient"]
        @timeit to "Deformation Gradient" datamanager = Deformation_Gradient.compute(datamanager, nodes, block_id)
    end
    if options["Bond Associated Deformation Gradient"]
        @timeit to "Deformation Gradient" datamanager = Bond_Deformation_Gradient.compute(datamanager, nodes, block_id)
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
        options["Deformed Bond Geometry"] = true
    end

    if options["Bond Associated Deformation Gradient"]
        datamanager.create_constant_bond_field("Bond Associated Deformation Gradient", Float64, "Matrix", dof)
        datamanager.create_constant_node_field("Weighted Volume", Float64, 1)
        datamanager.create_constant_bond_field("Lagrangian Gradient Weights", Float64, dof)
        options["Deformed Bond Geometry"] = true
    end
    if options["Shape Tensor"]
        datamanager.create_constant_node_field("Shape Tensor", Float64, "Matrix", dof)
        datamanager.create_constant_node_field("Inverse Shape Tensor", Float64, "Matrix", dof)
    end
    return datamanager
end

end