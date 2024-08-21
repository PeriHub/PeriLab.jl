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
include("../../Support/helpers.jl")
using .Helpers: find_active, get_active_update_nodes


function init_fields(datamanager::Module)
    dof = datamanager.get_dof()
    deformed_coorN, deformed_coorNP1 =
        datamanager.create_node_field("Deformed Coordinates", Float64, dof)
    deformed_coorN = copy(datamanager.get_field("Coordinates"))
    deformed_coorNP1 = copy(datamanager.get_field("Coordinates"))
    datamanager.create_node_field("Displacements", Float64, dof)
end

function synchronize(datamanager::Module, options::Dict, synchronise_field)
    if options["Bond Associated Deformation Gradient"]
        synchfield = Dict(
            "Deformation Gradient" => Dict(
                "upload_to_cores" => true,
                "dof" => datamanager.get_dof() * datamanager.get_dof(),
            ),
            "Weighted Volume" => Dict("upload_to_cores" => true, "dof" => 1),
        )
        synchronise_field(
            datamanager.get_comm(),
            synchfield,
            datamanager.get_overlap_map(),
            datamanager.get_field,
            "Deformation Gradient",
            "upload_to_cores",
        )
        synchronise_field(
            datamanager.get_comm(),
            synchfield,
            datamanager.get_overlap_map(),
            datamanager.get_field,
            "Weighted Volume",
            "upload_to_cores",
        )
    end
end

"""
    compute(datamanager::Module, block_nodes::Dict{Int64,Vector{Int64}})

Compute the pre-calculation.

# Arguments
- `datamanager`: Datamanager.
- `block_nodes::Dict{Int64,Vector{Int64}}`: List of block nodes.
# Returns
- `datamanager`: Datamanager.
"""
function compute(datamanager::Module, block_nodes::Dict{Int64,Vector{Int64}})
    active = datamanager.get_field("Active")
    update_list = datamanager.get_field("Update List")
    fem_option = datamanager.fem_active()
    models_options = datamanager.get_models_options()
    for pre_calculation_model in keys(models_options)
        if !(models_options[pre_calculation_model])
            continue
        end
        mod = datamanager.get_model_module(pre_calculation_model)
        for block in eachindex(block_nodes)
            nodes = block_nodes[block]
            active_nodes, update_nodes =
                get_active_update_nodes(active, update_list, block_nodes, block)
            if fem_option
                fe_nodes = datamanager.get_field("FE Nodes")
                update_nodes =
                    block_nodes[block][find_active(Vector{Bool}(.~fe_nodes[update_nodes]))]
            end
            datamanager = mod.compute(datamanager, update_nodes, block)
        end
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
        datamanager.set_model_module("Deformed Bond Geometry", Bond_Deformation)
    end

    if options["Deformation Gradient"]
        datamanager.create_constant_node_field(
            "Deformation Gradient",
            Float64,
            "Matrix",
            dof,
        )
        options["Shape Tensor"] = true
        options["Deformed Bond Geometry"] = true
        # order is important
        datamanager.set_model_module("Deformed Bond Geometry", Bond_Deformation)
        datamanager.set_model_module("Shape Tensor", Shape_Tensor)
        datamanager.set_model_module("Deformation Gradient", Deformation_Gradient)
    end

    if options["Bond Associated Deformation Gradient"]
        datamanager.create_constant_bond_field(
            "Bond Associated Deformation Gradient",
            Float64,
            "Matrix",
            dof,
        )
        datamanager.create_constant_node_field(
            "Deformation Gradient",
            Float64,
            "Matrix",
            dof,
        )
        datamanager.create_constant_bond_field("Lagrangian Gradient Weights", Float64, dof)
        datamanager.create_constant_node_field(
            "Weighted Deformation Gradient",
            Float64,
            "Matrix",
            dof,
        )
        options["Deformed Bond Geometry"] = true
        # order is important
        datamanager.set_model_module("Deformed Bond Geometry", Bond_Deformation)
        datamanager.set_model_module(
            "Bond Associated Deformation Gradient",
            Bond_Deformation_Gradient,
        )
    end
    if options["Shape Tensor"]
        datamanager.create_constant_node_field("Shape Tensor", Float64, "Matrix", dof)
        datamanager.create_constant_node_field(
            "Inverse Shape Tensor",
            Float64,
            "Matrix",
            dof,
        )
        # order is important
        datamanager.set_model_module("Deformed Bond Geometry", Bond_Deformation)
        datamanager.set_model_module("Shape Tensor", Shape_Tensor)
    end
    return datamanager
end

end
