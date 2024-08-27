# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Pre_calculation

include("../../Core/Module_inclusion/set_Modules.jl")
using .Set_modules
global module_list = Set_modules.find_module_files(@__DIR__, "pre_calculation_name")
Set_modules.include_files(module_list)

using TimerOutputs

export compute
export init_model
export synchronize
export init_fields

include("../../Support/helpers.jl")
using .Helpers: find_active, get_active_update_nodes


function init_fields(datamanager::Module)
    dof = datamanager.get_dof()
    deformed_coorN, deformed_coorNP1 =
        datamanager.create_node_field("Deformed Coordinates", Float64, dof)
    deformed_coorN = copy(datamanager.get_field("Coordinates"))
    deformed_coorNP1 = copy(datamanager.get_field("Coordinates"))
    datamanager.create_bond_field("Deformed Bond Geometry", Float64, dof)
    datamanager.create_bond_field("Deformed Bond Length", Float64, 1)
    datamanager.create_node_field("Displacements", Float64, dof)
    return datamanager
end

"""
    init_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}, block::Int64)

Initializes the model.

# Arguments
- `datamanager::Data_manager`: Datamanager
- `nodes::Union{SubArray,Vector{Int64}}`: The nodes.
- `block::Int64`: Block.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function init_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, block::Int64)
    dof = datamanager.get_dof()
    ## TODO options change
    #model_param = datamanager.get_properties(block, "Material Model") ??
    if options["Deformed Bond Geometry"]

    end

    if options["Deformation Gradient"]

        options["Shape Tensor"] = true
        options["Deformed Bond Geometry"] = true
        # order is important
        datamanager.set_model_module("Deformed Bond Geometry", Bond_Deformation)
        datamanager.set_model_module("Shape Tensor", Shape_Tensor)
        datamanager.set_model_module("Deformation Gradient", Deformation_Gradient)
    end

    if options["Bond Associated Deformation Gradient"]

        options["Deformed Bond Geometry"] = true
        # order is important
        datamanager.set_model_module("Deformed Bond Geometry", Bond_Deformation)
        datamanager.set_model_module(
            "Bond Associated Deformation Gradient",
            Bond_Deformation_Gradient,
        )
    end
    if options["Shape Tensor"]

        # order is important
        datamanager.set_model_module("Deformed Bond Geometry", Bond_Deformation)
        datamanager.set_model_module("Shape Tensor", Shape_Tensor)
    end
    return datamanager
end

"""
    fields_for_local_synchronization(model_param::Dict)

Finds all synchronization fields from the model class

# Arguments
- `model_param::Dict`: model parameter.
# Returns
- `synch_dict::Dict`: Synchronization Dictionary.
"""

function fields_for_local_synchronization((model_param::Dict))
    synch_dict = Dict()

    return synch_dict
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
function compute_model(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    model_param::Dict,
    time::Float64,
    dt::Float64,
    to::TimerOutput,
)
    models_options = datamanager.get_models_options()
    for (pre_calculation_model, active) in pairs(models_options)
        if !active
            continue
        end
        mod = datamanager.get_model_module(pre_calculation_model)
        datamanager = mod.compute(datamanager, nodes)
    end

    return datamanager
end


end
