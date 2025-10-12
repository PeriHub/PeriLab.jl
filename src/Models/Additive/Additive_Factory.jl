# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Additive
include("../../Core/Module_inclusion/set_Modules.jl")
global module_list = find_module_files(@__DIR__, "additive_name")
include_files(module_list)
using TimerOutputs
using .....Helpers: find_inverse_bond_id
export compute_model
export init_model
export init_fields
export fields_for_local_synchronization

"""
    init_fields(datamanager::Module)

Initialize additive model fields

# Arguments
- `datamanager::Data_Manager`: Datamanager.
# Returns
- `datamanager::Data_Manager`: Datamanager.
"""
function init_fields(datamanager::Module)
    if !datamanager.has_key("Activation_Time")
        @error "'Activation_Time' is missing. Please define an 'Activation_Time' for each point in the mesh file."
    end
    # must be specified, because it might be that no temperature model has been defined
    datamanager.create_node_field("Temperature", Float64, 1)
    datamanager.create_node_field("Heat Flow", Float64, 1)

    bond_damageN = datamanager.get_bond_damage("N")
    bond_damageNP1 = datamanager.get_bond_damage("NP1")
    nnodes = datamanager.get_nnodes()
    if !datamanager.has_key("Active")
        active = datamanager.create_constant_node_field("Active", Bool, 1, false)
        for iID in 1:nnodes
            bond_damageN[iID] .= 0
            bond_damageNP1[iID] .= 0
        end
    end
    nlist = datamanager.get_nlist()
    inverse_nlist = datamanager.set_inverse_nlist(find_inverse_bond_id(nlist))

    return datamanager
end

"""
    compute_model(datamanager::Module, nodes::AbstractVector{Int64}, model_param::Dict, block::Int64, time::Float64, dt::Float64,to::TimerOutput,)

Computes the addtive models

# Arguments
- `datamanager::Module`: The datamanager
- `nodes::AbstractVector{Int64}`: The nodes
- `model_param::Dict`: The model parameters
- `block::Int64`: The block
- `time::Float64`: The current time
- `dt::Float64`: The time step
# Returns
- `datamanager::Module`: The datamanager
"""
function compute_model(datamanager::Module,
                       nodes::AbstractVector{Int64},
                       model_param::Dict,
                       block::Int64,
                       time::Float64,
                       dt::Float64,
                       to::TimerOutput)
    mod = datamanager.get_model_module(model_param["Additive Model"])
    return mod.compute_model(datamanager, nodes, model_param, block, time, dt)
end

"""
    init_model(datamanager::Module, nodes::AbstractVector{Int64}, block::Int64)

Initialize the additive models.

# Arguments
- `datamanager::Module`: The data manager module where the additive model will be initialized.
- `nodes::AbstractVector{Int64}`: Nodes for the additive model.
- `block::Int64`: Block identifier for the additive model.

# Returns
- `datamanager`: The modified data manager module with the initialized additive model.

# Example
```julia
datamanager = init_model(my_data_manager, [1, 2, 3], 1)

"""
function init_model(datamanager::Module, nodes::AbstractVector{Int64},
                    block::Int64)
    model_param = datamanager.get_properties(block, "Additive Model")
    mod = create_module_specifics(model_param["Additive Model"],
                                  module_list,
                                  "additive_name")
    if isnothing(mod)
        @error "No additive model of name " * model_param["Additive Model"] * " exists."
        return nothing
    end
    datamanager.set_model_module(model_param["Additive Model"], mod)
    datamanager = mod.init_model(datamanager, nodes, model_param, block)
    return datamanager
end

"""
    fields_for_local_synchronization(datamanager, model, block)

Defines all synchronization fields for local synchronization

# Arguments
- `datamanager::Module`: datamanager.
- `model::String`: Model class.
- `block::Int64`: block ID
# Returns
- `datamanager::Module`: Datamanager.
"""
function fields_for_local_synchronization(datamanager, model, block)
    model_param = datamanager.get_properties(block, "Additive Model")
    mod = datamanager.get_model_module(model_param["Additive Model"])

    return mod.fields_for_local_synchronization(datamanager, model)
end

end
