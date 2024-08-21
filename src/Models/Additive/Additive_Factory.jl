# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Additive
include("../../Core/Module_inclusion/set_Modules.jl")
using .Set_modules
global module_list = Set_modules.find_module_files(@__DIR__, "additive_name")
Set_modules.include_files(module_list)
include("../../Support/helpers.jl")

using .Helpers: find_inverse_bond_id
export compute_additive_model
export init_additive_model
export init_additive_model_fields


"""
    init_additive_model_fields(datamanager::Module)

Initialize additive model fields

# Arguments
- `datamanager::Data_manager`: Datamanager.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function init_additive_model_fields(datamanager::Module)
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
        for iID = 1:nnodes
            bond_damageN[iID] .= 0
            bond_damageNP1[iID] .= 0
        end
    end
    nlist = datamanager.get_field("Neighborhoodlist")
    inverse_nlist = datamanager.set_inverse_nlist(find_inverse_bond_id(nlist))

    return datamanager
end

"""
    compute_additive_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, model_param::Dict, time::Float64, dt::Float64)

Computes the additive model

# Arguments
- `datamanager::Module`: The datamanager
- `nodes::Union{SubArray,Vector{Int64}}`: The nodes
- `model_param::Dict`: The model parameters
- `time::Float64`: The current time
- `dt::Float64`: The time step
# Returns
- `datamanager::Module`: The datamanager
"""
function compute_additive_model(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    model_param::Dict,
    time::Float64,
    dt::Float64,
)

    mod = datamanager.get_model_module(model_param["Additive Model"])
    return mod.compute_additive_model(datamanager, nodes, model_param, time, dt)

end


"""
    init_additive_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, additive_parameter::Dict, block::Int64)

Initialize an additive model within a given data manager.

# Arguments
- `datamanager::Module`: The data manager module where the additive model will be initialized.
- `nodes::Union{SubArray,Vector{Int64}}`: Nodes for the additive model.
- `additive_parameter::Dict`: Additive parameters
- `block::Int64`: Block identifier for the additive model.

# Returns
- `datamanager`: The modified data manager module with the initialized additive model.

# Example
```julia
datamanager = init_additive_model(my_data_manager, [1, 2, 3], 1)

"""
function init_additive_model(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    block::Int64,
)
    model_param = datamanager.get_properties(block, "Additive Model")
    mod = Set_modules.create_module_specifics(
        model_param["Additive Model"],
        module_list,
        "additive_name",
    )
    if isnothing(mod)
        @error "No additive model of name " * model_param["Additive Model"] * " exists."
        return nothing
    end
    datamanager.set_model_module(model_param["Additive Model"], mod)
    datamanager = mod.init_additive_model(datamanager, nodes, model_param, block)
    return datamanager
end

end
