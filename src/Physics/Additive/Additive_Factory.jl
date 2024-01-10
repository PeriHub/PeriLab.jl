# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Additive
include("../../Core/Module_inclusion/set_Modules.jl")
using .Set_modules
global module_list = Set_modules.find_module_files(@__DIR__, "additive_name")
Set_modules.include_files(module_list)
export compute_additive_model
export init_additive_model

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
function compute_additive_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, model_param::Dict, time::Float64, dt::Float64)

    specifics = Dict{String,String}("Call Function" => "compute_additive_model", "Name" => "additive_name")

    datamanager = Set_modules.create_module_specifics(model_param["Additive Model"], module_list, specifics, (datamanager, nodes, model_param, time, dt))
    if isnothing(datamanager)
        @error "No additive model of name " * model_param["Additive Model"] * " exists."
        return nothing
    end

    return datamanager
end


"""
    init_additive_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, block::Int64)

Initialize an additive model within a given data manager.

# Arguments
- `datamanager::Module`: The data manager module where the additive model will be initialized.
- `nodes::Union{SubArray,Vector{Int64}}`: Nodes for the additive model.
- `block::Int64`: Block identifier for the additive model.

# Returns
- `datamanager`: The modified data manager module with the initialized additive model.

# Example
```julia
datamanager = init_additive_model(my_data_manager, [1, 2, 3], 1)

"""
function init_additive_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, block::Int64)
    model_param = datamanager.get_properties(block, "Additive Model")
    specifics = Dict{String,String}("Call Function" => "init_additive_model", "Name" => "additive_name")

    datamanager = Set_modules.create_module_specifics(model_param["Additive Model"], module_list, specifics, (datamanager, nodes, model_param, time, dt))
    if isnothing(datamanager)
        @error "No additive model of name " * model_param["Additive Model"] * " exists."
        return nothing
    end

    return datamanager
end

end