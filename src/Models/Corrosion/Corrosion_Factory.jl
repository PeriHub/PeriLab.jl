# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Corrosion
include("../../Core/Module_inclusion/set_Modules.jl")
using .Set_modules
global module_list = Set_modules.find_module_files(@__DIR__, "corrosion_name")
Set_modules.include_files(module_list)
export compute_corrosion_model
export init_model
export init_fields



"""
init_fields(datamanager::Module)

Initialize concentration model fields

# Arguments
- `datamanager::Data_manager`: Datamanager.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function init_fields(datamanager::Module)
    datamanager.create_node_field("Concentration", Float64, 1)
    datamanager.create_node_field("Concentration Flux", Float64, 1)
    return datamanager
end

"""
    compute_corrosion_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, model_param::Dict, time::Float64, dt::Float64)

Computes the corrosion model

# Arguments
- `datamanager::Module`: The datamanager
- `nodes::Union{SubArray,Vector{Int64}}`: The nodes
- `model_param::Dict`: The model parameters
- `time::Float64`: The current time
- `dt::Float64`: The time step
# Returns
- `datamanager::Module`: The datamanager
"""
function compute_corrosion_model(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    model_param::Dict,
    time::Float64,
    dt::Float64,
)

    mod = datamanager.get_model_module(model_param["Corrosion Model"])
    return mod.compute_corrosion_model(datamanager, nodes, model_param, time, dt)

end


"""
    init_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, corrosion_parameter::Dict, block::Int64)

Initialize an corrosion model within a given data manager.

# Arguments
- `datamanager::Module`: The data manager module where the corrosion model will be initialized.
- `nodes::Union{SubArray,Vector{Int64}}`: Nodes for the corrosion model.
- `corrosion_parameter::Dict`: Corrosion parameters
- `block::Int64`: Block identifier for the corrosion model.

# Returns
- `datamanager`: The modified data manager module with the initialized corrosion model.

# Example
```julia
datamanager = init_model(my_data_manager, [1, 2, 3], 1)

"""
function init_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, block::Int64)
    model_param = datamanager.get_properties(block, "Corrosion Model")
    mod = Set_modules.create_module_specifics(
        model_param["Corrosion Model"],
        module_list,
        "corrosion_name",
    )
    if isnothing(mod)
        @error "No corrosion model of name " * model_param["Corrosion Model"] * " exists."
        return nothing
    end
    datamanager.set_model_module(model_param["Corrosion Model"], mod)
    datamanager = mod.init_model(datamanager, nodes, model_param, block)
    return datamanager
end

end
