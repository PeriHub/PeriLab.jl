# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Degradation
include("../../Core/Module_inclusion/set_Modules.jl")
using .Set_modules
global module_list = Set_modules.find_module_files(@__DIR__, "degradation_name")
Set_modules.include_files(module_list)
include("../../Support/Helpers.jl")
using TimerOutputs
using .Helpers: find_inverse_bond_id
export compute_model
export init_model
export init_fields
export fields_for_local_synchronization

"""
init_fields(datamanager::Module)

Initialize  model fields

# Arguments
- `datamanager::Data_manager`: Datamanager.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function init_fields(datamanager::Module)
    datamanager.create_node_field("Damage", Float64, 1)
    nlist = datamanager.get_nlist()
    inverse_nlist = datamanager.set_inverse_nlist(find_inverse_bond_id(nlist))
    return datamanager
end

"""
    compute_model(datamanager::Module, nodes::AbstractVector{Int64}, model_param::Dict, block::Int64, time::Float64, dt::Float64,to::TimerOutput,)

Computes the degradation models

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
    mod = datamanager.get_model_module(model_param["Degradation Model"])
    return mod.compute_model(datamanager, nodes, model_param, block, time, dt)
end

"""
    init_model(datamanager::Module, nodes::AbstractVector{Int64}, block::Int64)

Initialize a degradation models.

# Arguments
- `datamanager::Module`: The data manager module where the degradation model will be initialized.
- `nodes::AbstractVector{Int64}`: Nodes for the degradation model.
- `block::Int64`: Block identifier for the degradation model.

# Returns
- `datamanager`: The modified data manager module with the initialized degradation model.

# Example
```julia
datamanager = init_model(my_data_manager, [1, 2, 3], 1)

"""
function init_model(datamanager::Module, nodes::AbstractVector{Int64},
                    block::Int64)
    model_param = datamanager.get_properties(block, "Degradation Model")
    mod = Set_modules.create_module_specifics(model_param["Degradation Model"],
                                              module_list,
                                              "degradation_name")
    if isnothing(mod)
        @error "No degradation model of name " * model_param["Degradation Model"] *
               " exists."
        return nothing
    end
    datamanager.set_model_module(model_param["Degradation Model"], mod)
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
    model_param = datamanager.get_properties(block, "Degradation Model")
    mod = datamanager.get_model_module(model_param["Degradation Model"])

    return mod.fields_for_local_synchronization(datamanager, model)
end

end
