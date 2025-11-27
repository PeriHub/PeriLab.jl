# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Degradation

using ....Data_Manager
using ...Solver_Manager: find_module_files, create_module_specifics
global module_list = find_module_files(@__DIR__, "degradation_name")
for mod in module_list
    include(mod["File"])
end

using ....Helpers: find_inverse_bond_id
export compute_model
export init_model
export init_fields
export fields_for_local_synchronization

"""
init_fields()

Initialize  model fields
"""
function init_fields()
    Data_Manager.create_node_scalar_field("Damage", Float64)
    nlist = Data_Manager.get_nlist()
    inverse_nlist = Data_Manager.set_inverse_nlist(find_inverse_bond_id(nlist))
end

"""
    compute_model(nodes::AbstractVector{Int64}, model_param::Dict, block::Int64, time::Float64, dt::Float64)

Computes the degradation models

# Arguments
- `nodes::AbstractVector{Int64}`: The nodes
- `model_param::Dict`: The model parameters
- `block::Int64`: The block
- `time::Float64`: The current time
- `dt::Float64`: The time step
"""
function compute_model(nodes::AbstractVector{Int64},
                       model_param::Dict,
                       block::Int64,
                       time::Float64,
                       dt::Float64)
    mod = Data_Manager.get_model_module(model_param["Degradation Model"])
    return mod.compute_model(nodes, model_param, block, time, dt)
end

"""
    init_model(nodes::AbstractVector{Int64}, block::Int64)

Initialize a degradation models.

# Arguments
- `nodes::AbstractVector{Int64}`: Nodes for the degradation model.
- `block::Int64`: Block identifier for the degradation model.

# Example
```julia
init_model(my_data_manager, [1, 2, 3], 1)
```
"""
function init_model(nodes::AbstractVector{Int64},
                    block::Int64)
    model_param = Data_Manager.get_properties(block, "Degradation Model")
    mod = create_module_specifics(model_param["Degradation Model"],
                                  module_list,
                                  @__MODULE__,
                                  "degradation_name")
    if isnothing(mod)
        @error "No degradation model of name " * model_param["Degradation Model"] *
               " exists."
        return nothing
    end
    Data_Manager.set_model_module(model_param["Degradation Model"], mod)
    mod.init_model(nodes, model_param, block)
end

"""
    fields_for_local_synchronization(model, block)

Defines all synchronization fields for local synchronization

# Arguments
- `model::String`: Model class.
- `block::Int64`: block ID
"""
function fields_for_local_synchronization(model, block)
    model_param = Data_Manager.get_properties(block, "Degradation Model")
    mod = Data_Manager.get_model_module(model_param["Degradation Model"])

    return mod.fields_for_local_synchronization(model)
end

end
