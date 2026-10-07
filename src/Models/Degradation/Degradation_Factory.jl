# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Degradation

using ....Data_Manager
using ....PeriLabExceptions: @abort
using ....ModuleLoader: find_module_files
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

model_module(p) = parentmodule(typeof(p))

"""
    compute_model(nodes::AbstractVector{Int64}, p, block::Int64, time::Float64, dt::Float64)

Computes the degradation model of a block.

# Arguments
- `nodes::AbstractVector{Int64}`: The nodes
- `p`: The typed degradation model of the block
- `block::Int64`: The block
- `time::Float64`: The current time
- `dt::Float64`: The time step
"""
function compute_model(nodes::AbstractVector{Int64}, p, block::Int64, time::Float64,
                       dt::Float64)
    return model_module(p).compute_model(nodes, p, block, time, dt)
end

"""
    init_model(nodes::AbstractVector{Int64}, block::Int64)

Initialize the degradation model of a block (`Data_Manager.get_block_model("Degradation Model", block)`).

# Arguments
- `nodes::AbstractVector{Int64}`: Nodes for the degradation model.
- `block::Int64`: Block identifier for the degradation model.
"""
function init_model(nodes::AbstractVector{Int64}, block::Int64)
    p = Data_Manager.get_block_model("Degradation Model", block)
    model_module(p).init_model(nodes, p, block)
end

"""
    fields_for_local_synchronization(model, block)

Defines all synchronization fields for local synchronization

# Arguments
- `model::String`: Model class.
- `block::Int64`: block ID
"""
function fields_for_local_synchronization(model, block)
    p = Data_Manager.get_block_model("Degradation Model", block)
    return model_module(p).fields_for_local_synchronization(model)
end

end
