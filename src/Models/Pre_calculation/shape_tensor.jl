# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Shape_Tensor

using ......Data_Manager
using ......Helpers: find_active_nodes
using ......Geometry: compute_shape_tensors!
export init_model
export compute_model
export fields_for_local_synchronization

using ......ParameterSpec: @params, register_pre_calculation
"Switch only: this pre-calculation has no parameters."
@params struct ShapeTensorParams
end
__init__() = register_pre_calculation("Shape Tensor", ShapeTensorParams)

"""
    init_model(nodes, p, pre_calculation, block)

Inits the shape tensor calculation.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `p::ShapeTensorParams`: The model parameters.
- `pre_calculation`: `nothing` (pre-calculations have no shared parameters).
- `block::Int64`: The current block.

"""
function init_model(nodes::AbstractVector{Int64}, p::ShapeTensorParams, pre_calculation,
                    block::Int64)
    dof = Data_Manager.get_dof()
    Data_Manager.create_constant_node_tensor_field("Shape Tensor", Float64, dof)
    Data_Manager.create_constant_node_tensor_field("Inverse Shape Tensor", Float64, dof)

    # should be done here, because the shape tensor is needed for the matrix based correspondence models
    nlist = Data_Manager.get_nlist()
    volume = Data_Manager.get_field("Volume")
    omega = Data_Manager.get_field("Influence Function")
    bond_damage = Data_Manager.get_bond_damage("NP1")
    undeformed_bond = Data_Manager.get_field("Bond Geometry")
    shape_tensor = Data_Manager.get_field("Shape Tensor")
    inverse_shape_tensor = Data_Manager.get_field("Inverse Shape Tensor")

    compute_shape_tensors!(shape_tensor,
                           inverse_shape_tensor,
                           nodes,
                           nlist,
                           volume,
                           omega,
                           bond_damage,
                           undeformed_bond)
end

"""
    compute_model(nodes, p, pre_calculation, block, time, dt)

Compute the shape tensor.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `p::ShapeTensorParams`: The model parameters.
- `pre_calculation`: `nothing` (pre-calculations have no shared parameters).
- `block::Int64`: The current block.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
"""
function compute_model(nodes::AbstractVector{Int64}, p::ShapeTensorParams, pre_calculation,
                       block::Int64, time::Float64, dt::Float64)
    nlist::BondScalarState{Int64} = Data_Manager.get_nlist()
    volume::NodeScalarField{Float64} = Data_Manager.get_field("Volume")
    omega::BondScalarState{Float64} = Data_Manager.get_field("Influence Function")
    bond_damage::BondScalarState{Float64} = Data_Manager.get_bond_damage("NP1")
    undeformed_bond::BondVectorState{Float64} = Data_Manager.get_field("Bond Geometry")
    shape_tensor::NodeTensorField{Float64} = Data_Manager.get_field("Shape Tensor")
    inverse_shape_tensor::NodeTensorField{Float64} = Data_Manager.get_field("Inverse Shape Tensor")
    # update_list = Data_Manager.get_field("Update")
    # active_nodes = Data_Manager.get_field("Active Nodes")
    # active_nodes = find_active_nodes(update_list, active_nodes, nodes)

    compute_shape_tensors!(shape_tensor,
                           inverse_shape_tensor,
                           nodes,
                           nlist,
                           volume,
                           omega,
                           bond_damage,
                           undeformed_bond)
end

"""
    fields_for_local_synchronization(model::String)

Returns a user developer defined local synchronization. This happens before each model.



# Arguments

"""
function fields_for_local_synchronization(model::String)
    # download_from_cores = false
    # upload_to_cores = true
    # Data_Manager.set_local_synch(model, "Bond Forces", download_from_cores, upload_to_cores)
end

end
