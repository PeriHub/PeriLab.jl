# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Shape_Tensor
using DataStructures: OrderedDict

using ......Data_Manager
using ......Helpers: find_active_nodes
using ......Geometry: compute_shape_tensors!
export pre_calculation_name
export init_model
export compute
export fields_for_local_synchronization

"""
    pre_calculation_name()

Gives the pre_calculation name. It is needed for comparison with the yaml input deck.

# Arguments

# Returns
- `name::String`: The name of the Pre_Calculation.

Example:
```julia
println(pre_calculation_name())
"Shape Tensor"
```
"""
function pre_calculation_name()
    return "Shape Tensor"
end

"""
    init_model(nodes, parameter)

Inits the shape tensor calculation.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `parameter::Dict(String, Any)`: Dictionary with parameter.

"""
function init_model(nodes::AbstractVector{Int64},
                    parameter::Union{Dict,OrderedDict},
                    block::Int64)
    dof = Data_Manager.get_dof()
    Data_Manager.create_constant_node_field("Shape Tensor", Float64, dof,
                                            VectorOrMatrix = "Matrix")
    Data_Manager.create_constant_node_field("Inverse Shape Tensor", Float64, dof,
                                            VectorOrMatrix = "Matrix")
end

"""
    compute(nodes::AbstractVector{Int64})

Compute the shape tensor.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
"""
function compute(nodes::AbstractVector{Int64},
                 parameter::Union{Dict,OrderedDict},
                 block::Int64)
    nlist = Data_Manager.get_nlist()
    volume = Data_Manager.get_field("Volume")
    omega = Data_Manager.get_field("Influence Function")
    bond_damage = Data_Manager.get_bond_damage("NP1")
    undeformed_bond = Data_Manager.get_field("Bond Geometry")
    shape_tensor = Data_Manager.get_field("Shape Tensor")
    inverse_shape_tensor = Data_Manager.get_field("Inverse Shape Tensor")
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
