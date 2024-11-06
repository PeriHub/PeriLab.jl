# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Shape_Tensor
using DataStructures: OrderedDict
include("../../Support/helpers.jl")
using .Helpers: find_active_nodes
include("../../Support/geometry.jl")
using .Geometry: compute_shape_tensors!
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
    init_model(datamanager, nodes, parameter)

Inits the shape tensor calculation.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `parameter::Dict(String, Any)`: Dictionary with parameter.
# Returns
- `datamanager::Data_manager`: Datamanager.

"""
function init_model(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    parameter::Union{Dict,OrderedDict},
    block::Int64,
)
    dof = datamanager.get_dof()
    datamanager.create_constant_node_field("Shape Tensor", Float64, "Matrix", dof)
    datamanager.create_constant_node_field("Inverse Shape Tensor", Float64, "Matrix", dof)
    return datamanager
end

"""
    compute(datamanager::Module, nodes::Union{SubArray,Vector{Int64}})

Compute the shape tensor.

# Arguments
- `datamanager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
# Returns
- `datamanager`: Datamanager.
"""
function compute(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    parameter::Union{Dict,OrderedDict},
    block::Int64,
)

    nlist = datamanager.get_nlist()
    volume = datamanager.get_field("Volume")
    omega = datamanager.get_field("Influence Function")
    bond_damage = datamanager.get_bond_damage("NP1")
    undeformed_bond = datamanager.get_field("Bond Geometry")
    shape_tensor = datamanager.get_field("Shape Tensor")
    inverse_shape_tensor = datamanager.get_field("Inverse Shape Tensor")
    update_list = datamanager.get_field("Update")
    active_nodes = datamanager.get_field("Active Nodes")
    active_nodes = find_active_nodes(update_list, active_nodes, nodes)

    compute_shape_tensors!(
        shape_tensor,
        inverse_shape_tensor,
        active_nodes,
        nlist,
        volume,
        omega,
        bond_damage,
        undeformed_bond,
    )

    return datamanager
end



"""
    fields_for_local_synchronization(datamanager::Module, model::String)

Returns a user developer defined local synchronization. This happens before each model.



# Arguments

"""
function fields_for_local_synchronization(datamanager::Module, model::String)
    # download_from_cores = false
    # upload_to_cores = true
    # datamanager.set_local_synch(model, "Bond Forces", download_from_cores, upload_to_cores)
    return datamanager
end

end
