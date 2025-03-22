# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Deformation_Gradient
using DataStructures: OrderedDict
include("../../Support/Geometry.jl")
using .Geometry: compute_deformation_gradients!
export pre_calculation_name
export init_model
export compute

"""
    pre_calculation_name()

Gives the pre_calculation name. It is needed for comparison with the yaml input deck.

# Arguments

# Returns
- `name::String`: The name of the Pre_Calculation.

Example:
```julia
println(pre_calculation_name())
"Deformation Gradient"
```
"""
function pre_calculation_name()
    return "Deformation Gradient"
end

"""
    init_model(datamanager, nodes, parameter)

Inits the deformation gradient calculation.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: List of block nodes.
- `parameter::Dict(String, Any)`: Dictionary with parameter.
# Returns
- `datamanager::Data_manager`: Datamanager.

"""
function init_model(datamanager::Module,
                    nodes::Union{SubArray,Vector{Int64}},
                    parameter::Union{Dict,OrderedDict},
                    block::Int64)
    dof = datamanager.get_dof()
    datamanager.create_constant_node_field("Deformation Gradient", Float64, "Matrix", dof)
    return datamanager
end

"""
    compute(datamanager::Module, nodes::Union{SubArray,Vector{Int64}})

Compute the deformation gradient.

# Arguments
- `datamanager`: Datamanager.
- `nodes`: List of nodes.
# Returns
- `datamanager`: Datamanager.
"""
function compute(datamanager::Module,
                 nodes::Union{SubArray,Vector{Int64}},
                 parameter::Union{Dict,OrderedDict},
                 block::Int64)
    nlist = datamanager.get_nlist()
    volume = datamanager.get_field("Volume")
    omega = datamanager.get_field("Influence Function")
    bond_damage = datamanager.get_bond_damage("NP1")
    undeformed_bond = datamanager.get_field("Bond Geometry")
    deformed_bond = datamanager.get_field("Deformed Bond Geometry", "NP1")
    deformation_gradient = datamanager.get_field("Deformation Gradient")
    inverse_shape_tensor = datamanager.get_field("Inverse Shape Tensor")
    dof = datamanager.get_dof()
    compute_deformation_gradients!(deformation_gradient,
                                   nodes,
                                   dof,
                                   nlist,
                                   volume,
                                   omega,
                                   bond_damage,
                                   deformed_bond,
                                   undeformed_bond,
                                   inverse_shape_tensor)

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
