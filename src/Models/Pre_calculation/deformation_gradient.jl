# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Deformation_Gradient
using DataStructures: OrderedDict

using .......Data_Manager
using .......Geometry: compute_deformation_gradients!
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
    init_model(nodes, parameter)

Inits the deformation gradient calculation.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `parameter::Dict(String, Any)`: Dictionary with parameter.

"""
function init_model(nodes::AbstractVector{Int64},
                    parameter::Union{Dict,OrderedDict},
                    block::Int64)
    dof = Data_Manager.get_dof()
    Data_Manager.create_constant_node_field("Deformation Gradient", Float64, dof,
                                            VectorOrMatrix = "Matrix")
end

"""
    compute(nodes::AbstractVector{Int64})

Compute the deformation gradient.

# Arguments
- `nodes`: List of nodes.
"""
function compute(nodes::AbstractVector{Int64},
                 parameter::Union{Dict,OrderedDict},
                 block::Int64)
    nlist = Data_Manager.get_nlist()
    volume = Data_Manager.get_field("Volume")
    omega = Data_Manager.get_field("Influence Function")
    bond_damage = Data_Manager.get_bond_damage("NP1")
    undeformed_bond = Data_Manager.get_field("Bond Geometry")
    deformed_bond = Data_Manager.get_field("Deformed Bond Geometry", "NP1")
    deformation_gradient = Data_Manager.get_field("Deformation Gradient")
    inverse_shape_tensor = Data_Manager.get_field("Inverse Shape Tensor")
    dof = Data_Manager.get_dof()
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
