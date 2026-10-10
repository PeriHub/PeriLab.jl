# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Deformation_Gradient

using .......Data_Manager
using .......Geometry: compute_deformation_gradients!
export init_model
export compute_model

using ......ParameterSpec: @params, register_pre_calculation
"Switch only: this pre-calculation has no parameters."
@params struct DeformationGradientParams
end
__init__() = register_pre_calculation("Deformation Gradient", DeformationGradientParams)

"""
    init_model(nodes, p, pre_calculation, block)

Inits the deformation gradient calculation.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `p::DeformationGradientParams`: The model parameters.
- `pre_calculation`: `nothing` (pre-calculations have no shared parameters).
- `block::Int64`: The current block.

"""
function init_model(nodes::AbstractVector{Int64}, p::DeformationGradientParams,
                    pre_calculation, block::Int64)
    dof = Data_Manager.get_dof()
    Data_Manager.create_constant_node_tensor_field("Deformation Gradient", Float64, dof)
end

"""
    compute_model(nodes, p, pre_calculation, block, time, dt)

Compute the deformation gradient.

# Arguments
- `nodes`: List of nodes.
- `p::DeformationGradientParams`: The model parameters.
- `pre_calculation`: `nothing` (pre-calculations have no shared parameters).
- `block::Int64`: The current block.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
"""
function compute_model(nodes::AbstractVector{Int64}, p::DeformationGradientParams,
                       pre_calculation, block::Int64, time::Float64, dt::Float64)
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
