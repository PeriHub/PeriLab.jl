# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Correspondence

using TimerOutputs: @timeit

using .....Data_Manager
using .....PeriLabExceptions: @abort
using ....Zero_Energy_Control
using .....ModuleLoader: find_registered_modules
using ......ParameterSpec: model_parts, model_module, model_function
for file in find_registered_modules(@__DIR__, "register_material")
    include(file)
end
include("../Material_template/correspondence_template.jl")

using LinearAlgebra
using StaticArrays: SVector
using Rotations
include("./Bond_Associated_Correspondence.jl")
using .Bond_Associated_Correspondence
using ....Material_Basis: compute_Piola_Kirchhoff_stress!
using .......Helpers: invert, rotate, node_tensor, determinant, smat, matrix_diff!,
                      fast_mul!, mat_mul!
using .......Geometry: compute_strain!, compute_linear_strain!

export init_model
export compute_model
export fields_for_local_synchronization
"""
    init_model(nodes, material, block)

Initializes the correspondence model and its parts for a block.

# Arguments
- `nodes::AbstractVector{Int64}`: The block nodes.
- `material::BlockMaterial`: The typed block material (base, moduli, symmetry).
- `block::Int64`: The block.
"""
function init_model(nodes::AbstractVector{Int64}, material, block::Int64)
    # a missing Symmetry means isotropic (the Hooke matrix uses material.hooke_symmetry);
    # in 2D the Hooke matrix aborts without plane strain / plane stress
    dof = Data_Manager.get_dof()
    Data_Manager.create_node_tensor_field("Strain", Float64, dof)
    Data_Manager.create_constant_node_tensor_field("Strain Increment", Float64, dof)
    Data_Manager.create_node_tensor_field("Cauchy Stress", Float64, dof)
    Data_Manager.create_node_scalar_field("von Mises Stress", Float64)
    for part in model_parts(material.model)
        model_module(part).init_model(nodes, part, material, block)
    end
    if material.base.bond_associated
        return Bond_Associated_Correspondence.init_model(nodes, material)
    end
    Zero_Energy_Control.init_model(nodes, material, block)
end

"""
    fields_for_local_synchronization(model, material)

Registers the fields the correspondence parts synchronize locally.
"""
function fields_for_local_synchronization(model::String, material)
    for part in model_parts(material.model)
        model_module(part).fields_for_local_synchronization(model)
    end
    if material.base.bond_associated
        Bond_Associated_Correspondence.fields_for_local_synchronization(model)
    end
end

"""
    compute_model(nodes, material, block, time, dt)

Computes the correspondence forces of a block.

# Arguments
- `nodes::AbstractVector{Int64}`: The block nodes.
- `material::BlockMaterial`: The typed block material (base, moduli, symmetry).
- `block::Int64`: The block.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
"""
function compute_model(nodes::AbstractVector{Int64}, material, block::Int64, time::Float64,
                       dt::Float64)
    if material.base.bond_associated
        return Bond_Associated_Correspondence.compute_model(nodes, material, block, time,
                                                            dt)
    end
    return compute_correspondence_model(nodes, material, block, time, dt)
end

"""
    compute_correspondence_model(nodes, material, block, time, dt)

Computes deformation gradient, strain increment and stresses of the correspondence parts.

# Arguments
- `nodes::AbstractVector{Int64}`: The block nodes.
- `material::BlockMaterial`: The typed block material (base, moduli, symmetry).
- `block::Int64`: The block.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
"""
function compute_correspondence_model(nodes::AbstractVector{Int64},
                                      material,
                                      block::Int64,
                                      time::Float64,
                                      dt::Float64)
    rotation::Bool = Data_Manager.get_rotation()
    dof::Int64 = Data_Manager.get_dof()
    deformation_gradient::NodeTensorField{Float64} = Data_Manager.get_field("Deformation Gradient")
    bond_force::BondVectorField{Float64} = Data_Manager.get_field("Bond Forces")
    bond_damage::BondScalarState{Float64} = Data_Manager.get_bond_damage("NP1")
    undeformed_bond::BondVectorField{Float64} = Data_Manager.get_field("Bond Geometry")
    inverse_shape_tensor::NodeTensorField{Float64} = Data_Manager.get_field("Inverse Shape Tensor")
    omega::BondScalarState{Float64} = Data_Manager.get_field("Influence Function")
    strain_N::NodeTensorField{Float64} = Data_Manager.get_field("Strain", "N")
    strain_NP1::NodeTensorField{Float64} = Data_Manager.get_field("Strain", "NP1")
    stress_N::NodeTensorField{Float64} = Data_Manager.get_field("Cauchy Stress", "N")
    stress_NP1::NodeTensorField{Float64} = Data_Manager.get_field("Cauchy Stress", "NP1")
    strain_increment::NodeTensorField{Float64} = Data_Manager.get_field("Strain Increment")

    if material.base.linear_strain
        @timeit "compute linear strain" compute_linear_strain!(nodes, deformation_gradient,
                                                               strain_NP1)
    else
        @timeit "compute strain" compute_strain!(nodes, deformation_gradient, strain_NP1)
    end
    @timeit "compute matrix diff" matrix_diff!(strain_increment, nodes, strain_NP1,
                                               strain_N)
    #@abort ""
    if rotation
        @timeit "rotate forward" begin
            rotation_tensor::NodeTensorField{Float64} = Data_Manager.get_field("Rotation Tensor")
            rotate(nodes, stress_N, rotation_tensor, false)
            rotate(nodes, strain_increment,
                   rotation_tensor, false)
        end
    end

    @timeit "compute material" begin
        for part in model_parts(material.model)
            model_function(part, Val(:compute_stresses))(nodes, dof, part, material,
                                                         time, dt, strain_increment,
                                                         stress_N, stress_NP1)
        end
    end

    if rotation
        @timeit "rotate back" begin
            rotate(nodes, stress_NP1, rotation_tensor, true)
        end
    end

    @timeit "compute bond force" calculate_bond_force!(nodes,
                                                       dof,
                                                       deformation_gradient,
                                                       undeformed_bond,
                                                       bond_damage,
                                                       omega,
                                                       inverse_shape_tensor,
                                                       stress_NP1,
                                                       bond_force)

    force_filter::Bool = Data_Manager.has_key("Anisotropic Force Filter")

    @timeit "zero energy" Zero_Energy_Control.compute_zero_energy_control(nodes,
                                                                          material,
                                                                          block,
                                                                          time,
                                                                          dt)
    if force_filter
        @timeit "anisotropic bond damage" filter_failed_axes!(nodes, bond_force,
                                                              bond_damage)
    end
end

"""
    filter_failed_axes!(nodes, bond_force, bond_damage)

Anisotropic bond damage: removes the bond force components along the failed material
axes, f = (I - sum_k chi_k a_k a_k^T) f. The axes a_k are the columns of the rotation
tensor. The scalar bond damage still scales the force as before.
"""
function filter_failed_axes!(nodes::AbstractVector{Int64},
                             bond_force::BondVectorState{Float64,D},
                             bond_damage::BondScalarState{Float64}) where {D}
    axis_failure = Data_Manager.get_field("Anisotropic Axis Failure")::BondVectorState{Bool,
                                                                                       D}
    rotation_tensor = Data_Manager.get_field("Rotation Tensor")::NodeTensorField{Float64}
    @inbounds for iID in nodes
        failure = axis_failure[iID]
        any(any, failure) || continue
        forces = bond_force[iID]
        R = node_tensor(rotation_tensor, iID, Val(D))
        for jID in eachindex(forces)
            any(failure[jID]) || continue
            # material components of the force, those along failed axes removed
            forces[jID] = R * ((transpose(R) * forces[jID]) .* .!failure[jID])
        end
    end
end

"""
    calculate_bond_force(nodes::AbstractVector{Int64}, deformation_gradient::Array{Float64, 3}, undeformed_bond::Vector{Matrix{Float64}}, bond_damage::BondScalarState{Float64}, inverse_shape_tensor::Array{Float64, 3}, stress_NP1::Array{Float64, 3}, bond_force::Vector{Matrix{Float64}})

Calculate bond forces for specified nodes based on deformation gradients.

# Arguments
- `nodes::AbstractVector{Int64}`: List of block nodes.
- `deformation_gradient::SubArray`: Deformation gradient.
- `undeformed_bond::BondVectorState{Float64}`: Undeformed bond geometry.
- `bond_damage::BondScalarState{Float64}`: Bond damage.
- `omega::BondScalarState{Float64}`: Influence function.
- `inverse_shape_tensor::Array{Float64, 3}`: Inverse shape tensor.
- `stress_NP1::Array{Float64, 3}`: Stress at time step n+1.
- `bond_force::BondVectorState{Float64}`: Bond force.
# Returns
- `bond_force::BondVectorState{Float64}`: Bond force.
"""
function calculate_bond_force!(nodes::AbstractVector{Int64},
                               dof::Int64,
                               deformation_gradient::NodeTensorField{Float64},
                               undeformed_bond::BondVectorState{Float64},
                               bond_damage::BondScalarState{Float64},
                               omega::BondScalarState{Float64},
                               inverse_shape_tensor::NodeTensorField{Float64},
                               stress_NP1::NodeTensorField{Float64},
                               bond_force::BondVectorState{Float64})
    calculate_bond_force!(bond_force, nodes, deformation_gradient, undeformed_bond,
                          bond_damage, omega, inverse_shape_tensor, stress_NP1, Val(dof))
end

# f_ij = ω d P K⁻¹ ξ with the first Piola-Kirchhoff stress P = F S. "Cauchy Stress" holds
# S: the parts integrate C:ΔE with the Green-Lagrange strain E (rotation invariant); the
# "Rotation Tensor" is the fixed material orientation, not the polar rotation of F. P Fᵀ
# = F S Fᵀ is symmetric (angular momentum). For Linear Strain, F σ is consistent to first
# order. A singular F (node without bonds) is zero, so P = 0.
function calculate_bond_force!(bond_force::BondVectorState{Float64},
                               nodes::AbstractVector{Int64},
                               deformation_gradient::NodeTensorField{Float64},
                               undeformed_bond::BondVectorState{Float64},
                               bond_damage::BondScalarState{Float64},
                               omega::BondScalarState{Float64},
                               inverse_shape_tensor::NodeTensorField{Float64},
                               stress_NP1::NodeTensorField{Float64},
                               ::Val{D}) where {D}
    @inbounds for iID in nodes
        T = node_tensor(deformation_gradient, iID, Val(D)) *
            node_tensor(stress_NP1, iID, Val(D)) *
            node_tensor(inverse_shape_tensor, iID, Val(D))
        @fastmath for jID in eachindex(bond_damage[iID])
            bond_force[iID][jID] = (bond_damage[iID][jID] * omega[iID][jID]) *
                                   (T * SVector{D}(undeformed_bond[iID][jID]))
        end
    end
    return nothing
end

end
