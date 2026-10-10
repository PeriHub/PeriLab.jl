# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Correspondence

using TimerOutputs: @timeit

using .....Data_Manager
using .....PeriLabExceptions: @abort
using ....Zero_Energy_Control
using .....ModuleLoader: find_registered_modules
using ......ParameterSpec: model_parts, model_module
for file in find_registered_modules(@__DIR__, "register_material")
    include(file)
end
include("../Material_template/correspondence_template.jl")

using LinearAlgebra
a = using LoopVectorization
using Rotations
include("./Bond_Associated_Correspondence.jl")
using .Bond_Associated_Correspondence
using ....Material_Basis: compute_Piola_Kirchhoff_stress!
using .......Helpers: invert, rotate, determinant, smat, matrix_diff!, fast_mul!, mat_mul!
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
        return Bond_Associated_Correspondence.compute_model(nodes, material, block, time, dt)
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
    bond_force::BondVectorState{Float64} = Data_Manager.get_field("Bond Forces")
    bond_damage::BondScalarState{Float64} = Data_Manager.get_bond_damage("NP1")
    undeformed_bond::BondVectorState{Float64} = Data_Manager.get_field("Bond Geometry")
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
            model_module(part).compute_stresses(nodes, dof, part, material, time, dt,
                                                        strain_increment, stress_N,
                                                        stress_NP1)
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

    @timeit "zero energy" Zero_Energy_Control.compute_zero_energy_control(nodes,
                                                                          material,
                                                                          block,
                                                                          time,
                                                                          dt)
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
    if dof == 2
        return calculate_bond_force_2d!(bond_force,
                                        nodes,
                                        deformation_gradient,
                                        undeformed_bond,
                                        bond_damage,
                                        omega,
                                        inverse_shape_tensor,
                                        stress_NP1)
    elseif dof == 3
        return calculate_bond_force_3d!(bond_force,
                                        nodes,
                                        deformation_gradient,
                                        undeformed_bond,
                                        bond_damage,
                                        omega,
                                        inverse_shape_tensor,
                                        stress_NP1)
    end
end

function calculate_bond_force_2d!(bond_force::BondVectorState{Float64},
                                  nodes::AbstractVector{Int64},
                                  deformation_gradient::NodeTensorField{Float64},
                                  undeformed_bond::BondVectorState{Float64},
                                  bond_damage::BondScalarState{Float64},
                                  omega::BondScalarState{Float64},
                                  inverse_shape_tensor::NodeTensorField{Float64},
                                  stress_NP1::NodeTensorField{Float64})
    pk11, pk12, pk21, pk22 = 0.0, 0.0, 0.0, 0.0
    t11, t12, t21, t22 = 0.0, 0.0, 0.0, 0.0

    @inbounds for iID in nodes
        F11 = deformation_gradient[iID, 1, 1]
        F12 = deformation_gradient[iID, 1, 2]
        F21 = deformation_gradient[iID, 2, 1]
        F22 = deformation_gradient[iID, 2, 2]

        sigma11 = stress_NP1[iID, 1, 1]
        sigma12 = stress_NP1[iID, 1, 2]
        sigma21 = stress_NP1[iID, 2, 1]
        sigma22 = stress_NP1[iID, 2, 2]
        # optimized calculation of inverse of 2x2 matrix
        detF = F11 * F22 - F12 * F21
        if abs(detF) > 1e-15
            inv_detF = 1.0 / detF
            Finv11 = F22 * inv_detF
            Finv12 = -F12 * inv_detF
            Finv21 = -F21 * inv_detF
            Finv22 = F11 * inv_detF

            pk11 = sigma11 * Finv11 + sigma12 * Finv21
            pk12 = sigma11 * Finv12 + sigma12 * Finv22
            pk21 = sigma21 * Finv11 + sigma22 * Finv21
            pk22 = sigma21 * Finv12 + sigma22 * Finv22
        else
            pk11 = pk12 = pk21 = pk22 = 0.0
        end

        # temp = pk_stress * K^(-1) (inline)
        K11 = inverse_shape_tensor[iID, 1, 1]
        K12 = inverse_shape_tensor[iID, 1, 2]
        K21 = inverse_shape_tensor[iID, 2, 1]
        K22 = inverse_shape_tensor[iID, 2, 2]

        t11 = pk11 * K11 + pk12 * K21
        t12 = pk11 * K12 + pk12 * K22
        t21 = pk21 * K11 + pk22 * K21
        t22 = pk21 * K12 + pk22 * K22

        # Bond forces (inline)
        @fastmath for jID in eachindex(bond_damage[iID])
            bd = bond_damage[iID][jID] * omega[iID][jID]
            xi1 = undeformed_bond[iID][jID][1]
            xi2 = undeformed_bond[iID][jID][2]

            bond_force[iID][jID][1] = bd * (t11 * xi1 + t12 * xi2)
            bond_force[iID][jID][2] = bd * (t21 * xi1 + t22 * xi2)
        end
    end
    return nothing
end

function calculate_bond_force_3d!(bond_force::BondVectorState{Float64},
                                  nodes::AbstractVector{Int64},
                                  deformation_gradient::NodeTensorField{Float64},
                                  undeformed_bond::BondVectorState{Float64},
                                  bond_damage::BondScalarState{Float64},
                                  omega::BondScalarState{Float64},
                                  inverse_shape_tensor::NodeTensorField{Float64},
                                  stress_NP1::NodeTensorField{Float64})
    pk11, pk12, pk13 = 0.0, 0.0, 0.0
    pk21, pk22, pk23 = 0.0, 0.0, 0.0
    pk31, pk32, pk33 = 0.0, 0.0, 0.0

    t11, t12, t13 = 0.0, 0.0, 0.0
    t21, t22, t23 = 0.0, 0.0, 0.0
    t31, t32, t33 = 0.0, 0.0, 0.0

    @inbounds for iID in nodes
        F11 = deformation_gradient[iID, 1, 1]
        F12 = deformation_gradient[iID, 1, 2]
        F13 = deformation_gradient[iID, 1, 3]
        F21 = deformation_gradient[iID, 2, 1]
        F22 = deformation_gradient[iID, 2, 2]
        F23 = deformation_gradient[iID, 2, 3]
        F31 = deformation_gradient[iID, 3, 1]
        F32 = deformation_gradient[iID, 3, 2]
        F33 = deformation_gradient[iID, 3, 3]

        sigma11 = stress_NP1[iID, 1, 1]
        sigma12 = stress_NP1[iID, 1, 2]
        sigma13 = stress_NP1[iID, 1, 3]
        sigma21 = stress_NP1[iID, 2, 1]
        sigma22 = stress_NP1[iID, 2, 2]
        sigma23 = stress_NP1[iID, 2, 3]
        sigma31 = stress_NP1[iID, 3, 1]
        sigma32 = stress_NP1[iID, 3, 2]
        sigma33 = stress_NP1[iID, 3, 3]

        detF = F11 * (F22 * F33 - F23 * F32) -
               F12 * (F21 * F33 - F23 * F31) +
               F13 * (F21 * F32 - F22 * F31)

        if abs(detF) > 1e-15
            inv_detF = 1.0 / detF

            Finv11 = (F22 * F33 - F23 * F32) * inv_detF
            Finv12 = (F13 * F32 - F12 * F33) * inv_detF
            Finv13 = (F12 * F23 - F13 * F22) * inv_detF
            Finv21 = (F23 * F31 - F21 * F33) * inv_detF
            Finv22 = (F11 * F33 - F13 * F31) * inv_detF
            Finv23 = (F13 * F21 - F11 * F23) * inv_detF
            Finv31 = (F21 * F32 - F22 * F31) * inv_detF
            Finv32 = (F12 * F31 - F11 * F32) * inv_detF
            Finv33 = (F11 * F22 - F12 * F21) * inv_detF

            pk11 = sigma11 * Finv11 + sigma12 * Finv21 + sigma13 * Finv31
            pk12 = sigma11 * Finv12 + sigma12 * Finv22 + sigma13 * Finv32
            pk13 = sigma11 * Finv13 + sigma12 * Finv23 + sigma13 * Finv33

            pk21 = sigma21 * Finv11 + sigma22 * Finv21 + sigma23 * Finv31
            pk22 = sigma21 * Finv12 + sigma22 * Finv22 + sigma23 * Finv32
            pk23 = sigma21 * Finv13 + sigma22 * Finv23 + sigma23 * Finv33

            pk31 = sigma31 * Finv11 + sigma32 * Finv21 + sigma33 * Finv31
            pk32 = sigma31 * Finv12 + sigma32 * Finv22 + sigma33 * Finv32
            pk33 = sigma31 * Finv13 + sigma32 * Finv23 + sigma33 * Finv33
        else
            pk11 = pk12 = pk13 = 0.0
            pk21 = pk22 = pk23 = 0.0
            pk31 = pk32 = pk33 = 0.0
        end

        # temp = pk_stress * K^(-1) (inline)
        K11 = inverse_shape_tensor[iID, 1, 1]
        K12 = inverse_shape_tensor[iID, 1, 2]
        K13 = inverse_shape_tensor[iID, 1, 3]
        K21 = inverse_shape_tensor[iID, 2, 1]
        K22 = inverse_shape_tensor[iID, 2, 2]
        K23 = inverse_shape_tensor[iID, 2, 3]
        K31 = inverse_shape_tensor[iID, 3, 1]
        K32 = inverse_shape_tensor[iID, 3, 2]
        K33 = inverse_shape_tensor[iID, 3, 3]

        t11 = pk11 * K11 + pk12 * K21 + pk13 * K31
        t12 = pk11 * K12 + pk12 * K22 + pk13 * K32
        t13 = pk11 * K13 + pk12 * K23 + pk13 * K33

        t21 = pk21 * K11 + pk22 * K21 + pk23 * K31
        t22 = pk21 * K12 + pk22 * K22 + pk23 * K32
        t23 = pk21 * K13 + pk22 * K23 + pk23 * K33

        t31 = pk31 * K11 + pk32 * K21 + pk33 * K31
        t32 = pk31 * K12 + pk32 * K22 + pk33 * K32
        t33 = pk31 * K13 + pk32 * K23 + pk33 * K33

        # Bond forces (inline)
        @fastmath for jID in eachindex(bond_damage[iID])
            bd = bond_damage[iID][jID] * omega[iID][jID]
            xi1 = undeformed_bond[iID][jID][1]
            xi2 = undeformed_bond[iID][jID][2]
            xi3 = undeformed_bond[iID][jID][3]

            bond_force[iID][jID][1] = bd * (t11 * xi1 + t12 * xi2 + t13 * xi3)
            bond_force[iID][jID][2] = bd * (t21 * xi1 + t22 * xi2 + t23 * xi3)
            bond_force[iID][jID][3] = bd * (t31 * xi1 + t32 * xi2 + t33 * xi3)
        end
    end

    return nothing
end

end
