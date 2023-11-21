# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module global_zero_energy_control
include("../../../../Support/helpers.jl")
include("../../../../Support/geometry.jl")
include("../../material_basis.jl")
using TensorOperations
using .Geometry
export control_name
export compute_control

function control_name()
    return "Global"
end

function compute_control(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, material_parameter, time::Float64, dt::Float64)
    dof = datamanager.get_dof()
    deformation_gradient = datamanager.get_field("Deformation Gradient")
    bond_force = datamanager.create_constant_bond_field("Bond Forces", Float64, dof)
    bond_geometry = datamanager.get_field("Bond Geometry")
    bond_geometryNP1 = datamanager.get_field("Deformed Bond Geometry", "NP1")
    Kinv = datamanager.get_field("Inverse Shape Tensor")
    zStiff = datamanager.create_constant_node_field("Zero Energy Stiffness", Float64, "Matrix", dof)
    rotation::Bool, angles = datamanager.rotation_data()
    CVoigt = get_Hooke_matrix(material_parameter, material_parameter["Symmetry"], dof)
    zStiff = create_zero_energy_mode_stiffness(nodes, dof, CVoigt, angles, Kinv, zStiff, rotation)
    bond_force = get_zero_energy_mode_force(nodes, zStiff, deformation_gradient, bond_geometry, bond_geometryNP1, bond_force)
    return datamanager
end

function get_zero_energy_mode_force(nodes::Union{SubArray,Vector{Int64}}, zStiff::SubArray, deformation_gradient::SubArray, bond_geometry::SubArray, bond_geometryNP1::SubArray, bond_force::SubArray)
    for iID in nodes
        bond_force[iID][:, :] -= (bond_geometry[iID][:, 1:end-1] * deformation_gradient[iID, :, :] - bond_geometryNP1[iID][:, 1:end-1]) * zStiff[iID, :, :]
    end
    return bond_force
end

function create_zero_energy_mode_stiffness(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, CVoigt, angles, Kinv, zStiff, rotation::Bool)
    if rotation
        return create_zero_energy_mode_stiffness(nodes, dof, CVoigt, angles, Kinv, zStiff)
    end
    return create_zero_energy_mode_stiffness(nodes, dof, CVoigt, Kinv, zStiff)
end

function create_zero_energy_mode_stiffness(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, CVoigt, Kinv, zStiff)
    C = Array{Float64,4}(get_fourth_order(CVoigt, dof))
    for iID in nodes
        @tensor begin
            zStiff[iID, i, j] = C[i, j, k, l] * Kinv[iID, k, l]
        end
    end
    return zStiff
end

function create_zero_energy_mode_stiffness(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, CVoigt::Matrix{Float64}, angles::SubArray, Kinv::SubArray, zStiff::SubArray)
    C = Array{Float64,4}(get_fourth_order(CVoigt, dof))
    for iID in nodes
        C = rotate_fourth_order_tensor(angles[iID, :], C, dof, false)
        @tensor begin
            zStiff[iID, i, j] = C[i, j, k, l] * Kinv[iID, k, l]
        end
    end
    return zStiff
end

function rotate_fourth_order_tensor(angles::Union{Vector{Float64},Vector{Int64}}, C::Array{Float64,4}, dof::Int64, back::Bool)
    rot = Geometry.rotation_tensor(angles)
    R = rot[1:dof, 1:dof]
    if back
        @tensor begin
            C[m, n, o, p] = C[i, j, k, l] * R[m, i] * R[n, j] * R[o, k] * R[p, l]
        end
        return C
    end
    @tensor begin
        C[m, n, o, p] = C[i, j, k, l] * R[i, m] * R[j, n] * R[k, o] * R[l, p]
    end
    return C
end



end