# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module global_zero_energy_control
using TensorOperations
include("../../../../Support/tools.jl")
include("../../material_basis.jl")
using TensorOperations
export control_name
export compute_control

function control_name()
    return "Global"
end

function compute_control(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, material_parameter, time::Float32, dt::Float32)
    rotation::Bool = false
    if "Angles" in datamanager.get_all_field_keys()
        rotation = true
        angles = datamanager.get_field("Angles")
    else
        angles = nothing
    end
    dof = datamanager.get_dof()
    defGrad = datamanager.get_field("Deformation Gradient")
    bond_force = datamanager.create_constant_bond_field("Bond Forces", Float32, dof)
    bondGeom = datamanager.get_field("Deformed Bond Geometry", "N")
    bondGeomNP1 = datamanager.get_field("Deformed Bond Geometry", "NP1")
    Kinv = datamanager.get_field("Inverse Shape Tensor")
    zStiff = datamanager.create_constant_node_field("Zero Energy Stiffness", Float32, "Matrix", dof)
    CVoigt = get_Hooke_matrix(material_parameter, material_parameter["Symmetry"], dof)
    zStiff = create_zero_energy_mode_stiffness(nodes, dof, CVoigt, angles, Kinv, zStiff, rotation)
    bond_force = get_zero_energy_mode_force(nodes, zStiff, defGrad, bondGeom, bondGeomNP1, bond_force)
    return datamanager
end

function get_zero_energy_mode_force(nodes, zStiff, defGrad, bondGeom, bondGeomNP1, bond_force)
    for iID in nodes
        bond_force[iID][:, :] -= (bondGeom[iID][:, 1:end-1] * defGrad[iID, :, :] - bondGeomNP1[iID][:, 1:end-1]) * zStiff[iID, :, :]
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

function create_zero_energy_mode_stiffness(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, CVoigt, angles, Kinv, zStiff)
    C = Array{Float64,4}(get_fourth_order(CVoigt, dof))
    for iID in nodes
        C = rotate_stiffness_tensor(angles[iID, :, :], dof, C)
        @tensor begin
            zStiff[iID, i, j] = C[i, j, k, l] * Kinv[iID, k, l]
        end
    end
    return zStiff
end


function rotate_fourth_order_tensor(angles, C, dof::Int64, back::Bool)
    rot = Geometry.rotation_tensor(angles)
    R = rot[1:dof, dof]
    if back
        @tensor begin
            C[i, j, k, l] = C[i, j, k, l] * R[m, i] * R[n, j] * R[o, k] * R[p, l]
        end

    else
        @tensor begin
            C[i, j, k, l] = C[i, j, k, l] * R[i, m] * R[j, n] * R[k, o] * R[l, p]
        end
    end
    return C
end



end