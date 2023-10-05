# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module global_zero_energy_control
using TensorOperations

export get_name_control
export compute_control

function get_name_control()
    return "Global"
end

function compute_control(datamanager::Module, nodes::Vector{Int64}, material_parameter, time::Float32, dt::Float32)
    rotation::Bool = false
    if "Angles" in datamanager.get_all_field_keys()
        rotation = true
        angles = datamanager.get_field("Angles")
    end
    dof = datamanager.get_dof()
    defGradNP1 = datamanager.get_field("Deformation Gradient", "NP1")
    bond_force = datamanager.create_constant_bond_field("Bond Forces", Float32, dof)
    bondGeom = datamanager.get_field("Deformed Bond Geometry", "N")
    bondGeomNP1 = datamanager.get_field("Deformed Bond Geometry", "NP1")
    Kinv = datamanager.get_field("Inverse Shape Tensor")
    zStiff = datamanager.create_constant_node_field("Zero Energy Stiffness", Float32, "Matrix", dof)
    CVoigt = get_Hooke_matrix(material_parameter, material_parameter["Symmetry"], dof)
    zstiff = create_zero_energy_mode_stiffness(nodes, dof, CVoigt, angles, Kinv, zStiff, rotation)
    bond_force = get_zero_energy_mode_forde(nodes, zstiff, defGradNP1, bondGeom, bondGeomNP1, bond_force)
    return datamanager
end

function create_zero_energy_mode_stiffness(nodes::Vector{Int64}, dof::Int64, CVoigt, angles, Kinv, zStiff, rotation::Bool)
    if rotation
        return create_zero_energy_mode_stiffness(nodes, dof, CVoigt, angles, Kinv, zStiff)
    end
    return create_zero_energy_mode_stiffness(nodes, dof, CVoigt, Kinv, zStiff)
end

function create_zero_energy_mode_stiffness(nodes::Vector{Int64}, dof::Int64, CVoigt, Kinv, zStiff)
    C = get_fourth_order(CVoigt, dof)
    for iID in nodes
        @tensor begin
            zStiff[iID, i, j] = C[i, j, k, l] * Kinv[iID, k, l]
        end
    end
    return zStiff
end


function create_zero_energy_mode_stiffness(nodes::Vector{Int64}, dof::Int64, CVoigt, angles, Kinv, zStiff)
    C = get_fourth_order(CVoigt, dof)
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


function get_fourth_order(CVoigt, dof)
    C = zeros(dof, dof, dof, dof)
    for i in 1:dof
        for j in j:dof
            C[i, i, j, j] = CVoigt[i, j]#9
        end
    end

    if dof == 2
        C[1, 2, 1, 2] = C[2, 1, 2, 1] = C[1, 2, 2, 1] = C[1, 2, 2, 1] = CVoigt[3, 3]
        C[2, 1, 1, 1] = C[1, 2, 1, 1] = C[1, 1, 2, 1] = C[1, 1, 1, 2] = CVoigt[1, 3]
        C[1, 2, 2, 2] = C[2, 1, 2, 2] = C[2, 2, 1, 2] = C[2, 2, 2, 1] = CVoigt[2, 3]#12
    end
    #----------------------------------
    if dof == 3 # 81
        C[2, 3, 2, 3] = C[2, 3, 2, 3] = CVoigt[4, 4]#2
        C[2, 3, 1, 1] = C[3, 2, 1, 1] = C[1, 1, 2, 3] = C[1, 1, 3, 2] = CVoigt[1, 4]
        C[3, 2, 2, 2] = C[2, 3, 2, 2] = C[2, 2, 3, 2] = C[2, 2, 2, 3] = CVoigt[2, 4]
        C[3, 2, 3, 3] = C[2, 3, 3, 3] = C[3, 3, 3, 2] = C[3, 3, 2, 3] = CVoigt[3, 4]#12
        C[1, 3, 1, 3] = C[3, 1, 3, 1] = CVoigt[5, 5]#2
        C[3, 1, 1, 1] = C[1, 3, 1, 1] = C[1, 1, 3, 1] = C[1, 1, 1, 3] = CVoigt[1, 5]
        C[1, 3, 2, 2] = C[3, 1, 2, 2] = C[2, 2, 1, 3] = C[2, 2, 3, 1] = CVoigt[2, 5]
        C[1, 3, 3, 3] = C[3, 1, 3, 3] = C[3, 3, 1, 3] = C[3, 3, 3, 1] = CVoigt[3, 5]#12
        C[1, 2, 1, 2] = C[2, 1, 2, 1] = CVoigt[6, 6]#2
        C[2, 1, 1, 1] = C[1, 2, 1, 1] = C[1, 1, 2, 1] = C[1, 1, 1, 2] = CVoigt[1, 6]
        C[1, 2, 2, 2] = C[2, 1, 2, 2] = C[2, 2, 1, 2] = C[2, 2, 2, 1] = CVoigt[2, 6]
        C[1, 2, 3, 3] = C[2, 1, 3, 3] = C[3, 3, 1, 2] = C[3, 3, 2, 1] = CVoigt[3, 6]#12
        #----------------------------------
        C[1, 3, 1, 2] = C[1, 2, 1, 3] = C[3, 1, 1, 2] = C[2, 1, 1, 3] = C[1, 2, 1, 3] = CVoigt[5, 6]
        C[3, 1, 2, 1] = C[2, 1, 3, 1] = C[2, 1, 1, 3] = C[3, 1, 1, 2] = C[2, 1, 3, 1] = CVoigt[5, 6]

        C[2, 3, 1, 2] = C[3, 2, 1, 2] = C[3, 2, 2, 1] = C[2, 3, 2, 1] = C[1, 2, 2, 3] = CVoigt[4, 6]
        C[3, 2, 2, 1] = C[2, 3, 2, 1] = C[2, 3, 1, 2] = C[3, 2, 1, 2] = C[2, 1, 3, 2] = CVoigt[4, 6]

        C[2, 3, 1, 3] = C[3, 2, 1, 3] = C[3, 2, 3, 1] = C[2, 3, 3, 1] = C[1, 3, 2, 3] = CVoigt[4, 5]
        C[3, 2, 3, 1] = C[2, 3, 3, 1] = C[2, 3, 1, 3] = C[3, 2, 1, 3] = C[3, 1, 3, 2] = CVoigt[4, 5]#30


    end
    return C
end
end