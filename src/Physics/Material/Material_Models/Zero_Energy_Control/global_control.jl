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

"""
    control_name()

Returns the name of the zero energy control
"""
function control_name()
    return "Global"
end

"""
    compute_control(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, material_parameter, time::Float64, dt::Float64)

Computes the zero energy control

# Arguments
- `datamanager::Module`: The datamanager
- `nodes::Union{SubArray,Vector{Int64}}`: The nodes
- `material_parameter`: The material parameter
- `time::Float64`: The current time
- `dt::Float64`: The current time step
# Returns
- `datamanager::Module`: The datamanager
"""
function compute_control(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, material_parameter, time::Float64, dt::Float64)
    dof = datamanager.get_dof()
    deformation_gradient = datamanager.get_field("Deformation Gradient")
    bond_force = datamanager.create_constant_bond_field("Bond Forces", Float64, dof)
    undeformed_bond = datamanager.get_field("Bond Geometry")
    deformed_bond = datamanager.get_field("Deformed Bond Geometry", "NP1")
    Kinv = datamanager.get_field("Inverse Shape Tensor")
    zStiff = datamanager.create_constant_node_field("Zero Energy Stiffness", Float64, "Matrix", dof)
    rotation::Bool, angles = datamanager.rotation_data()
    CVoigt = get_Hooke_matrix(material_parameter, material_parameter["Symmetry"], dof)
    zStiff = create_zero_energy_mode_stiffness(nodes, dof, CVoigt, angles, Kinv, zStiff, rotation)
    bond_force = get_zero_energy_mode_force(nodes, zStiff, deformation_gradient, undeformed_bond, deformed_bond, bond_force)
    return datamanager
end

"""
    get_zero_energy_mode_force(nodes::Union{SubArray,Vector{Int64}}, zStiff::SubArray, deformation_gradient::SubArray, undeformed_bond::SubArray, deformed_bond::SubArray, bond_force::SubArray)

Computes the zero energy mode force

# Arguments
- `nodes::Union{SubArray,Vector{Int64}}`: The nodes
- `zStiff::SubArray`: The zero energy stiffness
- `deformation_gradient::SubArray`: The deformation gradient
- `undeformed_bond::SubArray`: The bond geometry
- `deformed_bond::SubArray`: The bond geometry at the next time step
- `bond_force::SubArray`: The bond force
# Returns
- `bond_force::SubArray`: The bond force
"""
function get_zero_energy_mode_force(nodes::Union{SubArray,Vector{Int64}}, zStiff::SubArray, deformation_gradient::SubArray, undeformed_bond::SubArray, deformed_bond::SubArray, bond_force::SubArray)
    for iID in nodes
        bond_force[iID][:, :] -= (undeformed_bond[iID][:, 1:end-1] * deformation_gradient[iID, :, :] - deformed_bond[iID][:, 1:end-1]) * zStiff[iID, :, :]
    end
    return bond_force
end

"""
    create_zero_energy_mode_stiffness(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, CVoigt, angles, Kinv, zStiff, rotation::Bool)

Creates the zero energy mode stiffness

# Arguments
- `nodes::Union{SubArray,Vector{Int64}}`: The nodes
- `dof::Int64`: The degree of freedom
- `CVoigt::Matrix{Float64}`: The Voigt matrix
- `angles::Vector{Float64}`: The angles
- `Kinv::SubArray`: The inverse shape tensor
- `zStiff::SubArray`: The zero energy stiffness
- `rotation::Bool`: The rotation
# Returns
- `zStiff::SubArray`: The zero energy stiffness
"""
function create_zero_energy_mode_stiffness(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, CVoigt, angles, Kinv, zStiff, rotation::Bool)
    if rotation
        return create_zero_energy_mode_stiffness(nodes, dof, CVoigt, angles, Kinv, zStiff)
    end
    return create_zero_energy_mode_stiffness(nodes, dof, CVoigt, Kinv, zStiff)
end

"""
    create_zero_energy_mode_stiffness(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, CVoigt, Kinv, zStiff)

Creates the zero energy mode stiffness

# Arguments
- `nodes::Union{SubArray,Vector{Int64}}`: The nodes
- `dof::Int64`: The degree of freedom
- `CVoigt::Matrix{Float64}`: The Voigt matrix
- `Kinv::SubArray`: The inverse shape tensor
- `zStiff::SubArray`: The zero energy stiffness
# Returns
- `zStiff::SubArray`: The zero energy stiffness
"""
function create_zero_energy_mode_stiffness(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, CVoigt, Kinv, zStiff)
    C = Array{Float64,4}(get_fourth_order(CVoigt, dof))
    for iID in nodes
        @tensor begin
            zStiff[iID, i, j] = C[i, j, k, l] * Kinv[iID, k, l]
        end
    end
    return zStiff
end

"""
    create_zero_energy_mode_stiffness(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, CVoigt::Matrix{Float64}, angles::SubArray, Kinv::SubArray, zStiff::SubArray)

Creates the zero energy mode stiffness

# Arguments
- `nodes::Union{SubArray,Vector{Int64}}`: The nodes
- `dof::Int64`: The degree of freedom
- `CVoigt::Matrix{Float64}`: The Voigt matrix
- `angles::SubArray`: The angles
- `Kinv::SubArray`: The inverse shape tensor
- `zStiff::SubArray`: The zero energy stiffness
# Returns
- `zStiff::SubArray`: The zero energy stiffness
"""
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

"""
    rotate_fourth_order_tensor(angles::Union{Vector{Float64},Vector{Int64}}, C::Array{Float64,4}, dof::Int64, back::Bool)

Rotates the fourth order tensor

# Arguments
- `angles::Union{Vector{Float64},Vector{Int64}}`: The angles
- `C::Array{Float64,4}`: The fourth order tensor
- `dof::Int64`: The degree of freedom
- `back::Bool`: The back
# Returns
- `C::Array{Float64,4}`: The fourth order tensor
"""
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