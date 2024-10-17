# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Global_zero_energy_control
include("../../../../Support/helpers.jl")
using .Helpers: get_fourth_order
using StaticArrays: MMatrix, MVector
using LoopVectorization
using LinearAlgebra: mul!
include("../../../../Support/geometry.jl")
include("../../material_basis.jl")
using .Material_Basis: get_Hooke_matrix
using TensorOperations
using .Geometry
export control_name
export compute_control
export global_zero_energy_mode_stiffness

"""
    control_name()

Returns the name of the zero energy control
"""
function control_name()
    return "Global"
end

"""
    compute_control(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, material_parameter::Dict, time::Float64, dt::Float64)

Computes the zero energy control

# Arguments
- `datamanager::Module`: The datamanager
- `nodes::Union{SubArray,Vector{Int64}}`: The nodes
- `material_parameter::Dict`: The material parameter
- `time::Float64`: The current time
- `dt::Float64`: The current time step
# Returns
- `datamanager::Module`: The datamanager
"""
function compute_control(
    datamanager::Module,
    nodes::Union{SubArray,Vector{Int64}},
    material_parameter::Dict,
    time::Float64,
    dt::Float64,
)
    dof = datamanager.get_dof()
    deformation_gradient = datamanager.get_field("Deformation Gradient")
    bond_force = datamanager.create_constant_bond_field("Bond Forces", Float64, dof)
    undeformed_bond = datamanager.get_field("Bond Geometry")
    deformed_bond = datamanager.get_field("Deformed Bond Geometry", "NP1")
    Kinv = datamanager.get_field("Inverse Shape Tensor")
    zStiff = datamanager.create_constant_node_field(
        "Zero Energy Stiffness",
        Float64,
        "Matrix",
        dof,
    )
    rotation::Bool = datamanager.get_rotation()
    angles = datamanager.get_field("Angles")
    CVoigt = get_Hooke_matrix(material_parameter, material_parameter["Symmetry"], dof)
    if !haskey(material_parameter, "UMAT Material Name")
        zStiff = create_zero_energy_mode_stiffness(
            nodes,
            dof,
            CVoigt,
            angles,
            Kinv,
            zStiff,
            rotation,
        )
    end

    if dof == 2
        get_zero_energy_mode_force_2d!(
            nodes,
            zStiff,
            deformation_gradient,
            undeformed_bond,
            deformed_bond,
            bond_force,
        )
    elseif dof == 3
        get_zero_energy_mode_force_3d!(
            nodes,
            zStiff,
            deformation_gradient,
            undeformed_bond,
            deformed_bond,
            bond_force,
        )
    end
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
function get_zero_energy_mode_force_2d!(
    nodes::Union{SubArray,Vector{Int64}},
    zStiff,
    deformation_gradient,
    undeformed_bond,
    deformed_bond,
    bond_force,
)
    # Working code
    # bond_force[iID][:, :] = (-zStiff[iID, :, :] * (deformation_gradient[iID, :, :] * (undeformed_bond[iID])' - (deformed_bond[iID])'))'
    @inbounds @fastmath for iID in nodes
        @inbounds @fastmath @views for nID ∈ axes(undeformed_bond[iID], 1)
            df = MVector{2}(zeros(2))
            @inbounds @fastmath @views for m ∈ axes(deformation_gradient, 2)
                df_i = zero(eltype(deformation_gradient))
                @inbounds @fastmath @views for n ∈ axes(deformation_gradient, 3)
                    df_i += deformation_gradient[iID, m, n] * undeformed_bond[iID][nID, n]
                end
                df[m] = df_i - deformed_bond[iID][nID, m]
            end
            @inbounds @fastmath @views for m ∈ axes(zStiff, 2)
                bf_i = zero(eltype(zStiff))
                @inbounds @fastmath @views for n ∈ axes(zStiff, 3)
                    bf_i -= zStiff[iID, m, n] * df[n]
                end
                bond_force[iID][nID, m] = bf_i
            end
        end

    end
end
function get_zero_energy_mode_force_3d!(
    nodes::Union{SubArray,Vector{Int64}},
    zStiff,
    deformation_gradient,
    undeformed_bond,
    deformed_bond,
    bond_force,
)
    # Working code
    # bond_force[iID][:, :] = (-zStiff[iID, :, :] * (deformation_gradient[iID, :, :] * (undeformed_bond[iID])' - (deformed_bond[iID])'))'
    @inbounds @fastmath for iID in nodes
        @inbounds @fastmath @views for nID ∈ axes(undeformed_bond[iID], 1)
            df = MVector{3}(zeros(3))
            @inbounds @fastmath @views for m ∈ axes(deformation_gradient, 2)
                df_i = zero(eltype(deformation_gradient))
                @inbounds @fastmath @views for n ∈ axes(deformation_gradient, 3)
                    df_i += deformation_gradient[iID, m, n] * undeformed_bond[iID][nID, n]
                end
                df[m] = df_i - deformed_bond[iID][nID, m]
            end
            @inbounds @fastmath @views for m ∈ axes(zStiff, 2)
                bf_i = zero(eltype(zStiff))
                @inbounds @fastmath @views for n ∈ axes(zStiff, 3)
                    bf_i -= zStiff[iID, m, n] * df[n]
                end
                bond_force[iID][nID, m] = bf_i
            end
        end

    end
end

"""
    create_zero_energy_mode_stiffness(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, CVoigt::Union{MMatrix,Matrix{Float64}}, angles::Vector{Float64}, Kinv::Array{Float64, 3}, zStiff::Array{Float64, 3}, rotation::Bool)

Creates the zero energy mode stiffness

# Arguments
- `nodes::Union{SubArray,Vector{Int64}}`: The nodes
- `dof::Int64`: The degree of freedom
- `CVoigt::Union{MMatrix, Matrix{Float64}}`: The Voigt matrix
- `angles::Vector{Float64}`: The angles
- `Kinv::Array{Float64, 3}`: The inverse shape tensor
- `zStiff::Array{Float64, 3}`: The zero energy stiffness
- `rotation::Bool`: The rotation
# Returns
- `zStiff::SubArray`: The zero energy stiffness
"""
function create_zero_energy_mode_stiffness(
    nodes::Union{SubArray,Vector{Int64}},
    dof::Int64,
    CVoigt::Union{MMatrix,Matrix{Float64}},
    angles::Union{Matrix{Float64},Nothing},
    Kinv::Array{Float64,3},
    zStiff::Array{Float64,3},
    rotation::Bool,
)
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
- `CVoigt::Union{MMatrix, Matrix{Float64}}`: The Voigt matrix
- `Kinv::Array{Float64, 3}`: The inverse shape tensor
- `zStiff::Array{Float64, 3}`: The zero energy stiffness
# Returns
- `zStiff::SubArray`: The zero energy stiffness
"""
function create_zero_energy_mode_stiffness(
    nodes::Union{SubArray,Vector{Int64}},
    dof::Int64,
    CVoigt::Union{MMatrix,Matrix{Float64}},
    Kinv::Array{Float64,3},
    zStiff::Array{Float64,3},
)

    for iID in nodes
        global_zero_energy_mode_stiffness(iID, get_fourth_order(CVoigt, dof), Kinv, zStiff)
    end
    return zStiff
end

"""
global_zero_energy_mode_stiffness(id::Int64, C, Kinv, zstiff)

Creates the zero energy mode stiffness, based on the UMAT interface

# Arguments
- `ID::Int64` : ID of the node
- `C::Union{MMatrix, Matrix{Float64}}`: The fourth order elasticity tensor
- `Kinv::Array{Float64, 3}`: The inverse shape tensor
- `zStiff::Array{Float64, 3}`: The zero energy stiffness
# Returns
- `zStiff::SubArray`: The zero energy stiffness
"""




function global_zero_energy_mode_stiffness(
    ID::Int64,
    C,
    Kinv::Array{Float64,3},
    zStiff::Array{Float64,3},
)

    # Perform matrix multiplication for each i, j
    @views @inbounds @fastmath for i ∈ axes(zStiff, 2), j ∈ axes(zStiff, 3)
        zij = zero(eltype(zStiff))
        @views @inbounds @fastmath for k ∈ axes(Kinv, 2), l ∈ axes(Kinv, 3)
            @views zij += C[i, j, k, l] * Kinv[ID, k, l]
        end
        zStiff[ID, i, j] = zij
    end

end


"""
    create_zero_energy_mode_stiffness(nodes::Union{SubArray,Vector{Int64}}, dof::Int64, CVoigt::Union{MMatrix,Matrix{Float64}}, angles::Matrix{Float64}, Kinv::Array{Float64, 3}, zStiff::Array{Float64, 3})

Creates the zero energy mode stiffness

# Arguments
- `nodes::Union{SubArray,Vector{Int64}}`: The nodes
- `dof::Int64`: The degree of freedom
- `CVoigt::Union{MMatrix, Matrix{Float64}}`: The Voigt matrix
- `angles::Matrix{Float64}`: The angles
- `Kinv::Array{Float64, 3}`: The inverse shape tensor
- `zStiff::Array{Float64, 3}`: The zero energy stiffness
# Returns
- `zStiff::SubArray`: The zero energy stiffness
"""
function create_zero_energy_mode_stiffness(
    nodes::Union{SubArray,Vector{Int64}},
    dof::Int64,
    CVoigt::Union{MMatrix,Matrix{Float64}},
    angles::Matrix{Float64},
    Kinv::Array{Float64,3},
    zStiff::Array{Float64,3},
)
    C = Array{Float64,4}(get_fourth_order(CVoigt, dof))
    for iID in nodes
        C = rotate_fourth_order_tensor(angles[iID, :], C, dof, false)
        global_zero_energy_mode_stiffness(iID, C, Kinv, zStiff)
    end

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
function rotate_fourth_order_tensor(
    angles::Union{Vector{Float64},Vector{Int64}},
    C::Array{Float64,4},
    dof::Int64,
    back::Bool,
)
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
