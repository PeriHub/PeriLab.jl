# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Global_Zero_Energy_Control

using StaticArrays: MMatrix, MVector, SMatrix, SVector
using LoopVectorization
using ....Data_Manager
using ....Helpers: get_fourth_order
using ...Material_Basis: hooke_matrix
using ....Geometry: rotation_tensor
using .....ParameterSpec: model_parts

export init_model
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

# the legacy dict marked UMAT materials by the key "UMAT Material Name"
function is_umat(material)
    any(part -> nameof(typeof(part)) === :CorrespondenceUMATParams,
        model_parts(material.model))
end

"""
    init_model(nodes, material)

Initializes the global zero energy control fields.

# Arguments
- `nodes::AbstractVector{Int64}`: The block nodes.
- `material::BlockMaterial`: The typed block material (base, moduli, symmetry).
"""
function init_model(nodes::AbstractVector{Int64}, material)
    dof::Int64 = Data_Manager.get_dof()
    Data_Manager.create_constant_node_tensor_field("Zero Energy Stiffness", Float64, dof)
    "Material Gradient" in Data_Manager.get_all_field_keys() && return
    hooke::NodeTensorField{Float64} = Data_Manager.create_constant_node_tensor_field("Material Gradient",
                                                                                     Float64,
                                                                                     Int64((dof *
                                                                                            (dof +
                                                                                             1)) /
                                                                                           2))
    for iID in nodes
        @views hooke[iID, :, :] = hooke_matrix(material, dof, iID)
    end
end

"""
	compute_control(nodes::AbstractVector{Int64}, material, time::Float64, dt::Float64)

Computes the zero energy control

# Arguments
- `nodes::AbstractVector{Int64}`: The nodes
- `material::BlockMaterial`: The typed block material
- `time::Float64`: The current time
- `dt::Float64`: The current time step


Global - J. Wan et al., "Improved method for zero-energy mode suppression in peridynamic correspondence model in Acta Mechanica Sinica https://doi.org/10.1007/s10409-019-00873-y

"""

function _compute_control!(nodes::AbstractVector{Int64},
                           apply_zero_energy::Bool,
                           time::Float64,
                           dt::Float64)
    # function barrier on the dof: the fields below get concrete types
    if Data_Manager.get_dof() == 2
        _compute_control!(nodes, apply_zero_energy, Val(2))
    else
        _compute_control!(nodes, apply_zero_energy, Val(3))
    end
end

function _compute_control!(nodes::AbstractVector{Int64}, apply_zero_energy::Bool,
                           ::Val{D}) where {D}
    deformation_gradient = Data_Manager.get_field("Deformation Gradient")::NodeTensorField{Float64}
    bond_force = Data_Manager.get_field("Bond Forces")::BondVectorState{Float64,D}
    undeformed_bond = Data_Manager.get_field("Bond Geometry")::BondVectorState{Float64,D}
    deformed_bond = Data_Manager.get_field("Deformed Bond Geometry",
                                           "NP1")::BondVectorState{Float64,D}
    omega = Data_Manager.get_field("Influence Function")::BondScalarState{Float64}
    Kinv = Data_Manager.get_field("Inverse Shape Tensor")::NodeTensorField{Float64}
    hooke_matrix = Data_Manager.get_field("Material Gradient")::NodeTensorField{Float64}
    zStiff = Data_Manager.get_field("Zero Energy Stiffness")::NodeTensorField{Float64}

    if apply_zero_energy
        if Data_Manager.get_rotation()
            angles = Data_Manager.get_field("Angles")::NodeVectorField{Float64}
            create_zero_energy_mode_stiffness!(nodes, D, hooke_matrix, angles, Kinv,
                                               zStiff)
        else
            zero_energy_mode_stiffness!(nodes, hooke_matrix, Kinv, zStiff, Val(D))
        end
    end
    zero_energy_mode_forces!(nodes, zStiff, deformation_gradient, undeformed_bond,
                             deformed_bond, omega, bond_force)
end

function compute_control(nodes::AbstractVector{Int64}, material, time::Float64, dt::Float64)
    _compute_control!(nodes,
                      !is_umat(material),
                      time,
                      dt)
end

"""
    zero_energy_mode_forces!(nodes, zStiff, deformation_gradient, undeformed_bond, deformed_bond, omega, bond_force)

Adds the zero energy mode forces `-ω Z (F X - x)` to the bond forces.

# Arguments
- `nodes::AbstractVector{Int64}`: The nodes
- `zStiff::NodeTensorField{Float64}`: The zero energy stiffness
- `deformation_gradient::NodeTensorField{Float64}`: The deformation gradient
- `undeformed_bond::BondVectorState{Float64,D}`: The bond geometry
- `deformed_bond::BondVectorState{Float64,D}`: The bond geometry at the next time step
- `omega::BondScalarState{Float64}`: The influence function
- `bond_force::BondVectorState{Float64,D}`: The bond force
"""
function zero_energy_mode_forces!(nodes::AbstractVector{Int64},
                                  zStiff::NodeTensorField{Float64},
                                  deformation_gradient::NodeTensorField{Float64},
                                  undeformed_bond::BondVectorState{Float64,D},
                                  deformed_bond::BondVectorState{Float64,D},
                                  omega::BondScalarState{Float64},
                                  bond_force::BondVectorState{Float64,D}) where {D}
    @inbounds for iID in nodes
        F = node_tensor(deformation_gradient, iID, Val(D))
        Z = node_tensor(zStiff, iID, Val(D))
        X = undeformed_bond[iID]
        x = deformed_bond[iID]
        w = omega[iID]
        f = bond_force[iID]
        @fastmath for jID in eachindex(X)
            f[jID] -= w[jID] * (Z * (F * X[jID] - x[jID]))
        end
    end
end

# the D×D tensor of node `iID` of a node tensor field
@inline function node_tensor(A::NodeTensorField{Float64}, iID::Int64, ::Val{D}) where {D}
    return SMatrix{D,D,Float64,D * D}(ntuple(k -> @inbounds(A[iID, (k - 1) % D + 1,
                                                              (k - 1) ÷ D + 1]),
                                             Val(D * D)))
end

"""
	create_zero_energy_mode_stiffness!(nodes::AbstractVector{Int64}, dof::Int64, CVoigt, Kinv, zStiff)

Creates the zero energy mode stiffness

# Arguments
- `nodes::AbstractVector{Int64}`: The nodes
- `dof::Int64`: The degree of freedom
- `CVoigt::Union{MMatrix, Matrix{Float64}}`: The Voigt matrix
- `Kinv::Array{Float64, 3}`: The inverse shape tensor
- `zStiff::Array{Float64, 3}`: The zero energy stiffness
# Returns
- `zStiff::SubArray`: The zero energy stiffness
"""

function create_zero_energy_mode_stiffness!(nodes::AbstractVector{Int64},
                                            dof::Int64,
                                            CVoigt::AbstractArray{Float64,3},
                                            Kinv::Array{Float64,3},
                                            zStiff::Array{Float64,3})
    if dof == 2
        zero_energy_mode_stiffness!(nodes, CVoigt, Kinv, zStiff, Val(2))
    else
        zero_energy_mode_stiffness!(nodes, CVoigt, Kinv, zStiff, Val(3))
    end
end

# Z_ij = C_ijkl Kinv_kl, read directly from the Voigt matrix (no fourth order tensor);
# like `get_fourth_order`, 3D uses the upper triangle of the Voigt matrix
function zero_energy_mode_stiffness!(nodes::AbstractVector{Int64},
                                     CVoigt::AbstractArray{Float64,3},
                                     Kinv::AbstractArray{Float64,3},
                                     zStiff::AbstractArray{Float64,3}, ::Val{D}) where {D}
    @inbounds for iID in nodes, j in 1:D, i in 1:D
        a = voigt_index(i, j, Val(D))
        zij = 0.0
        for l in 1:D, k in 1:D
            zij += voigt_entry(CVoigt, iID, a, voigt_index(k, l, Val(D)), Val(D)) *
                   Kinv[iID, k, l]
        end
        zStiff[iID, i, j] = zij
    end
end
@inline voigt_index(i, j, ::Val{2}) = i == j ? i : 3
@inline voigt_index(i, j, ::Val{3}) = i == j ? i : 9 - i - j   # 23 → 4, 13 → 5, 12 → 6
@inline voigt_entry(C, n, a, b, ::Val{2}) = @inbounds C[n, a, b]
@inline voigt_entry(C, n, a, b, ::Val{3}) = @inbounds C[n, min(a, b), max(a, b)]

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

function global_zero_energy_mode_stiffness(ID::Int64,
                                           C::Array{Float64,4},
                                           Kinv::NodeTensorField{Float64},
                                           zStiff::NodeTensorField{Float64})

    # Perform matrix multiplication for each i, j
    @views @inbounds @fastmath for i in axes(zStiff, 2), j in axes(zStiff, 3)
        zij = zero(eltype(zStiff))
        @views @inbounds @fastmath for k in axes(Kinv, 2), l in axes(Kinv, 3)
            @views zij += C[i, j, k, l] * Kinv[ID, k, l]
        end
        zStiff[ID, i, j] = zij
    end
end

"""
	create_zero_energy_mode_stiffness(nodes::AbstractVector{Int64}, dof::Int64, CVoigt::Union{MMatrix,Matrix{Float64}}, angles::Matrix{Float64}, Kinv::Array{Float64, 3}, zStiff::Array{Float64, 3})

Creates the zero energy mode stiffness

# Arguments
- `nodes::AbstractVector{Int64}`: The nodes
- `dof::Int64`: The degree of freedom
- `CVoigt::Union{MMatrix, Matrix{Float64}}`: The Voigt matrix
- `angles::Matrix{Float64}`: The angles
- `Kinv::Array{Float64, 3}`: The inverse shape tensor
- `zStiff::Array{Float64, 3}`: The zero energy stiffness
# Returns
- `zStiff::SubArray`: The zero energy stiffness
"""
function create_zero_energy_mode_stiffness!(nodes::AbstractVector{Int64},
                                            dof::Int64,
                                            CVoigt::AbstractArray{Float64,3},
                                            angles::AbstractArray{Float64},
                                            Kinv::AbstractArray{Float64,3},
                                            zStiff::AbstractArray{Float64,3})
    for iID in nodes
        @views C = get_fourth_order(CVoigt[iID, :, :], dof)
        @views C = rotate_fourth_order_tensor(angles[iID, :], C, dof, false)
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
- `C::Array{Float64,4}`: The rotated fourth order tensor (a new array; `C` is not changed,
  it may be the shared buffer of `get_fourth_order`)
"""
function rotate_fourth_order_tensor(angles::AbstractVector{T},
                                    C::Array{Float64,4},
                                    dof::Int64,
                                    back::Bool) where {T<:Union{Int64,Float64}}
    rot = rotation_tensor(angles, dof)
    @views R = rot[1:dof, 1:dof]
    return back ? fourth_order_rotation(R', C) : fourth_order_rotation(R, C)
end

function fourth_order_rotation(R, C::Array{Float64,4})
    C_rot = similar(C)
    @inbounds @fastmath for m in axes(C, 1)
        @inbounds @fastmath for n in axes(C, 1)
            @inbounds @fastmath for o in axes(C, 1)
                @inbounds @fastmath for p in axes(C, 1)
                    Cmnop = zero(eltype(C))
                    @inbounds @fastmath for i in axes(C, 1)
                        @inbounds @fastmath for j in axes(C, 1)
                            @inbounds @fastmath for k in axes(C, 1)
                                @inbounds @fastmath for l in axes(C, 1)
                                    Cmnop += C[i, j, k, l] *
                                             R[m, i] *
                                             R[n, j] *
                                             R[o, k] *
                                             R[p, l]
                                end
                            end
                        end
                    end
                    C_rot[m, n, o, p] = Cmnop
                end
            end
        end
    end
    return C_rot
end

end
