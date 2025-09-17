# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Global_zero_energy_control
include("../../../../Support/Helpers.jl")
using .Helpers: get_fourth_order, mul!
using StaticArrays: MMatrix, MVector
using LoopVectorization
include("../../../../Support/Geometry.jl")
include("../../Material_Basis.jl")
using .Material_Basis: get_Hooke_matrix
using .Geometry: rotation_tensor
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
    compute_control(datamanager::Module, nodes::AbstractVector{Int64}, material_parameter::Dict, time::Float64, dt::Float64)

Computes the zero energy control

# Arguments
- `datamanager::Module`: The datamanager
- `nodes::AbstractVector{Int64}`: The nodes
- `material_parameter::Dict`: The material parameter
- `time::Float64`: The current time
- `dt::Float64`: The current time step
# Returns
- `datamanager::Module`: The datamanager
"""
function compute_control(datamanager::Module,
                         nodes::AbstractVector{Int64},
                         material_parameter::Dict{String,Any},
                         time::Float64,
                         dt::Float64)
    dof = datamanager.get_dof()
    deformation_gradient = datamanager.get_field("Deformation Gradient")::Array{Float64,3}
    bond_force = datamanager.get_field("Bond Forces")::Vector{Vector{Vector{Float64}}}
    undeformed_bond = datamanager.get_field("Bond Geometry")::Vector{Vector{Vector{Float64}}}
    deformed_bond = datamanager.get_field("Deformed Bond Geometry",
                                          "NP1")::Vector{Vector{Vector{Float64}}}
    Kinv = datamanager.get_field("Inverse Shape Tensor")::Array{Float64,3}
    zStiff = datamanager.create_constant_node_field("Zero Energy Stiffness",
                                                    Float64,
                                                    dof,
                                                    VectorOrMatrix = "Matrix")::Array{Float64,
                                                                                      3}
    rotation = datamanager.get_rotation()::Bool

    symmetry = material_parameter["Symmetry"]::String
    CVoigt = get_Hooke_matrix(datamanager,
                              material_parameter,
                              symmetry,
                              dof)
    if !haskey(material_parameter, "UMAT Material Name")
        if rotation
            angles = datamanager.get_field("Angles")
            create_zero_energy_mode_stiffness!(nodes, Val(dof), CVoigt, angles, Kinv,
                                               zStiff)
        else
            create_zero_energy_mode_stiffness!(nodes, Val(dof), CVoigt, Kinv, zStiff)
        end
    end

    if dof == 2
        get_zero_energy_mode_force_2d!(nodes,
                                       zStiff,
                                       deformation_gradient,
                                       undeformed_bond,
                                       deformed_bond,
                                       bond_force)
    elseif dof == 3
        get_zero_energy_mode_force_3d!(nodes,
                                       zStiff,
                                       deformation_gradient,
                                       undeformed_bond,
                                       deformed_bond,
                                       bond_force)
    end
    return datamanager::Module
end

"""
    get_zero_energy_mode_force(nodes::AbstractVector{Int64}, zStiff::SubArray, deformation_gradient::SubArray, undeformed_bond::SubArray, deformed_bond::SubArray, bond_force::SubArray)

Computes the zero energy mode force

# Arguments
- `nodes::AbstractVector{Int64}`: The nodes
- `zStiff::SubArray`: The zero energy stiffness
- `deformation_gradient::SubArray`: The deformation gradient
- `undeformed_bond::SubArray`: The bond geometry
- `deformed_bond::SubArray`: The bond geometry at the next time step
- `bond_force::SubArray`: The bond force
# Returns
- `bond_force::SubArray`: The bond force
"""
function get_zero_energy_mode_force_2d!(nodes::AbstractVector{Int64},
                                        zStiff::Array{Float64,3},
                                        deformation_gradient::Array{Float64,3},
                                        undeformed_bond::Vector{Vector{Vector{Float64}}},
                                        deformed_bond::Vector{Vector{Vector{Float64}}},
                                        bond_force::Vector{Vector{Vector{Float64}}})
    df = MVector{2}(zeros(Float64, 2))
    @inbounds @fastmath for iID in nodes
        @inbounds @fastmath @views for nID in axes(undeformed_bond[iID], 1)
            @inbounds @fastmath @views for m in axes(deformation_gradient, 2)
                df_i = zero(eltype(deformation_gradient))
                @inbounds @fastmath @views for n in axes(deformation_gradient, 3)
                    df_i += deformation_gradient[iID, m, n] * undeformed_bond[iID][nID][n]
                end
                df[m] = df_i - deformed_bond[iID][nID][m]
            end
            @views bond_force_computation_2d!(zStiff[iID, :, :], df,
                                              bond_force[iID][nID])
        end
    end
end
function bond_force_computation_2d!(zStiff::AbstractArray{Float64}, df::MVector{2},
                                    bond_force::Vector{Float64})
    @inbounds @fastmath @views for m in axes(zStiff, 1)
        bf_i = zero(eltype(zStiff))
        @inbounds @fastmath @views for n in axes(zStiff, 2)
            bf_i -= zStiff[m, n] * df[n]
        end
        bond_force[m] += bf_i
    end
end
function bond_force_computation_3d!(zStiff::AbstractArray{Float64}, df::MVector{3},
                                    bond_force::Vector{Float64})
    @inbounds @fastmath @views for m in axes(zStiff, 1)
        bf_i = zero(eltype(zStiff))
        @inbounds @fastmath @views for n in axes(zStiff, 2)
            bf_i -= zStiff[m, n] * df[n]
        end
        bond_force[m] += bf_i
    end
end

function get_zero_energy_mode_force_3d!(nodes::AbstractVector{Int64},
                                        zStiff::AbstractArray{Float64},
                                        deformation_gradient::Array{Float64,3},
                                        undeformed_bond::Vector{Vector{Vector{Float64}}},
                                        deformed_bond::Vector{Vector{Vector{Float64}}},
                                        bond_force::Vector{Vector{Vector{Float64}}})
    df = MVector{3}(zeros(Float64, 3))
    @inbounds @fastmath for iID in nodes
        @inbounds @fastmath @views for nID in axes(undeformed_bond[iID], 1)
            @inbounds @fastmath @views for m in axes(deformation_gradient, 2)
                df_i = zero(eltype(deformation_gradient))
                @inbounds @fastmath @views for n in axes(deformation_gradient, 3)
                    df_i += deformation_gradient[iID, m, n] * undeformed_bond[iID][nID][n]
                end
                df[m] = df_i - deformed_bond[iID][nID][m]
            end
            @views bond_force_computation_3d!(zStiff[iID, :, :], df,
                                              bond_force[iID][nID])
        end
    end
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
                                            ::Val{DOF},
                                            CVoigt::MMatrix{N,N,Float64,N2},
                                            Kinv::Array{Float64,3},
                                            zStiff::Array{Float64,3}) where {DOF,N,N2}
    C = get_fourth_order(CVoigt, Val(DOF))  # construct once, if it's always same!
    for iID in nodes
        global_zero_energy_mode_stiffness(iID, C, Kinv, zStiff)
    end
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

function global_zero_energy_mode_stiffness(ID::Int64,
                                           C::AbstractArray{Float64,4},
                                           Kinv::Array{Float64,3},
                                           zStiff::Array{Float64,3})

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
                                            ::Val{DOF},
                                            CVoigt::MMatrix{N,N,Float64,N2},
                                            angles::AbstractArray{Float64},
                                            Kinv::AbstractArray{Float64,3},
                                            zStiff::AbstractArray{Float64,3}) where {DOF,N,
                                                                                     N2}
    C = get_fourth_order(CVoigt, Val{DOF})
    for iID in nodes
        C = rotate_fourth_order_tensor(angles[iID, :], C, Val{DOF}, false)
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
function rotate_fourth_order_tensor(angles::Vector{T},
                                    C::Array{Float64,4},
                                    dof::Int64,
                                    back::Bool) where {T<:Union{Int64,Float64}}
    rot = rotation_tensor(angles, dof)
    @views R = rot[1:dof, 1:dof]
    if back
        fourth_order_rotation(R', C)
    else
        fourth_order_rotation(R, C)
    end
    return C
end

function fourth_order_rotation(R, C::Array{Float64,4})
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
                    C[m, n, o, p] = Cmnop
                end
            end
        end
    end
end

end
