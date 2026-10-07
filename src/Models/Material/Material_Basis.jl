# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
module Material_Basis

using LinearAlgebra
using LoopVectorization
using StaticArrays
using ......Helpers: get_MMatrix, determinant, invert, smat, interpol_data,
                     mat_mul!, matrix_to_voigt, voigt_to_matrix
using ......Data_Manager
using ......ParameterSpec: value
using ......PeriLabExceptions: @abort
export hooke_matrix
export distribute_forces!
export local_damping_due_to_damage
export flaw_function
export get_von_mises_yield_stress
export compute_deviatoric_and_spherical_stresses
export get_strain
export compute_Piola_Kirchhoff_stress!
export apply_pointwise_E
export compute_bond_based_constants
export init_local_damping_due_to_damage
export local_damping_due_to_damage
function local_damping_due_to_damage(nodes::AbstractVector{Int64},
                                     params, dt)
    damage = Data_Manager.get_field("Damage", "NP1")
    nlist = Data_Manager.get_nlist()
    density = Data_Manager.get_field("Density")

    deformed_bond_lengthN = Data_Manager.get_field("Deformed Bond Length", "N")
    deformed_bond_lengthNP1 = Data_Manager.get_field("Deformed Bond Length", "NP1")
    deformend_bond_geometry = Data_Manager.get_field("Deformed Bond Geometry", "NP1")
    local_damping = params.damping_coefficient
    force_densities = Data_Manager.get_field("Force Densities", "NP1")
    volume = Data_Manager.get_field("Volume")
    E = params.representative_youngs_modulus
    constant = Data_Manager.get_field("Bond Based Constant")
    t = zeros(Float64, Data_Manager.get_dof())

    for iID in nodes
        v0 = sqrt(E / density[iID])
        for (jID, nID) in enumerate(nlist[iID])
            avg_damage = 0.5 * (damage[iID] - damage[nID])

            t = local_damping *
                E *
                constant[iID] *
                avg_damage *
                (deformed_bond_lengthNP1[iID][jID] - deformed_bond_lengthN[iID][jID]) /
                (dt * v0) * deformend_bond_geometry[iID][jID] /
                deformed_bond_lengthNP1[iID][jID]

            force_densities[iID, :] += t * volume[nID]
            force_densities[nID, :] -= t * volume[iID]
        end
    end
end

function init_local_damping_due_to_damage(nodes::AbstractVector{Int64},
                                          symmetry::String,
                                          local_damping)
    @info "Local damping is active with damping coefficient $(local_damping.damping_coefficient)"
    constant = Data_Manager.create_constant_node_scalar_field("Bond Based Constant",
                                                              Float64)
    horizon = Data_Manager.get_field("Horizon")
    compute_bond_based_constants(nodes, symmetry, constant, horizon)
end

function compute_bond_based_constants(nodes::AbstractVector{Int64}, symmetry, constant,
                                      horizon)
    for iID in nodes
        if symmetry == "plane stress"
            constant[iID] = 9 / (pi * horizon[iID]^3) # https://doi.org/10.1016/j.apm.2024.01.015 under EQ (9)
        elseif symmetry == "plane strain"
            constant[iID] = 48 / (5 * pi * horizon[iID]^3) # https://doi.org/10.1016/j.apm.2024.01.015 under EQ (9)
        else
            constant[iID] = 12 / (pi * horizon[iID]^4) # https://doi.org/10.1016/j.apm.2024.01.015 under EQ (9)
        end
    end
end

const _HOOKE_BASE_FIELDS = Dict("Young's Modulus X" => :youngs_modulus_x,
                                "Young's Modulus Y" => :youngs_modulus_y,
                                "Young's Modulus Z" => :youngs_modulus_z,
                                "Poisson's Ratio XY" => :poissons_ratio_xy,
                                "Poisson's Ratio YZ" => :poissons_ratio_yz,
                                "Poisson's Ratio XZ" => :poissons_ratio_xz,
                                "Shear Modulus XY" => :shear_modulus_xy,
                                "Shear Modulus YZ" => :shear_modulus_yz,
                                "Shear Modulus XZ" => :shear_modulus_xz)

@inline _at(x::Real, ::Int64) = x
@inline _at(x::AbstractVector, id::Int64) = x[id]

# constant of a typed block material at node `id`
function _typed_constant(material, key::String, id::Int64)
    key == "Poisson's Ratio" && return _at(material.moduli.poissons_ratio, id)
    key == "Young's Modulus" && return _at(material.moduli.youngs_modulus, id)
    key == "Shear Modulus" && return _at(material.moduli.shear_modulus, id)
    startswith(key, "C") && return getfield(material.base, Symbol(lowercase(key)))
    return value(getfield(material.base, _HOOKE_BASE_FIELDS[key]), id)
end

"""
	hooke_matrix(material, dof, ID = 1)

Hooke matrix of a typed block material (`BlockMaterial`) at node `ID`. Dependent
constants must be bound (`Material.bind_material!`).
"""
function hooke_matrix(material, dof::Int64, ID::Int64 = 1)
    return _hooke_matrix((key, id) -> _typed_constant(material, key, id),
                         material.hooke_symmetry, dof, ID)
end

# formulas of the Hooke matrix; c(key, id) returns a constant at node id
function _hooke_matrix(c, symmetry::String, dof::Int64, ID::Int64)
    """https://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_plane_stress.cfm
        https://de.wikipedia.org/wiki/Transversale_Isotropie"""

    symmetry = lowercase(symmetry)
    if occursin("anisotropic", symmetry)
        aniso_matrix = get_MMatrix(36)
        for iID in 1:6
            for jID in iID:6
                value = c("C" * string(iID) * string(jID), 1)
                aniso_matrix[iID, jID] = value
                aniso_matrix[jID, iID] = value
            end
        end
        return get_2D_Hooke_matrix(aniso_matrix, symmetry, dof)
    elseif occursin("orthotropic", symmetry)
        aniso_matrix = get_MMatrix(36)

        E_x = c("Young's Modulus X", ID)
        E_y = c("Young's Modulus Y", ID)
        E_z = c("Young's Modulus Z", ID)
        nu_xy = c("Poisson's Ratio XY", ID)
        nu_yz = c("Poisson's Ratio YZ", ID)
        nu_xz = c("Poisson's Ratio XZ", ID)
        g_xy = c("Shear Modulus XY", ID)
        g_yz = c("Shear Modulus YZ", ID)
        g_xz = c("Shear Modulus XZ", ID)

        nu_yx = nu_xy * E_y / E_x
        nu_zy = nu_yz * E_z / E_y
        nu_zx = nu_xz * E_z / E_x

        delta = (1 - nu_xy * nu_yx - nu_yz * nu_zy - nu_zx * nu_xz -
                 2 * nu_xy * nu_yz * nu_zx) / (E_x * E_y * E_z)

        aniso_matrix[1, 1] = (1 - nu_yz * nu_zy) / (E_y * E_z * delta)
        aniso_matrix[2, 2] = (1 - nu_zx * nu_xz) / (E_z * E_x * delta)
        aniso_matrix[3, 3] = (1 - nu_xy * nu_yx) / (E_x * E_y * delta)

        aniso_matrix[1, 2] = (nu_yx + nu_zx * nu_yz) / (E_y * E_z * delta)
        aniso_matrix[2, 1] = (nu_xy + nu_xz * nu_zy) / (E_z * E_x * delta)

        aniso_matrix[1, 3] = (nu_zx + nu_yx * nu_zy) / (E_y * E_z * delta)
        aniso_matrix[3, 1] = (nu_xz + nu_xy * nu_yz) / (E_x * E_y * delta)

        aniso_matrix[2, 3] = (nu_zy + nu_zx * nu_xy) / (E_z * E_x * delta)
        aniso_matrix[3, 2] = (nu_yz + nu_xz * nu_yx) / (E_x * E_y * delta)

        aniso_matrix[4, 4] = g_yz
        aniso_matrix[5, 5] = g_xz
        aniso_matrix[6, 6] = g_xy

        return get_2D_Hooke_matrix(aniso_matrix, symmetry, dof)
    elseif occursin("transverse isotropic", symmetry)
        if dof == 3
            aniso_matrix = get_MMatrix(36)

            E_x = c("Young's Modulus X", ID)
            E_y = c("Young's Modulus Y", ID)
            nu_xy = c("Poisson's Ratio XY", ID)
            nu_yz = c("Poisson's Ratio YZ", ID)
            g_xy = c("Shear Modulus XY", ID)
            g_yz = c("Shear Modulus YZ", ID)

            nu_yx = nu_xy * E_y / E_x

            delta = ((nu_xy * nu_yx + nu_yz) /
                     ((1 - nu_yz - 2 * nu_xy * nu_yx) * (1 + nu_yz))) *
                    E_y

            aniso_matrix[1, 1] = ((1 - nu_yz) / (1 - nu_yz - 2 * nu_xy * nu_yx)) * E_x
            aniso_matrix[2, 2] = delta + 2 * g_yz

            aniso_matrix[1, 2] = 2 * nu_xy * (delta + g_yz)
            aniso_matrix[2, 1] = aniso_matrix[1, 2]

            aniso_matrix[1, 3] = aniso_matrix[1, 2]
            aniso_matrix[3, 1] = aniso_matrix[1, 2]

            aniso_matrix[2, 3] = delta
            aniso_matrix[3, 2] = delta

            aniso_matrix[3, 3] = delta + 2 * g_yz

            aniso_matrix[4, 4] = g_yz
            aniso_matrix[5, 5] = g_xy
            aniso_matrix[6, 6] = g_xy

            return aniso_matrix
        elseif occursin("plane strain", symmetry)
            aniso_matrix = get_MMatrix(9)

            E_x = c("Young's Modulus X", ID)
            E_y = c("Young's Modulus Y", ID)
            nu_xy = c("Poisson's Ratio XY", ID)
            nu_yz = c("Poisson's Ratio YZ", ID)
            g_xy = c("Shear Modulus XY", ID)

            nu_yx = nu_xy * E_y / E_x
            D = (1 + nu_yz) * (1 - nu_yz - 2 * nu_xy * nu_yx)

            aniso_matrix[1, 1] = ((1 - nu_yz * nu_yz) / D) * E_x
            aniso_matrix[2, 2] = ((1 - nu_yx * nu_xy) / D) * E_y

            aniso_matrix[1, 2] = ((nu_yx * (1 + nu_yz)) / D) * E_x
            aniso_matrix[2, 1] = ((nu_xy * (1 + nu_yz)) / D) * E_y

            aniso_matrix[3, 3] = g_xy

            return aniso_matrix
        elseif occursin("plane stress", symmetry)
            aniso_matrix = get_MMatrix(9)

            E_x = c("Young's Modulus X", ID)
            E_y = c("Young's Modulus Y", ID)
            nu_xy = c("Poisson's Ratio XY", ID)
            g_xy = c("Shear Modulus XY", ID)

            nu_yx = nu_xy * E_y / E_x

            aniso_matrix[1, 1] = E_x / (1 - nu_xy * nu_yx)
            aniso_matrix[2, 2] = E_y / (1 - nu_xy * nu_yx)

            aniso_matrix[1, 2] = (nu_yx * E_x) / (1 - nu_xy * nu_yx)
            aniso_matrix[2, 1] = (nu_xy * E_y) / (1 - nu_xy * nu_yx)

            aniso_matrix[3, 3] = g_xy

            return aniso_matrix
        else
            @abort "2D model defintion is missing; plane stress or plane strain"
        end
    end

    if occursin("isotropic", symmetry)
        nu = c("Poisson's Ratio", ID)
        E = c("Young's Modulus", ID)
        G = c("Shear Modulus", ID)
        temp = E / ((1 + nu) * (1 - 2 * nu))

        if dof == 3
            matrix = get_MMatrix(36)
            matrix[1, 1] = (1 - nu) * temp
            matrix[2, 2] = (1 - nu) * temp
            matrix[3, 3] = (1 - nu) * temp
            matrix[1, 2] = nu * temp
            matrix[2, 1] = nu * temp
            matrix[1, 3] = nu * temp
            matrix[3, 1] = nu * temp
            matrix[2, 3] = nu * temp
            matrix[3, 2] = nu * temp
            matrix[4, 4] = G
            matrix[5, 5] = G
            matrix[6, 6] = G
            return matrix
        elseif occursin("plane strain", symmetry)
            matrix = get_MMatrix(9)
            matrix[1, 1] = (1 - nu) * temp
            matrix[2, 2] = (1 - nu) * temp
            matrix[3, 3] = G
            matrix[1, 2] = nu * temp
            matrix[2, 1] = nu * temp
            return matrix
        elseif occursin("plane stress", symmetry)
            matrix = get_MMatrix(9)
            matrix[1, 1] = E / (1 - nu * nu)
            matrix[1, 2] = E * nu / (1 - nu * nu)
            matrix[2, 1] = E * nu / (1 - nu * nu)
            matrix[2, 2] = E / (1 - nu * nu)
            matrix[3, 3] = G
            return matrix
        else
            @abort "2D model defintion is missing; plane stress or plane strain"
            return nothing
        end
    else
        matrix = get_MMatrix(9)

        @warn "material model defintion is missing; assuming isotropic plane stress "
        nu = c("Poisson's Ratio", ID)
        E = c("Young's Modulus", ID)
        G = c("Shear Modulus", ID)
        matrix[1, 1] = E / (1 - nu * nu)
        matrix[1, 2] = E * nu / (1 - nu * nu)
        matrix[2, 1] = E * nu / (1 - nu * nu)
        matrix[2, 2] = E / (1 - nu * nu)
        matrix[3, 3] = G
        return matrix
    end
end

function get_2D_Hooke_matrix(aniso_matrix::MMatrix{T}, symmetry::String,
                             dof::Int64) where {T}
    if dof == 3
        return aniso_matrix
    elseif occursin("plane strain", symmetry)
        matrix = get_MMatrix(9)
        matrix[1:2, 1:2] = aniso_matrix[1:2, 1:2]
        matrix[3, 1:2] = aniso_matrix[6, 1:2]
        matrix[1:2, 3] = aniso_matrix[1:2, 6]
        matrix[3, 3] = aniso_matrix[6, 6]
        return matrix
    elseif occursin("plane stress", symmetry)
        inv_aniso = invert(aniso_matrix, "Hooke matrix not invertable")
        matrix = get_MMatrix(36)
        matrix[1:2, 1:2] = inv_aniso[1:2, 1:2]
        matrix[3, 1:2] = inv_aniso[6, 1:2]
        matrix[1:2, 3] = inv_aniso[1:2, 6]
        matrix[3, 3] = inv_aniso[6, 6]
        return invert(matrix, "Hooke matrix not invertable")
    else
        @abort "2D model defintion is missing; plane stress or plane strain"
        return nothing
    end
end

"""
	distribute_forces!(nodes::AbstractVector{Int64}, nlist::BondScalarState{Int64}, nlist_filtered_ids::BondScalarState{Int64}, bond_force::Vector{Matrix{Float64}}, volume::NodeScalarField{Float64}, bond_damage::BondScalarState{Float64}, displacements::Matrix{Float64}, bond_norm::Vector{Matrix{Float64}}, force_densities::Matrix{Float64})

Distribute the forces on the nodes

# Arguments
- `nodes::AbstractVector{Int64}`: The nodes.
- `nlist::BondScalarState{Int64}`: The neighbor list.
- `nlist_filtered_ids::BondScalarState{Int64},`:  The filtered neighbor list.
- `bond_force::Vector{Matrix{Float64}}`: The bond forces.
- `volume::NodeScalarField{Float64}`: The volumes.
- `bond_damage::BondScalarState{Float64}`: The bond damage.
- `displacements::Matrix{Float64}`: The displacements.
- `bond_norm::Vector{Matrix{Float64}}`: The pre defined bond normal.
- `force_densities::Matrix{Float64}`: The force densities.
# Returns
- `force_densities::Matrix{Float64}`: The force densities.
"""
function distribute_forces!(force_densities::Matrix{Float64},
                            nodes::AbstractVector{Int64},
                            nlist::BondScalarState{Int64},
                            nlist_filtered_ids::BondScalarState{Int64},
                            bond_force::BondVectorState{Float64},
                            volume::NodeScalarField{Float64},
                            bond_damage::BondScalarState{Float64},
                            displacements::Matrix{Float64},
                            bond_norm::BondVectorState{Float64})
    @inbounds @fastmath for iID in nodes
        bond_mod = copy(bond_norm[iID])
        if length(nlist_filtered_ids[iID]) > 0
            for neighborID in nlist_filtered_ids[iID]
                if dot((displacements[nlist[iID][neighborID], :] - displacements[iID, :]),
                       bond_norm[iID][neighborID]) > 0
                    bond_mod[neighborID] .= 0
                else
                    bond_mod[neighborID] .= abs.(bond_norm[iID][neighborID])
                end
            end
        end

        @views @inbounds @fastmath for jID in axes(nlist[iID], 1)
            @views @inbounds @fastmath for m in axes(force_densities[iID, :], 1)
                #temp = bond_damage[iID][jID] * bond_force[iID][jID, m]
                force_densities[iID,
                                m] += bond_damage[iID][jID] *
                                      bond_force[iID][jID][m] *
                                      volume[nlist[iID][jID]] *
                                      bond_mod[jID][m]
                force_densities[nlist[iID][jID],
                                m] -= bond_damage[iID][jID] *
                                      bond_force[iID][jID][m] *
                                      volume[iID] *
                                      bond_mod[jID][m]
            end
        end
    end
end

"""
	distribute_forces!(nodes::AbstractVector{Int64}, nlist::BondScalarState{Int64}, bond_force::Vector{Matrix{Float64}}, volume::NodeScalarField{Float64}, bond_damage::BondScalarState{Float64}, force_densities::Matrix{Float64})

Distribute the forces on the nodes

# Arguments
- `nodes::AbstractVector{Int64}`: The nodes.
- `nlist::BondScalarState{Int64}`: The neighbor list.
- `bond_force::Vector{Matrix{Float64}}`: The bond forces.
- `volume::NodeScalarField{Float64}`: The volumes.
- `bond_damage::BondScalarState{Float64}`: The bond damage.
- `force_densities::Matrix{Float64}`: The force densities.
# Returns
- `force_densities::Matrix{Float64}`: The force densities.
"""
function distribute_forces!(force_densities::Matrix{Float64},
                            nodes::AbstractVector{Int64},
                            nlist::BondScalarState{Int64},
                            bond_force::BondVectorState{Float64},
                            volume::NodeScalarField{Float64},
                            bond_damage::BondScalarState{Float64})::Nothing
    @inbounds @fastmath for iID in nodes
        # Cache arrays to help type inference
        @views neighbors = nlist[iID]
        @views neighbor_forces = bond_force[iID]
        @views neighbor_damage = bond_damage[iID]
        vol_i = volume[iID]

        @inbounds @fastmath for jID_idx in eachindex(neighbors)
            jID = neighbors[jID_idx]
            damage_factor = neighbor_damage[jID_idx]
            vol_j = volume[jID]
            forces_j = neighbor_forces[jID_idx]

            @inbounds @fastmath for m in eachindex(forces_j)
                force_contribution = damage_factor * forces_j[m]
                force_densities[iID, m] += force_contribution * vol_j
                force_densities[jID, m] -= force_contribution * vol_i
            end
        end
    end
    return nothing
end

flaw_function(::Nothing, coor::AbstractVector{<:Real}, stress::Union{Int64,Float64}) = Float64(stress)

# flaw function of a block material (`FlawFunctionParams`)
function flaw_function(flaw, coor::AbstractVector{<:Real},
                       stress::T)::Float64 where {T<:Union{Int64,Float64}}
    flaw.active || return stress
    if flaw.flaw_size === nothing || flaw.flaw_magnitude === nothing ||
       flaw.flaw_location_x === nothing || flaw.flaw_location_y === nothing
        @abort "An active Flaw Function needs Flaw Size, Flaw Magnitude, Flaw Location X and Flaw Location Y."
    end
    flaw_size::Float64 = flaw.flaw_size
    flaw_magnitude::Float64 = flaw.flaw_magnitude
    if !(0 < flaw_magnitude <= 1)
        @abort "Flaw Magnitude should be between 0 and 1"
    end
    if flaw_size <= 0
        @abort "Flaw Size must be positive."
    end
    dx = Float64(coor[1]) - flaw.flaw_location_x
    dy = Float64(coor[2]) - flaw.flaw_location_y
    distance_squared = dx * dx + dy * dy
    if length(coor) == 3
        dz = Float64(coor[3]) - something(flaw.flaw_location_z, 0.0)
        distance_squared += dz * dz
    end
    return stress *
           (1 - flaw_magnitude * exp(-distance_squared / (flaw_size * flaw_size)))
end

"""
	get_von_mises_yield_stress(deviatoric_stress::AbstractMatrix{Float64})

# Arguments
- `deviatoric_stress_NP1::AbstractMatrix{Float64}`: Deviatoric stress
# returns
- `von_Mises_stress::Float64`: Von Mises stress

"""

function get_von_mises_yield_stress(deviatoric_stress::AbstractMatrix{Float64})
    temp = zero(eltype(deviatoric_stress))
    @views @inbounds @fastmath for i in axes(deviatoric_stress, 1)
        for j in axes(deviatoric_stress, 2)
            temp += deviatoric_stress[i, j] * deviatoric_stress[i, j]
        end
    end

    return sqrt(3.0 / 2.0 * temp)
end

function compute_deviatoric_and_spherical_stresses(stress,
                                                   spherical_stress,
                                                   deviatoric_stress,
                                                   dof)
    @views @inbounds @fastmath for i in axes(stress, 1)
        spherical_stress += stress[i, i]
    end
    spherical_stress /= dof

    @views @inbounds @fastmath for i in axes(stress, 1)
        for j in axes(stress, 2)
            deviatoric_stress[i, j] = stress[i, j]
        end
        deviatoric_stress[i, i] -= spherical_stress
    end
end

"""
	get_strain(stress_NP1::Matrix{Float64}, hooke_matrix::Matrix{Float64})

# Arguments
- `stress_NP1::Matrix{Float64}`: Stress.
- `hooke_matrix::Matrix{Float64}`: Hooke matrix
# returns
- `strain::Matrix{Float64}`: Strain
"""
function get_strain(stress_NP1::Matrix{Float64},
                    hooke_matrix::AbstractMatrix{Float64})
    strain_voigt = hooke_matrix' * matrix_to_voigt(stress_NP1)
    # hooke_matrix is a compliance matrix here (see calculate_strain), so this product
    # is engineering shear strain (gamma = 2*epsilon) in the Voigt shear slots, the
    # standard Voigt compliance convention. voigt_to_matrix expects tensor strain
    # (undoubled), so halve the shear entries back before the reverse conversion.
    scale = length(strain_voigt) == 3 ? (1.0, 1.0, 0.5) : (1.0, 1.0, 1.0, 0.5, 0.5, 0.5)
    return voigt_to_matrix(strain_voigt .* scale)
end

function compute_Piola_Kirchhoff_stress!(pk_stress::AbstractMatrix{Float64},
                                         stress::AbstractMatrix{Float64},
                                         deformation_gradient::AbstractMatrix{Float64})
    #50% less memory

    mat_mul!(pk_stress,
             smat(stress),
             invert(transpose(deformation_gradient),
                    "Deformation gradient is singular and cannot be inverted."))
    pk_stress .*= determinant(deformation_gradient)
    #return determinant(deformation_gradient) .* smat(stress) * invert(deformation_gradient,
    #              "Deformation gradient is singular and cannot be inverted.")
end

function apply_pointwise_E(nodes::AbstractVector{Int64}, E::Union{Int64,Float64},
                           bond_force::BondVectorState{Float64})
    @inbounds @fastmath for i in nodes
        @views @inbounds @fastmath for bf in bond_force[i]
            bf .*= E
        end
    end
end

function apply_pointwise_E(nodes::AbstractVector{Int64},
                           E::Union{SubArray,Vector{Float64},Vector{Int64}},
                           bond_force::BondVectorState{Float64})
    @inbounds @fastmath for i in nodes
        @views @inbounds @fastmath for bf in bond_force[i]
            bf .*= E[i]
        end
    end
end

function apply_pointwise_E(nodes::AbstractVector{Int64},
                           bond_force::BondVectorState{Float64}, dependent_field)
    warning_flag = true
    @inbounds @fastmath for i in nodes
        E_int = interpol_data(dependent_field[i],
                              damage_parameter["Young's Modulus"]["Data"],
                              warning_flag)
        @views @inbounds @fastmath for bf in bond_force[i]
            bf .*= E_int
        end
    end
end

end
