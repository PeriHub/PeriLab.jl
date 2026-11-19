# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
module Material_Basis

using LinearAlgebra
using StaticArrays
using ......Helpers: get_MMatrix, determinant, invert, smat, with_dof,
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
    # function barrier on the dof: the fields below get concrete types
    with_dof(D -> _local_damping_due_to_damage!(nodes, params, dt, D))
end

function _local_damping_due_to_damage!(nodes::AbstractVector{Int64}, params, dt,
                                       ::Val{D}) where {D}
    damage = Data_Manager.get_field("Damage", "NP1")::NodeScalarField{Float64}
    nlist = Data_Manager.get_nlist()::BondScalarState{Int64}
    density = Data_Manager.get_field("Density")::NodeScalarField{Float64}

    deformed_bond_lengthN = Data_Manager.get_field("Deformed Bond Length",
                                                   "N")::BondScalarState{Float64}
    deformed_bond_lengthNP1 = Data_Manager.get_field("Deformed Bond Length",
                                                     "NP1")::BondScalarState{Float64}
    deformend_bond_geometry = Data_Manager.get_field("Deformed Bond Geometry",
                                                     "NP1")::BondVectorState{Float64,D}
    local_damping::Float64 = params.damping_coefficient
    force_densities = Data_Manager.get_field("Force Densities",
                                             "NP1")::NodeVectorField{Float64}
    volume = Data_Manager.get_field("Volume")::NodeScalarField{Float64}
    E::Float64 = params.representative_youngs_modulus
    constant = Data_Manager.get_field("Bond Based Constant")::NodeScalarField{Float64}

    @inbounds for iID in nodes
        v0 = sqrt(E / density[iID])
        neighbors = nlist[iID]
        lengthN = deformed_bond_lengthN[iID]
        lengthNP1 = deformed_bond_lengthNP1[iID]
        bonds = deformend_bond_geometry[iID]
        for jID in eachindex(neighbors)
            nID = neighbors[jID]
            avg_damage = 0.5 * (damage[iID] + damage[nID])

            t = local_damping *
                E *
                constant[iID] *
                avg_damage *
                (lengthNP1[jID] - lengthN[jID]) /
                (dt * v0) * bonds[jID] /
                lengthNP1[jID]

            for d in 1:D
                force_densities[iID, d] += t[d] * volume[nID]
                force_densities[nID, d] -= t[d] * volume[iID]
            end
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
        constant[iID] = bond_based_constant(symmetry, horizon[iID])
    end
end

# https://doi.org/10.1016/j.apm.2024.01.015 under EQ (9)
function bond_based_constant(symmetry, horizon::Float64)
    symmetry == "plane stress" && return 9 / (pi * horizon^3)
    symmetry == "plane strain" && return 48 / (5 * pi * horizon^3)
    return 12 / (pi * horizon^4)
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
    contact_forces!(force_densities, nodes, nlist, nlist_filtered_ids, volume, bond_geometry, deformed_bond, bond_length, bond_norm, contact_stiffness)

Adds the forces of the contact bonds of a bond filter with `Allow Contact`. A contact
bond is broken for the material; it only pushes its nodes apart along the filter normal
`n` while they approach each other: f = c (Δu⋅n) / |ξ| n for Δu⋅n < 0 (bond-based style,
no friction).

# Arguments
- `nlist_filtered_ids::NestedNeighborhoodList`: The contact bonds (bond indices) of each node.
- `bond_norm`: The filter normal of a contact bond (oriented along the bond, see `apply_bond_filters`).
- `contact_stiffness`: The stiffness c of each contact bond, aligned to `nlist_filtered_ids`.
"""
function contact_forces!(force_densities::Matrix{Float64},
                         nodes::AbstractVector{Int64},
                         nlist::BondScalarState{Int64},
                         nlist_filtered_ids::NestedNeighborhoodList,
                         volume::NodeScalarField{Float64},
                         bond_geometry::BondVectorState{Float64,D},
                         deformed_bond::BondVectorState{Float64,D},
                         bond_length::BondScalarState{Float64},
                         bond_norm::BondVectorState{Float64,D},
                         contact_stiffness::Vector{Vector{Float64}}) where {D}
    @inbounds for iID in nodes
        neighbors = nlist[iID]
        for (k, jID) in enumerate(nlist_filtered_ids[iID])
            n = normalize(bond_norm[iID][jID])
            # relative displacement of the neighbor along n
            closing = dot(deformed_bond[iID][jID] - bond_geometry[iID][jID], n)
            closing < 0 || continue
            f = (contact_stiffness[iID][k] * closing / bond_length[iID][jID]) * n
            neighborID = neighbors[jID]
            for m in 1:D
                force_densities[iID, m] += f[m] * volume[neighborID]
                force_densities[neighborID, m] -= f[m] * volume[iID]
            end
        end
    end
end

"""
    youngs_modulus_along(material, n, dof, ID)

Young's modulus of a typed block material at node `ID` along the unit vector `n`:
`E` for an isotropic material, else from the compliance, 1/E_n = s' C⁻¹ s with the
Voigt vector s of the uniaxial stress n⊗n.
"""
function youngs_modulus_along(material, n::AbstractVector{Float64}, dof::Int64, ID::Int64)
    symmetry = lowercase(material.hooke_symmetry)
    any(s -> occursin(s, symmetry),
        ("anisotropic", "orthotropic", "transverse isotropic")) ||
        return _typed_constant(material, "Young's Modulus", ID)
    C = hooke_matrix(material, dof, ID)
    s = dof == 2 ? SVector(n[1]^2, n[2]^2, n[1] * n[2]) :
        SVector(n[1]^2, n[2]^2, n[3]^2, n[2] * n[3], n[1] * n[3], n[1] * n[2])
    return 1 / dot(s, C \ s)
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

function flaw_function(::Nothing, coor::AbstractVector{<:Real},
                       stress::Union{Int64,Float64})
    Float64(stress)
end

# flaw function of a block material (`FlawFunctionParams`); the input check makes
# sure an active one has a positive size and magnitude and a location
function flaw_function(flaw, coor::AbstractVector{<:Real},
                       stress::T)::Float64 where {T<:Union{Int64,Float64}}
    flaw.active || return stress
    flaw_size::Float64 = flaw.flaw_size
    flaw_magnitude::Float64 = flaw.flaw_magnitude
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
    return spherical_stress
end

"""
	get_strain(stress_NP1::Matrix{Float64}, hooke_matrix::Matrix{Float64})

# Arguments
- `stress_NP1::Matrix{Float64}`: Stress.
- `hooke_matrix::Matrix{Float64}`: Hooke matrix
# returns
- `strain::Matrix{Float64}`: Strain
"""
function get_strain(stress_NP1::AbstractMatrix{Float64},
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
        bf = bond_force[i]
        @inbounds @fastmath for j in eachindex(bf)
            bf[j] = bf[j] * E
        end
    end
end

function apply_pointwise_E(nodes::AbstractVector{Int64},
                           E::Union{SubArray,Vector{Float64},Vector{Int64}},
                           bond_force::BondVectorState{Float64})
    @inbounds @fastmath for i in nodes
        bf = bond_force[i]
        @inbounds @fastmath for j in eachindex(bf)
            bf[j] = bf[j] * E[i]
        end
    end
end

end
