# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
module Material_Basis

using LinearAlgebra
using LoopVectorization
using StaticArrays
include("../../Support/Helpers.jl")
using .Helpers: get_MMatrix, determinant, invert, smat, interpol_data, get_dependent_value
export get_value
export get_all_elastic_moduli
export get_Hooke_matrix
export distribute_forces!
export local_damping_due_to_damage
export flaw_function
export matrix_to_voigt
export voigt_to_matrix
export matrix_to_vector
export vector_to_matrix
export check_symmetry
export get_symmetry
export get_von_mises_yield_stress
export compute_deviatoric_and_spherical_stresses
export get_strain
export compute_Piola_Kirchhoff_stress
export apply_pointwise_E
export compute_bond_based_constants
export init_local_damping_due_to_damage
export local_damping_due_to_damage
function local_damping_due_to_damage(datamanager::Module, nodes, params, dt)
    damage = datamanager.get_field("Damage", "NP1")
    nlist = datamanager.get_nlist()
    density = datamanager.get_field("Density")

    deformed_bond_lengthN = datamanager.get_field("Deformed Bond Length", "N")
    deformed_bond_lengthNP1 = datamanager.get_field("Deformed Bond Length", "NP1")
    deformend_bond_geometry = datamanager.get_field("Deformed Bond Geometry", "NP1")
    local_damping = params["Damping coefficient"]
    force_densities = datamanager.get_field("Force Densities", "NP1")
    volume = datamanager.get_field("Volume")
    E = params["Representative Young's modulus"]
    constant = datamanager.get_field("Bond Based Constant")
    t = zeros(Float64, datamanager.get_dof())

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

function init_local_damping_due_to_damage(datamanager::Module,
                                          nodes,
                                          material_parameter,
                                          damage_parameter)
    if !haskey(damage_parameter["Local Damping"], "Representative Young's modulus")
        @error "Representative Young's modulus is missing."
        return nothing
    end
    if !haskey(damage_parameter["Local Damping"], "Damping coefficient")
        @error "Damping coefficient is missing."
        return nothing
    end
    @info "Local damping is active with damping coefficient $(damage_parameter["Local Damping"]["Damping coefficient"])"
    constant = datamanager.create_constant_node_field("Bond Based Constant", Float64, 1)
    horizon = datamanager.get_field("Horizon")
    symmetry::String = get_symmetry(material_parameter)
    compute_bond_based_constants(nodes, symmetry, constant, horizon)
end

function compute_bond_based_constants(nodes, symmetry, constant, horizon)
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

function get_value(datamanager::Module,
                   parameter::Union{Dict{Any,Any},Dict{String,Any}},
                   any_field_allocated::Bool,
                   key::String,
                   field_allocated::Bool)
    if field_allocated
        return datamanager.get_field(replace(key, " " => "_"))
    end
    if any_field_allocated
        if haskey(parameter, key)
            return datamanager.create_constant_node_field(replace(key, " " => "_"),
                                                          Float64,
                                                          1,
                                                          parameter[key])
        else
            return datamanager.create_constant_node_field(replace(key, " " => "_"),
                                                          Float64,
                                                          1)
        end
    elseif haskey(parameter, key)
        return parameter[key]
    end

    return Float64(0.0)
end

"""
    get_all_elastic_moduli(datamanager::Module, parameter::Union{Dict{Any,Any},Dict{String,Any}})

Returns the elastic moduli of the material.

# Arguments
- `parameter::Union{Dict{Any,Any},Dict{String,Any}}`: The material parameter.
"""
function get_all_elastic_moduli(datamanager::Module,
                                parameter::Union{Dict{Any,Any},Dict{String,Any}})
    state_factor_defined = haskey(parameter, "State Factor ID")

    if haskey(parameter, "Computed") &&
       !(state_factor_defined && datamanager.has_key("State Variables"))
        if parameter["Computed"]
            return nothing
        end
    end

    bond_based = occursin("Bond-based", parameter["Material Model"])
    if bond_based
        bond_based = !occursin("Unified Bond-based", parameter["Material Model"])
    end
    bulk_field = datamanager.has_key("Bulk_Modulus")
    youngs_field = datamanager.has_key("Young's_Modulus")
    poissons_field = datamanager.has_key("Poisson's_Ratio")
    shear_field = datamanager.has_key("Shear_Modulus")

    any_field_allocated = bulk_field | youngs_field | poissons_field | shear_field |
                          state_factor_defined

    bulk = haskey(parameter, "Bulk Modulus") | bulk_field
    youngs = haskey(parameter, "Young's Modulus") | youngs_field
    shear = haskey(parameter, "Shear Modulus") | shear_field
    poissons = haskey(parameter, "Poisson's Ratio") | poissons_field

    K = get_value(datamanager, parameter, any_field_allocated, "Bulk Modulus", bulk_field)
    E = get_value(datamanager,
                  parameter,
                  any_field_allocated,
                  "Young's Modulus",
                  youngs_field)
    G = get_value(datamanager, parameter, any_field_allocated, "Shear Modulus", shear_field)

    nu = get_value(datamanager,
                   parameter,
                   any_field_allocated,
                   "Poisson's Ratio",
                   poissons_field)

    if bond_based
        nu_fixed = datamanager.get_dof() == 2 ? 1 / 3 : 1 / 4
        if nu != 0.0 && nu != nu_fixed
            @warn "Chosen Bond-based model only supports a fixed Poisson's ratio of " *
                  string(nu_fixed)
        end
        nu = nu_fixed
        poissons = true
    end
    if haskey(parameter, "Symmetry")
        if occursin("anisotropic", lowercase(parameter["Symmetry"]))
            for iID in 1:6
                for jID in iID:6
                    if !("C" * string(iID) * string(jID) in keys(parameter))
                        @error "C" * string(iID) * string(jID) * " not defined"
                        return nothing
                    end
                end
            end
            return
        elseif occursin("orthotropic", lowercase(parameter["Symmetry"]))
            E_x = haskey(parameter, "Young's Modulus X")
            E_y = haskey(parameter, "Young's Modulus Y")
            E_z = haskey(parameter, "Young's Modulus Z")
            nu_xy = haskey(parameter, "Poisson's Ratio XY")
            nu_yz = haskey(parameter, "Poisson's Ratio YZ")
            nu_xz = haskey(parameter, "Poisson's Ratio XZ")
            g_xy = haskey(parameter, "Shear Modulus XY")
            g_yz = haskey(parameter, "Shear Modulus YZ")
            g_zx = haskey(parameter, "Shear Modulus XZ")
            if !E_x || !E_y || !E_z || !nu_xy || !nu_yz || !nu_xz || !g_xy || !g_yz || !g_zx
                @error "Orthotropic material requires Young's Modulus X, Y, Z, Poisson's Ratio XY, YZ, XZ, Shear Modulus XY, YZ, XZ"
                return nothing
            end
            return
        end
    else
        @warn "Material symmetry is not defined, assuming isotropic material"
        parameter["Symmetry"] = "isotropic"
    end

    # tbd non isotropic material check
    if bulk + youngs + shear + poissons < 2
        @error "Minimum of two parameters are needed for isotropic material"
        return nothing
    elseif bulk + youngs + shear + poissons > 2
        @warn "Only two parameters are needed for isotropic material, ignoring additional parameters"
    end

    if bulk && poissons
        E = 3 .* K .* (1 .- 2 .* nu)
        G = 3 .* K .* (1 .- 2 .* nu) ./ (2 .+ 2 .* nu)
    end
    if shear && poissons
        E = 2 .* G .* (1 .+ nu)
        K = 2 .* G .* (1 .+ nu) ./ (3 .- 6 .* nu)
    end
    if bulk && shear
        E = 9 .* K .* G ./ (3 .* K .+ G)
        nu = (3 .* K .- 2 .* G) ./ (6 .* K .+ 2 .* G)
    end
    if youngs && shear
        K = E .* G ./ (9 .* G .- 3 .* E)
        nu = E ./ (2 .* G) .- 1
    end

    if youngs && bulk
        G = 3 .* K .* E ./ (9 .* K .- E)
        nu = (3 .* K .- E) ./ (6 .* K)
    end
    if youngs && poissons
        K = E ./ (3 .- 6 .* nu)
        G = E ./ (2 .+ 2 .* nu)
    end

    if state_factor_defined && datamanager.has_key("State Variables")
        state_factor = datamanager.get_field("State Variables")[:,
                                                                parameter["State Factor ID"]]
        K .*= state_factor
        E .*= state_factor
        G .*= state_factor
    end

    parameter["Bulk Modulus"] = K
    parameter["Young's Modulus"] = E
    parameter["Shear Modulus"] = G
    parameter["Poisson's Ratio"] = nu
    parameter["Computed"] = true
    if any_field_allocated
        datamanager.get_field("Bulk_Modulus") .= K
        datamanager.get_field("Young's_Modulus") .= E
        datamanager.get_field("Shear_Modulus") .= G
        datamanager.get_field("Poisson's_Ratio") .= nu
    end
end

"""
    get_Hooke_matrix(datamanager::Module, parameter::Dict, symmetry::String, dof::Int64, ID::Int64=1)

Returns the Hooke matrix of the material.

# Arguments
- `datamanager::Module`: The data manager.
- `parameter::Union{Dict{Any,Any},Dict{String,Any}}`: The material parameter.
- `symmetry::String`: The symmetry of the material.
- `dof::Int64`: The degree of freedom.
- `ID::Int64=1`: ID of the point. Needed for point wise defined material properties.
# Returns
- `matrix::Matrix{Float64}`: The Hooke matrix.
"""
function get_Hooke_matrix(datamanager::Module,
                          parameter::Dict,
                          symmetry::String,
                          dof::Int64,
                          ID::Int64 = 1)
    """https://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_plane_stress.cfm"""

    symmetry = lowercase(symmetry)
    if occursin("anisotropic", symmetry)
        aniso_matrix = get_MMatrix(36)
        for iID in 1:6
            for jID in iID:6
                value = get_dependent_value(datamanager,
                                            "C" * string(iID) * string(jID),
                                            parameter)
                aniso_matrix[iID, jID] = value
                aniso_matrix[jID, iID] = value
            end
        end
        return get_2D_Hooke_matrix(aniso_matrix, symmetry, dof)
    elseif occursin("orthotropic", symmetry)
        aniso_matrix = get_MMatrix(36)

        E_x = get_dependent_value(datamanager, "Young's Modulus X", parameter, ID)
        E_y = get_dependent_value(datamanager, "Young's Modulus Y", parameter, ID)
        E_z = get_dependent_value(datamanager, "Young's Modulus Z", parameter, ID)
        nu_xy = get_dependent_value(datamanager, "Poisson's Ratio XY", parameter, ID)
        nu_yz = get_dependent_value(datamanager, "Poisson's Ratio YZ", parameter, ID)
        nu_xz = get_dependent_value(datamanager, "Poisson's Ratio XZ", parameter, ID)
        g_xy = get_dependent_value(datamanager, "Shear Modulus XY", parameter, ID)
        g_yz = get_dependent_value(datamanager, "Shear Modulus YZ", parameter, ID)
        g_xz = get_dependent_value(datamanager, "Shear Modulus XZ", parameter, ID)

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

        aniso_matrix[4, 4] = 2 * g_yz
        aniso_matrix[5, 5] = 2 * g_xz
        aniso_matrix[6, 6] = 2 * g_xy

        return get_2D_Hooke_matrix(aniso_matrix, symmetry, dof)
    end

    iID = ID
    if parameter["Poisson's Ratio"] isa Float64
        iID = 1
    end
    if occursin("isotropic", symmetry)
        nu = parameter["Poisson's Ratio"][iID]
        E = parameter["Young's Modulus"][iID]
        G = parameter["Shear Modulus"][iID]
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
            matrix[4, 4] = (1 - 2 * nu) * temp
            matrix[5, 5] = (1 - 2 * nu) * temp
            matrix[6, 6] = (1 - 2 * nu) * temp
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
            @error "2D model defintion is missing; plane stress or plane strain "
            return nothing
        end
    else
        matrix = get_MMatrix(9)

        @warn "material model defintion is missing; assuming isotropic plane stress "
        nu = parameter["Poisson's Ratio"][iID]
        E = parameter["Young's Modulus"][iID]
        G = parameter["Shear Modulus"][iID]
        matrix[1, 1] = E / (1 - nu * nu)
        matrix[1, 2] = E * nu / (1 - nu * nu)
        matrix[2, 1] = E * nu / (1 - nu * nu)
        matrix[2, 2] = E / (1 - nu * nu)
        matrix[3, 3] = G
        return matrix
    end
end

function get_2D_Hooke_matrix(aniso_matrix::MMatrix, symmetry::String, dof::Int64)
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
        @error "2D model defintion is missing; plane stress or plane strain "
        return nothing
    end
end

"""
    distribute_forces!(nodes::AbstractVector{Int64}, nlist::Vector{Vector{Int64}}, nlist_filtered_ids::Vector{Vector{Int64}}, bond_force::Vector{Matrix{Float64}}, volume::Vector{Float64}, bond_damage::Vector{Vector{Float64}}, displacements::Matrix{Float64}, bond_norm::Vector{Matrix{Float64}}, force_densities::Matrix{Float64})

Distribute the forces on the nodes

# Arguments
- `nodes::AbstractVector{Int64}`: The nodes.
- `nlist::Vector{Vector{Int64}}`: The neighbor list.
- `nlist_filtered_ids::Vector{Vector{Int64}},`:  The filtered neighbor list.
- `bond_force::Vector{Matrix{Float64}}`: The bond forces.
- `volume::Vector{Float64}`: The volumes.
- `bond_damage::Vector{Vector{Float64}}`: The bond damage.
- `displacements::Matrix{Float64}`: The displacements.
- `bond_norm::Vector{Matrix{Float64}}`: The pre defined bond normal.
- `force_densities::Matrix{Float64}`: The force densities.
# Returns
- `force_densities::Matrix{Float64}`: The force densities.
"""
function distribute_forces!(force_densities::Matrix{Float64},
                            nodes::AbstractVector{Int64},
                            nlist::Vector{Vector{Int64}},
                            nlist_filtered_ids::Vector{Vector{Int64}},
                            bond_force::Vector{Vector{Vector{Float64}}},
                            volume::Union{Vector{Float64},Vector{Int64}},
                            bond_damage::Vector{Vector{Float64}},
                            displacements::Matrix{Float64},
                            bond_norm::Vector{Vector{Vector{Float64}}})
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
    distribute_forces!(nodes::AbstractVector{Int64}, nlist::Vector{Vector{Int64}}, bond_force::Vector{Matrix{Float64}}, volume::Vector{Float64}, bond_damage::Vector{Vector{Float64}}, force_densities::Matrix{Float64})

Distribute the forces on the nodes

# Arguments
- `nodes::AbstractVector{Int64}`: The nodes.
- `nlist::Vector{Vector{Int64}}`: The neighbor list.
- `bond_force::Vector{Matrix{Float64}}`: The bond forces.
- `volume::Vector{Float64}`: The volumes.
- `bond_damage::Vector{Vector{Float64}}`: The bond damage.
- `force_densities::Matrix{Float64}`: The force densities.
# Returns
- `force_densities::Matrix{Float64}`: The force densities.
"""
function distribute_forces!(force_densities::Matrix{Float64},
                            nodes::AbstractVector{Int64},
                            nlist::Vector{Vector{Int64}},
                            bond_force::Vector{Vector{Vector{Float64}}},
                            volume::Union{Vector{Float64},Vector{Int64}},
                            bond_damage::Vector{Vector{Float64}})
    @inbounds @fastmath for iID in nodes
        @views @inbounds @fastmath for jID in axes(nlist[iID], 1)
            @views @inbounds @fastmath for m in axes(force_densities[iID, :], 1)
                #temp = bond_damage[iID][jID] * bond_force[iID][jID, m]
                force_densities[iID,
                                m] += bond_damage[iID][jID] *
                                      bond_force[iID][jID][m] *
                                      volume[nlist[iID][jID]]
                force_densities[nlist[iID][jID],
                                m] -= bond_damage[iID][jID] *
                                      bond_force[iID][jID][m] * volume[iID]
            end
        end
    end
end

"""
    matrix_to_voigt(matrix)

Convert a 2x2 or 3x3 matrix to Voigt notation (6x1 vector)

# Arguments
- `matrix::Matrix{Float64}`: The matrix.
# Returns
- `voigt::Vector{Float64}`: The Voigt notation.
"""
function matrix_to_voigt(matrix)
    if length(matrix) == 4
        return [matrix[1, 1]; matrix[2, 2]; 0.5 * (matrix[1, 2] + matrix[2, 1])]
    elseif length(matrix) == 9
        return [matrix[1, 1]
                matrix[2, 2]
                matrix[3, 3]
                0.5 * (matrix[2, 3] + matrix[3, 2])
                0.5 * (matrix[1, 3] + matrix[3, 1])
                0.5 * (matrix[1, 2] + matrix[2, 1])]
    else
        @error "Unsupported matrix size for matrix_to_voigt"
        return nothing
    end
end

"""
    voigt_to_matrix(voigt::Union{Vector{Float64},Vector{Int64}})

Convert a Voigt notation (6x1 or 3x1 vector) to a 2x2 or 3x3 matrix

# Arguments
- `voigt::Vector{Float64}`: The Voigt notation.
# Returns
- `matrix::Matrix{Float64}`: The matrix.
"""
function voigt_to_matrix(voigt::Union{MVector,SVector,Vector})
    if length(voigt) == 3
        return @SMatrix [voigt[1] voigt[3]; voigt[3] voigt[2]]
    elseif length(voigt) == 6
        return @SMatrix [voigt[1] voigt[6] voigt[5]
                         voigt[6] voigt[2] voigt[4]
                         voigt[5] voigt[4] voigt[3]]
    else
        @error "Unsupported matrix size for voigt_to_matrix"
        return nothing
    end
end

"""
    matrix_to_vector(matrix)

Convert a 3x3 matrix to a 6x1 vector

# Arguments
- `matrix::Matrix{Float64}`: The matrix.
# Returns
- `vector::Vector{Float64}`: The vector.
"""
function matrix_to_vector(matrix)
    if length(matrix) == 4
        return [matrix[1, 1]; matrix[2, 2]; 0.0; matrix[1, 2]; matrix[2, 1]]
    elseif length(matrix) == 9
        return [matrix[1, 1]
                matrix[2, 2]
                matrix[3, 3]
                matrix[1, 2]
                matrix[2, 3]
                matrix[3, 1]
                matrix[2, 1]
                matrix[3, 2]
                matrix[1, 3]]
    end
end

"""
    vector_to_matrix(matrix)

Convert a 6x1 vector to a 3x3 matrix

# Arguments
- `vector::Vector{Float64}`: The vector.
# Returns
- `matrix::Matrix{Float64}`: The matrix.
"""
function vector_to_matrix(vector)
    if length(vector) == 5
        return @SMatrix [vector[1] vector[3]
                         vector[4] vector[2]]
    elseif length(vector) == 9
        return @SMatrix [vector[1] vector[4] vector[9]
                         vector[7] vector[2] vector[5]
                         vector[6] vector[8] vector[3]]
    else
        @error "Unsupported vector size for vector_to_matrix"
        return nothing
    end
end

"""
    check_symmetry(prop::Dict, dof::Int64)

Check if the symmetry information is present in the material dictionary.

# Arguments
- `prop::Dict`: A dictionary containing material information.
- `dof::Int64`: The number of degrees of freedom.
# Returns
- `true`: If the symmetry information is present.
"""
function check_symmetry(prop::Dict, dof::Int64)
    if haskey(prop, "Symmetry")
        symmetry = prop["Symmetry"]
        if dof == 2
            if occursin("plane strain", symmetry) || occursin("plane stress", symmetry)
                return true
            else
                @error "Model definition is missing; plane stress or plane strain has to be defined for 2D"
                return
            end
        end
        return true
    end
end

"""
    flaw_function(params::Dict, coor::Union{Vector{Int64},Vector{Float64}}, stress::Float64)

Allows the modification of the yield stress at a specific position. This is typically used as starting point for plastic deformation.

# Arguments
- `params::Dict`: A dictionary containing material information.
- `coor::Union{Vector{Int64},Vector{Float64}, SubArray}`: Coordinate of the current point.
- `stress::Float64`: stresses to be modified.

# Returns
- `stress`::Float64: the modified stresses.
"""
function flaw_function(params::Dict,
                       coor::Union{Vector{Int64},Vector{Float64},SubArray{Float64}},
                       stress::Float64)
    if !haskey(params, "Flaw Function")
        return stress
    end
    if !haskey(params["Flaw Function"], "Active")
        @error "Flaw Function needs an entry ''Active''."
        return nothing
    end
    if !haskey(params["Flaw Function"], "Function")
        @error "Flaw Function needs an entry ''Function''."
        return nothing
    end
    if !params["Flaw Function"]["Active"]
        return stress
    end
    flaw_location::Vector{Float64} = zeros(length(coor))
    if params["Flaw Function"]["Function"] == "Pre-defined"
        flaw_size = params["Flaw Function"]["Flaw Size"]
        flaw_magnitude = params["Flaw Function"]["Flaw Magnitude"]
        if !(0 < flaw_magnitude <= 1)
            @error "Flaw Magnitude should be between 0 and 1"
            return nothing
        end
        flaw_location[1] = params["Flaw Function"]["Flaw Location X"]
        flaw_location[2] = params["Flaw Function"]["Flaw Location Y"]
        if haskey(params["Flaw Function"], "Flaw location Z") && length(coor) == 3
            flaw_location[3] = params["Flaw Function"]["Flaw Location Z"]
        end
        modified_stress = stress * (1 -
                           flaw_magnitude *
                           exp(-norm(coor - flaw_location) * norm(coor - flaw_location) /
                               flaw_size /
                               flaw_size))

    else
        @warn "Not very user friendly right now"
        #global x = coor[1]
        #global y = coor[2]
        #
        #modified_stress = stress * (1 - eval(Meta.parse(params["Flaw Function"]["Function"])))
    end
    return modified_stress
end

"""
    get_symmetry(material::Dict)

Return the symmetry information from the given material dictionary.

# Arguments
- `material::Dict`: A dictionary containing material information.

# Returns
- If the key "Symmetry" is present in the dictionary, the corresponding value is returned.
- If the key is not present, the default value "3D" is returned.

# Example
```julia
material_dict = Dict("Symmetry" => "Cubic", "Color" => "Red")
symmetry = get_sym(material_dict)
"""
function get_symmetry(material::Dict)
    if !haskey(material, "Symmetry")
        return "3D"
    end
    if occursin("plane strain", lowercase(material["Symmetry"]))
        return "plane strain"
    end
    if occursin("plane stress", lowercase(material["Symmetry"]))
        return "plane stress"
    end
    return "3D"
end

"""
    get_von_mises_yield_stress(von_Mises_stress::Float64, dof::Int64, stress_NP1::Matrix{Float64})

# Arguments
- `von_Mises_stress::Float64`: Von Mises stress
- `dof::Int64`: Degree of freedom.
- `stress_NP1::Matrix{Float64}`: Stress.
# returns
- `spherical_stress_NP1::Float64`: Spherical stress
- `deviatoric_stress_NP1::Matrix{Float64}`: Deviatoric stress
"""

function get_von_mises_yield_stress(von_Mises_stress::Float64,
                                    spherical_stress,
                                    deviatoric_stress)
    temp = zero(eltype(deviatoric_stress))
    @views @inbounds @fastmath for i in axes(deviatoric_stress, 1)
        for j in axes(deviatoric_stress, 2)
            temp += deviatoric_stress[i, j] * deviatoric_stress[i, j]
        end
    end

    von_Mises_stress = sqrt(3.0 / 2.0 * temp)

    return von_Mises_stress
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
                    hooke_matrix::Union{Matrix{Float64},MMatrix,SMatrix})
    return voigt_to_matrix(hooke_matrix' * matrix_to_voigt(stress_NP1))
end

function compute_Piola_Kirchhoff_stress(stress::Union{Matrix{Float64},SubArray{Float64}},
                                        deformation_gradient::Union{Matrix{Float64},
                                                                    SubArray{Float64}})
    #50% less memory
    return determinant(deformation_gradient) .* smat(stress) * invert(deformation_gradient,
                  "Deformation gradient is singular and cannot be inverted.")
end

function apply_pointwise_E(nodes, E::Union{Int64,Float64}, bond_force)
    @inbounds @fastmath for i in nodes
        @views @inbounds @fastmath for j in axes(bond_force, 2)
            bond_force[i, j] *= E
        end
    end
end

function apply_pointwise_E(nodes,
                           E::Union{SubArray,Vector{Float64},Vector{Int64}},
                           bond_force)
    @inbounds @fastmath for i in nodes
        @views @inbounds @fastmath for j in axes(bond_force, 2)
            bond_force[i, j] *= E[i]
        end
    end
end

function apply_pointwise_E(nodes, bond_force, dependent_field)
    warning_flag = true
    @inbounds @fastmath for i in nodes
        E_int = interpol_data(dependent_field[i],
                              damage_parameter["Young's Modulus"]["Data"],
                              warning_flag)
        @views @inbounds @fastmath for j in axes(bond_force, 2)
            bond_force[i, j] *= E_int
        end
    end
end

end
