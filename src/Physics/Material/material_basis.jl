# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
using LinearAlgebra
using StaticArrays

"""
    get_all_elastic_moduli(parameter::Union{Dict{Any,Any},Dict{String,Any}})

Returns the elastic moduli of the material.

# Arguments
- `parameter::Union{Dict{Any,Any},Dict{String,Any}}`: The material parameter.
"""
function get_all_elastic_moduli(parameter::Union{Dict{Any,Any},Dict{String,Any}})
    if haskey(parameter, "Computed")
        if parameter["Computed"]
            return nothing
        end
    end
    K::Float64 = 0
    E::Float64 = 0
    nu::Float64 = 0
    G::Float64 = 0
    bulk = false
    Youngs = false
    shear = false
    Poissons = false
    if haskey(parameter, "Bulk Modulus")
        K = parameter["Bulk Modulus"]
        bulk = true
    end
    if haskey(parameter, "Young's Modulus")
        E = parameter["Young's Modulus"]
        Youngs = true
    end
    if haskey(parameter, "Shear Modulus")
        G = parameter["Shear Modulus"]
        shear = true
    end
    if haskey(parameter, "Poisson's Ratio")
        nu = parameter["Poisson's Ratio"]
        Poissons = true
    end
    if bulk && Poissons
        E = 3 * K * (1 - 2 * nu)
        G = 3 * K * (1 - nu / (2 + 2 * nu))
    end
    if shear && Poissons
        E = 2 * G * (1 + nu)
        K = 2 * G * (1 + nu) / (3 - 6 * nu)
    end
    if bulk && shear
        E = 9 * K * G / (3 * K + G)
        nu = (3 * K - 2 * G) / (6 * K + 2 * G)
    end
    if Youngs && shear
        K = E * G / (9 * G - 3 * E)
        nu = E / (2 * G) - 1
    end

    if Youngs && bulk
        G = 3 * K * E / (9 * K - E)
        nu = (3 * K - E) / (6 * K)
    end
    if Youngs && Poissons
        K = E / (3 - 6 * nu)
        G = E / (2 + 2 * nu)
    end
    # tbd non isotropic material check
    if bulk + Youngs + shear + Poissons < 2
        @error "Minimum of two parameters are needed for isotropic material"
    end

    parameter["Bulk Modulus"] = K
    parameter["Young's Modulus"] = E
    parameter["Shear Modulus"] = G
    parameter["Poisson's Ratio"] = nu
    parameter["Computed"] = true
end

"""
    get_Hooke_matrix(parameter, symmetry, dof)

Returns the Hooke matrix of the material.

# Arguments
- `parameter::Union{Dict{Any,Any},Dict{String,Any}}`: The material parameter.
- `symmetry::String`: The symmetry of the material.
- `dof::Int64`: The degree of freedom.
# Returns
- `matrix::Matrix{Float64}`: The Hooke matrix.
"""
function get_Hooke_matrix(parameter, symmetry, dof)
    """https://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_plane_stress.cfm"""

    if occursin("anisotropic", symmetry)
        anisoMatrix = @MMatrix zeros(Float64, 6, 6)
        for iID in 1:6
            for jID in iID:6
                if "C" * string(iID) * string(jID) in keys(parameter)
                    value = parameter["C"*string(iID)*string(jID)]
                else
                    value = 0
                end
                anisoMatrix[iID, jID] = value
                anisoMatrix[jID, iID] = value
            end
        end
        if dof == 3
            return anisoMatrix
        elseif occursin("plane strain", symmetry)
            matrix = @MMatrix zeros(Float64, dof + 1, dof + 1)
            matrix[1:2, 1:2] = anisoMatrix[1:2, 1:2]
            matrix[3, 1:2] = anisoMatrix[6, 1:2]
            matrix[1:2, 3] = anisoMatrix[1:2, 6]
            matrix[3, 3] = anisoMatrix[6, 6]
            return matrix
        elseif occursin("plane stress", symmetry)
            matrix = @MMatrix zeros(Float64, dof + 1, dof + 1)
            invAniso = inv(anisoMatrix)
            matrix[1:2, 1:2] = invAniso[1:2, 1:2]
            matrix[3, 1:2] = invAniso[6, 1:2]
            matrix[1:2, 3] = invAniso[1:2, 6]
            matrix[3, 3] = invAniso[6, 6]
            return inv(matrix)
        else
            @error "2D model defintion is missing; plain stress or plain strain "
        end
    end
    if occursin("isotropic", symmetry)

        nu = parameter["Poisson's Ratio"]
        E = parameter["Young's Modulus"]
        G = parameter["Shear Modulus"]
        temp = E / ((1 + nu) * (1 - 2 * nu))

        if dof == 3
            matrix = @MMatrix zeros(Float64, 2 * dof, 2 * dof)
            matrix[1, 1] = (1 - nu) * temp
            matrix[2, 2] = (1 - nu) * temp
            matrix[3, 3] = (1 - nu) * temp
            matrix[1, 2] = nu * temp
            matrix[1, 3] = nu * temp
            matrix[2, 1] = nu * temp
            matrix[1, 2] = nu * temp
            matrix[3, 2] = nu * temp
            matrix[2, 3] = nu * temp
            matrix[4, 4] = G
            matrix[5, 5] = G
            matrix[6, 6] = G
            return matrix
        elseif occursin("plane strain", symmetry)
            matrix = @MMatrix zeros(Float64, dof + 1, dof + 1)
            matrix[1, 1] = (1 - nu) * temp
            matrix[2, 2] = (1 - nu) * temp
            matrix[3, 3] = G
            matrix[1, 2] = nu * temp
            matrix[2, 1] = nu * temp
            return matrix
        elseif occursin("plane stress", symmetry)
            matrix = @MMatrix zeros(Float64, dof + 1, dof + 1)
            matrix[1, 1] = E / (1 - nu * nu)
            matrix[1, 2] = E * nu / (1 - nu * nu)
            matrix[2, 1] = E * nu / (1 - nu * nu)
            matrix[2, 2] = E / (1 - nu * nu)
            matrix[3, 3] = G
            return matrix
        else
            @error "2D model defintion is missing; plain stress or plain strain "
        end
    else
        matrix = @MMatrix zeros(Float64, dof + 1, dof + 1)
        if haskey(parameter, "Poisson's Ratio") && haskey(parameter, "Young's Modulus")
            @warn "material model defintion is missing; assuming isotropic plain stress "
            nu = parameter["Poisson's Ratio"]
            E = parameter["Young's Modulus"]
            G = parameter["Shear Modulus"]
            matrix[1, 1] = E / (1 - nu * nu)
            matrix[1, 2] = E * nu / (1 - nu * nu)
            matrix[2, 1] = E * nu / (1 - nu * nu)
            matrix[2, 2] = E / (1 - nu * nu)
            matrix[3, 3] = G
            return matrix
        end
        @error "no valid definition"
        return nothing

    end
    return matrix
end

"""
    distribute_forces(nodes::Union{SubArray,Vector{Int64}}, nlist::SubArray, nlist_filtered_ids::SubArray, bond_force::SubArray, volume::SubArray, bond_damage::SubArray, displacements::SubArray, bond_norm::SubArray, force_densities::SubArray)

Distribute the forces on the nodes

# Arguments
- `nodes::Union{SubArray,Vector{Int64}}`: The nodes.
- `nlist::SubArray`: The neighbor list.
- `nlist_filtered_ids::SubArray`:  The filtered neighbor list.
- `bond_force::SubArray`: The bond forces.
- `volume::SubArray`: The volumes.
- `bond_damage::SubArray`: The bond damage.
- `displacements::SubArray`: The displacements.
- `bond_norm::SubArray`: The pre defined bond normal.
- `force_densities::SubArray`: The force densities.
# Returns
- `force_densities::SubArray`: The force densities.
"""
function distribute_forces(nodes::Union{SubArray,Vector{Int64}}, nlist::SubArray, nlist_filtered_ids::SubArray, bond_force::SubArray, volume::SubArray, bond_damage::SubArray, displacements::SubArray, bond_norm::SubArray, force_densities::SubArray)

    for iID in nodes
        bond_mod = copy(bond_norm[iID])
        if length(nlist_filtered_ids[iID]) > 0
            for neighborID in nlist_filtered_ids[iID]
                if dot((displacements[nlist[iID][neighborID], :] - displacements[iID, :]), bond_norm[iID][neighborID, :]) > 0
                    bond_mod[neighborID, :] .= 0
                else
                    bond_mod[neighborID, :] = bond_norm[iID][neighborID, :]
                end
            end
        end

        force_densities[iID, :] .+= transpose(sum(bond_damage[iID][:] .* bond_force[iID][:, :] .* bond_mod .* volume[nlist[iID][:]], dims=1))

        # force_densities[nlist[iID][:], :][:] .-= bond_damage[iID][:] .* bond_force[iID][:] .* bond_mod[:] .* volume[iID]
        force_densities[nlist[iID][:], :] .-= bond_damage[iID][:] .* bond_force[iID][:, :] .* bond_mod[:, :] .* volume[iID]
    end
    return force_densities
end

"""
    distribute_forces(nodes::Union{SubArray,Vector{Int64}}, nlist::SubArray, bond_force::SubArray, volume::SubArray, bond_damage::SubArray, force_densities::SubArray)

Distribute the forces on the nodes

# Arguments
- `nodes::Union{SubArray,Vector{Int64}}`: The nodes.
- `nlist::SubArray`: The neighbor list.
- `bond_force::SubArray`: The bond forces.
- `volume::SubArray`: The volumes.
- `bond_damage::SubArray`: The bond damage.
- `force_densities::SubArray`: The force densities.
# Returns
- `force_densities::SubArray`: The force densities.
"""
function distribute_forces(nodes::Union{SubArray,Vector{Int64}}, nlist::SubArray, bond_force::SubArray, volume::SubArray, bond_damage::SubArray, force_densities::SubArray)
    for iID in nodes
        force_densities[iID, :] .+= transpose(sum(bond_damage[iID][:] .* bond_force[iID][:, :] .* volume[nlist[iID][:]], dims=1))

        force_densities[nlist[iID][:], :] .-= bond_damage[iID][:] .* bond_force[iID][:, :] .* volume[iID]
    end
    return force_densities
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
    if size(matrix) == (2, 2)
        return [matrix[1, 1]; matrix[2, 2]; 0.5 * (matrix[1, 2] + matrix[2, 1])]
    elseif size(matrix) == (3, 3)
        return [matrix[1, 1]; matrix[2, 2]; matrix[3, 3]; 0.5 * (matrix[2, 3] + matrix[3, 2]); 0.5 * (matrix[1, 3] + matrix[3, 1]); 0.5 * (matrix[1, 2] + matrix[2, 1])]
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
        return [voigt[1] voigt[3]; voigt[3] voigt[2]]
    elseif length(voigt) == 6
        return [voigt[1] voigt[6] voigt[5]
            voigt[6] voigt[2] voigt[4]
            voigt[5] voigt[4] voigt[3]]
    else
        @error "Unsupported matrix size for voigt_to_matrix"
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
                @error "Model definition is missing; plain stress or plain strain has to be defined for 2D"
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
- `coor::Union{Vector{Int64},Vector{Float64}}`: Coordinate of the current point.
- `stress::Float64`: stresses to be modified.

# Returns
- `stress`::Float64: the modified stresses.
"""
function flaw_function(params::Dict, coor::Union{Vector{Int64},Vector{Float64}}, stress::Float64)
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
        if (1 < flaw_magnitude) || (flaw_magnitude < 0)
            @error "Flaw Magnitude should be between 0 and 1"
            return nothing
        end
        flaw_location[1] = params["Flaw Function"]["Flaw Location X"]
        flaw_location[2] = params["Flaw Function"]["Flaw Location Y"]
        if haskey(params["Flaw Function"], "Flaw location Z") && length(coor) == 3
            flaw_location[3] = params["Flaw Function"]["Flaw Location Z"]
        end
        modified_stress = stress * (1 - flaw_magnitude * exp(-norm(coor - flaw_location) * norm(coor - flaw_location) / flaw_size / flaw_size))

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
    get_symmmetry(material::Dict)

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
function get_symmmetry(material::Dict)
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
    get_von_mises_stress(von_Mises_stress::Float64, dof::Int64, stress_NP1::Matrix{Float64})

# Arguments
- `von_Mises_stress::Float64`: Von Mises stress
- `dof::Int64`: Degree of freedom.
- `stress_NP1::Matrix{Float64}`: Stress.
# returns
- `spherical_stress_NP1::Float64`: Spherical stress
- `deviatoric_stress_NP1::Matrix{Float64}`: Deviatoric stress
"""
function get_von_mises_stress(von_Mises_stress::Float64, dof::Int64, stress_NP1::Matrix{Float64})

    spherical_stress_NP1 = sum(stress_NP1[i, i] for i in 1:dof) / 3
    deviatoric_stress_NP1 = stress_NP1[:, :] - spherical_stress_NP1 .* I(dof)

    von_Mises_stress = sqrt(3.0 / 2.0 * sum(deviatoric_stress_NP1[:, :] .* deviatoric_stress_NP1[:, :]))

    return von_Mises_stress, spherical_stress_NP1, deviatoric_stress_NP1
end