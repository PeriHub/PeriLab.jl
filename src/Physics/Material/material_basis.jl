# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

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

    if occursin("isotropic", symmetry)
        matrix = zeros(Float64, 2 * dof, 2 * dof)
        nu = parameter["Poisson's Ratio"]
        E = parameter["Young's Modulus"]
        G = parameter["Shear Modulus"]
        temp = E / ((1 + nu) * (1 - 2 * nu))

        if dof == 3
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
            matrix = zeros(Float64, dof + 1, dof + 1)
            matrix[1, 1] = (1 - nu) * temp
            matrix[2, 2] = (1 - nu) * temp
            matrix[3, 3] = G
            matrix[1, 2] = nu * temp
            matrix[2, 1] = nu * temp
            return matrix
        elseif occursin("plane stress", symmetry)
            matrix = zeros(Float64, dof + 1, dof + 1)
            matrix[1, 1] = E / (1 - nu * nu)
            matrix[1, 2] = E * nu / (1 - nu * nu)
            matrix[2, 1] = E * nu / (1 - nu * nu)
            matrix[2, 2] = E / (1 - nu * nu)
            matrix[3, 3] = G
            return matrix
        else
            @error "2D model defintion is missing; plain stress or plain strain "
        end
        if occursin("anisotropic", symmetry)
            anisoMatrix = zeros(Float64, 6, 6)
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
                matrix = zeros(Float64, dof + 1, dof + 1)
                matrix[1:2, 1:2] = anisoMatrix[1:2, 1:2]
                matrix[3, 1:2] = anisoMatrix[6, 1:2]
                matrix[1:2, 3] = anisoMatrix[1:2, 6]
                matrix[3, 3] = anisoMatrix[6, 6]
                return matrix
            elseif occursin("plane stress", symmetry)
                matrix = zeros(Float64, dof + 1, dof + 1)
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
    else
        matrix = zeros(Float64, dof + 1, dof + 1)
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
        for (jID, neighborID) in enumerate(nlist[iID])
            force_densities[iID, :] += bond_damage[iID][jID] .* bond_force[iID][jID, :] .* volume[neighborID]
            force_densities[neighborID, :] -= bond_damage[iID][jID] .* bond_force[iID][jID, :] .* volume[iID]
        end
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
    end
end

"""
    voigt_to_matrix(voigt)

    Convert a Voigt notation (6x1 or 3x1 vector) to a 2x2 or 3x3 matrix

    # Arguments
    - `voigt::Vector{Float64}`: The Voigt notation.
    # Returns
    - `matrix::Matrix{Float64}`: The matrix.
"""
function voigt_to_matrix(voigt)
    if length(voigt) == 3
        return [voigt[1] voigt[3]; voigt[3] voigt[2]]
    elseif length(voigt) == 6
        return [voigt[1] voigt[6] voigt[5]
            voigt[6] voigt[2] voigt[4]
            voigt[5] voigt[4] voigt[3]]
    else
        @error "Unsupported matrix size for voigt_to_matrix"
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
