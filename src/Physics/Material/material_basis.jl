# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

function get_all_elastic_moduli(parameter::Dict{Any,Any})
    if parameter == Dict()
        return parameter
    end
    K::Float32 = 0
    E::Float32 = 0
    nu::Float32 = 0
    G::Float32 = 0
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
        nu = (3 * K - 2 * G) / (6 * K - 2 * G)
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

    if bulk + Youngs + shear + Poissons > 2
        @warn "More than two parameters for isotropic material. Please check consistency"
    end

    parameter["Bulk Modulus"] = K
    parameter["Young's Modulus"] = E
    parameter["Shear Modulus"] = G
    parameter["Poisson's Ratio"] = nu
    return parameter
end

function get_Hook_matrix(parameter, symmetry, dof)
    """https://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_plane_stress.cfm"""

    if occursin("isotropic", symmetry)
        matrix = zeros(Float32, 2 * dof, 2 * dof)
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
            matrix[1, 1] = nu * temp
            matrix[3, 2] = nu * temp
            matrix[2, 3] = nu * temp
            matrix[4, 4] = G
            matrix[5, 5] = G
            matrix[6, 6] = G
            return matrix
        elseif occursin("plain_strain", symmetry)
            matrix = zeros(Float32, dof + 1, dof + 1)
            matrix[1, 1] = E * (1 - nu) * temp
            matrix[2, 2] = E * (1 - nu) * temp
            matrix[3, 3] = G
            matrix[1, 2] = E * nu * temp
            matrix[2, 1] = E * nu * temp
            return matrix
        elseif occursin("plain_stress", symmetry)
            matrix = zeros(Float32, dof + 1, dof + 1)
            matrix[1, 1] = E / (1 - nu * nu)
            matrix[1, 2] = E * nu / (1 - nu * nu)
            matrix[2, 1] = E * nu / (1 - nu * nu)
            matrix[2, 2] = E / (1 - nu * nu)
            matrix[3, 3] = G
            return matrix
        else
            @error "2D model defintion is missing; plain_stress or plain_strain "
        end
        if occursin("anisotropic", symmetry)
            anisoMatrix = zeros(Float32, 6, 6)
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
            elseif occursin("plain_strain", symmetry)
                matrix = zeros(Float32, dof + 1, dof + 1)
                matrix[1:2, 1:2] = anisoMatrix[1:2, 1:2]
                matrix[3, 1:2] = anisoMatrix[6, 1:2]
                matrix[1:2, 3] = anisoMatrix[1:2, 6]
                matrix[3, 3] = anisoMatrix[6, 6]
                return matrix
            elseif occursin("plain_stress", symmetry)
                matrix = zeros(Float32, dof + 1, dof + 1)
                invAniso = inv(anisoMatrix)
                matrix[1:2, 1:2] = invAniso[1:2, 1:2]
                matrix[3, 1:2] = invAniso[6, 1:2]
                matrix[1:2, 3] = invAniso[1:2, 6]
                matrix[3, 3] = invAniso[6, 6]
                return inv(matrix)
            else
                @error "2D model defintion is missing; plain_stress or plain_strain "
            end
        end
    end
    return matrix
end

function distribute_forces(nnodes, nlist, bond_force, volume, force_densities)
    for iID in nnodes
        for (jID, neighborID) in enumerate(nlist[iID])
            force_densities[iID, :] += bond_force[iID][jID, :] .* volume[neighborID]
            force_densities[neighborID, :] -= bond_force[iID][jID, :] .* volume[iID]
        end
    end
    return force_densities
end