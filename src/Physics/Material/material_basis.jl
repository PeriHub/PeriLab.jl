

function get_all_elastic_moduli(parameter)
    bulkM::Float32 = 0
    E::Float32 = 0
    nu::Float32 = 0
    ShearM::Float32 = 0
    bulk = false
    Youngs = false
    shear = false
    Poissons = false
    if "Bulk Modulus" in keys(parameter)
        K = parameter["Bulk Modulus"]
        bulk = true
    end
    if "Young's Modulus" in keys(parameter)
        E = parameter["Young's Modulus"]
        Youngs = true
    end
    if "Shear Modulus" in keys(parameter)
        nu = parameter["Shear Modulus"]
        shear = true
    end
    if "Poisson's Ratio" in keys(parameter)
        G = parameter["Poisson's Ratio"]
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
    if bulk + Youngs + shear + Poissons < 2
        @error "Minimum of two parameters are needed for isotropic material"
    end
    parameter["Bulk Modulus"] = K
    parameter["Young's Modulus"] = E
    parameter["Shear Modulus"] = G
    parameter["Poisson's Ratio"] = nu
    return parameter
end
function get_Hook_matrix(parameter, dof)
    matrix = zeros(Float32, dof, dof)

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