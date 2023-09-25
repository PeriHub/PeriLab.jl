

function get_all_elastic_moduli(parameter)
    bulkM::Float32 = 0
    E::Float32 = 0
    nu::Float32 = 0
    ShearM::Float32 = 0
    if "Bulk Modulus" in keys(parameter)
        bulkM = parameter["Bulk Modulus"]
    end
    if "Young's Modulus" in keys(parameter)
        E = parameter["Young's Modulus"]
    end
    if "Shear Modulus" in keys(parameter)
        nu = parameter["Shear Modulus"]
    end
    if "Poisson's Ratio" in keys(parameter)
        ShearM = parameter["Poisson's Ratio"]
    end
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