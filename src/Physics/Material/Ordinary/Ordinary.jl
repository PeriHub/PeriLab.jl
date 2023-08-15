function compute_weighted_volume(nnodes, nlist, bond_geometry, bond_damage, omega, volume)
    weighted_volume = zeros(nnodes)
    for iID in 1:nnodes
        nbonds = size(bond_geometry[iID])[1]
        for jID in 1:nbonds
            weighted_volume[iID] += omega[iID][jID] * bond_damage[iID][jID] * bond_geometry[iID][jID, end] * bond_geometry[iID][jID, end] * volume[nlist[iID][jID]]
        end

    end
    return weighted_volume
end

function compute_dilatation(nnodes, bond_geometry, bond_stretch, bond_damage, volume, weighted_volume, omega)
    theta = zeros(nnodes)
    for iID in 1:nnodes
        nbonds = size(bond_geometry[iID])[1]
        for jID in 1:nbonds
            e = bond_stretch[iID][jID, end] - bond_geometry[iID][jID, end]
            theta[iID] += 3.0 * omega[iID][jID] * bond_damage[iID][jID] * bond_geometry[iID][jID, end] * e * volume[nlist[iID][jID]] / weighted_volume[iID]
        end
    end
    return theta
end

function elastic(nnodes, bond_geometry, bond_stretch, bond_damage, theta, weighted_volume, omega)

    for iID in 1:nnodes
        alpha = 15.0 * MU / weighted_volume[iID]
        beta = 3.0 * K / weighted_volume[iID]
        nbonds = size(bond_geometry[iID])[1]
        for jID in 1:nbonds
            c1 = omega * theta[iID] * (beta - alpha / 3.0)
            bond_force[iID][jID] = bond_damage[iID][jID] * omega[iID][jID] * (c1 * bond_geometry[iID][jID, end] + alpha * bond_stretch[iID][jID, end])
        end
    end
    return bond_force
end