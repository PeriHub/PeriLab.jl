function compute_weighted_volume(nnodes, bondlength, nlist, bonddamage, omega, volume, weightedVol)

    for iID in 1:nnodes
        nneighbors = length(nlist[iID])
        weightedVol[iID] = 0.0
        for jID in 1:nneighbors
            localID = nlist[iID, jID]
            weightedVol[iID] += omega[localID] * bonddamage[localID, jID] * bondlength[iID, jID, end] * bondlength[iID, jID, end] * volume[localID]
        end

    end
    return m
end

function compute_dilatation(nnodes, bond_geometry, bond_stretch, bond_damage, volume, weighted_volume, omega, dof)
    theta = zeros(nnodes)
    for iID in 1:nnodes
        nbonds = length(bond_geometry)
        for jID in 1:nbonds
            e = bond_stretch[iID][jID, dof+1] - bond_geometry[iID][jID, dof+1]
            theta[iID] += 3.0 * omega[iID][jID] * bond_damage[iID][jID] * bond_geometry[iID][jID, dof+1] * e * volume[iID] / weighted_volume[iID]
        end
    end
    return theta
end

function elastic(nnodes, bond_geometry, bond_stretch, bond_damage, theta, weighted_volume, omega, dof)

    for iID in 1:nnodes
        alpha = 15.0 * MU / weighted_volume[iID]
        beta = 3.0 * K / weighted_volume[iID]
        nbonds = length(bond_geometry)
        for jID in 1:nbonds

            c1 = omega * theta[iID] * (beta - alpha / 3.0)
            bond_force[iID][jID] = bond_damage[iID][jID] * omega[iID][jID] * (c1 * bond_geometry[iID][jID, dof+1] + alpha * bond_stretch[iID][jID, dof+1])
        end
    end
    return bond_force
end