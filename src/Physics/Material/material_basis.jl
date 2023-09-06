

function distribute_forces(nnodes, nneighbors, nlist, bond_force, volume, forces)
    for iID in 1:nnodes
        for jID in nneighbors[iID]
            forces[iID, :] += bond_force[iID][jID, :] * volume[nlist[iID][jID]]
            forces[nlist[iID][jID], :] -= bond_force[iID][jID, :] * volume[iID]
        end
    end
    return forces
end