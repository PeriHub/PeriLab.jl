

function distribute_forces(nnodes, nneighbors, nlist, bond_force, volume)
    for iID in 1:nnodes
        for jID in nneighbors[iID]
            forces[iID] += bond_force[iID][jID, :] * volume[nlist[jID]]
            forces[nlist[jID]] -= bond_force[iID][jID, :] * volume[iID]
        end
    end
end