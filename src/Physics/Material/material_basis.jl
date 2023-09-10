

function distribute_forces(nnodes, nneighbors, nlist, bond_force, volume, force_densities)
    for iID in nnodes
        for jID in nneighbors[iID]
            force_densities[iID, :] += bond_force[iID][jID, :] .* volume[nlist[iID][jID]]
            force_densities[nlist[iID][jID], :] -= bond_force[iID][jID, :] .* volume[iID]
        end
    end
    return force_densities
end