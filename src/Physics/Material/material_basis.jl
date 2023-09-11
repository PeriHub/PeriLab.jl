

function distribute_forces(nnodes, nlist, bond_force, volume, force_densities)
    for iID in nnodes
        for (jID, neighborID) in enumerate(nlist[iID])
            force_densities[iID, :] += bond_force[iID][jID, :] .* volume[neighborID]
            force_densities[neighborID, :] -= bond_force[iID][jID, :] .* volume[iID]
        end
    end
    return force_densities
end