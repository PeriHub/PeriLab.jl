include("../material_basis.jl")
module Ordinary

export compute_force
export material_name

function material_name()
    return "PD Solid Elastic"
end

function compute_weighted_volume(nodes, nneighbors, nlist, bond_geometry, bond_damage, omega, volume)
    """
    taken from Peridigm
    """
    # not optimal, because of many zeros, but simpler, because it avoids reorganization. Part of potential optimization
    weighted_volume = zeros(Float32, maximum(nodes))

    for iID in nodes
        for jID in 1:nneighbors[iID]
            weighted_volume[iID] += omega[iID] * bond_damage[iID][jID] * bond_geometry[iID][jID, end] * bond_geometry[iID][jID, end] * volume[nlist[iID][jID]]
        end
    end
    return weighted_volume
end

function compute_dilatation(nodes, nneighbors, nlist, bond_geometry, deformed_bond, bond_damage, volume, weighted_volume, omega)
    # not optimal, because of many zeros, but simpler, because it avoids reorganization. Part of potential optimization
    theta = zeros(Float32, maximum(nodes))
    for iID in nodes
        theta[iID] = 3.0 * omega[iID] * sum(
                         bond_damage[iID][jID] * bond_geometry[iID][jID, end] *
                         (deformed_bond[iID][jID, end] - bond_geometry[iID][jID, end]) *
                         volume[nlist[iID][jID]] / weighted_volume[iID]
                         for jID in 1:nneighbors[iID])
    end
    return theta
end

end