# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Ordinary

export compute_dilatation
export compute_weighted_volume

function compute_weighted_volume(nodes::SubArray, nneighbors, nlist, bond_geometry, bond_damage, omega, volume)
    """
    taken from Peridigm
    """
    # not optimal, because of many zeros, but simpler, because it avoids reorganization. Part of potential optimization
    weighted_volume = zeros(Float32, maximum(nodes))

    for iID in nodes
        for jID in 1:nneighbors[iID]
            weighted_volume[iID] += omega[iID][jID] * bond_damage[iID][jID] * bond_geometry[iID][jID, end] * bond_geometry[iID][jID, end] * volume[nlist[iID][jID]]
        end
    end
    return weighted_volume
end

function compute_dilatation(nodes::SubArray, nneighbors, nlist, bond_geometry, deformed_bond, bond_damage, volume, weighted_volume, omega)
    # not optimal, because of many zeros, but simpler, because it avoids reorganization. Part of potential optimization
    theta = zeros(Float32, maximum(nodes))
    for iID in nodes
        theta[iID] = 3.0 * sum(omega[iID][jID] *
                               bond_damage[iID][jID] * bond_geometry[iID][jID, end] *
                               (deformed_bond[iID][jID, end] - bond_geometry[iID][jID, end]) *
                               volume[nlist[iID][jID]] / weighted_volume[iID]
                               for jID in 1:nneighbors[iID])
    end
    return theta
end

end