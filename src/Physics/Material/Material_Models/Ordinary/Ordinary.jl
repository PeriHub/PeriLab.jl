# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Ordinary

export compute_dilatation
export compute_weighted_volume
"""
taken from Peridigm -> but adding the bond_damage; this is missing in Peridigm, but should be there
"""
function compute_weighted_volume(nodes::Union{SubArray,Vector{Int64}}, nneighbors::SubArray, nlist::SubArray, bond_geometry::SubArray, bond_damage::SubArray, omega::SubArray, volume::SubArray)

    if length(nodes) == 0
        return Float64[]
    end
    weighted_volume = zeros(Float64, maximum(nodes))

    for iID in nodes
        # in Peridigm the weighted volume is for some reason independend from damages
        weighted_volume[iID] = sum(omega[iID][:] .* bond_damage[iID][:] .* bond_geometry[iID][:, end] .* bond_damage[iID][:] .* bond_geometry[iID][:, end] .* volume[nlist[iID][:]])

    end
    return weighted_volume
end

function compute_dilatation(nodes::Union{SubArray,Vector{Int64}}, nneighbors::SubArray, nlist::SubArray, bond_geometry::SubArray, deformed_bond::SubArray, bond_damage::SubArray, volume::SubArray, weighted_volume::Vector{Float64}, omega::SubArray)
    # not optimal, because of many zeros, but simpler, because it avoids reorganization. Part of potential optimization
    if length(nodes) == 0
        return Float64[]
    end
    theta = zeros(Float64, maximum(nodes))
    for iID in nodes
        if weighted_volume[iID] == 0
            @warn "Weighted volume is zero for local point ID: $iID"
            continue
        end
        theta[iID] = 3.0 * sum(omega[iID][jID] *
                               bond_damage[iID][jID] * bond_geometry[iID][jID, end] *
                               (deformed_bond[iID][jID, end] - bond_geometry[iID][jID, end]) *
                               volume[nlist[iID][jID]] / weighted_volume[iID]
                               for jID in 1:nneighbors[iID])
    end
    return theta
end

end