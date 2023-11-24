# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Ordinary

export compute_dilatation
export compute_weighted_volume
"""
taken from Peridigm -> but adding the bond_damage; this is missing in Peridigm, but should be there
"""
function compute_weighted_volume(nodes::Union{SubArray,Vector{Int64}}, nlist::SubArray, undeformed_bond::SubArray, bond_damage::SubArray, omega::SubArray, volume::SubArray)

    if length(nodes) == 0
        return Float64[]
    end
    weighted_volume::Vector{Float64} = zeros(Float64, maximum(nodes))

    for iID in nodes
        # in Peridigm the weighted volume is for some reason independend from damages
        weighted_volume[iID] = sum(omega[iID][:] .* bond_damage[iID][:] .* undeformed_bond[iID][:, end] .* bond_damage[iID][:] .* undeformed_bond[iID][:, end] .* volume[nlist[iID][:]])

    end
    return weighted_volume
end


"""
    compute_dilatation(nodes::Union{SubArray,Vector{Int64}}, nneighbors::SubArray, nlist::SubArray, undeformed_bond::SubArray, deformed_bond::SubArray, bond_damage::SubArray, volume::SubArray, weighted_volume::Vector{Float64}, omega::SubArray)

Calculate the dilatation for each node.

# Arguments
- `nodes::Union{SubArray,Vector{Int64}}`: Nodes.
- `nneighbors::SubArray`: Number of neighbors.
- `nlist::SubArray`: Neighbor list.
- `undeformed_bond::SubArray`: Bond geometry.
- `deformed_bond::SubArray`: Deformed bond geometry.
- `bond_damage::SubArray`: Bond damage.
- `volume::SubArray`: Volume.
- `weighted_volume::Vector{Float64}`: Weighted volume.
- `omega::SubArray`: Influence function.
# Returns
- `theta::Vector{Float64}`: Dilatation.
"""
function compute_dilatation(nodes::Union{SubArray,Vector{Int64}}, nneighbors::SubArray, nlist::SubArray, undeformed_bond::SubArray, deformed_bond::SubArray, bond_damage::SubArray, volume::SubArray, weighted_volume::Vector{Float64}, omega::SubArray)
    # not optimal, because of many zeros, but simpler, because it avoids reorganization. Part of potential optimization
    if length(nodes) == 0
        return Float64[]
    end
    theta::Vector{Float64} = zeros(Float64, maximum(nodes))
    for iID in nodes
        if weighted_volume[iID] == 0
            @warn "Weighted volume is zero for local point ID: $iID"
            continue
        end
        theta[iID] = 3.0 * sum(omega[iID][jID] *
                               bond_damage[iID][jID] * undeformed_bond[iID][jID, end] *
                               (deformed_bond[iID][jID, end] - undeformed_bond[iID][jID, end]) *
                               volume[nlist[iID][jID]] / weighted_volume[iID]
                               for jID in 1:nneighbors[iID])
    end
    return theta
end

end