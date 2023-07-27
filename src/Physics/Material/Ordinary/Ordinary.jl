function compute_weighted_volume(nnodes, bondlength, nlist, bonddamage, omega, volume, weightedVol)

    for iID in 1:nnodes
        nneighbors = length(nlist[iID])
        weightedVol[iID] = 0.0
        for jID in 1:nneighbors
            localID = nlist[iID, jID]
            weightedVol[iID] += omega[localID] * bonddamage[localID, jID] * bondlength[iID, jID, end] * bondlength[iID, jID, end] * volune[localID]
        end

    end
    return m
end
