
function bond_geometry(data)
    coordinates = data["coordinates"]
    neigborlist = data["Neighborlist"]
    bondDamage = data["Bond Damage"]
    bondGeometry = data["Bond Geometry"] #dx,dy,dz,len
    invK = data["Inverse Shape Tensor"]

    nnodes = data.nnodes
    dof = data.dof
    for iID in 1:nnodes
        nneighbors = length(neigborlist[iID])
        for jID in 1:nneighbors
            pythagoras = 0.0
            for i in 1:dof
                bondGeometry[i] = coordinates[(iID-1)*dof+i] - coordinates[(neigborlist[jID]-1)*dof+i]
                pythagoras += bondGeometry[i] * bondGeometry[i]
            end
            bondGeometry[dof+1] = sqrt(pythagoras)
        end

    end
    return data
end
