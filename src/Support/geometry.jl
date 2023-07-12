
module Geometry
export bond_geometry
function bond_geometry(nnodes, dof, nlist, coor, bondgeom)
    for iID in 1:nnodes
        nneigbors = length(nlist[iID])
        for jID in 1:nneigbors
            bondgeom[iID, jID, 1:dof] = coor[nlist[iID, jID], :] - coor[iID, :]
            bondgeom[iID, jID, dof+1] = norm(bondgeom[iID, jID, 1:dof])
        end
    end
    return bondgeom
end
end