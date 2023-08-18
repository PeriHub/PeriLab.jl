module Geometry
using LinearAlgebra
export bond_geometry
function bond_geometry(nnodes, dof, nlist, coor, bondgeom)
    for iID in 1:nnodes
        for jID in eachindex(nlist[iID])
            bondgeom[iID][jID, 1:dof] = coor[nlist[iID][jID], :] - coor[iID, :]
            bondgeom[iID][jID, dof+1] = norm(bondgeom[iID][jID, 1:dof])
        end
    end
    return bondgeom
end
end