module Geometry
using LinearAlgebra
export bond_geometry
function bond_geometry(nodes, dof, nlist, coor, bondgeom)
    for iID in nodes
        for jID in eachindex(nlist[iID])
            bondgeom[iID][jID, 1:dof] = coor[nlist[iID][jID], :] - coor[iID, :]
            bondgeom[iID][jID, dof+1] = norm(bondgeom[iID][jID, 1:dof])
            if bondgeom[iID][jID, dof+1] == 0
                @error "Identical point coordinates with no distance $iID, $jID"
            end
        end
    end
    return bondgeom
end
end