import Geometry
function current_stretch(datamanger)
    nnodes = datamanger.get_nnodes()
    dof = datamanger.get_dof()
    nlist = datamanger.get_nlist()
    coor = datamanger.get_field("Deformed Coordinates")
    bondgeom = datamanger.get_field("Deformed Bond Geometry")
    bondgeom = bond_geometry(nnodes, dof, nlist, coor, bondgeom)
end