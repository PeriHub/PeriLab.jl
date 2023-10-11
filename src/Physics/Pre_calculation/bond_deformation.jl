# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Bond_Deformation
include("../../Support/geometry.jl")
using .Geometry
export compute

function compute(datamanager, nodes, time)
    dof = datamanager.get_dof()
    nlist = datamanager.get_nlist()
    defCoor = datamanager.get_field("Deformed Coordinates", "NP1")
    bond_defNP1 = datamanager.get_field("Deformed Bond Geometry", "NP1")
    bond_defNP1 = Geometry.bond_geometry(nodes, dof, nlist, defCoor, bond_defNP1)
    #if time == 0
    #    bondGeom = datamanager.get_field("Bond Geometry")
    #    bond_defN = datamanager.get_field("Deformed Bond Geometry", "N")
    #
    #    for iID in nodes
    #        for jID in eachindex(bond_defN[iID][:, 1])
    #            bond_defN[iID][jID, 1:end] = bondGeom[iID][jID, 1:end]
    #        end
    #    end
    #end
    return datamanager
end


end