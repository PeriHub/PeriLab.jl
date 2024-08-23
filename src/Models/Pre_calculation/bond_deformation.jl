# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Bond_Deformation
include("../../Support/geometry.jl")
using .Geometry
export compute

"""
    compute(datamanager::Module, nodes::Union{SubArray,Vector{Int64}})

Compute the bond deformation.

# Arguments
- `datamanager`: Datamanager.
- `nodes`: List of nodes.
- `time`: Time.
# Returns
- `datamanager`: Datamanager.
"""
function compute(datamanager::Module, nodes::Union{SubArray,Vector{Int64}})
    nlist = datamanager.get_nlist()
    deformed_coor = datamanager.get_field("Deformed Coordinates", "NP1")
    deformed_bond = datamanager.get_field("Deformed Bond Geometry", "NP1")
    deformed_bond_length = datamanager.get_field("Deformed Bond Length", "NP1")
    deformed_bond, deformed_bond_length = Geometry.bond_geometry(
        nodes,
        nlist,
        deformed_coor,
        deformed_bond,
        deformed_bond_length,
    )
    return datamanager
end


end
