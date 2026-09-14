# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
"""
    Text_Mesh

Text mesh importer. Reads a text file and converts the element
connectivities into the peridynamic mesh representation (center of each
element as a node, element volume, block id). Node sets are read separately
from the parameters via `get_node_sets`.
"""
module Text_Mesh

using ....PeriLabExceptions: @abort
using DataFrames
using ...IO: csv_reader

export mesh_import_name
export read_mesh

"""
    mesh_import_name()

Gives the mesh importer name. It is needed for comparison with the yaml input deck.

# Returns
- `name::String`: The name of the mesh importer.
"""
function mesh_import_name()
    return "Text File"
end

"""
    read_mesh(params::Dict, filename::String)

Reads a text mesh file and returns the mesh data as a DataFrame.

# Arguments
- `params::Dict`: The parameters.
- `filename::String`: The path to the text mesh file.
# Returns
- `mesh::DataFrame`: The mesh data as a DataFrame.
"""
function read_mesh(params::Dict, filename::String)
    return csv_reader(filename)
end

end # module Text_Mesh
