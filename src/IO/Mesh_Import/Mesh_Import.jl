# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
"""
    Mesh_Import

Mesh import factory. Loads the individual mesh importer modules (Abaqus, Gmsh,
Gcode, ...) which are discoverable via the ModuleLoader and dispatches the
`read_mesh` call to the one whose name matches the `Discretization.Type` from
the input deck.

This mirrors the pattern used by the other feature modules (Additive, Damage,
FEM, ...) so new importers can be added by dropping a new module into this
folder -- it just needs a `mesh_import_name()` and a `read_mesh(params,
filename)` function.
"""
module Mesh_Import

using ...Data_Manager
using ...PeriLabExceptions: @abort
using ...ModuleLoader: find_module_files, create_module_specifics

global module_list = find_module_files(@__DIR__, "mesh_import_name")
for mod in module_list
    include(mod["File"])
end

using .Mesh_Volume

export read_mesh

"""
    read_mesh(params::Dict, path::String, size::Int64, silent::Bool)

Reads the mesh file of the configured discretization type and returns the mesh
data as a DataFrame (and, for some importers, additional node sets).

Dispatches to the importer whose `mesh_import_name()` matches
`params["Discretization"]["Type"]`.

# Arguments
- `params::Dict`: The parameters.
- `path::String`: The path to the mesh file.
- `size::Int64`: The number of ranks.
- `silent::Bool`: Whether to run in silent mode.
# Returns
- `mesh::DataFrame`: The mesh data as a DataFrame.
- `nsets::Dict`: The node sets (if any).
"""
function read_mesh(params::Dict, path::String)
    mesh_name = get_mesh_name(params)
    type = params["Discretization"]["Type"]

    importer = create_module_specifics(type,
                                       module_list,
                                       @__MODULE__,
                                       "mesh_import_name")
    if isnothing(importer)
        @abort "No mesh importer for type '$type' exists."
        return nothing
    end

    @info "Read mesh file $(joinpath(path, mesh_name))"
    return importer.read_mesh(params, joinpath(path, mesh_name))
end

end # module Mesh_Import
