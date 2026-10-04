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
folder -- it just needs a `mesh_import_name()` and a `read_mesh(input::PeriLabInput,
filename)` function.
"""
module Mesh_Import

using ...Data_Manager
using ...PeriLabExceptions: @abort
using ...ModuleLoader: find_module_files, create_module_specifics
using ...InputDeck: PeriLabInput

include("Mesh_Volume.jl")
using .Mesh_Volume

global module_list = find_module_files(@__DIR__, "mesh_import_name")
for mod in module_list
    include(mod["File"])
end

export read_mesh

"""
    read_mesh(input::PeriLabInput, path::String)

Reads the mesh file of the configured discretization type and returns the mesh
data as a DataFrame (and, for some importers, additional node sets).

Dispatches to the importer whose `mesh_import_name()` matches
`Discretization` `Type`.

# Arguments
- `input::PeriLabInput`: The typed input deck.
- `path::String`: The path to the mesh file.
# Returns
- `mesh::DataFrame`: The mesh data as a DataFrame.
- `nsets::Dict`: The node sets (if any).
"""
function read_mesh(input::PeriLabInput, path::String)
    mesh_name = input.sections.discretization.input_mesh_file
    mesh_path = joinpath(path, mesh_name)
    if !isfile(mesh_path)
        @abort "Mesh file $mesh_path does not exist"
        return
    end

    type = input.sections.discretization.type

    importer = create_module_specifics(type,
                                       module_list,
                                       @__MODULE__,
                                       "mesh_import_name")
    if isnothing(importer)
        @abort "No mesh importer for type '$type' exists."
        return nothing
    end

    @info "Read mesh file $mesh_path using importer $(importer.mesh_import_name())"
    return importer.read_mesh(input, mesh_path)
end

end # module Mesh_Import
