# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
"""
    Gmsh_Mesh

Gmsh mesh importer. Reads a Gmsh `.msh` file and converts the element
connectivities into the peridynamic mesh representation (center of each
element as a node, element volume, block id). Node sets are read separately
from the parameters via `get_node_sets`.
"""
module Gmsh_Mesh

using ...PeriLabExceptions: @abort
using DataFrames
using Gmsh: gmsh
using .Mesh_Volume: tetrahedron_volume, area_of_polygon

export mesh_import_name
export read_mesh

"""
    mesh_import_name()

Gives the mesh importer name. It is needed for comparison with the yaml input deck.

# Returns
- `name::String`: The name of the mesh importer.
"""
function mesh_import_name()
    return "Gmsh"
end

"""
    read_mesh(params::Dict, filename::String)

Reads a Gmsh mesh file and returns the mesh data as a DataFrame.

# Arguments
- `params::Dict`: The parameters.
- `filename::String`: The path to the Gmsh `.msh` file.
# Returns
- `mesh::DataFrame`: The mesh data as a DataFrame.
"""
function read_mesh(params::Dict, filename::String)
    gmsh.initialize()

    gmsh.open(filename)

    dof = 3 # only 3d supported currntly

    if gmsh.model.mesh.getElements(3)[2] == []
        dof = 2
    end

    num_elements = length(gmsh.model.mesh.getElements(dof)[2][1])

    mesh_df = ifelse(dof == 2,
                     DataFrame(x = Array{Float64,1}(undef, num_elements),
                               y = Array{Float64,1}(undef, num_elements),
                               volume = Array{Float64,1}(undef, num_elements),
                               block_id = Array{Int64,1}(undef, num_elements)),
                     DataFrame(x = Array{Float64,1}(undef, num_elements),
                               y = Array{Float64,1}(undef, num_elements),
                               z = Array{Float64,1}(undef, num_elements),
                               volume = Array{Float64,1}(undef, num_elements),
                               block_id = Array{Int64,1}(undef, num_elements)))

    ids = gmsh.model.getPhysicalGroups(3)
    no_groups = false
    if ids == []
        ids = [1]
        no_groups = true
    end
    node_id = 1
    block_id = 0
    block_names = []
    for id in ids
        if no_groups
            block_id = 1
        else
            block_id = Int64(id[2])
        end
        name = gmsh.model.getPhysicalName(dof, block_id)
        push!(block_names, name)
        if no_groups
            element_tags = gmsh.model.mesh.gmsh.model.mesh.getElements(dof)[2][1]
        else
            element_tags = gmsh.model.mesh.gmsh.model.mesh.getElements(dof, block_id)[2][1]
        end
        for element_tag in element_tags
            element = gmsh.model.mesh.getElement(element_tag)
            node_tags = element[2]
            nodes = []
            for node_tag in node_tags
                node = gmsh.model.mesh.getNode(node_tag)[1]
                push!(nodes, node)
            end
            center = sum(nodes) / length(nodes)
            if dof == 2
                volume = area_of_polygon(nodes)
                mesh_df[node_id, :] = [center[1], center[2], volume, block_id]
            else
                volume = tetrahedron_volume(nodes)
                mesh_df[node_id,
                :] = [
                    center[1],
                    center[2],
                    center[3],
                    volume,
                    block_id
                ]
            end
            node_id += 1
        end
    end
    @info "Found $(block_id) block(s)"
    @info "Blocks: $block_names"

    return mesh_df
end

end # module Gmsh_Mesh
