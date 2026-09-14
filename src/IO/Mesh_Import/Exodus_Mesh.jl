# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
"""
    Exodus_Mesh

Exodus mesh importer. Reads an Exodus file and converts the element
connectivities into the peridynamic mesh representation (center of each
element as a node, element volume, block id). Node sets are read separately
from the parameters via `get_node_sets`.
"""
module Exodus_Mesh

using ....PeriLabExceptions: @abort
using DataFrames
using ..Mesh_Volume: tetrahedron_volume

export mesh_import_name
export read_mesh

"""
    mesh_import_name()

Gives the mesh importer name. It is needed for comparison with the yaml input deck.

# Returns
- `name::String`: The name of the mesh importer.
"""
function mesh_import_name()
    return "Exodus"
end

"""
    read_mesh(params::Dict, filename::String)

Reads an Exodus mesh file and returns the mesh data as a DataFrame.

# Arguments
- `params::Dict`: The parameters.
- `filename::String`: The path to the Exodus mesh file.
# Returns
- `mesh::DataFrame`: The mesh data as a DataFrame.
"""
function read_mesh(params::Dict, filename::String)
    exo = ExodusDatabase(filename, "r")

    coords = read_coordinates(exo)
    mesh_df = DataFrame(x = Float64[],
                        y = Float64[],
                        z = Float64[],
                        volume = Float64[],
                        block_id = Int64[])
    block_ids = read_ids(exo, Block)

    nodal_var_names = read_names(exo, NodalVariable)
    elem_var_names = read_names(exo, ElementVariable)
    num_nodal_var = length(nodal_var_names)
    num_elem_var = length(elem_var_names)
    if num_nodal_var > 0
        @info "Found $(num_nodal_var) nodal variables: $(nodal_var_names)"
    end
    if num_elem_var > 0
        @info "Found $(num_elem_var) element variables: $(elem_var_names)"
    end

    nodals = Dict()
    for nodal_var_name in nodal_var_names
        nodals[nodal_var_name] = read_values(exo, NodalVariable, 1, 1, nodal_var_name)
        nodals[nodal_var_name * "_sorted"] = []
    end

    for (iID, block_id) in enumerate(block_ids)
        block = read_block(exo, block_id)
        block_id_map = Exodus.read_block_connectivity(exo,
                                                      block_id,
                                                      block.num_nodes_per_elem *
                                                      block.num_elem)
        if block.elem_type == "TETRA"
            for i in 1:(block.num_elem)
                indices = (block.num_nodes_per_elem * (i - 1) + 1):(block.num_nodes_per_elem * i)
                node_ids = block_id_map[indices]
                vertices = coords[:, node_ids]
                center = sum(vertices, dims = 2) / size(vertices)[2]
                volume = tetrahedron_volume(vertices)
                push!(mesh_df,
                      (x = center[1],
                       y = center[2],
                       z = center[3],
                       volume = volume,
                       block_id = Int64(block_id)))
            end
        elseif block.elem_type == "HEX8"
            for i in 1:(block.num_elem)
                indices = (block.num_nodes_per_elem * (i - 1) + 1):(block.num_nodes_per_elem * i)
                node_ids = block_id_map[indices]
                vertices = coords[:, node_ids]
                center = sum(vertices, dims = 2) / size(vertices)[2]
                volume = hex8_volume(vertices)
                push!(mesh_df,
                      (x = center[1],
                       y = center[2],
                       z = center[3],
                       volume = volume,
                       block_id = Int64(block_id)))
            end
        elseif block.elem_type == "SPHERE"
            volume_nodal_name = nothing
            for name in nodal_var_names
                if lowercase(name) == "volume"
                    volume_nodal_name = name
                    break
                end
            end
            if volume_nodal_name === nothing
                @abort "Volume is missing. Please define a 'Volume' for each point in the mesh file."
            end

            for i in 1:(block.num_elem)
                node_ids = block_id_map[i]
                vertices = coords[:, node_ids]
                volume = nodals[volume_nodal_name][node_ids]
                push!(mesh_df,
                      (x = vertices[1],
                       y = vertices[2],
                       z = vertices[3],
                       volume = volume,
                       block_id = Int64(block_id)))
                for nodal_var_name in nodal_var_names
                    append!(nodals[nodal_var_name * "_sorted"],
                            nodals[nodal_var_name][node_ids])
                end
            end
        else
            @abort "Element type $(block.elem_type) not supported"
        end
    end
    for nodal_var_name in nodal_var_names
        mesh_df[!, nodal_var_name] = nodals[nodal_var_name * "_sorted"]
    end

    close(exo)
    coords = nothing
    block_ids = nothing

    return mesh_df
end

end # module Exodus_Mesh
