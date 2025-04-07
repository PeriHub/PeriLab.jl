# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

export get_all_blocks
export get_all_positions
export get_contact_relevant_global_ids
export get_free_surface_connections
export get_free_surface_nodes
export set_all_blocks
export set_all_positions
export set_contact_relevant_global_ids
export set_free_surface_connections
export set_free_surface_nodes
function get_all_blocks()
    return data["All Blocks"]
end

function get_all_positions()
    return data["All positions"]
end

function get_free_surfaces()
    return data["Free Surfaces"]
end

function get_contact_relevant_global_ids()
    return data["Contact IDs"]
end

function get_free_surface_nodes()
    return data["Free Surface Nodes"]
end

function get_free_surface_connections()
    return data["Free Surface Nodes Connections"]
end

function set_all_blocks(all_blocks)
    data["All Blocks"] = all_blocks
end

function set_all_positions(all_coordinate)
    data["All positions"] = all_coordinate
end

function set_contact_relevant_global_ids(ids::Vector{Int64})
    data["Contact IDs"] = ids
end

function set_free_surfaces(free_surfaces)
    data["Free Surfaces"] = free_surfaces
end

function set_free_surface_nodes(free_surface_nodes)
    merge(data["Free Surface Nodes"], free_surface_nodes)
end

function set_free_surface_connections(free_surface_nodes_connections)
    merge(data["Free Surface Nodes Connections"], free_surface_nodes_connections)
end
