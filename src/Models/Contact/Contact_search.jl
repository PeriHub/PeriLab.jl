# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Contact_search
include("../../Support/Helpers.jl")
using .Helpers: get_nearest_neighbors, nearest_point_id, get_block_nodes, compute_geometry,
                point_is_inside,
                get_surface_information, compute_distance_and_normals

export init_contact_search
export compute_contact_pairs

## for one contact pair

function init_contact_search(datamanager, contact_params, cm)

    # global list
    # do something with block nodes

    return datamanager
end

function compute_contact_pairs(datamanager::Module, cm::String, contact_params::Dict)
    all_positions = datamanager.get_all_positions()
    #-------------
    near_points, list_empty = find_potential_contact_pairs(datamanager, contact_params)
    if list_empty
        return
    end
    potential_contact_dict = create_potential_contact_dict(near_points, datamanager,
                                                           contact_params)
    contact_dict = datamanager.get_contact_dict(cm)

    # node ids to create the geometry for checking if master ids are inside
    slave_block_nodes = datamanager.get_contact_block_ids(contact_params["Slave Block ID"])

    poly = compute_geometry(all_positions[slave_block_nodes, :])
    contact_nodes = datamanager.get_field("Contact Nodes")
    mapping = datamanager.get_exchange_id_to_local_id()

    for (master_node, near_ids) in pairs(potential_contact_dict)
        try
            if !(all_positions[master_node, :] in poly)
                continue
            end
        catch
            # if poly is destroyed, due to damages
            continue
        end

        contact_dict[master_node] = Dict("Slaves" => [], "Normals" => [],
                                         "Distances" => [])
        for id in near_ids
            contact_nodes[mapping[id]] = 2
            contact_nodes[mapping[master_node]] = 3
            distance,
            normal = compute_distance_and_normals(all_positions[master_node, :],
                                                  all_positions[id, :])
            # TODO preallocation, e.g. max number of contact nodes -> size depended
            if contact_params["Search Radius"] < abs(distance)
                continue
            end
            contact_nodes[mapping[id]] = 4
            contact_nodes[mapping[master_node]] = 5
            append!(contact_dict[master_node]["Slaves"], id)
            append!(contact_dict[master_node]["Normals"], [normal])
            append!(contact_dict[master_node]["Distances"], distance)
        end
    end
    datamanager.set_contact_dict(cm, contact_dict)
end

function create_potential_contact_dict(near_points, datamanager, contact_params)
    # exchange vector ids of the free contact blocks
    # computed in function compute_and_set_free_surface_nodes
    master_nodes = datamanager.get_free_surface_nodes(contact_params["Master Block ID"])
    slave_nodes = datamanager.get_free_surface_nodes(contact_params["Slave Block ID"])

    contact_dict = Dict{Int64,Vector{Int64}}()
    for (pID, neighbors) in enumerate(near_points)
        if length(neighbors) > 0
            contact_dict[master_nodes[pID]] = slave_nodes[neighbors]
        end
    end
    return contact_dict
end

"""
    find_potential_contact_pairs(datamanager::Module, contact_params::Dict)

Finds a list of potential master slave pairs which are next to each other. Only the free surface nodes of the contact blocks are tested. The process is done equally at each computational core.

# Arguments
- `datamanager::module`: datamanager.
- `contact_params::Dict`: dictionary with contact relevant information.

# Returns
- pairs of potential contact partner in exchange vector ids.
"""## TODO test
function find_potential_contact_pairs(datamanager::Module, contact_params::Dict)
    all_positions = datamanager.get_all_positions()
    dof = datamanager.get_dof()
    # ids are exchange vector ids (''all position'' ids)
    master_nodes = datamanager.get_free_surface_nodes(contact_params["Master Block ID"])
    slave_nodes = datamanager.get_free_surface_nodes(contact_params["Slave Block ID"])

    nmaster = length(master_nodes)
    near_points = fill(Vector{Int64}([]), nmaster)
    return get_nearest_neighbors(1:nmaster,
                                 dof,
                                 all_positions[master_nodes, :],
                                 all_positions[slave_nodes, :],
                                 contact_params["Search Radius"],
                                 near_points,
                                 true)
end

end
