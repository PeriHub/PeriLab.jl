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
    return datamanager
end

function compute_contact_pairs(datamanager::Module, cg::String, contact_params::Dict,
                               max_cont::Int64)
    all_positions = datamanager.get_all_positions()
    #-------------
    near_points, list_empty = find_potential_contact_pairs(datamanager, contact_params)
    if list_empty
        return
    end
    potential_contact_dict = create_potential_contact_dict(near_points, datamanager,
                                                           contact_params)
    contact_dict = datamanager.get_contact_dict(cg)

    # node ids to create the geometry for checking if master ids are inside
    slave_block_nodes = datamanager.get_contact_block_ids(contact_params["Slave Block ID"])

    poly = compute_geometry(all_positions[slave_block_nodes, :])
    contact_nodes = datamanager.get_field("Contact Nodes")
    mapping = datamanager.get_exchange_id_to_local_id()
    contact_nodes .= 1

    for (master_node, near_ids) in pairs(potential_contact_dict)
        try
            if !(all_positions[master_node, :] in poly)
                continue
            end
        catch
            # if poly is destroyed, due to damages
            continue
        end
        dof = datamanager.get_dof()
        nslave::Int64 = 0
        contact_dict[master_node] = Dict("Slaves" => zeros(Int64, max_cont),
                                         "Normals" => zeros(Float64, max_cont, dof),
                                         "Distances" => zeros(Float64, max_cont),
                                         "nSlaves" => 0)
        for id in near_ids
            if contact_nodes[mapping[id]] == 1
                #contact_nodes[mapping[id]] = 2
            end
            if contact_nodes[mapping[master_node]] == 1
                #contact_nodes[mapping[master_node]] = 3
            end
            distance,
            normal = compute_distance_and_normals(all_positions[master_node, :],
                                                  all_positions[id, :])
            # TODO preallocation, e.g. max number of contact nodes -> size depended
            # @info "$(contact_params["Contact Radius"]), $distance, $(all_positions[master_node, :]),  $(all_positions[id, :]))"
            if contact_params["Contact Radius"] < abs(distance)
                continue
            end
            nslave += 1
            if nslave > max_cont
                @warn "Maximum number of potential contact pairs $max_cont is reached. The rest of the nearest neighbors are skipped."
                break
            end
            # for debugging
            if id in keys(mapping)
                contact_nodes[mapping[id]] = 4
            end
            if master_node in keys(mapping)
                contact_nodes[mapping[master_node]] = 5
            end

            contact_dict[master_node]["nSlaves"] = nslave
            contact_dict[master_node]["Slaves"][nslave] = id
            contact_dict[master_node]["Normals"][nslave, :] = normal
            contact_dict[master_node]["Distances"][nslave] = distance
        end
    end
    datamanager.set_contact_dict(cg, contact_dict)
end

function create_potential_contact_dict(near_points, datamanager, contact_params)
    # exchange vector ids of the free contact blocks
    # computed in function compute_and_set_free_contact_nodes
    master_nodes = datamanager.get_free_contact_nodes(contact_params["Master Block ID"])
    slave_nodes = datamanager.get_free_contact_nodes(contact_params["Slave Block ID"])

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
    master_nodes = datamanager.get_free_contact_nodes(contact_params["Master Block ID"])
    slave_nodes = datamanager.get_free_contact_nodes(contact_params["Slave Block ID"])

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
