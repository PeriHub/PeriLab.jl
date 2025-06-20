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
    datamanger.set_search_step(cm, 0)
    return datamanager
end

function global_contact_search(datamanager, contact_params)
    all_positions = datamanager.get_all_positions()
    #-------------
    dof = datamanager.get_dof()
    # ids are exchange vector ids (''all position'' ids)
    master_nodes = datamanager.get_free_contact_nodes(contact_params["Master Block ID"])
    slave_nodes = datamanager.get_free_contact_nodes(contact_params["Slave Block ID"])

    near_points,
    no_pairs = find_potential_contact_pairs(dof,
                                            all_positions[master_nodes, :],
                                            all_positions[slave_nodes, :],
                                            contact_params["Search Radius"])
    if no_pairs
        return no_pairs, nothing, nothing
    end

    return no_pairs, filter_global_contact_nodes(near_points, master_nodes, slave_nodes)
end

function compute_contact_pairs(datamanager::Module, cg::String, contact_params::Dict,
                               max_cont::Int64)
    n = datamanger.get_search_step()
    if datamanger.get_search_step() == 0
        synch
        no_pairs, global_master,
        global_slave = global_contact_search(datamanager,
                                             contact_params)
        datamanager.set_global_search_master_nodes(global_master)
        datamanager.set_global_search_slave_nodes(global_slave)
        datamanager.set_no_pairs_flag(cg, no_pairs)
        #
        if no_pairs # must be stored
            return
        end
        #local_contact_search(datamanager, contact_params,global_master,global_slave)
    end
    if datamanager.get_no_pairs_flag(cg) # must be stored
        return
    end
    #synch
    # local_contact_search(datamanager, contact_params,datamanager.get_global_search_master_nodes(global_master),datamanager.set_global_search_slave_nodes(global_slave))
end

function local_contact_search(datamanager, contact_params, master_nodes, slave_nodes)
    all_positions = datamanager.get_all_positions()
    #-------------
    dof = datamanager.get_dof()
    # ids are exchange vector ids (''all position'' ids)
    contact_points,
    no_pairs = find_potential_contact_pairs(dof,
                                            all_positions[master_nodes, :],
                                            all_positions[slave_nodes, :],
                                            contact_params["Contact Radius"])
    if no_pairs
        return
    end
    for (iID, neighbors) in enumerate(contact_points)
        if length(neighbors) == 0
            continue
        end
        master_id = master_nodes[iID]
        if master_id in keys(mapping)
            contact_nodes[mapping[master_id]] = 5
        end
        for jID in neighbors
            slave_id = slave_nodes[jID]
            distance,
            normal = compute_distance_and_normals(all_positions[master_id, :],
                                                  all_positions[slave_id, :])

            # for debugging
            if slave_id in keys(mapping)
                contact_nodes[mapping[slave_id]] = 4
            end

            contact_dict[master_node]["nSlaves"] = nslave
            contact_dict[master_node]["Slaves"][nslave] = id
            contact_dict[master_node]["Normals"][nslave, :] = normal
            contact_dict[master_node]["Distances"][nslave] = distance
        end
    end
    return
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
function find_potential_contact_pairs(dof, points_1, points_2, search_radius)
    nmaster = length(points_1)
    near_points = fill(Vector{Int64}([]), nmaster)

    return get_nearest_neighbors(1:nmaster,
                                 dof,
                                 points_1,
                                 points_2,
                                 search_radius,
                                 near_points,
                                 true)
end

function local_contact_search()
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
        dof = datamanager.get_dof()
        nslave::Int64 = 0
        contact_dict[master_node] = Dict("Slaves" => zeros(Int64, max_cont),
                                         "Normals" => zeros(Float64, max_cont, dof),
                                         "Distances" => zeros(Float64, max_cont),
                                         "nSlaves" => 0)
        for id in near_ids
            if contact_nodes[mapping[id]] == 1
                contact_nodes[mapping[id]] = 2
            end
            if contact_nodes[mapping[master_node]] == 1
                contact_nodes[mapping[master_node]] = 3
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
"""
create_potential_contact_dict(near_points, master_nodes, slave_nodes)

 Clear the near_points from the empty lists.

"""
function filter_global_contact_nodes(near_points, master_nodes, slave_nodes)
    contact_global_masters = []
    contact_global_slaves = []

    for (pID, neighbors) in enumerate(near_points)
        if length(neighbors) > 0
            append!(contact_global_masters, master_nodes[pID])
            append!(contact_global_slaves, slave_nodes[neighbors])
        end
    end
    return unique(contact_global_masters), unique(contact_global_slaves)
end

end
