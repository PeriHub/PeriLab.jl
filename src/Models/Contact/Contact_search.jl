# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Contact_search
include("../../Support/Helpers.jl")
using .Helpers: get_nearest_neighbors, nearest_point_id, get_block_nodes, compute_geometry,
                point_is_inside, get_surface_information, compute_distance_and_normals

export init_contact_search
export compute_contact_pairs

## for one contact pair

function init_contact_search(datamanager, contact_params, cm)
    datamanager.set_search_step(cm, 0)
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
        return no_pairs, Int64[], Int64[]
    end

    global_master,
    global_slave = filter_global_contact_nodes(near_points, master_nodes, slave_nodes)

    return no_pairs, global_master, global_slave
end

function compute_contact_pairs(datamanager::Module, cg::String, contact_params::Dict)
    if datamanager.get_search_step(cg) == 0
        datamanager.set_global_search_master_nodes(cg, Vector{Int64}([]))
        datamanager.set_global_search_slave_nodes(cg, Vector{Int64}([]))
        no_pairs, global_master,
        global_slave = global_contact_search(datamanager,
                                             contact_params)

        datamanager.set_global_search_master_nodes(cg, global_master)
        datamanager.set_global_search_slave_nodes(cg, global_slave)
        datamanager.set_no_pairs_flag(cg, no_pairs)

        datamanager.add_synchronization_list(global_master)
        datamanager.add_synchronization_list(global_slave)
    end
    if datamanager.get_no_pairs_flag(cg) # must be stored
        return
    end
    mapping = datamanager.get_exchange_id_to_local_id()

    contact_dict = local_contact_search(datamanager, contact_params,
                                        datamanager.get_global_search_master_nodes(cg),
                                        datamanager.get_global_search_slave_nodes(cg))
    datamanager.set_contact_dict(cg, contact_dict)
end

function local_contact_search(datamanager, contact_params, master_nodes, slave_nodes)
    all_positions = datamanager.get_all_positions()
    #-------------
    dof = datamanager.get_dof()
    contact_dict = Dict{Int64,Dict{String,Any}}()
    # ids are exchange vector ids (''all position'' ids)
    contact_points,
    no_pairs = find_potential_contact_pairs(dof,
                                            all_positions[master_nodes, :],
                                            all_positions[slave_nodes, :],
                                            contact_params["Contact Radius"])
    if no_pairs
        return contact_dict
    end

    mapping = datamanager.get_exchange_id_to_local_id()
    contact_nodes = datamanager.get_field("Contact Nodes")
    for (iID, neighbors) in enumerate(contact_points)
        if isempty(neighbors)
            continue
        end
        master_id = master_nodes[iID]
        if !isnothing(get(mapping, master_id, nothing))
            contact_nodes[mapping[master_id]] = 5
        end
        nslave::Int64 = length(neighbors)

        contact_dict[master_id] = Dict{String,Any}("nSlaves" => nslave,
                                                   "Slaves" => zeros(Int64, nslave),
                                                   "Normals" => zeros(Float64, nslave, dof),
                                                   "Distances" => zeros(Float64, nslave))
        @views for (jID, slave_id) in enumerate(slave_nodes[neighbors])
            distance,
            normal = compute_distance_and_normals(all_positions[master_id, :],
                                                  all_positions[slave_id, :])

            # for debugging
            if !isnothing(get(mapping, slave_id, nothing))
                contact_nodes[mapping[slave_id]] = 8
            end
            contact_dict[master_id]["nSlaves"] = nslave
            contact_dict[master_id]["Slaves"][jID] = slave_id
            contact_dict[master_id]["Normals"][jID, :] = normal
            contact_dict[master_id]["Distances"][jID] = distance
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
function find_potential_contact_pairs(dof, points_1, points_2, search_radius)
    nmaster = size(points_1)[1]
    near_points = fill(Vector{Int64}([]), nmaster)

    return get_nearest_neighbors(1:nmaster,
                                 dof,
                                 points_1,
                                 points_2,
                                 search_radius,
                                 near_points,
                                 true)
end

"""
filter_global_contact_nodes(near_points, master_nodes, slave_nodes)

 Clear the near_points from the empty lists and store the unique ids from the master and slave node lists.
 These nodes are used for the local search algorithm, because these are the ones in the search radius.

 return unique(contact_global_masters), unique(contact_global_slaves)
"""
function filter_global_contact_nodes(near_points, master_nodes, slave_nodes)
    contact_global_masters = Int64[]
    contact_global_slaves = Int64[]

    for (pID, neighbors) in enumerate(near_points)
        if !isempty(neighbors)
            push!(contact_global_masters, master_nodes[pID])
            append!(contact_global_slaves, slave_nodes[neighbors])
        end
    end
    return contact_global_masters, sort(unique(contact_global_slaves))
end

end
