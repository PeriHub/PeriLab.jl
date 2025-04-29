# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Contact_search
include("../../Support/Helpers.jl")
using .Helpers: get_nearest_neighbors, nearest_point_id, get_block_nodes, compute_geometry,
                point_is_inside, compute_surface_nodes_and_connections,
                get_surface_information, compute_distance_to_surfaces
using CDDLib, Polyhedra
export init_contact_search
export compute_contact_pairs

## for one contact pair

function init_contact_search(datamanager, contact_params, cm)

    # global list
    @info "Contact search: Create global position arrays"

    # local search field ids
    all_positions = datamanager.get_all_positions()
    # find all outer surfaces.
    master_id = contact_params["Master"]
    slave_id = contact_params["Slave"]
    # Polyhedra for the master block; can be 2D or 3D
    @info "Contact pair Master block $master_id - Slave block $slave_id"
    surface_ids_master = datamanager.get_contact_block_ids(master_id)
    free_surfaces = datamanager.get_free_contact_surfaces(master_id)
    poly_master = compute_geometry(all_positions[surface_ids_master, :])
    connection_master = compute_surface_nodes_and_connections(all_positions[surface_ids_master,
                                                                            :],
                                                              poly_master,
                                                              free_surfaces)

    surface_ids_slave = datamanager.get_contact_block_ids(slave_id)
    free_surfaces = datamanager.get_free_contact_surfaces(slave_id)
    # Polyhedra for the slave block; can be 2D or 3D
    poly_slave = compute_geometry(all_positions[surface_ids_slave, :])
    connection_slave = compute_surface_nodes_and_connections(all_positions[surface_ids_slave,
                                                                           :],
                                                             poly_slave,
                                                             free_surfaces)
    # it is merged and therefore works here
    datamanager.set_free_surface_connections(connection_master)
    datamanager.set_free_surface_connections(connection_slave)
    datamanager.set_free_surface_nodes(master_id,
                                       surface_ids_master[sort(collect(keys(connection_master)))])
    datamanager.set_free_surface_nodes(slave_id,
                                       surface_ids_slave[sort(collect(keys(connection_slave)))])
    # do something with block nodes

    return datamanager
end

"""
    which_surface()

checks which surface has to be used by checking the number of nodes in the nearest point list. The highest number at one surface, defines the surface.

"""

function which_surface()
end

"""
    function get_surface_connectivity(points::Union{Vector{Vector{Float64}}, Vector{Vector{Int64}}}, surface_nodes,  nlist, poly)

Identifies the surfaces which are not connected with other surfaces. These nodes are used for the near neighbor search. It is checked if the neighbors are outside the block. This leads to the issue that at boundary edges points lying at the surface are not considered.

!!! info "limitations"
    This might have limitations if non plane surfaces are used.


# Arguments
- `points::Union{Vector{Vector{Float64}}, Vector{Vector{Int64}}}`: All point coordinates.
- `surface_nodes::Dict{Int64, Int64}`: dictionary mapping contact block information to the surface nodes of that block.
- `nlist::Vector{Vector{Int64}}`: Neighborhoodlist
- `poly`: Polyhedra object 2D or 3D
# Returns
- `connections::Dict{Int64,Vector{Int64}}`: connections of the surface node with outer surface
"""

function filter_surface_connectivity(points::Union{Matrix{Float64},Matrix{Int64}}, nlist,
                                     connections, poly)
    msg = true
    for iID in eachindex(points[:, 1])
        if !check_neighbor_position(poly, points, nlist[iID], msg)
            delete!(connections, iID)
        end
    end
    return connections
end

function compute_contact_pairs(datamanager::Module, cm::String, contact_params::Dict)
    all_positions = datamanager.get_all_positions()
    #-------------
    near_points = find_potential_contact_pairs(datamanager, contact_params)

    potential_contact_dict = create_potential_contact_dict(near_points, datamanager,
                                                           contact_params)
    contact_dict = datamanager.get_contact_dict(cm)
    connectivity = datamanager.get_free_surface_connections()
    slave_block_nodes = datamanager.get_contact_block_ids(contact_params["Slave"])
    poly = polyhedron(vrep(all_positions[slave_block_nodes, :]), CDDLib.Library())
    normals, offsets = get_surface_information(poly)
    println(normals)
    @error ""
    for (master_node, near_ids) in pairs(potential_contact_dict)
        # if !(all_positions[master_node, :] in poly) # test if master point lies in slave block
        #     continue
        # end
        id = near_ids[nearest_point_id(all_positions[master_node, :],
                                       all_positions[near_ids, :])]

        distances = compute_distance_to_surfaces(all_positions[master_node, :], normals,
                                                 offsets,
                                                 connectivity[id])
        # TODO preallocation, e.g. max number of contact nodes -> size depended
        append!(contact_dict["Pairs: Master-Slave"], [(master_node, id)])
        append!(contact_dict["Normals"], [normals[connectivity[id], :]])
        append!(contact_dict["Distances"], [distances])
    end
    datamanager.set_contact_dict(cm, contact_dict)
end

function create_potential_contact_dict(near_points, datamanager, contact_params)
    # exchange vector ids
    master_nodes = datamanager.get_free_surface_nodes(contact_params["Master"])
    slave_nodes = datamanager.get_free_surface_nodes(contact_params["Slave"])

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
"""

## TODO test
function find_potential_contact_pairs(datamanager::Module, contact_params::Dict)
    all_positions = datamanager.get_all_positions()
    dof = datamanager.get_dof()
    # ids are exchange vector ids (''all position'' ids)
    master_nodes = datamanager.get_free_surface_nodes(contact_params["Master"])
    slave_nodes = datamanager.get_free_surface_nodes(contact_params["Slave"])

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

function get_surface_normals(points::Union{Vector{Vector{Float64}},Vector{Vector{Int64}}})
    return MixedMatHRep(hrep(polyhedron(vrep(points), CDDLib.Library())))
end

function get_surface_normals(poly)
    return MixedMatHRep(hrep(poly)).A
end
end
