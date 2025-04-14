# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Contact_search
include("../../Support/Helpers.jl")
using .Helpers: get_nearest_neighbors, nearest_point_id, get_block_nodes, compute_geometry,
                point_is_inside, compute_surface_nodes_and_connections,
                get_surface_information
using CDDLib, Polyhedra
export init_contact_search
export compute_contact_pairs

## for one contact pair

function init_contact_search(datamanager, contact_params, cm)

    # global list
    @info "Contact search: Create global position array"
    # local search field ids
    all_positions = datamanager.get_all_positions()
    # find all outer surfaces.
    master_id = contact_params["Master"]
    # Polyhedra for the master block; can be 2D or 3D

    surface_ids = datamanager.get_contact_block_ids(master_id)
    free_surfaces = datamanager.get_free_contact_surfaces(master_id)
    poly_master = compute_geometry(all_positions[surface_ids, :])
    connection_master = compute_surface_nodes_and_connections(all_positions[surface_ids, :],
                                                              poly_master,
                                                              free_surfaces)
    slave_id = contact_params["Slave"]
    surface_ids = datamanager.get_contact_block_ids(slave_id)
    free_surfaces = datamanager.get_free_contact_surfaces(slave_id)
    # Polyhedra for the slave block; can be 2D or 3D
    poly_slave = compute_geometry(all_positions[surface_ids, :])
    connection_slave = compute_surface_nodes_and_connections(all_positions[surface_ids, :],
                                                             poly_slave,
                                                             free_surfaces)
    # global ids -> muss gefiltered werden hier -> global all nodesID -> local contact ID ableiten
    # ich brauche  -> reduced global all nodesID und local nodesID

    # checken
    datamanager.set_free_surface_connections(Dict(master_id => connection_master))
    datamanager.set_free_surface_connections(Dict(slave_id => connection_slave))
    datamanager.set_free_surface_nodes(Dict(master_id => sort(collect(keys(connection_master)))))
    datamanager.set_free_surface_nodes(Dict(slave_id => sort(collect(keys(connection_slave)))))

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

function filter_all_positions(datamanager)
    master = data
end

function create_local_exchange_ids(filtered_all_positions)
    datamanager.get_local_nodes(filtered_all_positions)
end

function compute_contact_pairs(datamanager, contact_params)
    all_positions = datamanager.get_all_positions()
    #-------------

    near_points = find_potential_contact_pairs(datamanager, contact_params)
    println(near_points)
    #connectivity = datamanager.get_contact_connectivity(contact_params["Slave"])
    #backend = CDDLib.Library()
    #points = Vector{Vector{Float64}}(eachrow(slave_nodes))
    #poly = polyhedron(vrep(points), backend)
    #normals = get_surface_normals(points)
    #for (master_node, slave_nodes) in pairs(near_points)
    #    if !(master in poly)
    #        continue
    #    end
    #    id = nearest_point_id(all_positions[master_node, :], all_positions[slave_nodes, :])
    #    contact_dict[master_node] = id
    #    contact_normal[master_node] = normals[connectivity[id]]
    #    # shortest
    #end
    datamanager.set_contact_pairs_and_normals(contact_dict, contact_normal)    #connectivity -> glob to local
    # distance
end

"""
    find_potential_contact_pairs(datamanager::Module, contact_params::Dict)

Finds a list of potential master slave pairs which are next to each other. Only the free surface nodes of the contact blocks are tested. The process is done equally at each computational core.

# Arguments
- `datamanager::module`: datamanager.
- `contact_params::Dict`: dictionary with contact relevant information.

# Returns
- pairs of potential contact partner.
"""

function find_potential_contact_pairs(datamanager::Module, contact_params::Dict)
    all_positions = datamanager.get_all_positions()

    master_nodes = datamanager.get_free_surface_nodes(contact_params["Master"])
    slave_nodes = datamanager.get_free_surface_nodes(contact_params["Slave"])

    nmaster = length(master_nodes)
    near_points = fill(Vector{Int64}([]), length(search_radius))

    return get_nearest_neighbors(1:nmaster,
                                 dof::Int64,
                                 all_positions[master_nodes, :],
                                 all_positions[slave_nodes, :],
                                 contact_params["Search Radius"],
                                 near_points,
                                 true)
end

# convex_hull gibt mir die außenpunkte
# lokale convex hulls bilden, um master slave abzufangen -> inside?
# überlappung?
# alternativ bondbased -> abstand der punkte die drinne liegen zu ihren Partnern;

# für kleine Deformationen?
#E₁ = Ellipsoid(zeros(2), [1 0; 0 2.])
## den ganzen block als huelle?
## dann testen
## dann den nächsten Punkt finden; Radius ist horizont
## dann funktion aus den ueberlappenden horizonten machen -> bool aus hull

function get_surface_normals(points::Union{Vector{Vector{Float64}},Vector{Vector{Int64}}})
    return MixedMatHRep(hrep(polyhedron(vrep(points), CDDLib.Library())))
end

function get_surface_normals(poly)
    return MixedMatHRep(hrep(poly)).A
end

function synch_all_positions(datamanager)
    all_positions = datamanager.get_all_positions()
    deformed_coor = datamanager.get_field("Deformed Coordinates", "NP1")
    mapping = datamanager.local_to_global_contact()
    for iID in 1:datamanager.get_nnodes()
        all_positions[mapping[iID], :] .= deformed_coor[iID, :]
    end
    datamanager.set_all_positions(all_positions)
end

end
