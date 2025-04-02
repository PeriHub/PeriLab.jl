# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Contact_search
include("../../Support/Helpers.jl")
using .Helpers: get_nearest_neighbors, nearest_point_id, get_block_nodes, compute_geometry,
                point_is_inside, compute_surface_nodes, get_surface_information
using CDDLib, Polyhedra
export init_contact_search
export compute_contact_pairs

## for one contact pair

function init_contact_search(datamanager, contact_params, cm, surface_nodes)

    # global list
    @info "Contact search: Create global position array"
    all_positions = datamanager.get_all_positions()
    block_ids = datamanager.get_all_blocks()

    contact_block_nodes = get_block_nodes(block_ids, length(block_ids))
    # Polyhedra for the master block; can be 2D or 3D
    poly_master = compute_geometry(all_positions[contact_block_nodes[contact_params["Master"]],
                                                 :])
    # Polyhedra for the slave block; can be 2D or 3D
    poly_slave = compute_geometry(all_positions[contact_block_nodes[contact_params["Slave"]],
                                                :])
    # filter

    contact_block_nodes[contact_params["Master"]] = compute_surface_nodes(all_positions[contact_block_nodes[contact_params["Master"]],
                                                                                        :],
                                                                          poly_master)
    contact_block_nodes[contact_params["Slave"]] = compute_surface_nodes(all_positions[blocontact_block_nodesck_nodes[contact_params["Slave"]],
                                                                                       :],
                                                                         poly_slave)
    # global ids
    surface_nodes[cm] = unique(vcat(contact_block_nodes[contact_params["Master"]],
                                    contact_block_nodes[contact_params["Slave"]]))
    compute_contact_points_surface_connectivity(block_nodes[contact_params["Master"]],
                                                datamanager,
                                                contact_params["Master"], all_positions,
                                                poly_master)
    compute_contact_points_surface_connectivity(block_nodes[contact_params["Slave"]],
                                                datamanager,
                                                contact_params["Slave"], all_positions,
                                                poly_master)
    return datamanager
end

function compute_contact_points_surface_connectivity(datamanager::Module,
                                                     block_id::Int64,
                                                     block_nodes::Dict{Int64,Vector{Int64}},
                                                     all_positions,
                                                     poly)
    nlist = datamanager.get_nlist()
    @views connectivity = get_surface_connectivity(all_positions, block_nodes, nlist, poly)
    datamanager.set_contact_connectivity(block_id, connectivity)
    datamanager.set_contact_nodes(block_id, collect(keys(connectivity)))
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
function get_surface_connectivity(points::Union{Vector{Vector{Float64}},
                                                Vector{Vector{Int64}}},
                                  surface_nodes::Dict{Int64,Int64},
                                  nlist,
                                  poly)
    # H-Rep abrufen (Facetten als a⋅x ≤ b)

    normals, offsets = get_surface_information(poly)
    connections = Dict{Int64,Int64}()
    msg = true
    for (iID, point) in surface_nodes
        for surface_id in eachindex(offsets)
            # check to which surface the point is connected && check if all neighbors are inside the block.
            if isapprox(dot(point, normals[surface_id, :]) - offsets[surface_id], 0;
                        atol = 1e-6) &&
               !check_neighbor_position(poly, points, nlist[iID], msg)
                break
            end
            connections[iID] = surface_id
        end
    end
    return connections
end

function filter_all_positions(datamanager)
    master = data
end

function create_local_contact_nodes()
    #specifies the ids at core from the global contact nodes
    datamanager.get_local_nodes(global_nodes)
    datamanager.set_local_contact_nodes("Master")

    datamanager.set_local_contact_nodes("Slave")
end

function create_local_exchange_ids(filtered_all_positions)
    datamanager.get_local_nodes(filtered_all_positions)
end

function compute_contact_pairs(datamanager, contact_params)
    #synch_all_positions(datamanager)
    # Single core solution
    deformed_coor = datamanager.get_field("Deformed Coordinates", "NP1")

    for iID in 1:datamanager.get_nnodes()
        all_positions .= deformed_coor
    end
    datamanager.set_all_positions(all_positions)
    #-------------

    master_nodes = datamanager.get_contact_nodes(contact_params["Master"])
    slave_nodes = datamanager.get_contact_nodes(contact_params["Slave"])
    near_points = find_potential_contact_pairs(datamanager, contact_params)
    connectivity = datamanager.get_contact_connectivity(contact_params["Slave"])
    backend = CDDLib.Library()
    points = Vector{Vector{Float64}}(eachrow(slave_nodes))
    poly = polyhedron(vrep(points), backend)
    normals = get_surface_normals(points)
    for (master_node, slave_nodes) in pairs(near_points)
        if !(master in poly)
            continue
        end
        id = nearest_point_id(all_positions[master_node, :], all_positions[slave_nodes, :])
        contact_dict[master_node] = id
        contact_normal[master_node] = normals[connectivity[id]]
        # shortest
    end
    datamanager.set_contact_pairs_and_normals(contact_dict, contact_normal)    #connectivity -> glob to local
    # distance
end

function find_potential_contact_pairs(datamanager, contact_params)
    # global IDs
    all_positions = datamanager.get_all_positions()
    # optimization of number of positions. Reduction to surface nodes.

    master_nodes = datamanager.get_contact_nodes(contact_params["Master"])
    slave_nodes = datamanager.get_contact_nodes(contact_params["Slave"])
    #global_master_nodes = datamanager.loc_to_glob(local_master_nodes)
    # global_slave_nodes = datamanager.loc_to_glob(local_slave_nodes)
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
