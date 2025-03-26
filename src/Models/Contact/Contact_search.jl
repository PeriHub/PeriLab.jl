# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Contact_search
include("../../Support/Helpers.jl")
using .Helpers: get_nearest_neighbors, nearest_point_id
using LazySets: convex_hull

export init_contact
export compute_contact_pairs

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
    backend = CDDLib.Library()
    poly = polyhedron(vrep(points), backend)
    MixedMatHRep(hr)
    hrep_form = MixedMatHRep(hrep(poly))
    return hrep_form.A
end

function init_contact(datamanager, contact_params, block_nodes)
    ## block_nodes for all nodes!!!
    if !haskey(contact_params, "Master")
        @error "Contact model needs a ''Master''"
        return nothing
    end
    if !haskey(contact_params, "Slave")
        @error "Contact model needs a ''Slave''"
        return nothing
    end
    if contact_params["Master"] == contact_params["Slave"]
        @error "Contact master and slave are equal. Self contact is not implemented yet."
        return nothing
    end
    block_id = datamanager.get_field("Block_Id")
    if !(contact_params["Master"] in block_id)
        @error "Block defintion in master does not exist."
        return nothing
    end
    if !(contact_params["Slave"] in block_id)
        @error "Block defintion in slave does not exist."
        return nothing
    end
    if !haskey(contact_params, "Search Radius")
        @error "Contact model needs a ''Search Radius''"
        return nothing
    end
    # global list
    all_positions = datamanager.get_all_positions()
    datamanager = datamanager.get_all_blocks()
    define_contact_points_and_connectivity(datamanager, contact_params["Master"],
                                           block_nodes)
    define_contact_points_and_connectivity(datamanager, contact_params["Slave"],
                                           block_nodes)

    datamanager.set_all_positions()
    datamanager.blocks()
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

function define_contact_points_and_connectivity(datamanager::Module, block_id::Int64,
                                                block_nodes::Dict{Int64,Vector{Int64}})
    nlist = datamanager.get_nlist()
    points = Vector{Vector{Float64}}(eachrow(block_nodes[block_id]))
    @views connectivity = get_surface_connectivity(points, nlist)
    datamanager.set_contact_connectivity(block_id, connectivity)
    datamanager.set_contact_nodes(block_id, collect(keys(connectivity)))
end

"""
function get_surface_connectivity(points::Union{Vector{Vector{Float64}},Vector{Vector{Int64}}})

This function finds the point, surface cobination. However, no edge or corner nodes are included, because there are no unique normals. Only block outer surfaces will be found

!!! info "limitations"
    This might have limitations if non plane surfaces are used.

"""
function get_surface_connectivity(points::Union{Vector{Vector{Float64}},
                                                Vector{Vector{Int64}}}, nlist)
    backend = CDDLib.Library()
    poly = polyhedron(vrep(points), backend)
    # H-Rep abrufen (Facetten als a⋅x ≤ b)
    hrep_form = MixedMatHRep(hrep(poly))
    normals = hrep_form.A
    offsets = hrep_form.b
    connections = Dict{Int64,Int64}()
    msg = true
    for (iID, point) in enumerate(points)
        for id in eachindex(normals[:, 1])
            if abs(dot(normals[id, :], point) - offsets[id]) < 1e-6
                if !check_neighbor_position(poly, points, nlist[iID], msg)
                    break
                end
                if haskey(connections, iID)
                    delete!(connections, iID)
                    break
                end
                connections[iID] = id
            end
        end
    end
    return connections
end

function check_neighbor_position(poly, points, nlist::Vector{Int64}, msg::Bool)
    for nID in nlist
        if !(points[nID] in poly)
            if msg
                msg = false
                @warn "Make sure that your contact block is large enough. If it is surrounded by non contact blocks some of the points near the edgdes are ignored."
            end
            return false
        end
    end
    return true
end
end
