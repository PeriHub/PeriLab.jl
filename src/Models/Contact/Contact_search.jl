# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Contact_search
include("../../Support/Helpers.jl")
using .Helpers: get_nearest_neighbors
using LazySets: convex_hull

function init_contact(datamanager, contact_params)
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
    coor = datamanager.get_field("Coordinate")
    @views connectivity = get_surface_connectivity(Vector{Vector{Float64}}(eachrow(coor)))
    datamanager.set_connectivity(connectivity)
    datamanager.set_slave_nodes(collect(keys(connectivity)))
end

function compute_find_contact(datamanager, contact_params)
    find_potential_contact_pairs
    check_if_inside

    #connectivity -> glob to local
    # distance
end

function find_potential_contact_pairs(datamanager, contact_params)
    # global IDs
    all_positions = datamanager.get_all_contact_positions()
    master_nodes = datamanager.get_contact_master()
    slave_nodes = datamanager.get_contact_master()
    nmaster = length(master_nodes)
    near_points = fill(Vector{Int64}([]), length(search_radius))

    near_points = get_nearest_neighbors(1:nmaster,
                                        dof::Int64,
                                        all_positions[master_nodes, :],
                                        all_positions[slave_nodes, :],
                                        contact_params["Search Radius"],
                                        near_points,
                                        true)
end

function check_if_inside(masters, slaves, contact_dict)
    #slaves => neigbhboorhoodlist
    points_vec = [Vector{Float64}(row) for row in eachrow(slaves)]
    poly = polyhedron(vrep(points_vec))
    for master in masters
        if master in poly
            contact_dict[master] = slave # als vector?
        end
    end
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

"""
function get_surface_connectivity(points::Union{Vector{Vector{Float64}},Vector{Vector{Int64}}})

This function finds the point, surface cobination. However, no edge or corner nodes are included, because there are no unique normals.

!!! info "limitations"
    This might have limitations if non plane surfaces are used.

"""
function get_surface_connectivity(points::Union{Vector{Vector{Float64}},
                                                Vector{Vector{Int64}}})
    backend = CDDLib.Library()
    poly = polyhedron(vrep(points), backend)
    # H-Rep abrufen (Facetten als a⋅x ≤ b)
    hrep_form = MixedMatHRep(hrep(poly))
    normals = hrep_form.A
    offsets = hrep_form.b
    connections = Dict{Int64,Int64}()
    for (pID, point) in enumerate(points)
        for id in eachindex(normals[:, 1])
            if abs(dot(normals[id, :], point) - offsets[id]) < 1e-6
                if haskey(connections, pID)
                    delete!(connections, pID)
                    break
                end
                connections[pID] = id
            end
        end
    end
    return connections
end

end
