# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Contact_Search

using LinearAlgebra: dot
using StaticArrays: SVector
using PointNeighbors: GridNeighborhoodSearch, initialize_grid!, update!, foreach_neighbor
using ....Data_Manager
using ....Helpers: compute_distance_and_normals

export ContactPairs, ContactSearch, contact_pairs, compute_contact_pairs!

"""
    ContactPairs

The contact pairs of a contact group in the current step, as CSR over the masters: the
master `masters[k]` has the slaves `slaves[offsets[k]:(offsets[k + 1] - 1)]`; per pair the
`distances` and the unit `normals` (`dof` values per pair, from master to slave). All ids
are exchange ids (rows of `Data_Manager.get_all_positions()`). The vectors are reused in
every step; mutable, so a model gets it from the search state without a copy.
"""
mutable struct ContactPairs
    dof::Int64
    masters::Vector{Int64}
    offsets::Vector{Int64}
    slaves::Vector{Int64}
    distances::Vector{Float64}
    normals::Vector{Float64}
end

ContactPairs(dof::Int64) = ContactPairs(dof, Int64[], [1], Int64[], Float64[], Float64[])

function Base.empty!(p::ContactPairs)
    empty!(p.masters)
    resize!(p.offsets, 1)
    p.offsets[1] = 1
    empty!(p.slaves)
    empty!(p.distances)
    empty!(p.normals)
    return p
end

"""
    ContactSearch{D,S}

Search state of a contact group (`D` = dof, `S` = type of the grid). Every
`frequency` steps the global search collects the candidate pairs within `search_radius`
(at least `contact_radius`) from the free contact nodes of the master and the slave block
(a Verlet list, CSR as in `ContactPairs`). In every step only the candidate pairs are
checked against `contact_radius`. The grid of the slave nodes is kept and updated.
"""
mutable struct ContactSearch{D,S}
    contact_radius::Float64
    search_radius::Float64
    frequency::Int64
    # steps since the last global search; 0: the next step searches globally
    step::Int64
    master_nodes::Vector{Int64}
    slave_nodes::Vector{Int64}
    # coordinates of the master and slave nodes (dof x n), filled for the global search
    master_coords::Matrix{Float64}
    slave_coords::Matrix{Float64}
    grid::S
    grid_initialized::Bool
    candidate_masters::Vector{Int64}
    candidate_offsets::Vector{Int64}
    candidate_slaves::Vector{Int64}
    # sorted exchange ids of all candidate nodes, synchronized between the global searches
    candidate_nodes::Vector{Int64}
    pairs::ContactPairs
end

function ContactSearch(dof::Int64, master_nodes::Vector{Int64}, slave_nodes::Vector{Int64},
                       contact_radius::Float64, search_radius::Float64, frequency::Int64)
    # candidates within the search radius are a superset of the contact pairs only if the
    # search radius is not smaller than the contact radius
    radius = max(search_radius, contact_radius)
    # between two global searches only the candidates are checked: nodes that approach
    # each other in this time are missed without a margin
    if frequency > 1 && search_radius <= contact_radius
        @warn "Global Search Frequency $frequency > 1 with Search Radius $search_radius <= Contact Radius $contact_radius: contact pairs that come closer between two global searches are missed. Choose a Search Radius larger than the Contact Radius."
    end
    grid = GridNeighborhoodSearch{dof}(search_radius = radius,
                                       n_points = length(slave_nodes))
    return ContactSearch{dof,typeof(grid)}(contact_radius, radius, frequency, 0,
                                           master_nodes, slave_nodes,
                                           zeros(dof, length(master_nodes)),
                                           zeros(dof, length(slave_nodes)), grid, false,
                                           Int64[], [1], Int64[], Int64[],
                                           ContactPairs(dof))
end

"""
    contact_pairs(cg::String)

The `ContactPairs` of contact group `cg` in the current step.
"""
contact_pairs(cg::String) = (Data_Manager.get_contact_search(cg)::ContactSearch).pairs

"""
    global_search_due(s::ContactSearch)

`true` if the next `compute_contact_pairs!` of `s` searches globally.
"""
global_search_due(s::ContactSearch) = s.step == 0

"""
    compute_contact_pairs!(s::ContactSearch)

The contact pairs of the current step (`s.pairs`); every `s.frequency` steps after a
global search for the candidates.
"""
function compute_contact_pairs!(s::ContactSearch{D}) where {D}
    positions = Data_Manager.get_all_positions()::Matrix{Float64}
    s.step == 0 && global_search!(s, positions)
    s.step = s.step + 1 == s.frequency ? 0 : s.step + 1
    find_pairs!(s, positions)
    return nothing
end

function fill_coordinates!(coords::Matrix{Float64}, nodes::Vector{Int64},
                           positions::Matrix{Float64}, ::Val{D}) where {D}
    @inbounds for (k, id) in enumerate(nodes), d in 1:D
        coords[d, k] = positions[id, d]
    end
    return coords
end

"""
    global_search!(s::ContactSearch, positions::Matrix{Float64})

Collects the candidate pairs: the slaves within the search radius of every master.
"""
function global_search!(s::ContactSearch{D}, positions::Matrix{Float64}) where {D}
    empty!(s.candidate_masters)
    resize!(s.candidate_offsets, 1)
    s.candidate_offsets[1] = 1
    empty!(s.candidate_slaves)
    empty!(s.candidate_nodes)
    (isempty(s.master_nodes) || isempty(s.slave_nodes)) && return nothing
    fill_coordinates!(s.master_coords, s.master_nodes, positions, Val(D))
    fill_coordinates!(s.slave_coords, s.slave_nodes, positions, Val(D))
    # the grid holds the slave nodes; it is built once and then only updated
    if s.grid_initialized
        update!(s.grid, s.master_coords, s.slave_coords)
    else
        initialize_grid!(s.grid, s.slave_coords)
        s.grid_initialized = true
    end
    for k in eachindex(s.master_nodes)
        nfound = length(s.candidate_slaves)
        foreach_neighbor(s.master_coords, s.slave_coords, s.grid, k;
                         search_radius = s.search_radius) do _, j, _, _
            push!(s.candidate_slaves, s.slave_nodes[j])
        end
        length(s.candidate_slaves) == nfound && continue
        push!(s.candidate_masters, s.master_nodes[k])
        push!(s.candidate_offsets, length(s.candidate_slaves) + 1)
    end
    append!(s.candidate_nodes, s.candidate_masters)
    append!(s.candidate_nodes, s.candidate_slaves)
    unique!(sort!(s.candidate_nodes))
    return nothing
end

"""
    find_pairs!(s::ContactSearch, positions::Matrix{Float64})

The candidate pairs within the contact radius: `s.pairs` with the distances and the
normals; flags the nodes in the `Contact Nodes` field (master 5, slave 8).
"""
function find_pairs!(s::ContactSearch{D}, positions::Matrix{Float64}) where {D}
    p = empty!(s.pairs)
    r2 = s.contact_radius^2
    mapping = Data_Manager.get_exchange_id_to_local_id()::Dict{Int64,Int64}
    contact_nodes = Data_Manager.get_field("Contact Nodes")::Vector{Int64}
    @inbounds for k in eachindex(s.candidate_masters)
        master = s.candidate_masters[k]
        xm = SVector{D,Float64}(ntuple(d -> positions[master, d], Val(D)))
        nfound = length(p.slaves)
        for c in s.candidate_offsets[k]:(s.candidate_offsets[k + 1] - 1)
            slave = s.candidate_slaves[c]
            xs = SVector{D,Float64}(ntuple(d -> positions[slave, d], Val(D)))
            # the criterion of the grid search (`foreach_neighbor`)
            diff = xm - xs
            dot(diff, diff) <= r2 || continue
            distance, normal = compute_distance_and_normals(xm, xs)
            push!(p.slaves, slave)
            push!(p.distances, distance)
            append!(p.normals, normal)
            local_slave = get(mapping, slave, 0)
            local_slave > 0 && (contact_nodes[local_slave] = 8)
        end
        length(p.slaves) == nfound && continue
        push!(p.masters, master)
        push!(p.offsets, length(p.slaves) + 1)
        local_master = get(mapping, master, 0)
        local_master > 0 && (contact_nodes[local_master] = 5)
    end
    return p
end

end
