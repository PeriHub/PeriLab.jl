# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Contact

using TimerOutputs: @timeit

using .....Data_Manager
using .....PeriLabExceptions: @abort
using .....ModuleLoader: find_registered_modules
using ......ParameterSpec: model_module
using .....InputDeck: ContactInput, contact_blocks, contact_search_frequency
# before the contact models: they read the contact pairs (`Contact_Search.contact_pairs`)
include("Contact_search.jl")
using .Contact_Search: ContactSearch, compute_contact_pairs!, global_search_due
for file in find_registered_modules(@__DIR__, "register_contact")
    include(file)
end

using LinearAlgebra
using .....MPI_Communication: find_and_set_core_value_sum
using .....Helpers: remove_ids, get_block_nodes, compute_free_surface_nodes, find_indices,
                    compute_geometry, get_surface_information, each_call
export init_contact_model
export compute_contact_model

"""
    init_contact_model(contact)

Initializes the contact model.

# Arguments
- `contact::ContactInput`: The typed `Contact` section.
"""
function init_contact_model(contact::ContactInput)
    @info "Init Contact Model"
    Data_Manager.set_contact_properties(contact)

    check_valid_contact_model(contact, Data_Manager.get_all_blocks())
    contact_block_list = contact_blocks(contact)
    # get all the contact block surface global ids and reduce the exchange positions to this points.
    # all following functions deal with the local contact position id.
    only_surface = contact.globals.only_surface_contact_nodes

    global_contact_ids = identify_contact_block_nodes(contact_block_list,
                                                      only_surface)
    contact_nodes = Data_Manager.create_constant_node_scalar_field("Contact Nodes", Int64)
    block_list = Data_Manager.get_all_blocks()

    mapping = contact_block_ids(global_contact_ids, block_list, contact_block_list)

    Data_Manager.set_contact_block_ids(mapping)
    # identify all surface which have no neighboring nodes
    if only_surface
        free_surfaces = identify_free_contact_surfaces(contact_block_list)
    end
    points = Data_Manager.get_all_positions()

    block_nodes = get_block_nodes(block_list, length(block_list)) # all ids
    contact_index = Dict(id => index for (index, id) in enumerate(global_contact_ids))

    for model in values(contact.models)
        for (cg, group) in pairs(model.base.contact_groups)
            slave_id = group.slave_block_id
            master_id = group.master_block_id
            @info "Contact pair Master block $master_id - Slave block $slave_id"
            if !only_surface
                # the positions are reduced to global_contact_ids below, so the free
                # contact nodes are given as indices into that list, as in
                # compute_and_set_free_contact_nodes
                Data_Manager.set_free_contact_nodes(master_id,
                                                    [contact_index[id]
                                                     for id in block_nodes[master_id]])
                Data_Manager.set_free_contact_nodes(slave_id,
                                                    [contact_index[id]
                                                     for id in block_nodes[slave_id]])
                continue
            end

            compute_and_set_free_contact_nodes(points,
                                               block_nodes[slave_id],
                                               free_surfaces[slave_id],
                                               global_contact_ids,
                                               slave_id)
            compute_and_set_free_contact_nodes(points,
                                               block_nodes[master_id],
                                               free_surfaces[master_id],
                                               global_contact_ids,
                                               master_id)
        end
    end
    # reduce the block_list
    Data_Manager.set_all_blocks(block_list[global_contact_ids])
    Data_Manager.set_all_positions(points[global_contact_ids, :])

    shared_vol = Data_Manager.create_constant_free_size_field("Shared Volumes", Float64,
                                                              (length(global_contact_ids),
                                                               1))
    shared_horizon = Data_Manager.create_constant_free_size_field("Shared Horizon", Float64,
                                                                  (length(global_contact_ids),
                                                                   1))
    create_local_contact_id_mapping(global_contact_ids)

    shared_vol[:] = synchronize_contact_points("Volume", "Constant",
                                               shared_vol,
                                               Data_Manager.get_local_contact_ids())

    shared_horizon[:] = synchronize_contact_points("Horizon", "Constant",
                                                   shared_horizon,
                                                   Data_Manager.get_local_contact_ids())

    init_contact_search(contact)

    @info "Set contact models"
    for model in values(contact.models)
        contact_module(model).init_contact_model(model.model, model)
    end

    @info "Finish Init Contact Model"
end

"""
    ContactCall{S,M,F}

A contact group of a run: its name, search state, contact model of the input deck
(`model.model` the typed parameters) and the model's `compute_contact_model`.
"""
struct ContactCall{S<:ContactSearch,M,F}
    name::String
    search::S
    model::M
    compute::F
end

"""
    ContactState{C}

The contact groups as a tuple of `ContactCall` (its type holds the models, so every
group is computed statically) and the deformed coordinates this rank contributes to the
contact positions: its controller nodes `local_ids` with their rows `exchange_ids`.
Between two global searches only the `candidate_rows` (the candidate nodes of all
contact groups, `candidate_local` their local id or 0) are synchronized, over `buffer`.
"""
struct ContactState{C<:Tuple}
    calls::C
    local_ids::Vector{Int64}
    exchange_ids::Vector{Int64}
    candidate_rows::Vector{Int64}
    candidate_local::Vector{Int64}
    buffer::Vector{Float64}
end

"""
    init_contact_search(contact::ContactInput)

Creates the search state of every contact group (`Contact_Search.ContactSearch`) and
which coordinates this rank synchronizes.
"""
function init_contact_search(contact::ContactInput)
    dof = Data_Manager.get_dof()
    calls = []
    for model in values(contact.models), (cg, group) in pairs(model.base.contact_groups)
        search = ContactSearch(dof,
                               Data_Manager.get_free_contact_nodes(group.master_block_id),
                               Data_Manager.get_free_contact_nodes(group.slave_block_id),
                               model.base.contact_radius, group.search_radius,
                               contact_search_frequency(group, contact.globals))
        Data_Manager.set_contact_search(cg, search)
        push!(calls,
              ContactCall(cg, search, model, contact_module(model).compute_contact_model))
    end
    order = sort!(collect(Data_Manager.get_local_contact_ids()::Dict{Int64,Int64});
                  by = last)
    Data_Manager.set_contact_state(ContactState(Tuple(calls), first.(order), last.(order),
                                                Int64[], Int64[], Float64[]))
end

function create_local_contact_id_mapping(global_contact_ids)
    mapping = Dict{Int64,Int64}()
    inv_mapping = Dict{Int64,Int64}()
    # all controler nodes, because only they are allowed to communicate with the exchange_id
    nnodes = Data_Manager.get_nnodes()
    for (id, pid) in enumerate(global_contact_ids)
        local_id = Data_Manager.get_local_nodes([pid])
        if isempty(local_id) || local_id[1] > nnodes
            continue
        end
        # local id -> exchange id
        mapping[local_id[1]] = id
        # exchange id -> local id
        inv_mapping[id] = local_id[1]
    end
    Data_Manager.set_local_contact_ids(mapping)
    Data_Manager.set_exchange_id_to_local_id(inv_mapping)
end

"""
    contact_block_ids(global_ids::Vector{Int64}, block_list, contact_blocks)

Finds the ids which are in the specific contact related blocks.


"""
function contact_block_ids(global_ids::Vector{Int64}, block_list, contact_blocks)
    mapping = Dict{Int64,Vector{Int64}}()
    for cb in contact_blocks
        mapping[cb] = []
    end
    for (id, glob_id) in enumerate(global_ids)
        if block_list[glob_id] in contact_blocks
            # append the block node to the block list
            append!(mapping[block_list[glob_id]], id)
        end
    end
    return mapping
end

# the module of a contact model: the one that declares its parameter struct
contact_module(model) = model_module(model.model)

"""
    compute_contact_model(contact, time, dt)

Searches the contact pairs and computes the forces of every contact model.

# Arguments
- `contact::ContactInput`: The typed `Contact` section.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
"""
function compute_contact_model(contact::ContactInput,
                               time::Float64,
                               dt::Float64)
    # function barrier: the contact state holds the typed contact groups
    return compute_contacts!(Data_Manager.get_contact_state())
end

function compute_contacts!(state::ContactState)
    # all positions for a global search, else only the candidate nodes
    global_step = any_global_search(state.calls)
    @timeit "synchronize contact positions" synchronize_contact_positions!(state,
                                                                           global_step)
    @timeit "Contact search" each_call(compute_contact!, state.calls)
    global_step && update_candidate_rows!(state)
    return nothing
end

any_global_search(::Tuple{}) = false
function any_global_search(calls::Tuple)
    return global_search_due(first(calls).search) || any_global_search(Base.tail(calls))
end

function compute_contact!(call::ContactCall)
    @timeit "compute_contact_pairs" compute_contact_pairs!(call.search)
    @timeit "compute_contact_model" call.compute(call.name, call.model.model, call.model,
                                                 compute_master_force_density,
                                                 compute_slave_force_density)
    return nothing
end

"""
    synchronize_contact_positions!(sync::ContactState, global_step::Bool)

Writes the deformed coordinates into the contact positions (`get_all_positions`): all
contact nodes for a global search, else only the candidate nodes. With MPI every rank adds
its controller nodes and the rows are summed over the ranks.
"""
function synchronize_contact_positions!(sync::ContactState, global_step::Bool)
    positions = Data_Manager.get_all_positions()::Matrix{Float64}
    coor = Data_Manager.get_field("Deformed Coordinates", "NP1")::NodeVectorField{Float64}
    dof = size(positions, 2)
    parallel = Data_Manager.get_max_rank() > 1
    if global_step
        fill!(positions, 0.0)
        @inbounds for (local_id, row) in zip(sync.local_ids, sync.exchange_ids),
                      d in 1:dof
            positions[row, d] = coor[local_id, d]
        end
        parallel && find_and_set_core_value_sum(Data_Manager.get_comm(), positions)
        return nothing
    end
    if !parallel
        @inbounds for (row, local_id) in zip(sync.candidate_rows, sync.candidate_local)
            local_id == 0 && continue
            for d in 1:dof
                positions[row, d] = coor[local_id, d]
            end
        end
        return nothing
    end
    buffer = sync.buffer
    @inbounds for (k, (row, local_id)) in enumerate(zip(sync.candidate_rows,
                                                        sync.candidate_local)),
                  d in 1:dof
        buffer[(k - 1) * dof + d] = local_id == 0 ? 0.0 : coor[local_id, d]
    end
    find_and_set_core_value_sum(Data_Manager.get_comm(), buffer)
    @inbounds for (k, row) in enumerate(sync.candidate_rows), d in 1:dof
        positions[row, d] = buffer[(k - 1) * dof + d]
    end
    return nothing
end

# the candidate nodes of all contact groups after a global search, synchronized until the
# next one
function update_candidate_rows!(sync::ContactState)
    rows = sync.candidate_rows
    empty!(rows)
    each_call(append_candidates!, sync.calls, rows)
    unique!(sort!(rows))
    mapping = Data_Manager.get_exchange_id_to_local_id()::Dict{Int64,Int64}
    resize!(sync.candidate_local, length(rows))
    for (k, row) in enumerate(rows)
        sync.candidate_local[k] = get(mapping, row, 0)
    end
    resize!(sync.buffer, length(rows) * size(Data_Manager.get_all_positions(), 2))
    return nothing
end

function append_candidates!(call::ContactCall, rows::Vector{Int64})
    append!(rows,
            call.search.candidate_nodes)
end

function compute_master_force_density(id_m, id_s, contact_force)
    compute_force_density(id_m, id_s, -contact_force)
end

function compute_slave_force_density(id_s, id_m, contact_force)
    compute_force_density(id_s, id_m, contact_force)
end

function compute_force_density(id_1::Int64, id_2::Int64,
                               contact_force::AbstractVector{Float64})
    mapping::Dict{Int64,Int64} = Data_Manager.get_exchange_id_to_local_id()

    if isnothing(get(mapping, id_1, nothing))
        return
    end
    shared_volume::Matrix{Float64} = Data_Manager.get_field("Shared Volumes")
    force_densities::NodeVectorField{Float64} = Data_Manager.get_field("Force Densities",
                                                                       "NP1")
    # element-wise: a row view of the force densities would be created for every pair
    local_id = mapping[id_1]
    volume = shared_volume[id_2]
    @inbounds for d in eachindex(contact_force)
        force_densities[local_id, d] += contact_force[d] * volume
    end
    #println(force_densities[mapping[id_1], :])
end
"""
    identify_contact_block_surface_nodes(contact_blocks)
The ids of the surfaces is identified. This is needed to reduce the number of ids needed to create the polyhedra and to reduce the exchange vector.

"""
function identify_contact_block_nodes(contact_blocks, only_surface)
    global_ids = Vector{Int64}([])
    points = Data_Manager.get_all_positions()
    block_ids = Data_Manager.get_all_blocks() # all ids
    block_nodes = get_block_nodes(block_ids, length(block_ids)) # all ids
    for block in contact_blocks
        if only_surface
            append!(global_ids, get_surface_points(points, block_nodes[block]))
        else
            append!(global_ids, block_nodes[block])
        end
    end
    return sort(unique(global_ids))
end

function get_surface_points(points, block_nodes)
    global_ids = Vector{Int64}([])
    poly = compute_geometry(points[block_nodes, :])
    normals, offsets = get_surface_information(poly)
    for pID in block_nodes
        for id in eachindex(offsets)
            if isapprox(dot(normals[id, :], points[pID, :]) - offsets[id], 0;
                        atol = 1e-6, rtol = 1e-5)
                append!(global_ids, pID)
            end
        end
    end
    return global_ids
end

function identify_free_contact_surfaces(contact_blocks::Vector{Int64})
    # all blocks, not only the contact blocks to find the free surfaces
    all_positions = Data_Manager.get_all_positions()
    block_ids = Data_Manager.get_all_blocks()
    # I need all block nodes to check contact blocks with there neighbor non contact blocks
    block_nodes = get_block_nodes(block_ids, length(block_ids))
    free_surfaces = Dict{Int64,Vector{Int64}}()
    # loop to init the surfaces of a block in correct order (1:nnumber)
    for iID in contact_blocks
        poly_i = compute_geometry(all_positions[block_nodes[iID], :])
        normals_i, offsets_i = get_surface_information(poly_i)
        free_surfaces[iID] = Vector{Int64}(collect(1:length(offsets_i)))
    end
    # loop to remove the surfaces which are not free
    for iID in contact_blocks
        # Polyhedra for the master block; can be 2D or 3D
        poly_master = compute_geometry(all_positions[block_nodes[iID], :])
        normals_i, offsets_i = get_surface_information(poly_master)
        for jID in Data_Manager.get_block_id_list()
            if iID == jID
                continue
            end
            # Polyhedra for the slave block; can be 2D or 3D
            poly_slave = compute_geometry(all_positions[block_nodes[jID], :])
            normals_j, offsets_j = get_surface_information(poly_slave)
            # gives lists where all surfaces are included, which are equal
            ids = get_double_surfs(normals_i, offsets_i, normals_j, offsets_j)
            remove_ids(free_surfaces[iID], ids)
        end
    end
    # store these free surfaces
    return free_surfaces
end
"""
    compute_and_set_free_contact_nodes(all_positions, id::Int64)

Find the node ids of the free surfaces of the contact blocks and there connected surfaces. The node id is given in the exchange vector id.

"""
function compute_and_set_free_contact_nodes(all_positions, ids,
                                            free_surfaces, global_ids::Vector{Int64},
                                            block_id::Int64)
    free_nodes = compute_free_surface_nodes(all_positions,
                                            ids, free_surfaces)
    # in this routine the complete not reduced exchange ids are used. Therefore, the reduced global_ids list must be used to give the free surface node list the correct ids.
    sf = zeros(Int64, length(free_nodes))
    for (id, fid) in enumerate(free_nodes)
        sf[id] = find_indices(global_ids, fid)[1]
    end
    Data_Manager.set_free_contact_nodes(block_id, sf)
end

"""
    synchronize_contact_points(what::String, step::String,
                                    synch_vector))

Synchronises the deformation relevant for contact. These are the deformed coordinates of the surfaces of the contact blocks. The deformed surface vector is set to zero and than filled at each core from data from the deformed coordinate field. To bring it together the vector is summed at the root and send back to all cores.

tbd
# Returns
- nothing
"""

function synchronize_contact_points(what::String, step::String,
                                    synch_vector, local_map)
    sy = Data_Manager.get_field(what, step)

    synch_vector .= 0
    for (local_id, exchange_id) in pairs(local_map)
        synch_vector[exchange_id, :] .= sy[local_id, :]
    end
    if Data_Manager.get_max_rank() == 1
        return synch_vector
    end
    comm = Data_Manager.get_comm()

    return find_and_set_core_value_sum(comm, synch_vector)
end

## TODO check if this method holds true for all surfaces, because some with same size in the same half plane might be identified as equal
function get_double_surfs(normals_i, offsets_i, normals_j, offsets_j)
    ids = Vector{Int64}([])
    for iID in eachindex(offsets_i)
        for jID in eachindex(offsets_j)
            if isapprox(offsets_i[iID], offsets_j[jID], atol = 1e-6, rtol = 1e-5) &&
               isapprox(normals_i[iID, :], normals_j[jID, :];
                        atol = 1e-6)
                append!(ids, iID)
                break # only one surface pair can exist
            end
        end
    end
    return ids
end
function check_valid_contact_model(contact::ContactInput, block_ids)
    # an inverse pair (1-2 in one group, 2-1 in another) is not allowed
    check_dict = Dict{Int64,Int64}()
    for model in values(contact.models), group in values(model.base.contact_groups)
        master = group.master_block_id
        slave = group.slave_block_id
        if !(master in block_ids)
            @abort "Block defintion in master does not exist."
            return
        end
        if !(slave in block_ids)
            @abort "Block defintion in slave does not exist."
            return
        end
        check_dict[master] = slave
        if get(check_dict, slave, nothing) == master
            @abort "Master and Slave should be defined in an inverse way, e.g. Master = 1, Slave = 2 in model 1 and Master = 2, Slave = 1 in model 2."
            return
        end
    end
    return nothing
end
end
