# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Contact

using ...Solver_Manager: find_module_files, create_module_specifics
global module_list = find_module_files(@__DIR__, "contact_model_name")
for mod in module_list
    include(mod["File"])
end

using TimerOutputs
using LinearAlgebra
using .....MPI_Communication: find_and_set_core_value_sum
include("Contact_search.jl")
using .Contact_Search: init_contact_search, compute_geometry, get_surface_information,
                       compute_contact_pairs
using .....Helpers: remove_ids, get_block_nodes, compute_free_surface_nodes, find_indices
export init_contact_model
export compute_contact_model

"""
    init_model(datamanager::Module, params)

Initializes the contact model.

# Arguments
- `datamanager::Data_Manager`: Datamanager
- `params:`: Contact parameter.
# Returns
- `datamanager::Data_Manager`: Datamanager.
"""
function init_contact_model(datamanager::Module, params)
    @info "Init Contact Model"
    datamanager.set_contact_properties(params)

    check_valid_contact_model(params, datamanager.get_all_blocks())
    contact_blocks = get_all_contact_blocks(params)
    # get all the contact block surface global ids and reduce the exchange positions to this points.
    # all following functions deal with the local contact position id.
    if !haskey(params, "Globals")
        params["Globals"] = Dict()
        if isnothing(get(params["Globals"], "Only Surface Contact Nodes", nothing))
            params["Globals"]["Global Search Frequency"] = 1
        elseif params["Globals"]["Global Search Frequency"] < 1
            @warn "Globals Search Frequency must be greater than zero and set to 1."
            params["Globals"]["Global Search Frequency"] = 1
        end
    end
    only_surface = get(params["Globals"], "Only Surface Contact Nodes", true)

    global_contact_ids = identify_contact_block_nodes(datamanager, contact_blocks,
                                                      only_surface)
    contact_nodes = datamanager.create_constant_node_field("Contact Nodes", Int64, 1)
    block_list = datamanager.get_all_blocks()

    mapping = contact_block_ids(global_contact_ids, block_list, contact_blocks)

    datamanager.set_contact_block_ids(mapping)
    # identify all surface which have no neighboring nodes
    if only_surface
        free_surfaces = identify_free_contact_surfaces(datamanager, contact_blocks)
    end
    points = datamanager.get_all_positions()

    block_nodes = get_block_nodes(block_list, length(block_list)) # all ids

    for contact_model in filter(k -> k != "Globals", keys(params))
        for (cg, contact_params) in pairs(params[contact_model]["Contact Groups"])
            if !haskey(contact_params, "Global Search Frequency")
                contact_params["Global Search Frequency"] = params["Globals"]["Global Search Frequency"]
            end
            datamanager.set_search_step(cg, 0)
            slave_id = contact_params["Slave Block ID"]
            master_id = contact_params["Master Block ID"]
            @info "Contact pair Master block $master_id - Slave block $slave_id"
            if !only_surface
                datamanager.set_free_contact_nodes(master_id, block_nodes[master_id])
                datamanager.set_free_contact_nodes(slave_id, block_nodes[slave_id])
                continue
            end

            compute_and_set_free_contact_nodes(datamanager, points,
                                               block_nodes[slave_id],
                                               free_surfaces[slave_id],
                                               global_contact_ids,
                                               slave_id)
            compute_and_set_free_contact_nodes(datamanager, points,
                                               block_nodes[master_id],
                                               free_surfaces[master_id],
                                               global_contact_ids,
                                               master_id)
        end
    end
    # reduce the block_list
    datamanager.set_all_blocks(block_list[global_contact_ids])
    datamanager.set_all_positions(points[global_contact_ids, :])

    shared_vol = datamanager.create_constant_free_size_field("Shared Volumes", Float64,
                                                             (length(global_contact_ids),
                                                              1))
    shared_horizon = datamanager.create_constant_free_size_field("Shared Horizon", Float64,
                                                                 (length(global_contact_ids),
                                                                  1))
    create_local_contact_id_mapping(datamanager, global_contact_ids)

    shared_vol[:] = synchronize_contact_points(datamanager, "Volume", "Constant",
                                               shared_vol,
                                               datamanager.get_local_contact_ids())

    shared_horizon[:] = synchronize_contact_points(datamanager, "Horizon", "Constant",
                                                   shared_horizon,
                                                   datamanager.get_local_contact_ids())

    @info "Set contact models"
    for (cm, contact_params) in pairs(params)
        if cm == "Globals"
            continue
        end
        init_contact_search(datamanager, contact_params, cm)

        mod = create_module_specifics(contact_params["Type"],
                                      module_list,
                                      "contact_model_name")
        if isnothing(mod)
            @error "No contact model of type " * contact_params["Type"] *
                   " exists."
            return nothing
        end
        datamanager.set_model_module(contact_params["Type"], mod)
        datamanager = mod.init_contact_model(datamanager, contact_params)
    end

    @info "Finish Init Contact Model"
    return datamanager
end

function create_local_contact_id_mapping(datamanager, global_contact_ids)
    mapping = Dict{Int64,Int64}()
    inv_mapping = Dict{Int64,Int64}()
    # all controler nodes, because only they are allowed to communicate with the exchange_id
    nnodes = datamanager.get_nnodes()
    for (id, pid) in enumerate(global_contact_ids)
        local_id = datamanager.get_local_nodes([pid])
        if isempty(local_id) || local_id[1] > nnodes
            continue
        end
        # local id -> exchange id
        mapping[local_id[1]] = id
        # exchange id -> local id
        inv_mapping[id] = local_id[1]
    end
    datamanager.set_local_contact_ids(mapping)
    datamanager.set_exchange_id_to_local_id(inv_mapping)
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

function get_all_contact_blocks(params)
    contact_blocks = Vector{Int64}([])
    for contact_model in filter(k -> k != "Globals", keys(params))
        for contact_groups in values(params[contact_model]["Contact Groups"])
            append!(contact_blocks, contact_groups["Master Block ID"])
            append!(contact_blocks, contact_groups["Slave Block ID"])
        end
    end
    return Vector{Int64}(sort(unique(contact_blocks)))
end

"""
    compute_model(datamanager::Module, nodes::AbstractVector{Int64}, model_param::Dict, block::Int64, time::Float64, dt::Float64, to::TimerOutput)

Compute the forces of the contact model.

# Arguments
- `datamanager::Data_Manager`: Datamanager.
- `nodes::AbstractVector{Int64}`: The nodes.
- `model_param::Dict`: The contact parameter.
- `block::Int64`: The current block.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
# Returns
- `datamanager::Data_Manager`: Datamanager.
"""
function compute_contact_model(datamanager::Module,
                               contact_params::Dict,
                               time::Float64,
                               dt::Float64,
                               to::TimerOutput)
    # computes and synchronizes the relevant positions
    all_positions = datamanager.get_all_positions()
    if datamanager.global_synchronization()
        datamanager.clear_synchronization_list()
        all_positions = synchronize_contact_points(datamanager,
                                                   "Deformed Coordinates",
                                                   "NP1",
                                                   all_positions,
                                                   datamanager.get_local_contact_ids())
    else
        # local synchronization occurs if Global Search Frequency > 1
        update_list = datamanager.synchronization_list()
        all_positions = synchronize_contact_points(datamanager,
                                                   "Deformed Coordinates",
                                                   "NP1",
                                                   all_positions, update_list)
    end
    datamanager.set_all_positions(all_positions)
    @timeit to "Contact search" begin
        for (cm, block_contact_params) in pairs(contact_params)
            if cm == "Globals"
                continue
            end

            mod = datamanager.get_model_module(block_contact_params["Type"])
            for (cg, block_contact_group) in pairs(block_contact_params["Contact Groups"])
                n = datamanager.get_search_step(cg) + 1

                # needed in search and in model evaluation
                block_contact_group["Contact Radius"] = block_contact_params["Contact Radius"]
                datamanager.set_contact_dict(cg, Dict())

                @timeit to "compute_contact_pairs" compute_contact_pairs(datamanager,
                                                                         cg,
                                                                         block_contact_group)
                @timeit to "compute_contact_model" datamanager=mod.compute_contact_model(datamanager,
                                                                                         cg,
                                                                                         block_contact_params,
                                                                                         compute_master_force_density,
                                                                                         compute_slave_force_density)
                if n == block_contact_group["Global Search Frequency"]
                    datamanager.set_search_step(cg, 0)
                end
            end
        end
        #compute_contact()
    end
    return datamanager
end

function compute_master_force_density(datamanager, id_m, id_s, contact_force)
    compute_force_density(datamanager, id_m, id_s, -contact_force)
end

function compute_slave_force_density(datamanager, id_s, id_m, contact_force)
    compute_force_density(datamanager, id_s, id_m, contact_force)
end

function compute_force_density(datamanager, id_1, id_2, contact_force)
    mapping = datamanager.get_exchange_id_to_local_id()

    if isnothing(get(mapping, id_1, nothing))
        return
    end
    shared_volume = datamanager.get_field("Shared Volumes")
    force_densities = datamanager.get_field("Force Densities", "NP1")

    force_densities[mapping[id_1], :] .+= contact_force .* shared_volume[id_2]
    #println(force_densities[mapping[id_1], :])
end
"""
    identify_contact_block_surface_nodes(datamanager, contact_blocks)
The ids of the surfaces is identified. This is needed to reduce the number of ids needed to create the polyhedra and to reduce the exchange vector.

"""
function identify_contact_block_nodes(datamanager, contact_blocks, only_surface)
    global_ids = Vector{Int64}([])
    points = datamanager.get_all_positions()
    block_ids = datamanager.get_all_blocks() # all ids
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

function identify_free_contact_surfaces(datamanager::Module, contact_blocks::Vector{Int64})
    # all blocks, not only the contact blocks to find the free surfaces
    all_positions = datamanager.get_all_positions()
    block_ids = datamanager.get_all_blocks()
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
        for jID in datamanager.get_block_id_list()
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
    compute_and_set_free_contact_nodes(datamanager::Module, all_positions, id::Int64)

Find the node ids of the free surfaces of the contact blocks and there connected surfaces. The node id is given in the exchange vector id.

"""
function compute_and_set_free_contact_nodes(datamanager::Module, all_positions, ids,
                                            free_surfaces, global_ids::Vector{Int64},
                                            block_id::Int64)
    free_nodes = compute_free_surface_nodes(all_positions,
                                            ids, free_surfaces)
    # in this routine the complete not reduced exchange ids are used. Therefore, the reduced global_ids list must be used to give the free surface node list the correct ids.
    sf = zeros(Int64, length(free_nodes))
    for (id, fid) in enumerate(free_nodes)
        sf[id] = find_indices(global_ids, fid)[1]
    end
    datamanager.set_free_contact_nodes(block_id, sf)
end

"""
    synchronize_contact_points(datamanager::Module, what::String, step::String,
                                    synch_vector))

Synchronises the deformation relevant for contact. These are the deformed coordinates of the surfaces of the contact blocks. The deformed surface vector is set to zero and than filled at each core from data from the deformed coordinate field. To bring it together the vector is summed at the root and send back to all cores.

# Arguments
- `datamanager::Modul`: datamanager
tbd
# Returns
- nothing
"""

function synchronize_contact_points(datamanager::Module, what::String, step::String,
                                    synch_vector, local_map)
    sy = datamanager.get_field(what, step)

    synch_vector .= 0
    for (local_id, exchange_id) in pairs(local_map)
        synch_vector[exchange_id, :] .= sy[local_id, :]
    end
    if datamanager.get_max_rank() == 1
        return synch_vector
    end
    comm = datamanager.get_comm()

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
function check_valid_contact_model(params, block_ids::Vector{Int64})
    # inverse master slave check
    # tuple liste bauen und dann die neuen invers checken
    check_dict = Dict{Int64,Int64}()
    for contact_model in filter(k -> k != "Globals", keys(params))
        for contact_groups in values(params[contact_model]["Contact Groups"])
            if !haskey(contact_groups, "Master Block ID")
                @error "Contact model needs a ''Master''"
                return false
            end
            if !haskey(contact_groups, "Slave Block ID")
                @error "Contact model needs a ''Slave''"
                return false
            end
            if contact_groups["Master Block ID"] == contact_groups["Slave Block ID"]
                @error "Contact master and slave are equal. Self contact is not implemented yet."
                return false
            end

            if !(contact_groups["Master Block ID"] in block_ids)
                @error "Block defintion in master does not exist."
                return false
            end
            if !(contact_groups["Slave Block ID"] in block_ids)
                @error "Block defintion in slave does not exist."
                return false
            end
            check_dict[contact_groups["Master Block ID"]] = contact_groups["Slave Block ID"]
            if haskey(check_dict, contact_groups["Slave Block ID"]) &&
               check_dict[contact_groups["Slave Block ID"]] ==
               contact_groups["Master Block ID"]
                @error "Master and Slave should be defined in an inverse way, e.g. Master = 1, Slave = 2 in model 1 and Master = 2, Slave = 1 in model 2."
                return false
            end
            if !haskey(contact_groups, "Search Radius")
                @error "Contact model needs a ''Search Radius''."
                return false
            end
            if contact_groups["Search Radius"] <= 0
                @error "''Search Radius'' must be greater than zero."
                return false
            end
        end
    end
    return true
end
end
