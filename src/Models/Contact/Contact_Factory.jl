# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

module Contact_Factory
using TimerOutputs
using LinearAlgebra
include("../../MPI_communication/MPI_communication.jl")
using .MPI_communication: find_and_set_core_value_sum
include("Contact_search.jl")
using .Contact_search: init_contact_search, compute_geometry, get_surface_information,
                       compute_contact_pairs
include("Penalty_model.jl")
using .Penalty_model
include("../../Support/Helpers.jl")
using .Helpers: remove_ids, get_block_nodes
include("../../Core/Module_inclusion/set_Modules.jl")
using .Set_modules: find_module_files, include_files
global module_list = find_module_files(@__DIR__, "contact_model_name")
include_files(module_list)
export init_contact_model
export compute_contact_model

"""
    init_model(datamanager::Module, params)

Initializes the contact model.

# Arguments
- `datamanager::Data_manager`: Datamanager
- `params:`: Contact parameter.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function init_contact_model(datamanager::Module, params)
    @info "Init Contact Model"
    datamanager.set_contact_properties(params)

    check_valid_contact_model(params, datamanager.get_all_blocks())
    contact_blocks = get_all_contact_blocks(params)
    # get all the contact block surface global ids and reduce the exchange positions to this points.
    # all following functions deal with the local contact position id.
    global_contact_ids = identify_contact_block_surface_nodes(datamanager, contact_blocks)

    block_list = datamanager.get_all_blocks()

    mapping = contact_block_ids(global_contact_ids, block_list, contact_blocks)
    # hier ist der Fehler; welche nummer brauche ich die lokale Liste oder die globale?

    datamanager.set_contact_block_ids(mapping)

    identify_outer_contact_surfaces(datamanager, contact_blocks)

    #mapping = contact_block_free_surface_ids(global_contact_ids, block_list, contact_blocks)
    #datamanager.set_contact_block_free_surface_ids(mapping)
    # reduce the block_list
    datamanager.set_all_blocks(block_list[global_contact_ids])
    points = datamanager.get_all_positions()

    datamanager.set_all_positions(points[global_contact_ids, :])

    shared_vol = datamanager.create_constant_free_size_field("Shared Volumes", Float64,
                                                             (length(global_contact_ids),
                                                              1))

    create_local_contact_id_mapping(datamanager, global_contact_ids)
    shared_vol = synchronize_contact_points(datamanager, "Volume", "Constant", shared_vol)
    for (cm, contact_params) in pairs(params)
        init_contact_search(datamanager, contact_params, cm)
        Penalty_model.init_contact_model(datamanager, contact_params)
    end

    @info "Finish Init Contact Model"
    return datamanager
end

function create_local_contact_id_mapping(datamanager, global_contact_ids)
    mapping = Dict{Int64,Int64}()
    inv_mapping = Dict{Int64,Int64}()
    for (id, pid) in enumerate(global_contact_ids)
        local_id = datamanager.get_local_nodes([pid])
        if local_id == []
            continue
        end
        mapping[local_id[1]] = id
        inv_mapping[id] = local_id[1]
    end
    datamanager.set_local_contact_ids(mapping)
    datamanager.set_exchange_id_to_local_id(inv_mapping)
end

function loc_to_contact_exchange_id(global_contact_ids::Vector{Int64})
    mapping = Dict{Int64,Int64}()
    for (id, num) in enumerate(global_contact_ids)
        mapping[id] = num
    end
    return mapping
end
# give the ids for the contact blocks in the exchange vector
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

function apply_contact_forces(datamanager)
    contact_force = datamanager.get_field("Contact Forces")
    force = datamanager.get_field("Force Densities", "NP1")
    local_ids = datamanager.get_local_contact_ids()
    for (id, local_id) in enumerate(values(local_ids))
        force[local_id, :] .+= contact_force[id, :]
    end
    contact_force .= 0.0
end

function get_all_contact_blocks(params)
    contact_blocks = []
    for contact_params in values(params)
        append!(contact_blocks, contact_params["Master"])
        append!(contact_blocks, contact_params["Slave"])
    end
    return sort(unique(contact_blocks))
end

"""
    compute_model(datamanager::Module, nodes::Union{SubArray,Vector{Int64}}, model_param::Dict, block::Int64, time::Float64, dt::Float64, to::TimerOutput)

Compute the forces of the contact model.

# Arguments
- `datamanager::Data_manager`: Datamanager.
- `nodes::Union{SubArray,Vector{Int64}}`: The nodes.
- `model_param::Dict`: The contact parameter.
- `block::Int64`: The current block.
- `time::Float64`: The current time.
- `dt::Float64`: The current time step.
# Returns
- `datamanager::Data_manager`: Datamanager.
"""
function compute_contact_model(datamanager::Module,
                               contact_params::Dict,
                               time::Float64,
                               dt::Float64,
                               to::TimerOutput)
    # computes and synchronizes the relevant positions
    all_positions = datamanager.get_all_positions()
    all_positions = synchronize_contact_points(datamanager, "Deformed Coordinates", "NP1",
                                               all_positions)
    datamanager.set_all_positions(all_positions)

    for (cm, block_contact_params) in pairs(contact_params)
        datamanager.set_contact_dict(cm,
                                     Dict("Pairs: Master-Slave" => Vector{Tuple{Int64,
                                                                                Int64}}([]),
                                          "Normals" => Vector{Array{Float64}}([]),
                                          "Distances" => Vector{Vector{Float64}}([])))
        compute_contact_pairs(datamanager, cm, block_contact_params)
        Penalty_model.compute_contact_model(datamanager, cm, block_contact_params,
                                            compute_master_force_density,
                                            compute_slave_force_density)
    end

    #compute_contact()
    return datamanager
end

function compute_master_force_density(datamanager, id_m, id_s, contact_force)
    compute_force_density(datamanager, id_m, id_s, contact_force)
end

function compute_slave_force_density(datamanager, id_s, id_m, contact_force)
    compute_force_density(datamanager, id_s, id_m, -contact_force)
end

function compute_force_density(datamanager, id_1, id_2, contact_force)
    mapping = datamanager.get_exchange_id_to_local_id()
    if !(id_1 in keys(mapping))
        return
    end
    shared_volume = datamanager.get_field("Shared Volumes")
    force_densities = datamanager.get_field("Force Densities", "NP1")
    force_densities[mapping[id_1], :] .+= 0.5 .* contact_force ./ shared_volume[id_2]
end

function identify_contact_block_surface_nodes(datamanager, contact_blocks)
    global_ids = Vector{Int64}([])
    points = datamanager.get_all_positions()
    block_ids = datamanager.get_all_blocks()
    block_nodes = get_block_nodes(block_ids, length(block_ids))
    for block in contact_blocks
        poly = compute_geometry(points[block_nodes[block], :])
        normals, offsets = get_surface_information(poly)
        for pID in eachindex(points[:, 1])
            if !(points[pID, :] in poly)
                # can happen, e.g. points are on hplane surfraces but not connected to the block
                continue
            end
            for id in eachindex(offsets)
                if isapprox(dot(normals[id, :], points[pID, :]) - offsets[id], 0;
                            atol = 1e-6)
                    append!(global_ids, pID)
                end
            end
        end
    end
    return unique(global_ids)
end

function identify_outer_contact_surfaces(datamanager, contact_blocks)
    # all blocks, not only the contact blocks to find the free surfaces
    all_positions = datamanager.get_all_positions()
    block_ids = datamanager.get_all_blocks()
    block_nodes = get_block_nodes(block_ids, length(block_ids))

    free_surfaces = Dict{Int64,Vector{Int64}}()
    for iID in 1:maximum(keys(block_nodes))
        poly_i = compute_geometry(all_positions[block_nodes[iID], :])
        normals_i, offsets_i = get_surface_information(poly_i)
        free_surfaces[iID] = Vector{Int64}(collect(1:length(offsets_i)))
    end

    for iID in 1:(maximum(block_ids) - 1)
        # Polyhedra for the master block; can be 2D or 3D
        poly_master = compute_geometry(all_positions[block_nodes[iID], :])
        normals_i, offsets_i = get_surface_information(poly_master)
        for jID in (iID + 1):maximum(block_ids)
            # Polyhedra for the slave block; can be 2D or 3D
            poly_slave = compute_geometry(all_positions[block_nodes[jID], :])
            normals_j, offsets_j = get_surface_information(poly_slave)
            # gives lists where all surfaces are included, which are equal
            ids, jds = get_double_surfs(normals_i, offsets_i, normals_j, offsets_j)
            remove_ids(free_surfaces[iID], ids)
            remove_ids(free_surfaces[jID], jds)
        end
    end
    # remove all blocks from the free surface list, which are not master or slave
    for block in keys(free_surfaces)
        if !(block in contact_blocks)
            delete!(free_surfaces, block)
        end
    end
    # store these free surfaces
    datamanager.set_free_contact_surfaces(free_surfaces)
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
                                    synch_vector)
    sy = datamanager.get_field(what, step)
    local_map = datamanager.get_local_contact_ids()

    for (local_id, exchange_id) in pairs(local_map)
        synch_vector[exchange_id, :] .= sy[local_id, :]
    end
    if datamanager.get_max_rank() == 1
        return synch_vector
    end
    comm = datamanager.get_comm()
    find_and_set_core_value_sum(comm, synch_vector)
    return synch_vector
end

function get_local_ids(datamanager::Module)
    local_map = Dict{Int64,Int64}()
    global_nodes = datamanager.get_contact_relevant_global_ids()
    for (id, global_node) in enumerate(global_nodes)
        local_node = datamanager.get_local_nodes([global_node])
        if !(local_node == [])
            local_map[local_node[1]] = id
        end
    end
    datamanager.set_contact_relevant_local_ids(local_map)
end
function get_double_surfs(normals_i, offsets_i, normals_j, offsets_j)
    ids = Vector{Int64}([])
    jds = Vector{Int64}([])
    for iID in eachindex(offsets_i)
        for jID in eachindex(offsets_j)
            if offsets_i[iID] == offsets_j[jID] && normals_i[iID, :] == normals_j[jID, :]
                append!(ids, iID)
                append!(jds, jID)
                break # only one surface pair can exist
            end
        end
    end
    return ids, jds
end

function check_valid_contact_model(params, block_ids::Vector{Int64})
    # inverse master slave check
    # tuple liste bauen und dann die neuen invers checken
    check_dict = Dict{Int64,Int64}()
    for contact_params in values(params)
        if !haskey(contact_params, "Master")
            @error "Contact model needs a ''Master''"
            return false
        end
        if !haskey(contact_params, "Slave")
            @error "Contact model needs a ''Slave''"
            return false
        end
        if contact_params["Master"] == contact_params["Slave"]
            @error "Contact master and slave are equal. Self contact is not implemented yet."
            return false
        end

        if !(contact_params["Master"] in block_ids)
            @error "Block defintion in master does not exist."
            return false
        end
        if !(contact_params["Slave"] in block_ids)
            @error "Block defintion in slave does not exist."
            return false
        end
        check_dict[contact_params["Master"]] = contact_params["Slave"]
        if haskey(check_dict, contact_params["Slave"]) &&
           check_dict[contact_params["Slave"]] == contact_params["Master"]
            @error "Master and Slave should be defined in an inverse way, e.g. Master = 1, Slave = 2 in model 1 and Master = 2, Slave = 1 in model 2."
            return false
        end
        if !haskey(contact_params, "Search Radius")
            @error "Contact model needs a ''Search Radius''."
            return false
        end
        if contact_params["Search Radius"] <= 0
            @error "''Search Radius'' must be greater than zero."
            return false
        end
    end
    return true
end
end
