# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

export get_all_blocks
export get_all_positions
export get_contact_overlap_map
export get_contact_block_ids
export get_contact_relevant_global_ids
export get_local_contact_ids
export get_exchange_id_to_local_id

export get_free_contact_nodes
export get_contact_properties
export set_all_blocks
export set_contact_block_ids
export set_all_positions
export set_contact_overlap_map
export set_contact_relevant_global_ids
export set_local_contact_ids
export set_exchange_id_to_local_id
export set_free_contact_nodes
export set_contact_properties
"""
    get_all_blocks()

Gives back a global list of all block, initially. Is reduced to the outer (not free) surfaces of the contact blocks.

"""
function get_all_blocks()
    return data["All Blocks"]
end
"""
    get_all_positions()

Gives back a global list of current positions, initially. Is reduced to the outer (not free) surfaces of the contact blocks.

"""
function get_all_positions()
    return data["All positions"]
end

"""
    get_contact_overlap_map()

Get the contact overlap map
"""
function get_contact_overlap_map()
    return data["contact_overlap_map"]
end

function get_contact_block_ids(block::Int64)
    return data["Contact block IDs"][block]
end

function get_contact_relevant_global_ids()
    return data["Global Contact IDs"]
end

function get_local_contact_ids()
    return data["Local Contact IDs"]
end

function get_free_contact_nodes(block::Int64)
    return data["Free Surface Nodes"][block]
end
function get_contact_properties()
    return data["Contact Properties"]
end

function set_all_blocks(all_blocks)
    data["All Blocks"] = all_blocks
end

function set_all_positions(all_coordinate)
    data["All positions"] = all_coordinate
end

function set_search_step(cg, n)
    data["Contact Search Step"][cg] = n
end

function get_search_step(cg)
    return data["Contact Search Step"][cg]
end

function clear_synchronization_list()
    data["Synchronization List"] = Dict{Int64,Int64}()
end

function add_synchronization_list(list::Vector{Int64})
    mapping = get_local_contact_ids()
    for iID in list
        if !isnothing(get(mapping, iID, nothing))
            data["Synchronization List"][mapping[iID]] = iID
        end
    end
end
function set_global_search_master_nodes(cg, global_master)
    data["Global Master Search Nodes"][cg] = global_master
end

function get_global_search_master_nodes(cg)
    return data["Global Master Search Nodes"][cg]
end
function set_global_search_slave_nodes(cg, global_slave)
    data["Global Slave Search Nodes"][cg] = global_slave
end

function get_global_search_slave_nodes(cg)
    return data["Global Slave Search Nodes"][cg]
end

function get_synchronization_list()
    return data["Synchronization List"]
end

function global_synchronization()
    return 0 in values(data["Contact Search Step"])
end

function set_no_pairs_flag(cg, n::Bool)
    data["Contact Search No Pairs"][cg] = n
end

function get_no_pairs_flag(cg)
    return data["Contact Search No Pairs"][cg]
end

"""
    set_contact_overlap_map(topo)

Sets the contact overlap map globally.

# Arguments
- `topo`: The overlap map.
"""
function set_contact_overlap_map(topo)
    data["contact_overlap_map"] = topo
end

function set_contact_block_ids(mapping::Dict{Int64,Vector{Int64}})
    data["Contact block IDs"] = mapping
end

function set_contact_relevant_global_ids(ids::Vector{Int64})
    data["Global Contact IDs"] = append!(data["Global Contact IDs"], ids)
end

function set_local_contact_ids(ids::Dict{Int64,Int64})
    data["Local Contact IDs"] = merge(data["Local Contact IDs"], ids)
end

function set_exchange_id_to_local_id(ids::Dict{Int64,Int64})
    data["Exchange id to local id"] = ids
end

function get_exchange_id_to_local_id()
    return data["Exchange id to local id"]
end

function set_free_contact_nodes(block::Int64, free_surface_nodes::Vector{Int64})
    data["Free Surface Nodes"][block] = free_surface_nodes
end

function set_contact_properties(params::Dict)
    data["Contact Properties"] = params
end

function get_contact_dict(id::String)
    return data["Contact Dictionary"][id]
end

function set_contact_dict(id::String, params::Dict)
    data["Contact Dictionary"][id] = params
end
