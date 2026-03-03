# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

export set_distribution
export set_glob_to_loc
export synch_manager
export SynchBuffer
export get_or_create_synch_buffer!

# Cache for pre-allocated MPI send/receive buffers, keyed by field name.
# Buffers are created on first use and reused in subsequent time steps.
const synch_buffer_cache = Dict{String,Any}()

"""
    SynchBuffer{T,N}

Pre-allocated send/receive buffers for a single MPI-synchronised field.
Avoids repeated memory allocation during time stepping by reusing contiguous
buffers that are compatible with `MPI.Send` and `MPI.Recv!`.

# Type parameters
- `T`: Element type (e.g. `Float64`)
- `N`: Array dimension (`1` for vectors, `2` for matrices)

# Fields
- `send_bufs`: Contiguous send buffers per `(icore, jcore)` pair
- `recv_bufs`: Contiguous receive buffers per `(icore, jcore)` pair
- `send_indices`: Controller node indices to copy into the send buffer
- `recv_indices`: Responder node indices to write the receive buffer into
"""
struct SynchBuffer{T,N}
    send_bufs::Dict{Tuple{Int,Int},Array{T,N}}
    recv_bufs::Dict{Tuple{Int,Int},Array{T,N}}
    send_indices::Dict{Tuple{Int,Int},Vector{Int}}
    recv_indices::Dict{Tuple{Int,Int},Vector{Int}}
end

"""
    SynchBuffer(overlapnodes, array::Array{T,N}) where {T,N}

Construct a [`SynchBuffer`](@ref) from the overlap map and a reference array.
Buffer shape and element type are inferred automatically from `array`.

# Arguments
- `overlapnodes`: Overlap map `Dict{Int, Dict{Int, Dict{String, Vector{Int}}}}` defining controller/responder node indices per core pair
- `array`: Reference array used to infer element type `T` and dimension `N`
"""
function SynchBuffer(overlapnodes, array::Array{T,N}) where {T,N}
    send_bufs = Dict{Tuple{Int,Int},Array{T,N}}()
    recv_bufs = Dict{Tuple{Int,Int},Array{T,N}}()
    send_indices = Dict{Tuple{Int,Int},Vector{Int}}()
    recv_indices = Dict{Tuple{Int,Int},Vector{Int}}()

    for (icore, jcore_dict) in overlapnodes
        for (jcore, node_dict) in jcore_dict
            ctrl = node_dict["Controller"]
            resp = node_dict["Responder"]
            if !isempty(ctrl)
                send_indices[(icore, jcore)] = ctrl
                send_bufs[(icore, jcore)] = N == 1 ? similar(array, length(ctrl)) :
                                            similar(array, length(ctrl), size(array, 2))
            end
            if !isempty(resp)
                recv_indices[(icore, jcore)] = resp
                recv_bufs[(icore, jcore)] = N == 1 ? similar(array, length(resp)) :
                                            similar(array, length(resp), size(array, 2))
            end
        end
    end
    return SynchBuffer{T,N}(send_bufs, recv_bufs, send_indices, recv_indices)
end

"""
    fill_send!(buf::SynchBuffer{T,1}, key, vector) where T

Copy controller node values from `vector` into the pre-allocated send buffer.
"""
function fill_send!(buf::SynchBuffer{T,1}, key, vector) where {T}
    buf.send_bufs[key] .= vector[buf.send_indices[key]]
end

"""
    fill_send!(buf::SynchBuffer{T,2}, key, vector) where T

Copy controller node rows from `vector` into the pre-allocated send buffer.
"""
function fill_send!(buf::SynchBuffer{T,2}, key, vector) where {T}
    buf.send_bufs[key] .= vector[buf.send_indices[key], :]
end

"""
    flush_recv!(buf::SynchBuffer{T,1}, key, vector) where T

Write the receive buffer back into the responder node positions of `vector`.
"""
function flush_recv!(buf::SynchBuffer{T,1}, key, vector) where {T}
    vector[buf.recv_indices[key]] .= buf.recv_bufs[key]
end

"""
    flush_recv!(buf::SynchBuffer{T,2}, key, vector) where T

Write the receive buffer rows back into the responder node positions of `vector`.
"""
function flush_recv!(buf::SynchBuffer{T,2}, key, vector) where {T}
    vector[buf.recv_indices[key], :] .= buf.recv_bufs[key]
end

"""
    get_or_create_synch_buffer!(cache, field_name, vector, overlapnodes)

Return a cached [`SynchBuffer`](@ref) for `field_name`, creating a new one if
none exists or if the cached buffer is no longer compatible with `vector`.

Compatibility is checked by element type `T`, dimension `N`, and — for matrices
— the number of columns. If incompatible, the buffer is rebuilt and the cache
is updated.

# Arguments
- `cache::Dict`: Buffer cache, typically `synch_buffer_cache`
- `field_name::String`: Unique field identifier used as cache key
- `vector`: Current field array used to check compatibility and build new buffers
- `overlapnodes`: Overlap map passed to [`SynchBuffer`](@ref) constructor if rebuild is needed
"""
function get_or_create_synch_buffer!(cache::Dict, field_name::String,
                                     vector::Array{T,N}, overlapnodes) where {T,N}
    if haskey(cache, field_name)
        buf = cache[field_name]
        if buf isa SynchBuffer{T,N} &&
           (N == 1 || isempty(buf.send_bufs) ||
            size(first(values(buf.send_bufs)), 2) == size(vector, 2))
            return buf
        end
        @debug "Rebuilding SynchBuffer for field '$field_name' (type or size changed)"
    end
    cache[field_name] = SynchBuffer(overlapnodes, vector)
    return cache[field_name]
end

"""
    set_distribution(values::Vector{Int64})

Sets the distribution globally.

# Arguments
- `values::Vector{Int64}`: The distribution.
"""
function set_distribution(values::Vector{Int64})
    data["distribution"] = values
end

"""
    set_glob_to_loc(dict::Dict{Int64,Int64})

Sets the global-to-local mapping dict globally.

# Arguments
- `dict::Dict{Int64,Int64}`: The dict representing the global-to-local mapping.
"""
function set_glob_to_loc(dict::Dict{Int64,Int64})
    data["glob_to_loc"] = dict
end

"""
    synch_manager(synchronise_field, direction::String)

Synchronises the fields.

# Arguments
- `synchronise_field`: The function to synchronise the field.
- `direction::String`: The direction of the synchronisation.
"""
function synch_manager(synchronise_field::Function, direction::String)
    synch_fields::Dict{String,Any} = get_synch_fields()
    for synch_field in keys(synch_fields)
        synchronise_field(get_comm(),
                          synch_fields,
                          get_overlap_map(),
                          get_field,
                          synch_field,
                          direction,
                          synch_buffer_cache)
    end
end
