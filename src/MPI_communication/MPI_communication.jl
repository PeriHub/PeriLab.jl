# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause
module MPI_Communication
import MPI
using ...PeriLabExceptions: @abort
export send_single_value_from_vector
export synch_responder_to_controller
export synch_controller_to_responder
export synch_controller_bonds_to_responder
export send_vector_from_root_to_core_i
export broadcast_value
export find_and_set_core_value_min
export find_and_set_core_value_sum
export find_and_set_core_value_avg
export gather_values

"""
TODO
Contact
send all information to first core and synch to all otherwise
optimization is possible by reducing it to slave and master. Therefore its only the surface.

Master are known and their core. local to global is know and the sending can occur.

"""

"""
    send_single_value_from_vector(comm::MPI.Comm, controller::Int64, values::Union{Int64,Vector{Float64},Vector{Int64},Vector{Bool}}, type::Type)

Sends a single value from a vector to a controller

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `controller::Int64`: The controller
- `values::Union{Int64,Vector{Float64},Vector{Int64},Vector{Bool}}`: The values
- `type::Type`: The type
# Returns
- `recv_msg::Union{Int64,Vector{Float64},Vector{Int64},Vector{Bool}}`: The received message
"""
function send_single_value_from_vector(comm::MPI.Comm,
                                       controller::Int64,
                                       values::Union{Int64,Vector{Float64},Vector{Int64},
                                                     Vector{Bool}},
                                       type::Type)
    ncores = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)
    requests = Vector{MPI.Request}()
    if type == String
        @abort "Wrong type - String in function send_single_value_from_vector"
        return nothing
    end
    recv_msg = zeros(type, 1, 1)
    if rank == controller
        send_msg = zeros(type, 1, 1)
        for i in 0:(ncores - 1)
            # +1 because the index of cores is zero based and julia matrices are one based
            send_msg[1] = values[i + 1]
            if i != controller
                push!(requests, MPI.Isend(send_msg, comm; dest = i, tag = 0))
                # @debug "Sending   $rank -> $i"
            else
                recv_msg[1] = send_msg[1]
            end
        end

    else
        MPI.Recv!(recv_msg, comm; source = controller, tag = 0)
        # @debug "Receiving $controller -> $rank"
    end
    if !isempty(requests)
        MPI.Waitall(requests)
    end
    MPI.Barrier(comm)
    return recv_msg[1]
end

# --------------------------------------------------------------------------------------
# Exchange of overlap data between the ranks
#
# Every exchange packs the entries of one field for each neighbouring rank into a single
# contiguous buffer, posts all receives and sends non-blocking and waits for all of them,
# then unpacks. The buffers are kept per element type and neighbouring rank, grown when a
# field needs more and reused by every call; since a call completes all its transfers
# before it returns, the fields can share them. After the first calls an exchange
# therefore allocates nothing of the size of the data.
# --------------------------------------------------------------------------------------

const SEND_BUFFERS = Dict{Tuple{DataType,Int},Vector}()
const RECV_BUFFERS = Dict{Tuple{DataType,Int},Vector}()
# created on first use, MPI handles must not be built at precompile time
const REQUESTS = Ref{Union{Nothing,MPI.MultiRequest}}(nothing)

"""
    exchange_buffer(buffers, T, jcore, len)

Leading `len` entries of the cached buffer for element type `T` and rank `jcore`, grown
if it is shorter.
"""
function exchange_buffer(buffers::Dict{Tuple{DataType,Int},Vector}, ::Type{T},
                         jcore::Int, len::Int) where {T}
    buffer = get!(Vector{T}, buffers, (T, jcore))::Vector{T}
    length(buffer) < len && resize!(buffer, len)
    return view(buffer, 1:len)
end

"""
    exchange_requests(n)

The cached request set, enlarged to at least `n` requests.
"""
function exchange_requests(n::Int)
    requests = REQUESTS[]
    if isnothing(requests) || length(requests) < n
        requests = MPI.MultiRequest(n)
        REQUESTS[] = requests
    end
    return requests
end

number_type(::Type{T}) where {T<:Number} = T
number_type(::Type{<:AbstractArray{S}}) where {S} = number_type(S)

# --- node fields: an array whose first dimension runs over the nodes ---------------------

node_width(field::AbstractArray) = size(field, 1) == 0 ? 0 : length(field) ÷ size(field, 1)

function pack_nodes!(buffer::AbstractVector, field::AbstractArray,
                     nodes::AbstractVector{Int64})
    n = size(field, 1)
    k = 0
    @inbounds for c in 1:node_width(field), i in nodes
        k += 1
        buffer[k] = field[i + (c - 1) * n]
    end
    return buffer
end

function unpack_nodes!(field::AbstractArray, buffer::AbstractVector,
                       nodes::AbstractVector{Int64}, add::Bool)
    # flags are not summed: a node that is true on its own rank stays true
    (add && eltype(field) == Bool) && return field
    n = size(field, 1)
    k = 0
    @inbounds for c in 1:node_width(field), i in nodes
        k += 1
        if add
            field[i + (c - 1) * n] += buffer[k]
        else
            field[i + (c - 1) * n] = buffer[k]
        end
    end
    return field
end

function node_count(field::AbstractArray, nodes::AbstractVector{Int64})
    node_width(field) *
    length(nodes)
end

# --- bond fields: one array of bond values per node -----------------------------------

bond_count(bonds::AbstractArray{<:Number}) = length(bonds)
bond_count(bonds::AbstractVector{<:AbstractArray}) = sum(length, bonds; init = 0)

function bond_count(field::AbstractVector, nodes::AbstractVector{Int64})
    total = 0
    @inbounds for i in nodes
        total += bond_count(field[i])
    end
    return total
end

function pack_bonds!(buffer::AbstractVector, k::Int, bonds::AbstractArray{<:Number})
    @inbounds for value in bonds
        k += 1
        buffer[k] = value
    end
    return k
end

function pack_bonds!(buffer::AbstractVector, k::Int, bonds::AbstractVector{<:AbstractArray})
    @inbounds for bond in bonds
        k = pack_bonds!(buffer, k, bond)
    end
    return k
end

function pack_bonds!(buffer::AbstractVector, field::AbstractVector,
                     nodes::AbstractVector{Int64})
    k = 0
    @inbounds for i in nodes
        k = pack_bonds!(buffer, k, field[i])
    end
    return buffer
end

function unpack_bonds!(bonds::AbstractArray{<:Number}, buffer::AbstractVector, k::Int)
    @inbounds for index in eachindex(bonds)
        k += 1
        bonds[index] = buffer[k]
    end
    return k
end

function unpack_bonds!(bonds::AbstractVector{<:AbstractArray}, buffer::AbstractVector,
                       k::Int)
    @inbounds for bond in bonds
        k = unpack_bonds!(bond, buffer, k)
    end
    return k
end

function unpack_bonds!(field::AbstractVector, buffer::AbstractVector,
                       nodes::AbstractVector{Int64}, ::Bool)
    k = 0
    @inbounds for i in nodes
        k = unpack_bonds!(field[i], buffer, k)
    end
    return field
end

"""
    exchange!(comm, overlapnodes, field, send_role, recv_role, add, count, pack!, unpack!)

Sends the entries of `field` at the nodes of role `send_role` to every neighbouring rank
and receives into the nodes of role `recv_role`, adding the received values if `add`,
replacing them otherwise. `count`, `pack!` and `unpack!` define the layout of the field.
"""
function exchange!(comm::MPI.Comm, overlapnodes, field, send_role::String,
                   recv_role::String, add::Bool, count::F1, pack!::F2,
                   unpack!::F3) where {F1,F2,F3}
    ncores = MPI.Comm_size(comm)
    ncores == 1 && return field
    rank = MPI.Comm_rank(comm)
    T = number_type(typeof(field))
    neighbours = overlapnodes[rank + 1]
    requests = exchange_requests(2 * ncores)

    r = 0
    for jcore in 1:ncores
        jcore == rank + 1 && continue
        nodes = neighbours[jcore][recv_role]
        isempty(nodes) && continue
        buffer = exchange_buffer(RECV_BUFFERS, T, jcore, count(field, nodes))
        r += 1
        MPI.Irecv!(buffer, comm, requests[r]; source = jcore - 1, tag = 0)
    end
    for jcore in 1:ncores
        jcore == rank + 1 && continue
        nodes = neighbours[jcore][send_role]
        isempty(nodes) && continue
        buffer = exchange_buffer(SEND_BUFFERS, T, jcore, count(field, nodes))
        pack!(buffer, field, nodes)
        r += 1
        MPI.Isend(buffer, comm, requests[r]; dest = jcore - 1, tag = 0)
    end
    MPI.Waitall(requests)

    for jcore in 1:ncores
        jcore == rank + 1 && continue
        nodes = neighbours[jcore][recv_role]
        isempty(nodes) && continue
        buffer = exchange_buffer(RECV_BUFFERS, T, jcore, count(field, nodes))
        unpack!(field, buffer, nodes, add)
    end
    return field
end

"""
    synch_responder_to_controller(comm::MPI.Comm, overlapnodes, vector, dof)

Adds the values of the responder copies of every node to its controller, e.g. the force
densities a rank computed for nodes it does not own. Boolean fields are not summed.

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `overlapnodes::AbstractDict`: The overlap nodes
- `vector::AbstractArray`: Node field, first dimension over the nodes
- `dof::Int`: The degrees of freedom; the layout is taken from `vector` itself
# Returns
- `vector::AbstractArray`: The field
"""
function synch_responder_to_controller(comm::MPI.Comm, overlapnodes, vector, dof)
    return exchange!(comm, overlapnodes, vector, "Responder", "Controller", true,
                     node_count, pack_nodes!, unpack_nodes!)
end

"""
    synch_controller_to_responder(comm::MPI.Comm, overlapnodes, vector, dof)

Copies the values of every controller node to its responder copies on the other ranks.

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `overlapnodes::AbstractDict`: The overlap nodes
- `vector::AbstractArray`: Node field, first dimension over the nodes
- `dof::Int`: The degrees of freedom; the layout is taken from `vector` itself
# Returns
- `vector::AbstractArray`: The field
"""
function synch_controller_to_responder(comm::MPI.Comm, overlapnodes, vector, dof)
    return exchange!(comm, overlapnodes, vector, "Controller", "Responder", false,
                     node_count, pack_nodes!, unpack_nodes!)
end

"""
    synch_controller_bonds_to_responder(comm::MPI.Comm, overlapnodes, array, dof)

Copies the bond values of every controller node to its responder copies on the other
ranks. Every node has the same bonds on all ranks, so the values are written into the
existing bond arrays of the responders.

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `overlapnodes::AbstractDict`: The overlap nodes
- `array::AbstractVector`: Bond field, one array of bond values per node
- `dof::Int`: The degrees of freedom; the layout is taken from `array` itself
# Returns
- `array::AbstractVector`: The field
"""
function synch_controller_bonds_to_responder(comm::MPI.Comm, overlapnodes, array, dof)
    return exchange!(comm, overlapnodes, array, "Controller", "Responder", false,
                     bond_count, pack_bonds!, unpack_bonds!)
end

"""
    send_vector_from_root_to_core_i(comm::MPI.Comm, send_msg, recv_msg, distribution)

Sends a vector from the root to the core i

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `send_msg::Union{Int64,Vector{Float64},Vector{Int64},Vector{Bool}}`: The send message
- `recv_msg::Union{Int64,Vector{Float64},Vector{Int64},Vector{Bool}}`: The receive message
- `distribution::Vector{Int64}`: The distribution
# Returns
- `recv_msg::Union{Int64,Vector{Float64},Vector{Int64},Vector{Bool}}`: The received message
"""
function send_vector_from_root_to_core_i(comm::MPI.Comm, send_msg, recv_msg, distribution)
    currentRank = MPI.Comm_rank(comm)
    # MPI.Barrier(comm)
    requests = Vector{MPI.Request}()
    if currentRank == 0
        for rank in 1:(MPI.Comm_size(comm) - 1)
            push!(requests,
                  MPI.Isend(send_msg[distribution[rank + 1]], comm; dest = rank, tag = 0))
            # @debug "Sending $currentRank -> $rank"
        end
        recv_msg .= send_msg[distribution[1]]
    else
        MPI.Recv!(recv_msg, comm; source = 0, tag = 0)
        # @debug "Receiving 0 -> $currentRank"
    end
    if !isempty(requests)
        MPI.Waitall(requests)
    end
    MPI.Barrier(comm)
    return recv_msg
end

"""
    broadcast_value(comm::MPI.Comm, send_msg)

Broadcast a value to all ranks

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `send_msg::Union{Int64,Vector{Float64},Vector{Int64},Vector{Bool}}`: The send message
# Returns
- `recv_msg::Union{Int64,Vector{Float64},Vector{Int64},Vector{Bool}}`: The received message
"""
function broadcast_value(comm::MPI.Comm,
                         send_msg::T) where {T<:Union{Float64,Int64,
                                                      Vector{Float64},
                                                      Vector{Vector{Float64}},
                                                      Vector{Int64},Vector{Vector{Int64}},
                                                      Matrix{Float64},
                                                      Matrix{Int64},
                                                      Dict,
                                                      Bool,
                                                      Nothing,
                                                      Any}}
    # recv_msg = MPI.Comm_rank(comm) == controller ? send_msg : nothing
    # recv_msg = MPI.bcast(send_msg, controller, comm)
    # return recv_msg
    return MPI.bcast(send_msg, 0, comm)
end

"""
    find_and_set_core_value_sum(comm::MPI.Comm, value::Union{Float64,Int64})

Find and set core value sum

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `value::Union{Float64,Int64}`: The value
# Returns
- `recv_msg::Union{Int64,Vector{Float64},Vector{Int64},Vector{Bool}}`: The received message
"""
function find_and_set_core_value_sum(comm::MPI.Comm,
                                     value::T) where {T<:Union{Float64,Int64,
                                                               Vector{Float64},
                                                               Vector{Int64},
                                                               Matrix{Float64},
                                                               Matrix{Int64}}}
    return MPI.Allreduce(value, MPI.SUM, comm)
end

"""
    find_and_set_core_value_max(comm::MPI.Comm, value::Union{Float64,Int64})

Find and set core value max

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `value::Union{Float64,Int64}`: The value
# Returns
- `recv_msg::Union{Int64,Vector{Float64},Vector{Int64},Vector{Bool}}`: The received message
"""
function find_and_set_core_value_max(comm::MPI.Comm,
                                     value::T) where {T<:Union{Float64,Int64}}
    return MPI.Allreduce(value, MPI.MAX, comm)
end

"""
    find_and_set_core_value_min(comm::MPI.Comm, value::Union{Float64,Int64})

Find and set core value sum

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `value::Union{Float64,Int64}`: The value
# Returns
- `recv_msg::Union{Int64,Vector{Float64},Vector{Int64},Vector{Bool}}`: The received message
"""
function find_and_set_core_value_min(comm::MPI.Comm,
                                     value::T) where {T<:Union{Float64,Int64,
                                                               Vector{Float64},
                                                               Vector{Int64},
                                                               Matrix{Float64},
                                                               Matrix{Int64}}}
    return MPI.Allreduce(value, MPI.MIN, comm)
end

"""
    find_and_set_core_value_avg(comm::MPI.Comm,
                                     value::T,
                                     nnodes::Int64) where {T<:Union{Float64,Int64}}

Find and set core value avg

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `value::Union{Float64,Int64}`: The value
# Returns
- `recv_msg::Float64`: The received a Float64 message
"""
function find_and_set_core_value_avg(comm::MPI.Comm,
                                     value::T,
                                     nnodes::Int64) where {T<:Union{Float64,Int64}}
    return MPI.Allreduce(value, MPI.SUM, comm) / MPI.Allreduce(nnodes, MPI.SUM, comm)
end

"""
    gather_values(comm::MPI.Comm, value::Any)

Gather values

# Arguments
- `comm::MPI.Comm`: The MPI communicator
- `value::Any`: The value
# Returns
- `recv_msg::Any`: The received message
"""
function gather_values(comm::MPI.Comm,
                       value::T) where {T<:Union{Float64,Int64,
                                                 Vector{Float64},
                                                 Vector{Int64},
                                                 Matrix{Float64},
                                                 Matrix{Int64},
                                                 Any}}
    return MPI.gather(value, comm; root = 0)
end

end
